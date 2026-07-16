#!/usr/bin/env python3
"""Export ViroForge collection metadata as static JSON for the web explorer.

The interactive educational app is driven by the *model* layer (per-genome
metadata plus per-collection abundances) and never needs the genome sequences,
which are ~500 MB of the 525 MB database. This utility runs the same
collection join the generator uses (scripts/generate_fastq_dataset.py), drops
the sequence column, and writes a small JSON bundle (~1 MB for all 20
collections) that ships as a static asset. Sliders then run client-side math
over it with no backend.

Two outputs, both carrying a provenance stamp (ViroForge version, git SHA,
export date, schema version) so the portal can pin to a versioned snapshot
rather than track ViroForge's moving main:

  collections.json      every collection with its genome table (no sequences)
  model_constants.json  the ViroForge model constants for the client-side port

Provenance and drift safety. ViroForge is under active development, so the
model constants are handled in two ways:
  - Source-derived: constants exposed as stable public config (VLP protocols,
    virion-size relationships, amplification and reverse-transcription defaults)
    are read straight from the imported modules, so a parameter change flows
    through automatically on the next export.
  - Drift-guarded: values that live inside method bodies of modules under
    active refactoring (contamination presets, per-type removal, RT ranges) are
    kept as documented literals rather than introspected out of churning code.
    provenance.source_sha256 records a hash of each source module, so a consumer
    can detect when a drift-guarded literal's source changed and re-verify.

Collection metadata reads use only the standard library, so that output works
even without the ViroForge environment. Full source derivation of the model
constants needs the ViroForge deps; run under them for a clean derivation:

  uv run --with numpy,pandas,biopython,scipy,pyyaml,rich \
      python scripts/export_web_data.py

Usage:
  python3 scripts/export_web_data.py [--db PATH] [--out-dir DIR] [--pretty]
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime
import hashlib
import inspect
import json
import re
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DB = REPO_ROOT / "viroforge" / "data" / "viral_genomes.db"
DEFAULT_OUT = REPO_ROOT / "data" / "web"

# Running the script by path puts scripts/ on sys.path, not the repo root, so
# `import viroforge` would miss. Put the repo root first so source derivation of
# the model constants works when the ViroForge deps are installed. Harmless for
# the stdlib-only collections export.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Bump when the output JSON shape changes in a way the portal must adapt to.
SCHEMA_VERSION = "1"

# Source modules whose content is hashed into provenance for drift detection.
SOURCE_FILES = {
    "vlp": REPO_ROOT / "viroforge" / "enrichment" / "vlp.py",
    "amplification": REPO_ROOT / "viroforge" / "amplification.py",
    "rna_virome": REPO_ROOT / "viroforge" / "workflows" / "rna_virome.py",
    "contamination": REPO_ROOT / "viroforge" / "core" / "contamination.py",
}

# Precision for float rounding, chosen to keep the JSON small without losing
# the small-abundance tail (values run down to ~1e-5).
GC_DIGITS = 5
ABUNDANCE_DIGITS = 8
PREVALENCE_DIGITS = 5


def _round(value: Any, digits: int) -> Any:
    return round(value, digits) if isinstance(value, (int, float)) else value


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #

def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO_ROOT), *args],
            capture_output=True, text=True, timeout=5,
        )
        return out.stdout.strip() or None
    except Exception:
        return None


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except Exception:
        return None


def _viroforge_version() -> str:
    try:
        import viroforge  # noqa: PLC0415 - optional, degrades to file parse

        return str(getattr(viroforge, "__version__", "unknown"))
    except Exception:
        try:
            text = (REPO_ROOT / "viroforge" / "__init__.py").read_text()
            m = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', text)
            return m.group(1) if m else "unknown"
        except Exception:
            return "unknown"


def provenance() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "viroforge_version": _viroforge_version(),
        "viroforge_git_sha": _git("rev-parse", "--short", "HEAD"),
        "viroforge_git_dirty": bool(_git("status", "--porcelain")),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "source_sha256": {k: _sha256(p) for k, p in SOURCE_FILES.items()},
        "note": (
            "Regenerate on each ViroForge release and copy into the portal. "
            "source_sha256 lets consumers detect when a drift-guarded literal's "
            "source module changed and needs re-verification."
        ),
    }


# --------------------------------------------------------------------------- #
# Collections (stdlib only)
# --------------------------------------------------------------------------- #

def export_collections(db_path: Path) -> dict[str, Any]:
    """Read every collection and its genome table (no sequences) from the DB."""
    if not db_path.exists():
        raise FileNotFoundError(f"ViroForge database not found: {db_path}")

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    try:
        collection_rows = conn.execute(
            """
            SELECT collection_id, collection_name, description, n_genomes,
                   selection_criteria, literature_references, version
            FROM body_site_collections
            ORDER BY collection_id
            """
        ).fetchall()

        collections: list[dict[str, Any]] = []
        total_genomes = 0

        for c in collection_rows:
            genome_rows = conn.execute(
                """
                SELECT
                    cg.genome_id,
                    cg.relative_abundance,
                    cg.abundance_rank,
                    cg.prevalence,
                    g.genome_name,
                    g.length,
                    g.gc_content,
                    g.genome_type,
                    g.genome_structure,
                    g.n_segments,
                    t.realm, t.kingdom, t.phylum, t.class, t.order_name,
                    t.family, t.subfamily, t.genus, t.species
                FROM collection_genomes cg
                JOIN genomes g ON cg.genome_id = g.genome_id
                LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
                WHERE cg.collection_id = ?
                ORDER BY cg.relative_abundance DESC
                """,
                (c["collection_id"],),
            ).fetchall()

            genomes = [
                {
                    "genome_id": r["genome_id"],
                    "name": r["genome_name"],
                    "length": r["length"],
                    "gc_content": _round(r["gc_content"], GC_DIGITS),
                    "genome_type": r["genome_type"],
                    "genome_structure": r["genome_structure"],
                    "n_segments": r["n_segments"],
                    "relative_abundance": _round(r["relative_abundance"], ABUNDANCE_DIGITS),
                    "abundance_rank": r["abundance_rank"],
                    "prevalence": _round(r["prevalence"], PREVALENCE_DIGITS),
                    "taxonomy": {
                        "realm": r["realm"],
                        "kingdom": r["kingdom"],
                        "phylum": r["phylum"],
                        "class": r["class"],
                        "order": r["order_name"],
                        "family": r["family"],
                        "subfamily": r["subfamily"],
                        "genus": r["genus"],
                        "species": r["species"],
                    },
                }
                for r in genome_rows
            ]
            total_genomes += len(genomes)

            collections.append(
                {
                    "id": c["collection_id"],
                    "name": c["collection_name"],
                    "description": c["description"],
                    "n_genomes_declared": c["n_genomes"],
                    "n_genomes_loaded": len(genomes),
                    "selection_criteria": c["selection_criteria"],
                    "literature_references": c["literature_references"],
                    "version": c["version"],
                    "genomes": genomes,
                }
            )
    finally:
        conn.close()

    return {
        "provenance": provenance(),
        "source": {
            "database": str(db_path),
            "n_collections": len(collections),
            "n_collection_genomes": total_genomes,
            "note": "Per-genome metadata only; genome sequences are omitted by design.",
        },
        "collections": collections,
    }


# --------------------------------------------------------------------------- #
# Model constants
# --------------------------------------------------------------------------- #

def _literal_constants() -> dict[str, Any]:
    """The drift-guarded literal blocks (always present).

    These live inside method bodies of modules under active development; they
    are transcribed here rather than introspected out of churning code. Any
    change is flagged by provenance.source_sha256 for the named module.
    """
    return {
        "virion_size": {
            "derived": False,
            "guard": "vlp",
            "formula": "diameter_nm = a * log10(length_bp) + b, clamped [15, 500]",
            "relationships": {
                "dsDNA": {"a": 35, "b": -70, "variance": 15},
                "ssDNA": {"a": 8, "b": -5, "variance": 5},
                "dsRNA": {"a": 25, "b": -50, "variance": 10},
                "ssRNA": {"a": 15, "b": -20, "variance": 8},
                "ssRNA-RT": {"a": 30, "b": -60, "variance": 12},
                "unknown": {"a": 25, "b": -45, "variance": 20},
            },
        },
        "vlp": {
            "derived": False,
            "guard": "vlp",
            "retention": "1 / (1 + exp(-steepness * (diameter_nm - pore_size_um * 1000)))",
            "protocols": {
                "tangential_flow_standard": {"pore_size_um": 0.2, "retention_curve_type": "sigmoid", "retention_curve_steepness": 0.008, "nuclease_efficiency": 0.98, "recovery_rate": 0.85, "contamination_reduction": 0.95},
                "syringe_filter_standard": {"pore_size_um": 0.2, "retention_curve_type": "sigmoid", "retention_curve_steepness": 0.010, "nuclease_efficiency": 0.90, "recovery_rate": 0.60, "contamination_reduction": 0.85},
                "ultracentrifugation": {"pore_size_um": None, "retention_curve_type": "step", "nuclease_efficiency": 0.95, "recovery_rate": 0.90, "contamination_reduction": 0.75},
                "norgen_kit": {"pore_size_um": None, "retention_curve_type": "sigmoid", "retention_curve_steepness": 0.015, "nuclease_efficiency": 0.92, "recovery_rate": 0.70, "contamination_reduction": 0.88},
                "no_vlp": {"pore_size_um": None, "retention_curve_type": "step", "nuclease_efficiency": 0.0, "recovery_rate": 1.0, "contamination_reduction": 0.0},
            },
        },
        "vlp_per_type_removal": {
            "derived": False,
            "guard": "vlp",
            "host_dna": "nuclease_efficiency, clipped [0.85, 0.99]",
            "rrna": "nuclease_efficiency * 0.92, clipped [0.80, 0.97]",
            "reagent_bacteria": "contamination_reduction * 0.98, clipped [0.70, 0.99]",
            "phix": "filtration retention of a 27 nm virus, clipped [0.60, 0.95]",
        },
        "contamination_dna_presets_pct": {
            "derived": False,
            "guard": "contamination",
            "keys": ["host_dna", "rrna", "reagent_bacteria", "phix"],
            "clean": [0.1, 0.5, 0.01, 0.1],
            "realistic": [2.0, 5.0, 0.5, 0.1],
            "heavy": [10.0, 15.0, 2.0, 0.1],
            "failed": [15.0, 20.0, 5.0, 0.1],
        },
        "contamination_rna_presets_pct": {
            "derived": False,
            "guard": "contamination",
            "keys": ["host_rna_before", "host_rna_after", "bacterial_rna", "reagent", "phix"],
            "clean": [90, 5, 2.0, 0.1, 0.1],
            "realistic": [90, 10, 5.0, 0.5, 0.1],
            "heavy": [95, 20, 10.0, 2.0, 0.1],
            "failed": [95, 90, 5.0, 1.0, 0.1],
        },
        "amplification": {
            "derived": False,
            "guard": "amplification",
            "length_efficiency": "exp(-0.015 * length_bias_strength * length_kb)",
            "gc_efficiency": "exp(-((abs(gc - optimal_gc) / gc_tolerance) ** 2) * gc_bias_strength)",
            "factor": "(length_efficiency * gc_efficiency) ** cycles",
        },
        "rna_workflow": {
            "derived": False,
            "guard": "rna_virome",
            "rt_by_virus_type": {"ssRNA+": [0.70, 0.90], "ssRNA-": [0.50, 0.70], "dsRNA": [0.40, 0.80], "ssRNA-RT": [0.75, 0.85]},
            "primer_modifiers": {"hexamer": 1.0, "octamer": 1.05, "oligo_dt": 0.9, "specific": 1.1},
            "ribo_depletion": {"ribo_zero": [0.90, 0.95], "ribominus": [0.85, 0.90]},
        },
        "coverage": {
            "derived": "structural",
            "expected_completeness": "1 - exp(-coverage)",
            "note": "Structural formula, not a tunable constant.",
        },
    }


def _overlay_derived(constants: dict[str, Any]) -> str:
    """Overwrite derivable blocks with values read from imported ViroForge.

    Returns a short status string. On any import failure the literal blocks
    (already in `constants`) stand, so this never raises.
    """
    try:
        from viroforge.enrichment.vlp import VLPProtocol, VirionSizeEstimator  # noqa: PLC0415
    except Exception as exc:
        return f"literals-only (viroforge import failed: {exc})"

    # Virion size + VLP protocols (stable public config).
    constants["virion_size"] = {
        "derived": True,
        "source": "viroforge.enrichment.vlp.VirionSizeEstimator.SIZE_RELATIONSHIPS",
        "formula": "diameter_nm = a * log10(length_bp) + b, clamped [15, 500]",
        "relationships": VirionSizeEstimator.SIZE_RELATIONSHIPS,
    }
    protocols: dict[str, Any] = {}
    for method in ("tangential_flow_standard", "syringe_filter_standard", "ultracentrifugation", "norgen_kit", "no_vlp"):
        cfg = getattr(VLPProtocol, method)()
        protocols[method] = {
            "name": cfg.name,
            "filtration_method": cfg.filtration_method,
            "pore_size_um": cfg.pore_size_um,
            "retention_curve_type": cfg.retention_curve_type,
            "retention_curve_steepness": cfg.retention_curve_steepness,
            "nuclease_efficiency": cfg.nuclease_efficiency,
            "recovery_rate": cfg.recovery_rate,
            "contamination_reduction": cfg.contamination_reduction,
        }
    constants["vlp"] = {
        "derived": True,
        "source": "viroforge.enrichment.vlp.VLPProtocol",
        "retention": "1 / (1 + exp(-steepness * (diameter_nm - pore_size_um * 1000)))",
        "protocols": protocols,
    }

    # Amplification defaults (constructor defaults are stable public config).
    try:
        from viroforge.amplification import RdABAmplification  # noqa: PLC0415

        sig = inspect.signature(RdABAmplification.__init__)
        defaults = {
            name: sig.parameters[name].default
            for name in ("cycles", "optimal_gc", "gc_tolerance", "length_bias_strength", "gc_bias_strength")
            if name in sig.parameters and sig.parameters[name].default is not inspect.Parameter.empty
        }
        constants["amplification"] = {
            **constants["amplification"],
            "derived": "partial",
            "source": "viroforge.amplification.RdABAmplification (constructor defaults)",
            "rdab_defaults": defaults,
            "note": "The 0.015 length-decay coefficient is a literal in _length_efficiency(); drift-guarded via provenance.source_sha256.amplification.",
        }
    except Exception:
        pass

    # Reverse-transcription defaults (dataclass field defaults).
    try:
        from viroforge.workflows.rna_virome import ReverseTranscription  # noqa: PLC0415

        rt_defaults: dict[str, Any] = {}
        if dataclasses.is_dataclass(ReverseTranscription):
            for f in dataclasses.fields(ReverseTranscription):
                if f.default is not dataclasses.MISSING and isinstance(f.default, (int, float, str, bool)):
                    rt_defaults[f.name] = f.default
        if rt_defaults:
            constants["rna_workflow"] = {
                **constants["rna_workflow"],
                "derived": "partial",
                "source": "viroforge.workflows.rna_virome.ReverseTranscription (field defaults)",
                "rt_defaults": rt_defaults,
                "note": "Per-type efficiency ranges and ribo-depletion factors are method-local; drift-guarded via provenance.source_sha256.rna_virome.",
            }
    except Exception:
        pass

    return "source-derived (VLP, virion size, amplification + RT defaults)"


def export_model_constants() -> tuple[dict[str, Any], str]:
    constants: dict[str, Any] = {
        "provenance": provenance(),
        "note": (
            "Model constants for the client-side mean-model replay. Blocks with "
            '"derived": true are read from ViroForge source; blocks with '
            '"derived": false are drift-guarded literals whose source module is '
            "hashed in provenance.source_sha256."
        ),
    }
    constants.update(_literal_constants())
    status = _overlay_derived(constants)
    return constants, status


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help=f"ViroForge SQLite DB (default: {DEFAULT_DB})")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT, help=f"output directory (default: {DEFAULT_OUT})")
    parser.add_argument("--pretty", action="store_true", help="pretty-print collections.json (larger file)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    indent = 2 if args.pretty else None
    separators = None if args.pretty else (",", ":")

    collections = export_collections(args.db)
    constants, status = export_model_constants()

    collections_path = args.out_dir / "collections.json"
    constants_path = args.out_dir / "model_constants.json"
    collections_path.write_text(json.dumps(collections, indent=indent, separators=separators))
    constants_path.write_text(json.dumps(constants, indent=2))

    prov = collections["provenance"]
    src = collections["source"]
    print(f"ViroForge {prov['viroforge_version']} @ {prov['viroforge_git_sha']}"
          f"{' (dirty)' if prov['viroforge_git_dirty'] else ''}  schema {prov['schema_version']}")
    print(f"collections.json     {src['n_collections']} collections, {src['n_collection_genomes']} genomes"
          f"  ({collections_path.stat().st_size / 1024:,.0f} KB) -> {collections_path}")
    print(f"model_constants.json ({constants_path.stat().st_size / 1024:,.0f} KB) -> {constants_path}")
    print(f"  model constants: {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
