#!/usr/bin/env python3
"""Build a verified family-level property map from the ICTV VMR.

The ViroForge DB's `genomes.genome_type` column is unreliable (RNA-virus genomes
are frequently mislabeled `dsDNA`, e.g. the entire arbovirus collection). This
script derives, straight from the authoritative ICTV VMR (data/ictv/VMR_current.xlsx,
sheet "VMR MSL40"), a per-family map of:

  - nucleic-acid / Baltimore type  (from the VMR "Genome" column)
  - host type: phage | eukaryotic | unknown  (from the VMR "Host source" column)

aggregated to the family level by majority vote, with a purity fraction so callers
can see how clean each family is. The output TSV is the trusted source the
composition evaluator uses instead of `genome_type`.

Usage:
  python3 scripts/build_composition_reference.py \
      --vmr data/ictv/VMR_current.xlsx \
      --out data/reference_profiles/family_properties.tsv
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

# Host-source domain -> host type. VMR host_source uses domain words for real
# hosts and "(S)" / "unknown" tags for environmental samples with no known host.
PHAGE_HOSTS = {"bacteria", "archaea"}
EUK_HOSTS = {"vertebrates", "invertebrates", "plants", "fungi", "protists", "algae"}


def _host_type(raw: object) -> str:
    """Classify a VMR Host source value into phage / eukaryotic / unknown.

    Values can be comma-separated combinations ("invertebrates, vertebrates") and
    carry a trailing "(S)" environmental-source tag which we strip. A value that
    mixes phage and eukaryotic hosts is rare; we treat any eukaryotic host as
    eukaryotic (a virus that infects a eukaryote is not a phage)."""
    if not isinstance(raw, str) or not raw.strip():
        return "unknown"
    tokens = [t.strip().replace("(s)", "").strip() for t in raw.lower().split(",")]
    tokens = [t for t in tokens if t]
    has_euk = any(t in EUK_HOSTS for t in tokens)
    has_phage = any(t in PHAGE_HOSTS for t in tokens)
    if has_euk:
        return "eukaryotic"
    if has_phage:
        return "phage"
    return "unknown"


def _na_class(baltimore: str) -> str:
    """Coarse DNA vs RNA from a Baltimore/Genome string."""
    b = baltimore.upper()
    if "RT" in b:
        return "RT"  # reverse-transcribing (spans DNA/RNA)
    if "RNA" in b:
        return "RNA"
    if "DNA" in b:
        return "DNA"
    return "unknown"


def _na_coarse(baltimore: str) -> str:
    """Strand+class type in the DB's 4-value vocabulary (dsDNA/ssDNA/ssRNA/dsRNA).

    Reverse-transcribing genomes are folded to their packaged form (dsDNA-RT ->
    dsDNA, ssRNA-RT -> ssRNA). Returns 'unknown' if neither strand+class matches."""
    for tag in ("dsDNA", "ssDNA", "dsRNA", "ssRNA"):
        if tag in baltimore:
            return tag
    return "unknown"


def build(vmr_path: Path) -> pd.DataFrame:
    df = pd.read_excel(vmr_path, sheet_name="VMR MSL40")
    df = df[df["Family"].notna()].copy()

    records = []
    for family, grp in df.groupby("Family"):
        n = len(grp)

        balt = grp["Genome"].dropna().astype(str)
        balt_top, balt_purity = _mode_purity(balt, n)
        na_class = _na_class(balt_top) if balt_top else "unknown"

        # coarse strand+class type, with its own (higher) purity: families whose
        # exact Baltimore string varies but whose ds/ss + DNA/RNA does not.
        coarse_top, coarse_purity = _mode_purity(balt.map(_na_coarse), n)

        hosts = grp["Host source"].map(_host_type)
        host_top, host_purity = _mode_purity(hosts, n)

        records.append(
            {
                "family": family,
                "na_type": balt_top or "unknown",
                "na_class": na_class,
                "na_purity": round(balt_purity, 3),
                "na_coarse": coarse_top,
                "na_coarse_purity": round(coarse_purity, 3),
                "host_type": host_top,
                "host_purity": round(host_purity, 3),
                "n_vmr_records": n,
            }
        )
    out = pd.DataFrame(records).sort_values("family").reset_index(drop=True)
    return out


def _mode_purity(values, total: int) -> tuple[str, float]:
    vals = [v for v in values if v and v != "unknown"]
    if not vals:
        return ("unknown", 0.0)
    counts = Counter(vals)
    top, cnt = counts.most_common(1)[0]
    return (top, cnt / total)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vmr", type=Path, default=Path("data/ictv/VMR_current.xlsx"))
    ap.add_argument("--out", type=Path, default=Path("data/reference_profiles/family_properties.tsv"))
    args = ap.parse_args()

    out = build(args.vmr)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"{len(out)} families -> {args.out}")
    print("\nhost_type distribution:")
    print(out["host_type"].value_counts().to_string())
    print("\nna_class distribution:")
    print(out["na_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
