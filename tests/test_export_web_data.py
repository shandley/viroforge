"""CI test for scripts/export_web_data.py.

The tracked fixture DB (viroforge/data/test_viral_genomes.db) has no collections,
so this test builds its own minimal SQLite fixture with the tables and columns the
export join reads, runs the export against it, and asserts the output schema. It
exercises the export mechanism, not the full production DB (which is too large to
track).
"""
from __future__ import annotations

import json
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "export_web_data.py"

HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _build_fixture(db_path: Path) -> None:
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE body_site_collections(
            collection_id INTEGER PRIMARY KEY, collection_name TEXT, description TEXT,
            n_genomes INTEGER, selection_criteria TEXT, literature_references TEXT, version INTEGER);
        CREATE TABLE genomes(
            genome_id TEXT PRIMARY KEY, genome_name TEXT, sequence TEXT, length INTEGER,
            gc_content REAL, genome_type TEXT, genome_structure TEXT, n_segments INTEGER);
        CREATE TABLE taxonomy(
            genome_id TEXT PRIMARY KEY, realm TEXT, kingdom TEXT, phylum TEXT, class TEXT,
            order_name TEXT, family TEXT, subfamily TEXT, genus TEXT, species TEXT);
        CREATE TABLE collection_genomes(
            collection_id INTEGER, genome_id TEXT, relative_abundance REAL,
            prevalence REAL, abundance_rank INTEGER);
        """
    )
    # two collections, three genomes each, abundances normalized to 1.0
    for cid in (1, 2):
        conn.execute(
            "INSERT INTO body_site_collections VALUES (?,?,?,?,?,?,?)",
            (cid, f"Test collection {cid}", "desc", 3, "criteria", "refs", 1),
        )
    fams = ["Suoliviridae", "Microviridae", "Anelloviridae"]
    for cid in (1, 2):
        abunds = [0.6, 0.3, 0.1]
        for i, (fam, ab) in enumerate(zip(fams, abunds)):
            gid = f"G{cid}_{i}"
            conn.execute(
                "INSERT INTO genomes VALUES (?,?,?,?,?,?,?,?)",
                (gid, f"phage {gid}", "ACGT", 40000, 0.4, "dsDNA", "linear", 1),
            )
            conn.execute(
                "INSERT INTO taxonomy VALUES (?,?,?,?,?,?,?,?,?,?)",
                (gid, "Duplodnaviria", None, None, "Caudoviricetes", "Crassvirales", fam, None, f"genus{i}", f"species{i}"),
            )
            conn.execute(
                "INSERT INTO collection_genomes VALUES (?,?,?,?,?)",
                (cid, gid, ab, 1.0, i + 1),
            )
    conn.commit()
    conn.close()


def test_export_schema(tmp_path):
    fixture = tmp_path / "fixture.db"
    _build_fixture(fixture)
    out = tmp_path / "out"

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--db", str(fixture), "--out-dir", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    coll_file = out / "collections.json"
    const_file = out / "model_constants.json"
    assert coll_file.exists(), "collections.json not written"
    assert const_file.exists(), "model_constants.json not written"

    data = json.loads(coll_file.read_text())

    # provenance
    prov = data["provenance"]
    assert prov["schema_version"] == "1"
    for k in ("viroforge_version", "viroforge_git_sha", "viroforge_git_dirty", "generated_utc"):
        assert k in prov, f"missing provenance field: {k}"
    sha = prov["source_sha256"]
    for mod in ("vlp", "amplification", "rna_virome", "contamination"):
        assert mod in sha, f"missing source_sha256 for {mod}"
        assert HEX64.match(sha[mod]), f"source_sha256[{mod}] is not a 64-hex string"

    # source summary
    assert data["source"]["n_collections"] == 2
    assert data["source"]["n_collection_genomes"] == 6

    # collections + genome fields
    assert len(data["collections"]) == 2
    for col in data["collections"]:
        for k in ("id", "name", "n_genomes_loaded", "genomes"):
            assert k in col, f"collection missing {k}"
        assert len(col["genomes"]) == 3
        # abundances sum to 1.0 (guards the normalization invariant)
        s = sum(g["relative_abundance"] for g in col["genomes"])
        assert abs(s - 1.0) < 1e-6, f"collection {col['id']} abundances sum to {s}, not 1.0"
        for g in col["genomes"]:
            for k in ("genome_id", "name", "relative_abundance", "taxonomy"):
                assert k in g, f"genome missing {k}"
            assert "family" in g["taxonomy"]
        # sequences must never be exported
        assert "sequence" not in col["genomes"][0]
