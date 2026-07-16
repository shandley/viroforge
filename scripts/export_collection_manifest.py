#!/usr/bin/env python3
"""Export a provenance manifest of collection membership from the database.

The 500 MB SQLite database is not tracked in git, so this manifest is the
committed, diffable record of exactly which genomes are in each of the 20
collections. Regenerate it after any change to the collections (for example a
fresh `viroforge setup-db`) and commit the diff, so the ground truth stays
auditable and reproducible.

Usage:
    python scripts/export_collection_manifest.py            # write manifest
    python scripts/export_collection_manifest.py --check    # verify, non-zero on drift
"""

from __future__ import annotations

import argparse
import hashlib
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = PROJECT_ROOT / "viroforge" / "data" / "viral_genomes.db"
TSV_PATH = PROJECT_ROOT / "data" / "collection_membership.tsv"
SUMMARY_PATH = PROJECT_ROOT / "data" / "collection_membership.sha256.txt"

TSV_HEADER = "collection_id\tcollection_name\tabundance_rank\tgenome_id\tgenome_name\trelative_abundance"


def _fetch_rows(conn: sqlite3.Connection) -> list[tuple]:
    """Return membership rows sorted deterministically for stable diffs."""
    cur = conn.execute(
        """
        SELECT cg.collection_id, c.collection_name, cg.abundance_rank,
               cg.genome_id, g.genome_name, cg.relative_abundance
        FROM collection_genomes cg
        JOIN body_site_collections c ON c.collection_id = cg.collection_id
        JOIN genomes g ON g.genome_id = cg.genome_id
        ORDER BY cg.collection_id, cg.genome_id
        """
    )
    return cur.fetchall()


def _render(rows: list[tuple]) -> tuple[str, str]:
    """Build the TSV body and the checksum summary from membership rows."""
    tsv_lines = [TSV_HEADER]
    for cid, cname, rank, gid, gname, ab in rows:
        tsv_lines.append(f"{cid}\t{cname}\t{rank}\t{gid}\t{gname}\t{ab}")
    tsv_body = "\n".join(tsv_lines) + "\n"

    # Per-collection checksum over sorted genome_ids only, so it is stable
    # against abundance-value formatting and reflects membership alone.
    per_collection: dict[int, list[str]] = {}
    names: dict[int, str] = {}
    for cid, cname, _rank, gid, _gname, _ab in rows:
        per_collection.setdefault(cid, []).append(gid)
        names[cid] = cname

    summary_lines = [
        "# ViroForge collection membership manifest (provenance record)",
        "# sha256 per collection is over its sorted genome_id list (membership only).",
        "# Regenerate with: python scripts/export_collection_manifest.py",
        "#",
        f"# collections: {len(per_collection)}   total_associations: {len(rows)}",
        "",
        "collection_id\tn_genomes\tsha256\tcollection_name",
    ]
    global_hash = hashlib.sha256()
    for cid in sorted(per_collection):
        gids = sorted(per_collection[cid])
        h = hashlib.sha256("\n".join(gids).encode()).hexdigest()
        summary_lines.append(f"{cid}\t{len(gids)}\t{h}\t{names[cid]}")
        global_hash.update(f"{cid}:{h}\n".encode())
    summary_lines.append("")
    summary_lines.append(f"global\t{len(rows)}\t{global_hash.hexdigest()}\tALL")
    summary = "\n".join(summary_lines) + "\n"
    return tsv_body, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify the committed manifest matches the database; exit 1 on drift.",
    )
    parser.add_argument("--database", default=str(DB_PATH), help="Path to the database.")
    args = parser.parse_args()

    db = Path(args.database)
    if not db.exists():
        print(f"ERROR: database not found: {db}", file=sys.stderr)
        return 2

    with sqlite3.connect(db) as conn:
        rows = _fetch_rows(conn)
    tsv_body, summary = _render(rows)

    if args.check:
        current_tsv = TSV_PATH.read_text() if TSV_PATH.exists() else ""
        if current_tsv == tsv_body:
            print(f"OK: manifest matches database ({len(rows)} associations).")
            return 0
        print("DRIFT: committed manifest does not match the database.", file=sys.stderr)
        print("Run: python scripts/export_collection_manifest.py", file=sys.stderr)
        return 1

    TSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    TSV_PATH.write_text(tsv_body)
    SUMMARY_PATH.write_text(summary)
    print(f"Wrote {TSV_PATH} ({len(rows)} associations)")
    print(f"Wrote {SUMMARY_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
