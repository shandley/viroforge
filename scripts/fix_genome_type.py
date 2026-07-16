#!/usr/bin/env python3
"""Correct genomes.genome_type from the verified ICTV family property map.

The stored genome_type mislabels RNA viruses as dsDNA (a silent default; e.g. the
entire arbovirus collection). This relabels each genome's coarse nucleic-acid type
(dsDNA / ssDNA / ssRNA / dsRNA, the existing 4-value vocabulary) from its family's
ICTV Baltimore class, but ONLY where the family is confidently typed
(na_purity >= threshold) and the stored value actually disagrees. Genomes with an
Unknown or low-purity family are left untouched.

Audit-by-default; pass --apply to write. genome_type is display/stats metadata
(cli browse + db_utils), not used for molecule-type logic, so this is low risk.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path


def load_map(path: Path, min_purity: float) -> dict[str, str]:
    """family -> coarse nucleic-acid type, using the coarse-level purity so that
    families with a certain ds/ss + DNA/RNA type but a varying exact Baltimore
    string (e.g. Phenuiviridae) still relabel."""
    m = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c = row.get("na_coarse", "")
            if c and c != "unknown" and float(row["na_coarse_purity"]) >= min_purity:
                m[row["family"]] = c
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    ap.add_argument("--properties", type=Path, default=Path("data/reference_profiles/family_properties.tsv"))
    ap.add_argument("--min-purity", type=float, default=0.9)
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    fam2type = load_map(args.properties, args.min_purity)
    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    rows = con.execute(
        """SELECT g.genome_id, g.genome_type AS old, t.family AS family
           FROM genomes g JOIN taxonomy t ON g.genome_id = t.genome_id
           WHERE t.family IS NOT NULL AND t.family != 'Unknown'"""
    ).fetchall()

    changes = []
    for r in rows:
        new = fam2type.get(r["family"])
        if new and new != r["old"]:
            changes.append((r["genome_id"], r["old"], new, r["family"]))

    # summarize by transition
    from collections import Counter
    trans = Counter((o, n) for _, o, n, _ in changes)
    print(f"{len(rows)} classified-family genomes; {len(changes)} need relabeling "
          f"(min_purity={args.min_purity}):")
    for (o, n), c in sorted(trans.items(), key=lambda kv: -kv[1]):
        print(f"  {o:8} -> {n:8} : {c}")

    # collection-14 arbovirus sanity sample
    print("\nsample (arbovirus collection 14):")
    cur = con.execute(
        """SELECT g.genome_id, g.genome_type, t.family FROM collection_genomes cg
           JOIN genomes g ON cg.genome_id=g.genome_id JOIN taxonomy t ON g.genome_id=t.genome_id
           WHERE cg.collection_id=14 LIMIT 4""")
    for gr in cur.fetchall():
        new = fam2type.get(gr["family"], gr["genome_type"])
        print(f"  {gr['genome_id']:14} {gr['family']:18} {gr['genome_type']} -> {new}")

    if not args.apply:
        print("\naudit only. re-run with --apply to write.")
        con.close()
        return

    con.executemany("UPDATE genomes SET genome_type = ? WHERE genome_id = ?",
                    [(n, gid) for gid, _, n, _ in changes])
    con.commit()
    print(f"\napplied {len(changes)} updates.")
    print("new genome_type distribution:")
    for r in con.execute("SELECT genome_type, COUNT(*) FROM genomes GROUP BY genome_type ORDER BY 2 DESC"):
        print(f"  {r[0]:8} {r[1]}")
    con.close()


if __name__ == "__main__":
    main()
