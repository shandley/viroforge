#!/usr/bin/env python3
"""Audit and fix per-collection relative-abundance normalization.

Each collection's genome `relative_abundance` values are meant to sum to 1.0.
Some collections (currently the RNA-workflow ones) sum to less than 1.0 because
post-hoc filtering scripts remove genomes without renormalizing the remainder.
This script audits every collection and, with --apply, rescales each collection's
abundances so they sum to 1.0. It is idempotent: running it twice is a no-op.

Default mode is a read-only audit. Pass --apply to write.

Usage:
  python3 scripts/renormalize_abundances.py --db viroforge/data/viral_genomes.db
  python3 scripts/renormalize_abundances.py --db viroforge/data/viral_genomes.db --apply
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

TOL = 1e-3


def audit(conn: sqlite3.Connection) -> list[tuple[int, str, float, int]]:
    rows = conn.execute(
        """
        SELECT bc.collection_id, bc.collection_name,
               SUM(cg.relative_abundance) AS s, COUNT(*) AS n
        FROM body_site_collections bc
        JOIN collection_genomes cg ON bc.collection_id = cg.collection_id
        GROUP BY bc.collection_id
        ORDER BY bc.collection_id
        """
    ).fetchall()
    return [(r[0], r[1], r[2] or 0.0, r[3]) for r in rows]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, required=True, help="ViroForge SQLite DB")
    ap.add_argument("--apply", action="store_true", help="write renormalized abundances (default: audit only)")
    args = ap.parse_args()
    if not args.db.exists():
        raise SystemExit(f"DB not found: {args.db}")

    conn = sqlite3.connect(str(args.db))
    try:
        rows = audit(conn)
        off = [r for r in rows if abs(r[2] - 1.0) > TOL]
        print(f"{len(rows)} collections; {len(off)} not summing to 1.0 (tol {TOL}):")
        for cid, name, s, n in off:
            print(f"  [{cid}] {name[:40]:40s} sum={s:.4f} n={n}")
        if not off:
            print("all collections normalized; nothing to do.")
            return
        if not args.apply:
            print("\naudit only. re-run with --apply to renormalize the collections above.")
            return
        for cid, name, s, n in off:
            if s <= 0:
                print(f"  [{cid}] sum={s}; cannot rescale, skipping")
                continue
            conn.execute(
                "UPDATE collection_genomes SET relative_abundance = relative_abundance / ? "
                "WHERE collection_id = ?",
                (s, cid),
            )
        conn.commit()
        # verify
        recheck = [r for r in audit(conn) if abs(r[2] - 1.0) > TOL]
        print(f"\napplied. collections still off after renormalization: {len(recheck)}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
