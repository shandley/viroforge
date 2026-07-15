#!/usr/bin/env python3
"""
One-time migration: renumber collection IDs from the legacy 9-28 scheme to 1-20.

Background
----------
Older databases assigned collection IDs by AUTOINCREMENT, which left the eight
body-site collections at IDs 9-16 and the additional collections at 17-28, with
IDs 1-8 either empty or holding orphaned collection_genomes rows (near-duplicate
associations with no body_site_collections metadata). Issue #5 tracks the gap.

The curation scripts now assign explicit IDs 1-20 (body sites 1-8, additional
collections 9-20). This script brings an existing legacy database in line with
that scheme without a full rebuild. The mapping is a uniform shift: new = old - 8
(e.g. Gut 9->1, Wastewater 17->9, Urinary 28->20). Orphaned collection_genomes
rows at IDs 1-8 are removed first so they do not merge into the shifted rows.

The script is idempotent: run against an already-renumbered database it is a
no-op. Make a backup before running (the database is not tracked in git).

Usage
-----
    python scripts/migrate_renumber_collections.py [--database PATH] [--dry-run]
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

DEFAULT_DB = Path("viroforge/data/viral_genomes.db")
SHIFT = 8


def collection_ids(conn: sqlite3.Connection) -> list[int]:
    return [r[0] for r in conn.execute(
        "SELECT collection_id FROM body_site_collections ORDER BY collection_id")]


def summarize(conn: sqlite3.Connection) -> list[tuple[int, str, int]]:
    return list(conn.execute(
        """SELECT c.collection_id, c.collection_name,
                  (SELECT COUNT(*) FROM collection_genomes g
                   WHERE g.collection_id = c.collection_id)
           FROM body_site_collections c ORDER BY c.collection_id"""))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--database", type=Path, default=DEFAULT_DB)
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would change without writing")
    args = ap.parse_args()

    if not args.database.exists():
        print(f"ERROR: database not found: {args.database}", file=sys.stderr)
        return 1

    conn = sqlite3.connect(args.database)
    ids = collection_ids(conn)
    if not ids:
        print("ERROR: body_site_collections is empty; nothing to migrate.", file=sys.stderr)
        return 1

    if max(ids) <= 20 and min(ids) >= 1:
        print(f"Database already uses the 1-20 scheme (IDs {min(ids)}-{max(ids)}). No-op.")
        return 0
    if min(ids) < 9:
        print(f"ERROR: unexpected collection IDs {ids}; refusing to migrate.", file=sys.stderr)
        return 1

    orphans = conn.execute(
        "SELECT COUNT(*) FROM collection_genomes WHERE collection_id BETWEEN 1 AND 8").fetchone()[0]

    print("Before:")
    for cid, name, n in summarize(conn):
        print(f"  {cid:>2}  {name}  ({n} genomes)")
    print(f"Orphan collection_genomes rows at IDs 1-8 to delete: {orphans}")
    print(f"Shift: new = old - {SHIFT}")

    if args.dry_run:
        print("\n[dry-run] no changes written.")
        return 0

    try:
        conn.execute("PRAGMA foreign_keys = OFF")
        conn.execute("BEGIN")
        conn.execute("DELETE FROM collection_genomes WHERE collection_id BETWEEN 1 AND 8")
        conn.execute("UPDATE body_site_collections SET collection_id = collection_id - ? "
                     "WHERE collection_id >= 9", (SHIFT,))
        conn.execute("UPDATE collection_genomes SET collection_id = collection_id - ? "
                     "WHERE collection_id >= 9", (SHIFT,))
        # Keep AUTOINCREMENT bookkeeping consistent with the new max ID.
        conn.execute("UPDATE sqlite_sequence SET seq = "
                     "(SELECT MAX(collection_id) FROM body_site_collections) "
                     "WHERE name = 'body_site_collections'")
        conn.execute("COMMIT")
    except Exception as exc:  # noqa: BLE001
        conn.execute("ROLLBACK")
        print(f"ERROR: migration failed, rolled back: {exc}", file=sys.stderr)
        return 1
    finally:
        conn.execute("PRAGMA foreign_keys = ON")

    violations = conn.execute("PRAGMA foreign_key_check").fetchall()
    if violations:
        print(f"ERROR: foreign key check found {len(violations)} violations after migration.",
              file=sys.stderr)
        return 1

    print("\nAfter:")
    for cid, name, n in summarize(conn):
        print(f"  {cid:>2}  {name}  ({n} genomes)")
    print("\nMigration complete; foreign key check clean.")
    conn.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
