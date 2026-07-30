#!/usr/bin/env python3
"""Populate per-collection contamination baselines in body_site_collections.

Reads the literature-informed defaults from
data/reference_profiles/contamination_defaults.tsv and writes them into the
default_host_pct / default_rrna_pct / default_reagent_pct / default_phix_pct /
host_organism columns, keyed by collection_id. Adds the columns first if an older
DB predates them. Idempotent; audit-by-default (pass --apply to write). Folded into
setup-db so rebuilds carry the defaults.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path

COLUMNS = {
    "default_host_pct": "REAL",
    "default_rrna_pct": "REAL",
    "default_reagent_pct": "REAL",
    "default_phix_pct": "REAL",
    # The sample's own microbiome, not scaled by --contamination-level.
    "default_bacterial_pct": "REAL",
    "default_fungal_pct": "REAL",
    "default_archaeal_pct": "REAL",
    "bacterial_community": "TEXT",
    "host_organism": "TEXT",
}


def load_defaults(tsv: Path) -> list[dict]:
    rows = []
    with open(tsv, newline="") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in reader:
            rows.append({
                "collection_id": int(r["collection_id"]),
                "default_host_pct": float(r["host_pct"]),
                "default_rrna_pct": float(r["rrna_pct"]),
                "default_reagent_pct": float(r["reagent_pct"]),
                "default_phix_pct": float(r["phix_pct"]),
                "default_bacterial_pct": float(r["bacterial_pct"]),
                "default_fungal_pct": float(r["fungal_pct"]),
                "default_archaeal_pct": float(r["archaeal_pct"]),
                "bacterial_community": r["bacterial_community"].strip(),
                "host_organism": r["host_organism"].strip(),
            })
    return rows


def ensure_columns(conn: sqlite3.Connection) -> list[str]:
    existing = {r[1] for r in conn.execute("PRAGMA table_info(body_site_collections)")}
    added = []
    for col, typ in COLUMNS.items():
        if col not in existing:
            conn.execute(f"ALTER TABLE body_site_collections ADD COLUMN {col} {typ}")
            added.append(col)
    return added


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    ap.add_argument("--defaults", type=Path,
                    default=Path("data/reference_profiles/contamination_defaults.tsv"))
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    defaults = load_defaults(args.defaults)
    conn = sqlite3.connect(str(args.db))
    try:
        valid_ids = {r[0] for r in conn.execute("SELECT collection_id FROM body_site_collections")}
        rows = [d for d in defaults if d["collection_id"] in valid_ids]
        missing = sorted({d["collection_id"] for d in defaults} - valid_ids)

        print(f"{len(rows)}/{len(defaults)} defaults map to existing collections"
              + (f"; not in DB: {missing}" if missing else ""))
        if not args.apply:
            added = [c for c in COLUMNS
                     if c not in {r[1] for r in conn.execute("PRAGMA table_info(body_site_collections)")}]
            print(f"columns to add: {added or 'none (already present)'}")
            print("\naudit only. re-run with --apply to write.")
            return

        added = ensure_columns(conn)
        if added:
            print(f"added columns: {added}")
        # Built from COLUMNS rather than spelled out, so adding a column in one
        # place cannot leave it silently NULL here.
        assignments = ", ".join(f"{c}=?" for c in COLUMNS)
        sql = f"UPDATE body_site_collections SET {assignments} WHERE collection_id=?"
        for d in rows:
            missing = [c for c in COLUMNS if c not in d]
            if missing:
                raise KeyError(
                    f"collection {d['collection_id']}: TSV supplied no value for "
                    f"{missing}. Add it to load_defaults()."
                )
            conn.execute(sql, tuple(d[c] for c in COLUMNS) + (d["collection_id"],))
        conn.commit()
        print(f"applied contamination defaults to {len(rows)} collections.")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
