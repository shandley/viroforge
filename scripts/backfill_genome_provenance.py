#!/usr/bin/env python3
"""Backfill genome_provenance column on existing ViroForge databases.

Run once after upgrading to add the genome_provenance column and populate it
by parsing genome_name. Safe to run multiple times (idempotent).

Usage:
    python scripts/backfill_genome_provenance.py [--db PATH]
"""

import argparse
import logging
import sqlite3
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

VALID_PROVENANCE = {"prophage", "provirus", "endogenous", "satellite", "viroid", "isolate"}


def detect_genome_provenance(genome_name: str) -> str:
    name_lower = genome_name.lower()
    if "prophage" in name_lower:
        return "prophage"
    if "endogenous" in name_lower:
        return "endogenous"
    if "proviral" in name_lower or "provirus" in name_lower:
        return "provirus"
    if "satellite" in name_lower:
        return "satellite"
    if "viroid" in name_lower:
        return "viroid"
    return "isolate"


def backfill(db_path: str) -> None:
    conn = sqlite3.connect(db_path)

    columns = [row[1] for row in conn.execute("PRAGMA table_info(genomes)").fetchall()]
    if "genome_provenance" not in columns:
        logger.info("Adding genome_provenance column...")
        conn.execute("ALTER TABLE genomes ADD COLUMN genome_provenance TEXT NOT NULL DEFAULT 'isolate'")
        conn.commit()
    else:
        logger.info("genome_provenance column already exists.")

    rows = conn.execute("SELECT genome_id, genome_name FROM genomes").fetchall()
    updates = []
    for gid, gname in rows:
        prov = detect_genome_provenance(gname)
        if prov != "isolate":
            updates.append((prov, gid))

    if updates:
        conn.executemany("UPDATE genomes SET genome_provenance = ? WHERE genome_id = ?", updates)
        conn.commit()

    index_exists = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='index' AND name='idx_provenance'"
    ).fetchone()
    if not index_exists:
        logger.info("Creating idx_provenance index...")
        conn.execute("CREATE INDEX idx_provenance ON genomes(genome_provenance)")
        conn.commit()

    counts = conn.execute(
        "SELECT genome_provenance, COUNT(*) FROM genomes GROUP BY genome_provenance ORDER BY COUNT(*) DESC"
    ).fetchall()
    total = sum(c for _, c in counts)
    logger.info(f"Provenance distribution ({total:,} genomes):")
    for prov, n in counts:
        logger.info(f"  {prov:12s}  {n:6,}  ({n / total * 100:.1f}%)")

    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Backfill genome_provenance column")
    parser.add_argument("--db", default="viroforge/data/viral_genomes.db",
                        help="Path to ViroForge database")
    args = parser.parse_args()

    if not Path(args.db).exists():
        logger.error(f"Database not found: {args.db}")
        sys.exit(1)

    backfill(args.db)
    logger.info("Done.")


if __name__ == "__main__":
    main()
