#!/usr/bin/env python3
"""
Migration script: Add collection-specific contamination defaults.

Adds default_host_dna_pct, default_rrna_pct, default_reagent_pct,
default_phix_pct, and host_organism columns to body_site_collections,
then populates them with sample-type-appropriate defaults.

These defaults are based on published literature:
- Roux et al. 2016 (ViromeQC) - contamination ranges by sample type
- Shkoporov et al. 2019 - gut virome contamination
- Thurber et al. 2009 - marine virome protocols
- Salter et al. 2014 - reagent contamination
"""

import sqlite3
import sys
from pathlib import Path

# Collection-specific contamination defaults (for "realistic" level)
# These represent typical contamination BEFORE VLP enrichment
COLLECTION_DEFAULTS = {
    # id: (host_dna_pct, rrna_pct, reagent_pct, phix_pct, host_organism)
    1:  (5.0,  3.0,  0.5, 0.1, 'human'),   # Gut Virome - Adult Healthy (Western Diet)
    2:  (10.0, 5.0,  0.5, 0.1, 'human'),   # Oral Virome - Saliva (Healthy)
    3:  (15.0, 2.0,  0.5, 0.1, 'human'),   # Skin Virome - Sebaceous Sites (Healthy)
    4:  (20.0, 5.0,  0.5, 0.1, 'human'),   # Respiratory Virome - Nasopharynx (Healthy)
    5:  (0.05, 5.0,  0.2, 0.1, 'none'),    # Marine Virome - Coastal Surface Water
    6:  (0.05, 8.0,  0.3, 0.1, 'none'),    # Soil Virome - Agricultural
    7:  (0.05, 6.0,  0.2, 0.1, 'none'),    # Freshwater Virome - Lake Surface Water
    8:  (3.0,  3.0,  0.5, 0.1, 'mouse'),   # Mouse Gut Virome - Laboratory (C57BL/6)
    9:  (1.0,  5.0,  0.3, 0.1, 'human'),   # Wastewater Virome - Urban Treatment Plant
    10: (8.0,  4.0,  0.5, 0.1, 'human'),   # IBD Gut Virome
    11: (10.0, 5.0,  0.5, 0.1, 'human'),   # HIV+ Gut Virome
    12: (25.0, 5.0,  0.5, 0.1, 'human'),   # Cystic Fibrosis Respiratory Virome
    13: (20.0, 8.0,  0.5, 0.1, 'human'),   # Human Respiratory RNA Virome
    14: (0.1,  3.0,  0.3, 0.1, 'none'),    # Arbovirus Environmental (Mosquito Virome)
    15: (5.0,  10.0, 0.5, 0.1, 'human'),   # Fecal RNA Virome
    16: (30.0, 3.0,  0.5, 0.1, 'human'),   # Vaginal Virome (Healthy)
    17: (40.0, 0.5,  0.3, 0.1, 'human'),   # Blood/Plasma Virome (Healthy)
    18: (20.0, 2.0,  0.3, 0.1, 'human'),   # Ocular Surface Virome (Healthy)
    19: (25.0, 5.0,  0.5, 0.1, 'human'),   # Lower Respiratory (Lung) Virome (Healthy)
    20: (15.0, 3.0,  0.3, 0.1, 'human'),   # Urinary Virome (Healthy)
}


def migrate(db_path: str):
    """Add contamination columns and populate defaults."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Check if columns already exist
    columns = [row[1] for row in cursor.execute("PRAGMA table_info(body_site_collections)")]

    new_columns = [
        ("default_host_dna_pct", "REAL"),
        ("default_rrna_pct", "REAL"),
        ("default_reagent_pct", "REAL"),
        ("default_phix_pct", "REAL DEFAULT 0.1"),
        ("host_organism", "TEXT DEFAULT 'human'"),
    ]

    for col_name, col_type in new_columns:
        if col_name not in columns:
            cursor.execute(f"ALTER TABLE body_site_collections ADD COLUMN {col_name} {col_type}")
            print(f"  Added column: {col_name}")

    # Populate defaults
    for cid, (host, rrna, reagent, phix, organism) in COLLECTION_DEFAULTS.items():
        cursor.execute("""
            UPDATE body_site_collections
            SET default_host_dna_pct = ?,
                default_rrna_pct = ?,
                default_reagent_pct = ?,
                default_phix_pct = ?,
                host_organism = ?
            WHERE collection_id = ?
        """, (host, rrna, reagent, phix, organism, cid))

        if cursor.rowcount > 0:
            print(f"  Collection {cid}: host={host}%, rrna={rrna}%, reagent={reagent}%, organism={organism}")

    conn.commit()

    # Verify
    rows = cursor.execute("""
        SELECT collection_id, collection_name, default_host_dna_pct, default_rrna_pct,
               default_reagent_pct, default_phix_pct, host_organism
        FROM body_site_collections
        ORDER BY collection_id
    """).fetchall()

    print(f"\n  Migrated {len(rows)} collections:")
    for r in rows:
        print(f"    {r[0]:2d}. {r[1][:45]:45s} host={r[2]}% rrna={r[3]}% reagent={r[4]}% phix={r[5]}% org={r[6]}")

    conn.close()
    print("\nMigration complete.")


if __name__ == "__main__":
    db_path = sys.argv[1] if len(sys.argv) > 1 else "viroforge/data/viral_genomes.db"
    print(f"Migrating: {db_path}")
    migrate(db_path)
