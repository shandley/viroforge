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
# These represent typical contamination BEFORE VLP enrichment.
#
# IMPORTANT: These are ESTIMATES based on available literature and expert
# knowledge. Exact read-level proportions vary significantly depending on
# DNA extraction method, library prep, sequencing depth, and biological
# variability between samples.
#
# What is well-supported by data:
# - Gut bacterial fraction (60-80%): HMP (Qin et al. 2010), MetaHIT (Lloyd-Price et al. 2019)
# - Host DNA fractions: Roux et al. 2016 (ViromeQC), Thurber et al. 2009
# - Marine/environmental: TARA Oceans
# - Gut fungal presence and dominant species: Nash et al. 2017, Sokol et al. 2017
#
# What is approximate/estimated:
# - Exact fungal read fractions (most mycobiome papers report relative abundance
#   within the fungal community, not absolute fraction of total sequencing reads)
# - Bacterial/fungal fractions for non-gut body sites (respiratory, urinary,
#   ocular, skin) — fewer quantitative shotgun metagenome studies available
# - Disease-specific elevations (IBD fungal, HIV+ fungal) — qualitative trends
#   are documented but exact read-level proportions are not well-quantified
#
# These defaults should be refined as more quantitative metagenomic data
# becomes available. Users can override with --bacterial-fraction and
# --fungal-fraction flags.
COLLECTION_DEFAULTS = {
    # id: (host_dna_pct, rrna_pct, reagent_pct, phix_pct, host_organism, bacterial_pct, fungal_pct)
    1:  (5.0,  3.0,  0.5, 0.1, 'human', 70.0, 1.0),   # Gut Virome - Adult Healthy
    2:  (10.0, 5.0,  0.5, 0.1, 'human', 60.0, 0.5),   # Oral Virome - Saliva
    3:  (15.0, 2.0,  0.5, 0.1, 'human', 50.0, 5.0),   # Skin Virome - Sebaceous (high Malassezia)
    4:  (20.0, 5.0,  0.5, 0.1, 'human', 45.0, 0.5),   # Respiratory Virome - Nasopharynx
    5:  (0.05, 5.0,  0.2, 0.1, 'none',  70.0, 0.1),   # Marine Virome - Coastal
    6:  (0.05, 8.0,  0.3, 0.1, 'none',  65.0, 3.0),   # Soil Virome - Agricultural
    7:  (0.05, 6.0,  0.2, 0.1, 'none',  60.0, 1.0),   # Freshwater Virome - Lake
    8:  (3.0,  3.0,  0.5, 0.1, 'mouse', 70.0, 0.5),   # Mouse Gut Virome
    9:  (1.0,  5.0,  0.3, 0.1, 'human', 55.0, 0.5),   # Wastewater Virome
    10: (8.0,  4.0,  0.5, 0.1, 'human', 65.0, 2.0),   # IBD Gut Virome (elevated Candida)
    11: (10.0, 5.0,  0.5, 0.1, 'human', 60.0, 3.0),   # HIV+ Gut Virome (elevated fungi)
    12: (25.0, 5.0,  0.5, 0.1, 'human', 40.0, 2.0),   # CF Respiratory Virome
    13: (20.0, 8.0,  0.5, 0.1, 'human', 35.0, 1.0),   # Human Respiratory RNA Virome
    14: (0.1,  3.0,  0.3, 0.1, 'none',  50.0, 0.5),   # Arbovirus Environmental (Mosquito)
    15: (5.0,  10.0, 0.5, 0.1, 'human', 65.0, 1.0),   # Fecal RNA Virome
    16: (30.0, 3.0,  0.5, 0.1, 'human', 40.0, 8.0),   # Vaginal Virome (high Candida)
    17: (40.0, 0.5,  0.3, 0.1, 'human', 5.0,  0.1),   # Blood/Plasma Virome (minimal microbes)
    18: (20.0, 2.0,  0.3, 0.1, 'human', 30.0, 1.0),   # Ocular Surface Virome
    19: (25.0, 5.0,  0.5, 0.1, 'human', 35.0, 2.0),   # Lower Respiratory (Lung) Virome
    20: (15.0, 3.0,  0.3, 0.1, 'human', 40.0, 1.0),   # Urinary Virome
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
        ("default_bacterial_pct", "REAL DEFAULT 0.0"),
        ("default_fungal_pct", "REAL DEFAULT 0.0"),
    ]

    for col_name, col_type in new_columns:
        if col_name not in columns:
            cursor.execute(f"ALTER TABLE body_site_collections ADD COLUMN {col_name} {col_type}")
            print(f"  Added column: {col_name}")

    # Populate defaults
    for cid, (host, rrna, reagent, phix, organism, bacterial, fungal) in COLLECTION_DEFAULTS.items():
        cursor.execute("""
            UPDATE body_site_collections
            SET default_host_dna_pct = ?,
                default_rrna_pct = ?,
                default_reagent_pct = ?,
                default_phix_pct = ?,
                host_organism = ?,
                default_bacterial_pct = ?,
                default_fungal_pct = ?
            WHERE collection_id = ?
        """, (host, rrna, reagent, phix, organism, bacterial, fungal, cid))

        if cursor.rowcount > 0:
            print(f"  Collection {cid}: host={host}% rrna={rrna}% bact={bacterial}% fungal={fungal}% org={organism}")

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
