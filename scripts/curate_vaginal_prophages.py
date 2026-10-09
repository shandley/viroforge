#!/usr/bin/env python3
"""
Curate vaginal prophage sequences for differential prophage analysis.

Downloads the Gardnerella phage vB_Gva_AB1 (MW387018.1) — the only
standalone vaginal prophage genome in GenBank — and adds it to the
ViroForge database. Also reclassifies existing Lactobacillus phages
that are known induced prophages (Lv-1, phi jlb1, phiadh, Lj771,
Lj928, Lj965) from 'isolate' to 'prophage'.

Then wires prophages into the two vaginal collections with different
abundance profiles to enable differential prophage analysis:

  Collection 16 (Healthy vaginal):
    - Lactobacillus prophages present at LOW abundance (quiescent)
    - No Gardnerella prophage

  Collection 24 (HPV-infected vaginal):
    - Gardnerella prophage vB_Gva_AB1 at HIGH abundance (induced)
    - Lactobacillus prophages at LOW abundance (reduced Lactobacillus)

This creates a detectable differential signal for prophage analysis
pipelines: HPV-infected samples will show Gardnerella prophage reads
that are absent from healthy samples.

Known Lactobacillus phages that are induced prophages:
  - Lv-1 (GCF_000882635.1): from L. jensenii vaginal isolate
  - phi jlb1 (GCF_001507495.1): from L. gasseri ATCC 33323
  - phiadh (GCF_000848805.1): from L. gasseri
  - Lj771 (GCF_000871405.1): from L. johnsonii (already 'prophage')
  - Lj928 (GCF_000843425.1): from L. johnsonii (already 'prophage')
  - Lj965 (GCF_000842565.1): from L. johnsonii (already 'prophage')
  - KC5a (GCF_000868065.1): from L. gasseri
  - A2 (GCF_000848025.1): from L. casei/paracasei

Gardnerella phage vB_Gva_AB1 (MW387018.1):
  - 50,268 bp, circular dsDNA
  - Induced from G. vaginalis strain AB1
  - Published: Putonti et al. 2021, Microbiology Resource Announcements
  - BioProject: PRJNA687336

Author: ViroForge Development Team
Date: 2026-09-01
"""

import sqlite3
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

DB_PATH = 'viroforge/data/viral_genomes.db'

# Lactobacillus phages to reclassify as prophage
RECLASSIFY_TO_PROPHAGE = [
    'GCF_000882635.1',  # Lv-1 (L. jensenii vaginal)
    'GCF_001507495.1',  # phi jlb1 (L. gasseri)
    'GCF_000868065.1',  # KC5a (L. gasseri)
    'GCF_000848025.1',  # A2 (L. casei)
]


def download_gardnerella_prophage():
    """Download Gardnerella phage vB_Gva_AB1 from NCBI."""
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        logger.error("Biopython required: pip install biopython")
        raise

    Entrez.email = 'viroforge@example.com'
    accession = 'MW387018.1'

    logger.info(f"Downloading {accession} (Gardnerella phage vB_Gva_AB1)...")
    handle = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    handle.close()

    logger.info(f"  Length: {len(record.seq)} bp")
    logger.info(f"  Description: {record.description}")

    return {
        'genome_id': accession,
        'genome_name': 'Gardnerella phage vB_Gva_AB1, complete genome',
        'length': len(record.seq),
        'gc_content': (record.seq.count('G') + record.seq.count('C')) / len(record.seq),
        'sequence': str(record.seq),
        'n_segments': 1,
        'genome_provenance': 'prophage',
    }


def add_genome_to_db(conn, genome_data):
    """Insert a new genome into the database if not already present."""
    cursor = conn.cursor()

    cursor.execute("SELECT genome_id FROM genomes WHERE genome_id = ?",
                   (genome_data['genome_id'],))
    if cursor.fetchone():
        logger.info(f"  {genome_data['genome_id']} already in database — updating provenance")
        cursor.execute("UPDATE genomes SET genome_provenance = ? WHERE genome_id = ?",
                       (genome_data['genome_provenance'], genome_data['genome_id']))
        return False

    from datetime import date
    cursor.execute("""
        INSERT INTO genomes
        (genome_id, genome_name, sequence, length, gc_content, genome_type,
         n_segments, genome_provenance, source_database, date_added)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        genome_data['genome_id'],
        genome_data['genome_name'],
        genome_data['sequence'],
        genome_data['length'],
        genome_data['gc_content'],
        'dsDNA',
        genome_data['n_segments'],
        genome_data['genome_provenance'],
        'GenBank',
        date.today().isoformat(),
    ))

    cursor.execute("""
        INSERT OR REPLACE INTO taxonomy
        (genome_id, realm, kingdom, phylum, class, order_name, family, genus, species)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        genome_data['genome_id'],
        'Duplodnaviria', 'Heunggongvirae', 'Uroviricota',
        'Caudoviricetes', 'Unknown', 'Unknown',
        'Unknown', 'Gardnerella phage vB_Gva_AB1',
    ))

    logger.info(f"  Inserted {genome_data['genome_id']}: {genome_data['genome_name']}")
    return True


def reclassify_prophages(conn):
    """Reclassify known induced prophages from 'isolate' to 'prophage'."""
    logger.info("\nReclassifying known induced prophages...")
    cursor = conn.cursor()

    for gid in RECLASSIFY_TO_PROPHAGE:
        cursor.execute("SELECT genome_name, genome_provenance FROM genomes WHERE genome_id = ?",
                       (gid,))
        row = cursor.fetchone()
        if row:
            name, prov = row
            if prov != 'prophage':
                cursor.execute("UPDATE genomes SET genome_provenance = 'prophage' WHERE genome_id = ?",
                               (gid,))
                logger.info(f"  {gid}: {name[:50]} — {prov} -> prophage")
            else:
                logger.info(f"  {gid}: already prophage")
        else:
            logger.warning(f"  {gid}: not found in database")


def wire_prophages_to_collections(conn):
    """Add prophage genomes to Collections 16 and 24."""
    logger.info("\nWiring prophages into vaginal collections...")
    cursor = conn.cursor()

    # Prophages for healthy vaginal (Collection 16)
    # Lactobacillus prophages at low abundance (quiescent state)
    healthy_prophages = {
        'GCF_000882635.1': 0.002,   # Lv-1 (L. jensenii vaginal)
        'GCF_001507495.1': 0.001,   # phi jlb1 (L. gasseri)
        'GCF_000848805.1': 0.001,   # phiadh (L. gasseri)
    }

    # Prophages for HPV-infected vaginal (Collection 24)
    # Gardnerella prophage at HIGH abundance (induced in dysbiosis)
    # Lactobacillus prophages at LOW abundance (reduced Lactobacillus)
    hpv_prophages = {
        'MW387018.1': 0.04,          # Gardnerella vB_Gva_AB1 (INDUCED)
        'GCF_000882635.1': 0.0005,   # Lv-1 (reduced with Lactobacillus decline)
        'GCF_001507495.1': 0.0003,   # phi jlb1
    }

    for coll_id, prophages, label in [
        (16, healthy_prophages, 'Healthy vaginal'),
        (24, hpv_prophages, 'HPV-infected vaginal'),
    ]:
        logger.info(f"\n  {label} (Collection {coll_id}):")

        # Get current max abundance_rank
        cursor.execute("SELECT MAX(abundance_rank) FROM collection_genomes WHERE collection_id = ?",
                       (coll_id,))
        max_rank = cursor.fetchone()[0] or 0

        # Get current total abundance for renormalization
        cursor.execute("SELECT SUM(relative_abundance) FROM collection_genomes WHERE collection_id = ?",
                       (coll_id,))
        current_total = cursor.fetchone()[0] or 1.0

        rank = max_rank
        for gid, target_abundance in prophages.items():
            # Check if already in collection
            cursor.execute(
                "SELECT genome_id FROM collection_genomes WHERE collection_id = ? AND genome_id = ?",
                (coll_id, gid))
            if cursor.fetchone():
                logger.info(f"    {gid} already in collection {coll_id}")
                continue

            # Verify genome exists
            cursor.execute("SELECT genome_name FROM genomes WHERE genome_id = ?", (gid,))
            row = cursor.fetchone()
            if not row:
                logger.warning(f"    {gid} not in genomes table — skipping")
                continue

            rank += 1
            cursor.execute("""
                INSERT INTO collection_genomes
                (collection_id, genome_id, relative_abundance, prevalence, abundance_rank)
                VALUES (?, ?, ?, ?, ?)
            """, (coll_id, gid, target_abundance, 1.0, rank))
            logger.info(f"    Added {row[0][:50]} at abundance {target_abundance}")

        # Renormalize abundances to sum to 1.0
        cursor.execute("SELECT SUM(relative_abundance) FROM collection_genomes WHERE collection_id = ?",
                       (coll_id,))
        new_total = cursor.fetchone()[0]
        if new_total and abs(new_total - 1.0) > 0.001:
            cursor.execute("""
                UPDATE collection_genomes
                SET relative_abundance = relative_abundance / ?
                WHERE collection_id = ?
            """, (new_total, coll_id))
            logger.info(f"    Renormalized {coll_id}: {new_total:.4f} -> 1.0")

        # Update n_genomes
        cursor.execute("SELECT COUNT(*) FROM collection_genomes WHERE collection_id = ?",
                       (coll_id,))
        n = cursor.fetchone()[0]
        cursor.execute("UPDATE body_site_collections SET n_genomes = ? WHERE collection_id = ?",
                       (n, coll_id))
        logger.info(f"    Collection {coll_id} now has {n} genomes")


def set_collection24_defaults(conn):
    """Set contamination defaults for Collection 24."""
    logger.info("\nSetting contamination defaults for Collection 24...")
    cursor = conn.cursor()

    cursor.execute("""
        UPDATE body_site_collections
        SET default_host_pct = 20.0,
            default_rrna_pct = 3.0,
            default_reagent_pct = 0.5,
            default_phix_pct = 0.1,
            default_bacterial_pct = 65.0,
            default_fungal_pct = 5.0,
            default_archaeal_pct = 0.0,
            bacterial_community = 'vaginal',
            host_organism = 'human'
        WHERE collection_id = 24
    """)
    logger.info("  Set: bacterial=65%, fungal=5%, host=20%, community=vaginal")


def verify(conn):
    """Verify the prophage setup."""
    logger.info("\n" + "=" * 80)
    logger.info("VERIFICATION")
    logger.info("=" * 80)

    cursor = conn.cursor()

    for coll_id, label in [(16, 'Healthy'), (24, 'HPV-infected')]:
        cursor.execute("""
            SELECT g.genome_name, g.genome_provenance, cg.relative_abundance
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            WHERE cg.collection_id = ?
              AND g.genome_provenance = 'prophage'
            ORDER BY cg.relative_abundance DESC
        """, (coll_id,))
        rows = cursor.fetchall()
        logger.info(f"\n  {label} (Collection {coll_id}) — {len(rows)} prophages:")
        for name, prov, abund in rows:
            logger.info(f"    {name[:55]:55s} {abund:.4f}")

    # Check Gardnerella prophage is only in Collection 24
    cursor.execute("""
        SELECT cg.collection_id
        FROM collection_genomes cg
        WHERE cg.genome_id = 'MW387018.1'
    """)
    in_collections = [r[0] for r in cursor.fetchall()]
    if 24 in in_collections and 16 not in in_collections:
        logger.info("\n  Gardnerella prophage is in Collection 24 only")
    else:
        logger.warning(f"\n  Gardnerella prophage in collections: {in_collections}")


def main():
    logger.info("=" * 80)
    logger.info("CURATING VAGINAL PROPHAGE SEQUENCES")
    logger.info("=" * 80)

    conn = sqlite3.connect(DB_PATH)

    try:
        # Step 1: Download and add Gardnerella prophage
        gard = download_gardnerella_prophage()
        add_genome_to_db(conn, gard)

        # Step 2: Reclassify known Lactobacillus prophages
        reclassify_prophages(conn)

        # Step 3: Wire prophages into collections
        wire_prophages_to_collections(conn)

        # Step 4: Set contamination defaults for Collection 24
        set_collection24_defaults(conn)

        conn.commit()

        # Step 5: Verify
        verify(conn)

        logger.info("\n" + "=" * 80)
        logger.info("VAGINAL PROPHAGE CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nDifferential prophage signal:")
        logger.info("  Healthy (16): Lactobacillus prophages only (quiescent)")
        logger.info("  HPV (24):     Gardnerella prophage INDUCED + reduced Lactobacillus prophages")
        logger.info("\nExpected tutorial result:")
        logger.info("  Prophage analysis should find Gardnerella vB_Gva_AB1 reads")
        logger.info("  enriched in HPV samples and absent from healthy samples.")

    finally:
        conn.close()


if __name__ == '__main__':
    main()
