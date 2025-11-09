#!/usr/bin/env python3
"""
Fix citation errors in database collection metadata.

This script corrects the literature_references field for Collections 17-20
with accurate journal names and updated citations.

Author: ViroForge Development Team
Date: 2025-11-09
"""

import sqlite3
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def fix_citations():
    """Update literature_references in database with corrected citations."""

    conn = sqlite3.connect('viroforge/data/viral_genomes.db')
    cursor = conn.cursor()

    # Collection 17: Wastewater (fix Crank year and journal names)
    logger.info("Updating Collection 17 (Wastewater) citations...")
    cursor.execute("""
        UPDATE body_site_collections
        SET literature_references = ?
        WHERE collection_id = 17
    """, (
        'Crits-Christoph et al. 2021 (mBio), Crank et al. 2022 (Sci Total Environ), '
        'Symonds et al. 2019 (Curr Opin Virol)',
    ))

    # Collection 18: IBD Gut (fix journal names, remove Duerkop)
    logger.info("Updating Collection 18 (IBD) citations...")
    cursor.execute("""
        UPDATE body_site_collections
        SET literature_references = ?
        WHERE collection_id = 18
    """, (
        'Norman et al. 2015 (Cell), Zuo et al. 2019 (Gut), '
        'Clooney et al. 2019 (Cell Host Microbe)',
    ))

    # Collection 19: HIV+ Gut (fix journal names)
    logger.info("Updating Collection 19 (HIV+) citations...")
    cursor.execute("""
        UPDATE body_site_collections
        SET literature_references = ?
        WHERE collection_id = 19
    """, (
        'Handley et al. 2012 (Cell), Monaco et al. 2016 (Cell Host Microbe), '
        'Vujkovic-Cvijin et al. 2013 (Sci Transl Med), '
        'Nganou-Makamdop et al. 2018 (Cell Host Microbe)',
    ))

    # Collection 20: CF Respiratory (fix journal names)
    logger.info("Updating Collection 20 (CF) citations...")
    cursor.execute("""
        UPDATE body_site_collections
        SET literature_references = ?
        WHERE collection_id = 20
    """, (
        'Lim et al. 2014 (J Clin Microbiol), Wat et al. 2008 (J Cyst Fibros), '
        'Esther Jr et al. 2014 (Pediatr Pulmonol), '
        'Cuthbertson et al. 2020 (J Cyst Fibros)',
    ))

    conn.commit()

    # Verify updates
    logger.info("\n" + "=" * 80)
    logger.info("VERIFICATION: Updated citations in database")
    logger.info("=" * 80)

    cursor.execute("""
        SELECT collection_id, collection_name, literature_references
        FROM body_site_collections
        WHERE collection_id BETWEEN 17 AND 20
        ORDER BY collection_id
    """)

    for row in cursor.fetchall():
        coll_id, coll_name, lit_refs = row
        logger.info(f"\nCollection {coll_id}: {coll_name}")
        logger.info(f"  References: {lit_refs}")

    conn.close()
    logger.info("\n✓ Database citations corrected!")


if __name__ == '__main__':
    fix_citations()
