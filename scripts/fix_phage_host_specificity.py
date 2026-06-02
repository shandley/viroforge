#!/usr/bin/env python3
"""
Fix phage host specificity in body-site collections.

Problem: Several collections contain phages that infect bacteria NOT found
at the body site (e.g., Vibrio/Ralstonia/Xanthomonas phages in oral/gut,
Coliphages in skin/respiratory). This happened because the curation scripts
selected phages by broad family (Microviridae, Inoviridae) rather than by
host bacterium.

This script replaces non-site-appropriate phages with phages that infect
bacteria actually found at each body site.

Usage:
    python scripts/fix_phage_host_specificity.py [--dry-run] [--apply]
"""

import sqlite3
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DB_PATH = Path(__file__).parent.parent / 'viroforge' / 'data' / 'viral_genomes.db'

# ============================================================
# Define which phages are WRONG for each body site
# (bacteria not found at that body site)
# ============================================================

NON_SITE_PHAGES = [
    # Marine/aquatic bacteria phages — wrong for ALL human body sites
    '%Vibrio phage%',
    '%Bdellovibrio%',
    # Plant bacteria phages — wrong for ALL human body sites
    '%Ralstonia phage%',
    '%Xanthomonas phage%',
    '%Erwinia phage%',
    # Insect bacteria phages
    '%Spiroplasma phage%',
    # Animal-associated
    '%Guinea pig Chlamydia%',
]

# Per-collection wrong phages (in addition to the global list above)
COLLECTION_WRONG_PHAGES = {
    # Gut: Stenotrophomonas is not a typical gut bacterium
    # Note: Coliphages/Enterobacteria phages are CORRECT for gut (E. coli is gut)
    1: [
        '%Stenotrophomonas phage%',
    ],
    # Oral: E. coli phages don't belong in oral cavity
    2: [
        '%Coliphage%', '%Enterobacteria phage%', '%Escherichia phage%',
        'Genome of phage G4%',
        '%Pseudomonas phage pf%', '%Stenotrophomonas phage%',
    ],
    # Skin: E. coli phages don't belong on skin
    3: [
        '%Coliphage%', '%Enterobacteria phage%', '%Escherichia phage%',
        'Genome of phage G4%',
    ],
    # Respiratory: E. coli phages don't belong in respiratory tract
    4: [
        '%Coliphage%', '%Enterobacteria phage%', '%Escherichia phage%',
        'Genome of phage G4%', 'Bacteriophage phiK%',
    ],
    # IBD Gut: Stenotrophomonas not typical gut bacterium
    # Note: Coliphages/Enterobacteria phages are CORRECT for gut
    10: [
        '%Pseudomonas phage Pf%', '%Stenotrophomonas phage%',
    ],
    # HIV+ Gut: remaining marine/plant phages
    11: [
        # Global NON_SITE_PHAGES handles Vibrio/Ralstonia/Xanthomonas
        # Coliphages are correct for gut
    ],
    # Vaginal: E. coli phages don't belong in vaginal
    16: [
        '%Coliphage%', '%Enterobacteria phage%', '%Escherichia phage%',
        'Genome of phage G4%',
    ],
}

# ============================================================
# Define body-site-appropriate replacement phages
# ============================================================

SITE_APPROPRIATE_PHAGES = {
    # Gut: Bacteroides, Faecalibacterium, Clostridium are dominant gut bacteria
    1: {
        'primary': "g.genome_name LIKE '%Bacteroides phage%' OR g.genome_name LIKE '%Faecalibacterium%phage%'",
        'secondary': "g.genome_name LIKE '%Clostridium phage%' OR g.genome_name LIKE '%Enterococcus phage%'",
        'fallback': "g.genome_name LIKE '%Klebsiella phage%'",
    },
    # Oral: dominated by Streptococcus in the oral cavity
    2: {
        'primary': "g.genome_name LIKE '%Streptococcus phage%'",
        'secondary': "g.genome_name LIKE '%Actinomyces phage%' OR g.genome_name LIKE '%Fusobacterium phage%' OR g.genome_name LIKE '%Neisseria phage%' OR g.genome_name LIKE '%Veillonella phage%'",
        'fallback': "g.genome_name LIKE '%Streptococcus phage%'",
    },
    # Skin: dominated by Cutibacterium (Propionibacterium) and Staphylococcus
    3: {
        'primary': "g.genome_name LIKE '%Propionibacterium phage%' OR g.genome_name LIKE '%Cutibacterium phage%'",
        'secondary': "g.genome_name LIKE '%Staphylococcus phage%'",
        'fallback': "g.genome_name LIKE '%Propionibacterium phage%' OR g.genome_name LIKE '%Staphylococcus phage%'",
    },
    # Respiratory: Streptococcus, Haemophilus, Moraxella
    4: {
        'primary': "g.genome_name LIKE '%Streptococcus phage%'",
        'secondary': "g.genome_name LIKE '%Haemophilus phage%' OR g.genome_name LIKE '%Moraxella phage%' OR g.genome_name LIKE '%Neisseria phage%'",
        'fallback': "g.genome_name LIKE '%Streptococcus phage%'",
    },
    # IBD Gut: Bacteroides, Faecalibacterium, E. coli phages are appropriate
    10: {
        'primary': "g.genome_name LIKE '%Bacteroides phage%' OR g.genome_name LIKE '%Faecalibacterium%phage%'",
        'secondary': "g.genome_name LIKE '%Escherichia phage%' OR g.genome_name LIKE '%Klebsiella phage%'",
        'fallback': "g.genome_name LIKE '%Enterococcus phage%'",
    },
    # HIV+ Gut: same as gut
    11: {
        'primary': "g.genome_name LIKE '%Bacteroides phage%' OR g.genome_name LIKE '%Faecalibacterium%phage%'",
        'secondary': "g.genome_name LIKE '%Clostridium phage%' OR g.genome_name LIKE '%Enterococcus phage%'",
        'fallback': "g.genome_name LIKE '%Klebsiella phage%'",
    },
    # Vaginal: Lactobacillus phages
    16: {
        'primary': "g.genome_name LIKE '%Lactobacillus phage%'",
        'secondary': "g.genome_name LIKE '%Lactococcus phage%'",
        'fallback': "g.genome_name LIKE '%Lactobacillus phage%'",
    },
}


def find_wrong_phages(conn, collection_id: int) -> List[Dict]:
    """Find non-site-appropriate phages in a collection."""
    # Build WHERE clause for wrong phages
    conditions = []
    for pattern in NON_SITE_PHAGES:
        conditions.append(f"g.genome_name LIKE '{pattern}'")

    # Add collection-specific wrong phages
    for pattern in COLLECTION_WRONG_PHAGES.get(collection_id, []):
        conditions.append(f"g.genome_name LIKE '{pattern}'")

    where = " OR ".join(conditions)

    cursor = conn.execute(f"""
        SELECT cg.genome_id, g.genome_name, cg.relative_abundance, cg.abundance_rank
        FROM collection_genomes cg
        JOIN genomes g ON cg.genome_id = g.genome_id
        WHERE cg.collection_id = ? AND ({where})
        ORDER BY g.genome_name
    """, (collection_id,))

    return [dict(row) for row in cursor.fetchall()]


def find_replacements(conn, collection_id: int, count: int, exclude_ids: Set[str]) -> List[Dict]:
    """Find body-site-appropriate phage replacements."""
    site_config = SITE_APPROPRIATE_PHAGES.get(collection_id)
    if not site_config:
        return []

    placeholders = ','.join('?' for _ in exclude_ids) if exclude_ids else "''"
    exclude_params = list(exclude_ids) if exclude_ids else []

    # Try primary, then secondary, then fallback
    for source in ['primary', 'secondary', 'fallback']:
        condition = site_config[source]
        query = f"""
            SELECT g.genome_id, g.genome_name
            FROM genomes g
            WHERE ({condition})
            AND g.genome_id NOT IN ({placeholders})
            AND g.sequence IS NOT NULL
            ORDER BY RANDOM()
            LIMIT ?
        """
        params = exclude_params + [count]
        cursor = conn.execute(query, params)
        results = [dict(row) for row in cursor.fetchall()]

        if results:
            # Add found IDs to exclude set for next rounds
            for r in results:
                exclude_ids.add(r['genome_id'])

            if len(results) >= count:
                return results[:count]
            else:
                # Got some but not enough, try next source
                remaining = count - len(results)
                more = find_replacements(conn, collection_id, remaining, exclude_ids)
                return results + more

    return []


def fix_collection(conn, collection_id: int, collection_name: str, dry_run: bool) -> Dict:
    """Fix phage host specificity for one collection."""
    wrong = find_wrong_phages(conn, collection_id)

    if not wrong:
        logger.info(f"  Collection {collection_id} ({collection_name}): No wrong phages")
        return {'removed': 0, 'added': 0}

    logger.info(f"  Collection {collection_id} ({collection_name}): "
                f"{len(wrong)} non-site-appropriate phages")

    # Get current genome IDs
    cursor = conn.execute(
        "SELECT genome_id FROM collection_genomes WHERE collection_id = ?",
        (collection_id,)
    )
    current_ids = {row['genome_id'] for row in cursor.fetchall()}

    # Find replacements
    wrong_ids = {w['genome_id'] for w in wrong}
    exclude_ids = current_ids - wrong_ids  # Don't pick phages already in collection
    replacements = find_replacements(conn, collection_id, len(wrong), exclude_ids)

    # Log changes
    for w in wrong:
        logger.info(f"    REMOVE: {w['genome_name']}")
    for r in replacements:
        logger.info(f"    ADD:    {r['genome_name']}")

    if len(replacements) < len(wrong):
        logger.warning(f"    Only found {len(replacements)} replacements for {len(wrong)} removals")

    # Apply
    if not dry_run:
        for w in wrong:
            conn.execute(
                "DELETE FROM collection_genomes WHERE collection_id = ? AND genome_id = ?",
                (collection_id, w['genome_id'])
            )

        for i, r in enumerate(replacements):
            # Use average abundance from removed phages
            avg_abundance = sum(w['relative_abundance'] for w in wrong) / len(wrong)
            avg_rank = int(sum(w['abundance_rank'] for w in wrong) / len(wrong))

            conn.execute(
                "INSERT OR IGNORE INTO collection_genomes "
                "(collection_id, genome_id, relative_abundance, abundance_rank) "
                "VALUES (?, ?, ?, ?)",
                (collection_id, r['genome_id'], avg_abundance, avg_rank)
            )

        # Update genome count
        cursor = conn.execute(
            "SELECT COUNT(*) as cnt FROM collection_genomes WHERE collection_id = ?",
            (collection_id,)
        )
        new_count = cursor.fetchone()['cnt']
        conn.execute(
            "UPDATE body_site_collections SET n_genomes = ? WHERE collection_id = ?",
            (new_count, collection_id)
        )
        conn.commit()

        logger.info(f"  Applied: removed {len(wrong)}, added {len(replacements)}, "
                    f"new total: {new_count}")

    return {'removed': len(wrong), 'added': len(replacements)}


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Fix phage host specificity')
    parser.add_argument('--database', default=str(DB_PATH))
    parser.add_argument('--dry-run', action='store_true', default=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()

    if args.apply:
        args.dry_run = False

    conn = sqlite3.connect(args.database)
    conn.row_factory = sqlite3.Row

    collections_to_fix = {
        1: 'Gut Virome',
        2: 'Oral Virome',
        3: 'Skin Virome',
        4: 'Respiratory Virome',
        8: 'Mouse Gut Virome',
        10: 'IBD Gut Virome',
        11: 'HIV+ Gut Virome',
        16: 'Vaginal Virome',
    }

    mode = "DRY RUN" if args.dry_run else "APPLYING FIXES"
    logger.info(f"{'='*70}")
    logger.info(f"Phage Host Specificity Fix ({mode})")
    logger.info(f"{'='*70}")

    total_removed = 0
    total_added = 0

    for cid, cname in sorted(collections_to_fix.items()):
        result = fix_collection(conn, cid, cname, args.dry_run)
        total_removed += result['removed']
        total_added += result['added']
        print()

    logger.info(f"{'='*70}")
    logger.info(f"Total: {total_removed} wrong phages removed, {total_added} site-appropriate added")

    if args.dry_run:
        print("\n*** DRY RUN — no changes made. Use --apply to fix. ***")

    conn.close()


if __name__ == '__main__':
    main()
