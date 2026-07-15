#!/usr/bin/env python3
"""
Populate the host_associations table by parsing host information from genome names.

Most phage genome names in RefSeq follow the pattern:
    "{Bacterium} phage {name}" → host_genus = "Bacterium"
    e.g., "Streptococcus phage Dp-1" → host = "Streptococcus"

This also handles:
    - "Enterobacteria phage X" → host_genus = "Escherichia" (common synonym)
    - "Coliphage X" → host_genus = "Escherichia"
    - "Cyanophage X" → host_genus = "Synechococcus" (or Prochlorococcus)
    - "MAG: Bacterium phage X" → strips MAG prefix

Usage:
    python scripts/populate_host_associations.py [--database PATH] [--dry-run]
"""

import sqlite3
import re
import logging
from pathlib import Path
from typing import Optional, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DB_PATH = Path(__file__).parent.parent / 'viroforge' / 'data' / 'viral_genomes.db'

# Special host name mappings (common names → standard genus)
HOST_MAPPINGS = {
    'Enterobacteria': 'Escherichia',
    'Enterobacterial': 'Escherichia',
    'Coliphage': 'Escherichia',
    'Cyanophage': 'Cyanobacteria',  # broad group
    'Mycobacteriophage': 'Mycobacterium',
    'Actinophage': 'Streptomyces',
}

# Body site associations for common host genera
HOST_BODY_SITES = {
    # Oral
    'Streptococcus': 'oral,respiratory,skin',
    'Actinomyces': 'oral',
    'Fusobacterium': 'oral',
    'Porphyromonas': 'oral',
    'Prevotella': 'oral,vaginal',
    'Treponema': 'oral',
    'Veillonella': 'oral',
    'Neisseria': 'oral,respiratory',
    'Aggregatibacter': 'oral',

    # Skin
    'Propionibacterium': 'skin',
    'Cutibacterium': 'skin',
    'Staphylococcus': 'skin,nasal,respiratory',
    'Corynebacterium': 'skin',
    'Micrococcus': 'skin',

    # Gut
    'Bacteroides': 'gut',
    'Faecalibacterium': 'gut',
    'Escherichia': 'gut',
    'Klebsiella': 'gut',
    'Enterococcus': 'gut',
    'Clostridium': 'gut',
    'Clostridioides': 'gut',
    'Lactobacillus': 'gut,vaginal',
    'Bifidobacterium': 'gut',
    'Salmonella': 'gut',
    'Shigella': 'gut',
    'Campylobacter': 'gut',
    'Helicobacter': 'gut',
    'Ruminococcus': 'gut',
    'Roseburia': 'gut',
    'Eubacterium': 'gut',
    'Parabacteroides': 'gut',

    # Respiratory
    'Pseudomonas': 'respiratory,environmental',
    'Haemophilus': 'respiratory',
    'Moraxella': 'respiratory',
    'Bordetella': 'respiratory',
    'Burkholderia': 'respiratory',
    'Acinetobacter': 'respiratory,skin',
    'Mycobacterium': 'respiratory,environmental',
    'Legionella': 'respiratory,environmental',

    # Vaginal
    'Gardnerella': 'vaginal',
    'Atopobium': 'vaginal',
    'Sneathia': 'vaginal',
    'Mobiluncus': 'vaginal',

    # Urinary
    'Proteus': 'urinary',
    'Ureaplasma': 'urinary,vaginal',
    'Mycoplasma': 'urinary,vaginal,respiratory',

    # Environmental
    'Vibrio': 'marine',
    'Synechococcus': 'marine,freshwater',
    'Prochlorococcus': 'marine',
    'Pelagibacter': 'marine',
    'Cellulophaga': 'marine',
    'Roseobacter': 'marine',
    'Ralstonia': 'plant_soil',
    'Xanthomonas': 'plant',
    'Erwinia': 'plant',
    'Agrobacterium': 'plant_soil',
    'Rhizobium': 'plant_soil',
    'Bradyrhizobium': 'soil',
    'Bacillus': 'soil,environmental',
    'Streptomyces': 'soil',
    'Arthrobacter': 'soil',
    'Rhodococcus': 'soil,environmental',
    'Caulobacter': 'freshwater',
    'Bdellovibrio': 'aquatic',
    'Spiroplasma': 'insect',
}


def parse_host_from_name(genome_name: str) -> Optional[Tuple[str, str]]:
    """
    Parse host genus from a phage genome name.

    Returns:
        Tuple of (host_genus, evidence) or None if not parseable.
    """
    name = genome_name.strip()

    # Strip "MAG: " prefix
    if name.startswith('MAG: '):
        name = name[5:]

    # Pattern 1: "Coliphage X" → Escherichia
    if name.startswith('Coliphage ') or name.startswith('coliphage '):
        return ('Escherichia', 'name_pattern:Coliphage')

    # Pattern 2: "Genome of phage G4 (coliphage)"
    if 'coliphage' in name.lower():
        return ('Escherichia', 'name_pattern:coliphage_mention')

    # Pattern 3: "Cyanophage X"
    if name.startswith('Cyanophage ') or name.startswith('cyanophage '):
        return ('Cyanobacteria', 'name_pattern:Cyanophage')

    # Pattern 4: "{Genus} phage {name}" — the main pattern
    # Match: word(s) before "phage"
    match = re.match(r'^(\w+)\s+phage\b', name, re.IGNORECASE)
    if match:
        host_word = match.group(1)

        # Check special mappings
        if host_word in HOST_MAPPINGS:
            return (HOST_MAPPINGS[host_word], f'name_pattern:mapped_{host_word}')

        # If it looks like a genus name (capitalized, not a common word)
        if host_word[0].isupper() and len(host_word) > 3:
            return (host_word, 'name_pattern:genus_phage')

    # Pattern 5: "Bacteriophage X" — too generic, skip
    if name.startswith('Bacteriophage '):
        return None

    # Pattern 6: "{Genus} virus {name}" for some phages
    match = re.match(r'^(\w+)\s+virus\b', name, re.IGNORECASE)
    if match:
        # Only for known bacterial hosts (not eukaryotic viruses)
        host_word = match.group(1)
        if host_word in HOST_BODY_SITES:
            return (host_word, 'name_pattern:genus_virus')

    return None


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Populate host_associations table')
    parser.add_argument('--database', default=str(DB_PATH))
    parser.add_argument('--dry-run', action='store_true', default=False)
    parser.add_argument('--force', action='store_true', default=False,
                        help='Clear and repopulate the table (default: only add missing genomes)')
    args = parser.parse_args()

    conn = sqlite3.connect(args.database)
    conn.row_factory = sqlite3.Row

    # Check current state
    cursor = conn.execute("SELECT COUNT(*) as cnt FROM host_associations")
    existing = cursor.fetchone()['cnt']
    logger.info(f"Current host_associations rows: {existing}")

    if existing > 0:
        if args.force and not args.dry_run:
            logger.info("--force: clearing existing host_associations rows")
            conn.execute("DELETE FROM host_associations")
            conn.commit()
        else:
            logger.info("Table already populated; adding only missing genomes "
                        "(use --force to repopulate from scratch).")

    # Get all genomes
    cursor = conn.execute("SELECT genome_id, genome_name FROM genomes")
    genomes = [dict(row) for row in cursor.fetchall()]
    logger.info(f"Total genomes: {len(genomes)}")

    # Get already-populated genome IDs
    cursor = conn.execute("SELECT DISTINCT genome_id FROM host_associations")
    already_done = {row['genome_id'] for row in cursor.fetchall()}

    parsed = 0
    skipped = 0
    failed = 0
    body_site_counts = {}

    for genome in genomes:
        if genome['genome_id'] in already_done:
            skipped += 1
            continue

        result = parse_host_from_name(genome['genome_name'])
        if result is None:
            failed += 1
            continue

        host_genus, evidence = result
        body_sites = HOST_BODY_SITES.get(host_genus, 'unknown')

        if not args.dry_run:
            conn.execute("""
                INSERT INTO host_associations
                (genome_id, host_genus, host_name, body_site, association_type, evidence)
                VALUES (?, ?, ?, ?, 'predicted', ?)
            """, (
                genome['genome_id'],
                host_genus,
                host_genus,
                body_sites,
                evidence,
            ))

        parsed += 1

        # Track body site distribution
        for site in body_sites.split(','):
            body_site_counts[site] = body_site_counts.get(site, 0) + 1

    if not args.dry_run:
        conn.commit()

    logger.info(f"{'='*60}")
    logger.info(f"Results:")
    logger.info(f"  Parsed and inserted: {parsed}")
    logger.info(f"  Already in table:    {skipped}")
    logger.info(f"  Could not parse:     {failed}")
    logger.info(f"  Total:               {parsed + skipped + failed}")
    logger.info(f"")
    logger.info(f"Body site distribution of parsed hosts:")
    for site, count in sorted(body_site_counts.items(), key=lambda x: -x[1]):
        logger.info(f"  {site:25s} {count:>5d}")

    conn.close()

    if args.dry_run:
        print("\n*** DRY RUN — no changes made. Remove --dry-run to populate. ***")


if __name__ == '__main__':
    main()
