#!/usr/bin/env python3
"""
Fix animal virus contamination in human-associated collections.

Problem: Many curation scripts use family-level taxonomy queries
(e.g., WHERE t.family = 'Adenoviridae') without filtering for human
viruses. This pulls animal viruses (bat, bovine, porcine, avian, etc.)
into human virome collections.

Additionally, human herpesviruses are under family='Unknown' in the
taxonomy table (not Orthoherpesviridae), so family-based queries
for herpesviruses return only animal herpesviruses.

This script:
1. Identifies animal/non-human eukaryotic viruses in each collection
2. Replaces them with correctly-selected human viruses
3. Reports all changes made

Usage:
    python scripts/fix_animal_virus_contamination.py [--dry-run] [--database PATH]
"""

import sqlite3
import argparse
import logging
import random
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DB_PATH = Path(__file__).parent.parent / 'viroforge' / 'data' / 'viral_genomes.db'

# ============================================================
# Human-specific genome name patterns
# ============================================================
# Each entry: (description, SQL WHERE clause fragment for genome_name)
# These patterns identify viruses known to infect humans.

HUMAN_VIRUS_PATTERNS = {
    'Adenoviridae': {
        'include': [
            "g.genome_name LIKE 'Human adenovirus%'",
            "g.genome_name LIKE 'Human mastadenovirus%'",
        ],
        'exclude': [],
    },
    'Herpesviridae': {
        # Human herpesviruses are under family='Unknown' in taxonomy!
        # Must query by genome_name, not family.
        'include': [
            "g.genome_name LIKE 'Human herpesvirus%'",
            "g.genome_name LIKE 'Human alphaherpesvirus%'",
            "g.genome_name LIKE 'Human betaherpesvirus%'",
            "g.genome_name LIKE 'Human gammaherpesvirus%'",
        ],
        'exclude': [],
    },
    'Papillomaviridae': {
        'include': [
            "g.genome_name LIKE 'Human papillomavirus%'",
            "g.genome_name LIKE 'Alphapapillomavirus%'",
            "g.genome_name LIKE 'Betapapillomavirus%'",
            "g.genome_name LIKE 'Gammapapillomavirus%'",
        ],
        'exclude': [
            "g.genome_name NOT LIKE '%Felis%'",
            "g.genome_name NOT LIKE '%Rattus%'",
            "g.genome_name NOT LIKE '%Bos%'",
            "g.genome_name NOT LIKE '%Canis%'",
            "g.genome_name NOT LIKE '%Equus%'",
            "g.genome_name NOT LIKE '%Trichechus%'",
            "g.genome_name NOT LIKE '%Rousettus%'",
            "g.genome_name NOT LIKE '%Alces%'",
            "g.genome_name NOT LIKE '%Mustela%'",
        ],
    },
    'Polyomaviridae': {
        'include': [
            "g.genome_name LIKE 'Human polyomavirus%'",
            "g.genome_name LIKE 'Merkel cell polyomavirus%'",
            "g.genome_name LIKE 'JC polyomavirus%'",
            "g.genome_name LIKE 'BK polyomavirus%'",
            "g.genome_name LIKE 'Trichodysplasia spinulosa%'",
            "g.genome_name LIKE 'MW polyomavirus%'",
            "g.genome_name LIKE 'STL polyomavirus%'",
        ],
        'exclude': [],
    },
    'Coronaviridae': {
        'include': [
            "g.genome_name LIKE 'Human coronavirus%'",
            "g.genome_name LIKE 'SARS-CoV-2%'",
            "g.genome_name LIKE 'Severe acute respiratory syndrome%'",
            "g.genome_name LIKE 'Middle East respiratory syndrome%'",
        ],
        'exclude': [],
    },
    'Orthomyxoviridae': {
        'include': [
            "g.genome_name LIKE 'Influenza A virus%'",
            "g.genome_name LIKE 'Influenza B virus%'",
            "g.genome_name LIKE 'Influenza C virus%'",
        ],
        'exclude': [
            # Exclude clearly avian/swine strains by host indicator
            "g.genome_name NOT LIKE '%/Goose/%'",
            "g.genome_name NOT LIKE '%/Duck/%'",
            "g.genome_name NOT LIKE '%/Chicken/%'",
            "g.genome_name NOT LIKE '%/Swine/%'",
            "g.genome_name NOT LIKE '%/Mallard/%'",
            "g.genome_name NOT LIKE '%/Turkey/%'",
            "g.genome_name NOT LIKE '%/Equine/%'",
        ],
    },
    'Anelloviridae': {
        'include': [
            # Human TTV: named "Torque teno virus N" (plain number)
            "g.genome_name LIKE 'Torque teno virus %'",
            "g.genome_name LIKE 'Torque teno mini virus%'",
            "g.genome_name LIKE 'Torque teno midi virus%'",
        ],
        'exclude': [
            "g.genome_name NOT LIKE '%canis%'",
            "g.genome_name NOT LIKE '%felis%'",
            "g.genome_name NOT LIKE '%tamarin%'",
            "g.genome_name NOT LIKE '%indri%'",
            "g.genome_name NOT LIKE '%simian%'",
            "g.genome_name NOT LIKE '%Simian%'",
            "g.genome_name NOT LIKE '%porcine%'",
            "g.genome_name NOT LIKE '%Porcine%'",
            "g.genome_name NOT LIKE '%bovine%'",
            "g.genome_name NOT LIKE '%Rodent%'",
            "g.genome_name NOT LIKE '%rodent%'",
            "g.genome_name NOT LIKE 'Chicken%'",
            "g.genome_name NOT LIKE 'Avian%'",
            "g.genome_name NOT LIKE 'Gyrovirus%'",
        ],
    },
    'Parvoviridae': {
        'include': [
            "g.genome_name LIKE 'Human bocavirus%'",
            "g.genome_name LIKE 'Human parvovirus%'",
            "g.genome_name LIKE 'Adeno-associated virus%'",
            "g.genome_name LIKE 'Primate erythroparvovirus%'",
            "g.genome_name LIKE 'Human erythrovirus%'",
        ],
        'exclude': [],
    },
    'Paramyxoviridae': {
        'include': [
            "g.genome_name LIKE 'Human parainfluenza%'",
            "g.genome_name LIKE 'Human respirovirus%'",
            "g.genome_name LIKE 'Human metapneumovirus%'",
            "g.genome_name LIKE 'Human orthopneumovirus%'",
            "g.genome_name LIKE 'Respiratory syncytial virus%'",
            "g.genome_name LIKE 'Human respiratory syncytial%'",
            "g.genome_name LIKE 'Measles%'",
            "g.genome_name LIKE 'Mumps%'",
        ],
        'exclude': [
            "g.genome_name NOT LIKE '%Bovine%'",
            "g.genome_name NOT LIKE '%bovine%'",
            "g.genome_name NOT LIKE '%Porcine%'",
            "g.genome_name NOT LIKE '%porcine%'",
            "g.genome_name NOT LIKE '%Avian%'",
            "g.genome_name NOT LIKE '%avian%'",
        ],
    },
    'Picornaviridae': {
        'include': [
            "g.genome_name LIKE 'Enterovirus %'",
            "g.genome_name LIKE 'Human enterovirus%'",
            "g.genome_name LIKE 'Human rhinovirus%'",
            "g.genome_name LIKE 'Rhinovirus %'",
            "g.genome_name LIKE 'Hepatitis A virus%'",
            "g.genome_name LIKE 'Aichivirus A%'",
            "g.genome_name LIKE 'Human parechovirus%'",
            "g.genome_name LIKE 'Salivirus%'",
            "g.genome_name LIKE 'Cosavirus%'",
            "g.genome_name LIKE 'Cardiovirus%'",
        ],
        'exclude': [
            "g.genome_name NOT LIKE '%Possum%'",
            "g.genome_name NOT LIKE '%possum%'",
            "g.genome_name NOT LIKE '%Bovine%'",
            "g.genome_name NOT LIKE '%bovine%'",
            "g.genome_name NOT LIKE '%Porcine%'",
            "g.genome_name NOT LIKE '%porcine%'",
            "g.genome_name NOT LIKE '%Simian%'",
            "g.genome_name NOT LIKE '%Aichivirus B%'",
            "g.genome_name NOT LIKE '%Aichivirus C%'",
        ],
    },
    'Astroviridae': {
        'include': [
            "g.genome_name LIKE 'Human astrovirus%'",
            "g.genome_name LIKE 'Mamastrovirus 1%'",
        ],
        'exclude': [],
    },
    'Caliciviridae': {
        'include': [
            "g.genome_name LIKE '%Norovirus%'",
            "g.genome_name LIKE '%Norwalk%'",
            "g.genome_name LIKE '%Sapovirus%'",
            "g.genome_name LIKE 'Human calicivirus%'",
        ],
        'exclude': [
            "g.genome_name NOT LIKE '%Murine%'",
            "g.genome_name NOT LIKE '%murine%'",
            "g.genome_name NOT LIKE '%Bovine%'",
            "g.genome_name NOT LIKE '%Porcine%'",
            "g.genome_name NOT LIKE '%Feline%'",
        ],
    },
    'Reoviridae_Sedoreoviridae': {
        'include': [
            "g.genome_name LIKE 'Rotavirus A%'",
            "g.genome_name LIKE 'Human rotavirus%'",
        ],
        'exclude': [],
    },
}

# Patterns to detect non-human eukaryotic viruses (for flagging)
ANIMAL_INDICATORS = [
    'Bat ', 'bat ', 'Bovine', 'bovine', 'Porcine', 'porcine',
    'Canine', 'canine', 'Feline', 'feline', 'Murine', 'murine',
    'Avian', 'avian', 'Simian', 'simian', 'Equine', 'equine',
    'Turkey ', 'Duck ', 'Chicken ', 'Goose', 'Rodent',
    'California sea lion', 'Tree shrew', 'Bottlenose dolphin',
    'Psittac', 'Ranid', 'Cacatuid', 'Macropodid', 'Phascolarctid',
    'Sphenicid', 'Vombatid', 'Colobine', 'Cercopithecus',
    'Pan troglodytes', 'Rousettus', 'Alces alces', 'Mustela',
    'Trichechus', 'Piliocolobus', 'Ursus americanus',
    'Roe deer', 'Mink', 'Ovine', 'Otarine', 'Possum',
    'Chinook salmon', 'Fathead minnow', 'Callinectes', 'Eriocheir',
    'Micromonas pusilla', 'Beluga', 'Ball python', 'Veiled chameleon',
    'Bellinger River', 'Infectious salmon',
    'Torque teno canis', 'Torque teno felis', 'Torque teno tamarin',
    'Torque teno indri', 'Gyrovirus', 'Chicken anemia',
]

# Phage families that are NOT eukaryotic viruses (keep these)
PHAGE_FAMILIES = {
    'Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae',
    'Inoviridae', 'Tectiviridae', 'Leviviridae', 'Cystoviridae',
    'Corticoviridae', 'Plasmaviridae', 'Fuselloviridae',
    'Steitzviridae', 'Fiersviridae', 'Autographiviridae',
    'Demerecviridae', 'Drexlerviridae', 'Herelleviridae',
    'Ackermannviridae', 'Chaseviridae', 'Schitoviridae',
    'Zobellviridae', 'Peduoviridae', 'Casjensviridae',
    'Intestiviridae', 'Suoliviridae', 'Salasmaviridae',
    'Straboviridae', 'Zierdtviridae',
}

# Plant virus indicators
PLANT_INDICATORS = [
    'Lettuce', 'Peanut', 'Sida ', 'Trachyspermum', 'Kalanchoe',
    'Cereal', 'Sclerotinia', 'Sorghum', 'Sweet potato', 'Pineapple',
    'Asian prunus', 'Groundnut', 'Potato', 'Citrus', 'Commelina',
    'Phaius', 'Plasmopara', 'Cassava',
]

# Lab/spike-in phage to exclude from all collections
LAB_PHAGES = {
    'GCF_000819615.1',  # PhiX174 (Illumina spike-in)
}

# Collections that are human-associated (need human virus filtering)
# Environmental collections (marine, soil, freshwater, wastewater, arbovirus)
# do not need this filtering.
HUMAN_COLLECTIONS = {
    1: 'Gut Virome',
    2: 'Oral Virome',
    3: 'Skin Virome',
    4: 'Respiratory Virome',
    8: 'Mouse Gut Virome',
    10: 'IBD Gut Virome',
    11: 'HIV+ Gut Virome',
    12: 'CF Respiratory Virome',
    13: 'Human Respiratory RNA Virome',
    15: 'Fecal RNA Virome',
    16: 'Vaginal Virome',
    17: 'Blood/Plasma Virome',
    18: 'Ocular Surface Virome',
    19: 'Lung Virome',
    20: 'Urinary Virome',
}


def get_collection_ids(conn: sqlite3.Connection) -> Dict[int, str]:
    """Get all collection IDs and names from database."""
    cursor = conn.execute(
        "SELECT collection_id, collection_name FROM body_site_collections ORDER BY collection_id"
    )
    return {row[0]: row[1] for row in cursor.fetchall()}


def find_suspicious_genomes(
    conn: sqlite3.Connection,
    collection_id: int
) -> List[Dict]:
    """Find non-human eukaryotic viruses in a human-associated collection."""
    cursor = conn.execute("""
        SELECT cg.genome_id, g.genome_name, t.family, t.genus,
               cg.relative_abundance, cg.abundance_rank
        FROM collection_genomes cg
        JOIN genomes g ON cg.genome_id = g.genome_id
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE cg.collection_id = ?
        ORDER BY g.genome_name
    """, (collection_id,))

    suspicious = []
    for row in cursor.fetchall():
        genome_id, genome_name, family, genus, abundance, rank = row

        # Skip phage families (they're expected in all collections)
        if family in PHAGE_FAMILIES:
            continue

        # Check for lab phages
        if genome_id in LAB_PHAGES:
            suspicious.append({
                'genome_id': genome_id,
                'genome_name': genome_name,
                'family': family,
                'genus': genus,
                'abundance': abundance,
                'rank': rank,
                'reason': 'Lab spike-in phage (PhiX174)',
            })
            continue

        # Check for animal virus indicators
        for indicator in ANIMAL_INDICATORS:
            if indicator in genome_name:
                suspicious.append({
                    'genome_id': genome_id,
                    'genome_name': genome_name,
                    'family': family,
                    'genus': genus,
                    'abundance': abundance,
                    'rank': rank,
                    'reason': f'Animal virus (matched: "{indicator}")',
                })
                break

        # Check for plant virus indicators
        for indicator in PLANT_INDICATORS:
            if indicator in genome_name:
                suspicious.append({
                    'genome_id': genome_id,
                    'genome_name': genome_name,
                    'family': family,
                    'genus': genus,
                    'abundance': abundance,
                    'rank': rank,
                    'reason': f'Plant virus (matched: "{indicator}")',
                })
                break

    return suspicious


    # Map taxonomy families in DB to our pattern keys
FAMILY_TO_PATTERN_KEY = {
    'Adenoviridae': 'Adenoviridae',
    'Orthoherpesviridae': 'Herpesviridae',
    'Alloherpesviridae': 'Herpesviridae',
    'Papillomaviridae': 'Papillomaviridae',
    'Polyomaviridae': 'Polyomaviridae',
    'Coronaviridae': 'Coronaviridae',
    'Orthomyxoviridae': 'Orthomyxoviridae',
    'Anelloviridae': 'Anelloviridae',
    'Parvoviridae': 'Parvoviridae',
    'Paramyxoviridae': 'Paramyxoviridae',
    'Picornaviridae': 'Picornaviridae',
    'Astroviridae': 'Astroviridae',
    'Caliciviridae': 'Caliciviridae',
    'Sedoreoviridae': 'Reoviridae_Sedoreoviridae',
    'Tobaniviridae': None,  # No human replacements in this family
    'Rhabdoviridae': None,  # Plant rhabdoviruses, no direct replacement
    'Betaflexiviridae': None,
    'Bromoviridae': None,
    'Caulimoviridae': None,
    'Potyviridae': None,
    'Alphaflexiviridae': None,
    'Botourmiaviridae': None,
    'Secoviridae': None,
    'Tospoviridae': None,
    'Unknown': None,  # Will try Herpesviridae patterns for herpes-like names
}


def find_human_replacements(
    conn: sqlite3.Connection,
    family: str,
    count: int,
    already_in_collection: Set[str],
) -> List[str]:
    """Find human-specific viruses from a family to use as replacements."""
    # Map DB family name to our pattern key
    pattern_key = FAMILY_TO_PATTERN_KEY.get(family, family)
    if pattern_key is None:
        return []
    patterns = HUMAN_VIRUS_PATTERNS.get(pattern_key)
    if not patterns:
        return []

    include_clauses = " OR ".join(patterns['include'])
    exclude_clauses = " AND ".join(patterns['exclude']) if patterns['exclude'] else "1=1"

    query = f"""
        SELECT g.genome_id, g.genome_name
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE ({include_clauses})
        AND {exclude_clauses}
        AND g.genome_id NOT IN ({','.join('?' for _ in already_in_collection)})
        ORDER BY RANDOM()
        LIMIT ?
    """
    params = list(already_in_collection) + [count]
    cursor = conn.execute(query, params)
    return [row[0] for row in cursor.fetchall()]


def fix_collection(
    conn: sqlite3.Connection,
    collection_id: int,
    collection_name: str,
    dry_run: bool = True
) -> Dict:
    """Fix a single collection by replacing animal viruses with human ones."""
    suspicious = find_suspicious_genomes(conn, collection_id)

    if not suspicious:
        logger.info(f"  Collection {collection_id} ({collection_name}): No issues found")
        return {'collection_id': collection_id, 'removed': [], 'added': [], 'kept': []}

    logger.info(f"  Collection {collection_id} ({collection_name}): "
                f"Found {len(suspicious)} suspicious genomes")

    # Get all current genome IDs in collection
    cursor = conn.execute(
        "SELECT genome_id FROM collection_genomes WHERE collection_id = ?",
        (collection_id,)
    )
    current_ids = {row[0] for row in cursor.fetchall()}

    # Group suspicious genomes by family for replacement
    by_family = {}
    for s in suspicious:
        family = s['family'] or 'Unknown'
        if family not in by_family:
            by_family[family] = []
        by_family[family].append(s)

    removed = []
    added = []
    kept = []

    for family, items in by_family.items():
        # Try to find human replacements
        ids_to_remove = {s['genome_id'] for s in items}
        remaining_ids = current_ids - ids_to_remove

        replacements = find_human_replacements(
            conn, family, len(items), remaining_ids
        )

        if replacements:
            for s in items:
                removed.append(s)
                logger.info(f"    REMOVE: {s['genome_name']} ({s['reason']})")

            for r_id in replacements:
                cursor = conn.execute(
                    "SELECT genome_name FROM genomes WHERE genome_id = ?", (r_id,)
                )
                r_name = cursor.fetchone()[0]
                added.append({'genome_id': r_id, 'genome_name': r_name, 'family': family})
                logger.info(f"    ADD:    {r_name}")

            current_ids = (current_ids - ids_to_remove) | set(replacements)
        else:
            # No human replacements available — remove plant/animal viruses
            # that clearly don't belong, without replacement
            for s in items:
                is_plant = s['reason'].startswith('Plant virus')
                is_lab = s['reason'].startswith('Lab spike-in')
                is_animal = s['reason'].startswith('Animal virus')

                if is_lab or is_plant or is_animal:
                    removed.append(s)
                    current_ids.discard(s['genome_id'])
                    logger.info(f"    REMOVE (no replacement): {s['genome_name']} ({s['reason']})")
                else:
                    kept.append(s)
                    logger.warning(
                        f"    KEEP (unclear): {s['genome_name']} "
                        f"({s['reason']}, family={family})"
                    )

    # Apply changes
    if not dry_run and (removed or added):
        for s in removed:
            conn.execute(
                "DELETE FROM collection_genomes WHERE collection_id = ? AND genome_id = ?",
                (collection_id, s['genome_id'])
            )

        for a in added:
            # Use average abundance from removed genomes of same family
            family_removed = [r for r in removed if r['family'] == a['family']]
            avg_abundance = sum(r['abundance'] for r in family_removed) / max(len(family_removed), 1)
            avg_rank = int(sum(r['rank'] for r in family_removed) / max(len(family_removed), 1))

            conn.execute(
                "INSERT INTO collection_genomes (collection_id, genome_id, relative_abundance, abundance_rank) "
                "VALUES (?, ?, ?, ?)",
                (collection_id, a['genome_id'], avg_abundance, avg_rank)
            )

        # Update genome count
        cursor = conn.execute(
            "SELECT COUNT(*) FROM collection_genomes WHERE collection_id = ?",
            (collection_id,)
        )
        new_count = cursor.fetchone()[0]
        conn.execute(
            "UPDATE body_site_collections SET n_genomes = ? WHERE collection_id = ?",
            (new_count, collection_id)
        )

        conn.commit()
        logger.info(f"  Applied: removed {len(removed)}, added {len(added)}, "
                     f"kept {len(kept)}, new total: {new_count}")

    return {
        'collection_id': collection_id,
        'collection_name': collection_name,
        'removed': removed,
        'added': added,
        'kept': kept,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Fix animal virus contamination in human-associated collections'
    )
    parser.add_argument('--database', type=str, default=str(DB_PATH),
                        help='Path to viral_genomes.db')
    parser.add_argument('--dry-run', action='store_true', default=True,
                        help='Report issues without making changes (default)')
    parser.add_argument('--apply', action='store_true',
                        help='Actually apply fixes to the database')
    parser.add_argument('--collection-id', type=int, default=None,
                        help='Fix a specific collection only')
    args = parser.parse_args()

    if args.apply:
        args.dry_run = False

    db_path = Path(args.database)
    if not db_path.exists():
        logger.error(f"Database not found: {db_path}")
        return

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    # Get all collections
    all_collections = get_collection_ids(conn)

    if args.collection_id:
        if args.collection_id not in all_collections:
            logger.error(f"Collection {args.collection_id} not found")
            return
        collections_to_check = {args.collection_id: all_collections[args.collection_id]}
    else:
        # Only check human-associated collections
        collections_to_check = {
            cid: name for cid, name in all_collections.items()
            if cid in HUMAN_COLLECTIONS
        }

    mode = "DRY RUN" if args.dry_run else "APPLYING FIXES"
    logger.info(f"{'='*70}")
    logger.info(f"Animal Virus Contamination Audit ({mode})")
    logger.info(f"{'='*70}")
    logger.info(f"Checking {len(collections_to_check)} human-associated collections\n")

    all_results = []
    total_removed = 0
    total_added = 0

    for cid, cname in sorted(collections_to_check.items()):
        result = fix_collection(conn, cid, cname, dry_run=args.dry_run)
        all_results.append(result)
        total_removed += len(result['removed'])
        total_added += len(result['added'])
        print()

    # Summary
    logger.info(f"{'='*70}")
    logger.info(f"SUMMARY")
    logger.info(f"{'='*70}")
    logger.info(f"Collections checked: {len(collections_to_check)}")
    logger.info(f"Total animal/plant viruses found: {total_removed}")
    logger.info(f"Total human replacements: {total_added}")

    print(f"\n{'Collection':<35s} {'Removed':>8s} {'Added':>8s} {'Kept':>8s}")
    print("-" * 65)
    for r in all_results:
        if r['removed'] or r['kept']:
            cname = r.get('collection_name', f"ID {r['collection_id']}")
            print(f"{cname:<35s} "
                  f"{len(r['removed']):>8d} {len(r['added']):>8d} {len(r['kept']):>8d}")

    if args.dry_run:
        print(f"\n*** DRY RUN — no changes made. Use --apply to fix. ***")

    conn.close()


if __name__ == '__main__':
    main()
