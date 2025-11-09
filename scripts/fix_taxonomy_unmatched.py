#!/usr/bin/env python3
"""
Fix taxonomy for unmatched genomes using fuzzy matching.

This script implements enhanced matching strategies to assign ICTV taxonomy
to genomes that failed initial matching, particularly for:
- Strain-specific names (e.g., Influenza A virus (A/California/07/2009))
- Isolate-specific names
- Case sensitivity issues
- Segmented virus genomes

Author: ViroForge Development Team
Date: 2025-11-09
"""

import argparse
import sqlite3
import json
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class TaxonomyFixer:
    """Fix taxonomy assignments using enhanced fuzzy matching."""

    def __init__(self, db_path: str, ictv_path: str):
        self.db_path = db_path
        self.ictv_path = ictv_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row

        # Load ICTV taxonomy
        logger.info(f"Loading ICTV taxonomy from {ictv_path}")
        with open(ictv_path, 'r') as f:
            self.ictv_records = json.load(f)
        logger.info(f"  Loaded {len(self.ictv_records)} ICTV records")

        # Build lookup structures
        self.ictv_lookup = self._build_ictv_lookup()

    def _build_ictv_lookup(self) -> Dict[str, Dict]:
        """Build multiple lookup structures for ICTV data."""
        lookup = {
            'species': {},      # Species name -> record
            'virus_name': {},   # Virus common name -> record
            'genus': {},        # Genus -> list of records
            'family': {}        # Family -> list of records
        }

        for record in self.ictv_records:
            species = (record.get('species') or '').lower().strip()
            virus_name = (record.get('virus_name') or '').lower().strip()
            genus = (record.get('genus') or '').lower().strip()
            family = (record.get('family') or '').lower().strip()

            if species:
                lookup['species'][species] = record

            if virus_name:
                lookup['virus_name'][virus_name] = record
                # Also without "virus" suffix
                virus_name_noVirus = virus_name.replace(' virus', '').strip()
                if virus_name_noVirus:
                    lookup['virus_name'][virus_name_noVirus] = record

            if genus:
                if genus not in lookup['genus']:
                    lookup['genus'][genus] = []
                lookup['genus'][genus].append(record)

            if family:
                if family not in lookup['family']:
                    lookup['family'][family] = []
                lookup['family'][family].append(record)

        logger.info(f"  Built lookup: {len(lookup['species'])} species, {len(lookup['virus_name'])} virus names")
        return lookup

    def normalize_virus_name(self, name: str) -> str:
        """
        Normalize virus name for better matching.

        Handles:
        - Strain/isolate information: (A/California/07/2009), isolate XYZ, strain ABC
        - Segment information: segment 1, RNA 1
        - Type keywords: "type 60" -> "60"
        - Trailing numbers: "virus - 4" -> "virus", "virus 1" -> "virus"
        - Subtype information: "subtype 1" -> ""
        - Case normalization
        """
        # Remove strain/isolate in parentheses
        cleaned = re.sub(r'\([^)]*\)', '', name)

        # Remove "isolate" and following word
        cleaned = re.sub(r'\s+isolate\s+\S+', '', cleaned, flags=re.IGNORECASE)

        # Remove "strain" and following word (including strain names after "strain")
        cleaned = re.sub(r'\s+strain\s+\S+', '', cleaned, flags=re.IGNORECASE)

        # Remove segment information
        cleaned = re.sub(r'\s+segment\s+\d+', '', cleaned, flags=re.IGNORECASE)
        cleaned = re.sub(r'\s+RNA\s+\d+', '', cleaned, flags=re.IGNORECASE)

        # Remove "type" keyword (e.g., "papillomavirus type 60" -> "papillomavirus 60")
        cleaned = re.sub(r'\s+type\s+', ' ', cleaned, flags=re.IGNORECASE)

        # Remove "subtype" keyword
        cleaned = re.sub(r'\s+subtype\s+\S+', '', cleaned, flags=re.IGNORECASE)

        # Remove trailing dash-number patterns: " - 4", " - 3"
        cleaned = re.sub(r'\s*-\s*\d+\s*$', '', cleaned)

        # Remove trailing single letter/number after space: "virus A", "virus 1"
        # But be conservative - only if it ends with single char/digit
        cleaned = re.sub(r'\s+[A-Z0-9]$', '', cleaned, flags=re.IGNORECASE)

        # Remove extra whitespace
        cleaned = ' '.join(cleaned.split())

        return cleaned.strip().lower()

    def clean_species_name(self, species: str) -> str:
        """Legacy method - calls normalize_virus_name for compatibility."""
        return self.normalize_virus_name(species)

    def match_by_pattern(self, species_name: str) -> Optional[Dict]:
        """
        Match by virus name patterns to infer family.

        Common patterns:
        - *herpesvirus* -> Orthoherpesviridae, Alloherpesviridae
        - *papillomavirus* -> Papillomaviridae
        - *adenovirus*, *mastadenovirus*, *adeno-associated* -> Adenoviridae or Parvoviridae
        - *parvovirus*, *densovirus* -> Parvoviridae
        - *polyomavirus* -> Polyomaviridae
        - *poxvirus*, *fibroma virus*, *myxoma* -> Poxviridae
        - *retrovirus*, *lentivirus* -> Retroviridae
        - *nucleopolyhedrovirus*, *granulovirus* -> Baculoviridae
        - *picornavirus*, *enterovirus*, *cardiovirus* -> Picornaviridae
        - *respirovirus*, *morbillivirus* -> Paramyxoviridae
        - *cytomegalovirus* -> Orthoherpesviridae
        """
        species_lower = species_name.lower()

        # Define pattern-to-family mappings
        # We'll search for any ICTV record from that family and use it as template
        family_patterns = {
            # Herpesviruses
            'herpesvirus': ['Orthoherpesviridae', 'Alloherpesviridae'],
            'cytomegalovirus': ['Orthoherpesviridae'],

            # Papillomaviruses and polyomaviruses
            'papillomavirus': ['Papillomaviridae'],
            'polyomavirus': ['Polyomaviridae'],

            # Adenoviruses and parvoviruses
            'adenovirus': ['Adenoviridae'],
            'mastadenovirus': ['Adenoviridae'],
            'adeno-associated': ['Parvoviridae'],  # AAV
            'parvovirus': ['Parvoviridae'],
            'densovirus': ['Parvoviridae'],

            # Poxviruses
            'poxvirus': ['Poxviridae'],
            'fibroma virus': ['Poxviridae'],
            'myxoma': ['Poxviridae'],

            # Retroviruses
            'retrovirus': ['Retroviridae'],
            'lentivirus': ['Retroviridae'],

            # Baculoviruses
            'nucleopolyhedrovirus': ['Baculoviridae'],
            'granulovirus': ['Baculoviridae'],

            # Picornaviruses
            'picornavirus': ['Picornaviridae'],
            'enterovirus': ['Picornaviridae'],
            'cardiovirus': ['Picornaviridae'],

            # Paramyxoviruses
            'respirovirus': ['Paramyxoviridae'],
            'morbillivirus': ['Paramyxoviridae'],
        }

        for pattern, families in family_patterns.items():
            if pattern in species_lower:
                # Find any ICTV record with this family
                for family_name in families:
                    family_key = family_name.lower()
                    if family_key in self.ictv_lookup['family']:
                        # Return first record from this family
                        # This gives us the family hierarchy
                        return self.ictv_lookup['family'][family_key][0]

        return None

    def match_species(self, species_name: str) -> Optional[Dict]:
        """
        Match a species name to ICTV using multiple strategies.

        Returns ICTV record if matched, None otherwise.
        """
        species_lower = species_name.lower().strip()

        # Strategy 1: Direct match on species name
        if species_lower in self.ictv_lookup['species']:
            return self.ictv_lookup['species'][species_lower]

        # Strategy 2: Direct match on virus_name
        if species_lower in self.ictv_lookup['virus_name']:
            return self.ictv_lookup['virus_name'][species_lower]

        # Strategy 3: Normalized match (remove strain/isolate/type keywords)
        normalized = self.normalize_virus_name(species_name)
        if normalized != species_lower:
            # Try normalized against species
            if normalized in self.ictv_lookup['species']:
                return self.ictv_lookup['species'][normalized]

            # Try normalized against virus_name
            if normalized in self.ictv_lookup['virus_name']:
                return self.ictv_lookup['virus_name'][normalized]

            # Also try without "virus" suffix
            normalized_noVirus = normalized.replace(' virus', '').strip()
            if normalized_noVirus in self.ictv_lookup['virus_name']:
                return self.ictv_lookup['virus_name'][normalized_noVirus]

        # Strategy 4: Partial match - extract base name before qualifiers
        base_name = species_lower.split('(')[0].split(' isolate')[0].split(' strain')[0].strip()
        if base_name in self.ictv_lookup['virus_name']:
            return self.ictv_lookup['virus_name'][base_name]

        # Also try normalized base name
        base_normalized = self.normalize_virus_name(base_name)
        if base_normalized in self.ictv_lookup['virus_name']:
            return self.ictv_lookup['virus_name'][base_normalized]

        # Strategy 5: Pattern-based family matching
        # This handles herpesviruses, papillomaviruses, etc.
        pattern_match = self.match_by_pattern(species_name)
        if pattern_match:
            return pattern_match

        return None

    def fix_unmatched_genomes(self) -> Dict[str, int]:
        """
        Fix taxonomy for genomes with family='Unknown'.

        Returns statistics about fixes.
        """
        logger.info("=" * 80)
        logger.info("FIXING UNMATCHED GENOMES")
        logger.info("=" * 80)

        # Get genomes with Unknown family
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT genome_id, species
            FROM taxonomy
            WHERE family = 'Unknown' OR family IS NULL
        """)
        unmatched = cursor.fetchall()

        logger.info(f"Found {len(unmatched)} genomes with Unknown family")

        stats = {
            'total_unmatched': len(unmatched),
            'fixed': 0,
            'still_unmatched': 0
        }

        fixes = []  # (genome_id, matched_record)
        still_unmatched = []

        for genome_id, species_name in unmatched:
            matched_record = self.match_species(species_name)

            # Only use match if it has a family (required by database)
            if matched_record and matched_record.get('family'):
                fixes.append((genome_id, matched_record))
                stats['fixed'] += 1
            else:
                still_unmatched.append((genome_id, species_name))
                stats['still_unmatched'] += 1

        logger.info(f"\nMatching results:")
        logger.info(f"  Fixed: {stats['fixed']} ({stats['fixed']/stats['total_unmatched']*100:.1f}%)")
        logger.info(f"  Still unmatched: {stats['still_unmatched']} ({stats['still_unmatched']/stats['total_unmatched']*100:.1f}%)")

        # Apply fixes to database
        if fixes:
            logger.info(f"\nUpdating {len(fixes)} taxonomy records...")
            for genome_id, record in fixes:
                cursor.execute("""
                    UPDATE taxonomy
                    SET realm = ?,
                        kingdom = ?,
                        phylum = ?,
                        class = ?,
                        order_name = ?,
                        family = ?,
                        subfamily = ?,
                        genus = ?
                    WHERE genome_id = ?
                """, (
                    record.get('realm'),
                    record.get('kingdom'),
                    record.get('phylum'),
                    record.get('class'),
                    record.get('order'),
                    record.get('family'),
                    record.get('subfamily'),
                    record.get('genus'),
                    genome_id
                ))

            self.conn.commit()
            logger.info("✓ Database updated")

        # Show examples of fixes
        if fixes:
            logger.info("\nExample fixes:")
            for i, (genome_id, record) in enumerate(fixes[:10], 1):
                # Get species name
                cursor.execute("SELECT species FROM taxonomy WHERE genome_id = ?", (genome_id,))
                species = cursor.fetchone()[0]
                logger.info(f"  {i}. {species[:60]:60s} -> {record.get('family')}")

        # Save still unmatched for review
        if still_unmatched:
            output_file = Path("data/ictv/still_unmatched_after_fix.tsv")
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                f.write("genome_id\tspecies_name\n")
                for genome_id, species_name in still_unmatched:
                    f.write(f"{genome_id}\t{species_name}\n")
            logger.info(f"\n✓ Saved still unmatched genomes: {output_file}")

        return stats

    def verify_fixes(self):
        """Verify that fixes were applied correctly."""
        logger.info("\n" + "=" * 80)
        logger.info("VERIFICATION")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        # Count by family
        cursor.execute("""
            SELECT family, COUNT(*) as count
            FROM taxonomy
            GROUP BY family
            ORDER BY count DESC
            LIMIT 15
        """)

        logger.info("\nTop 15 families after fix:")
        for row in cursor.fetchall():
            family, count = row
            logger.info(f"  {family:30s} {count:6d} genomes")

        # Check influenza specifically
        cursor.execute("""
            SELECT COUNT(*), family
            FROM taxonomy
            WHERE species LIKE '%influenza%'
            GROUP BY family
        """)

        logger.info("\nInfluenza virus family assignments:")
        for row in cursor.fetchall():
            count, family = row
            logger.info(f"  {family:30s} {count:6d} genomes")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main workflow."""
    parser = argparse.ArgumentParser(description="Fix taxonomy for unmatched genomes")
    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to genome database'
    )
    parser.add_argument(
        '--ictv',
        default='data/ictv/ictv_taxonomy.json',
        help='Path to ICTV taxonomy JSON'
    )

    args = parser.parse_args()

    fixer = TaxonomyFixer(args.database, args.ictv)

    try:
        # Fix unmatched genomes
        stats = fixer.fix_unmatched_genomes()

        # Verify fixes
        fixer.verify_fixes()

        logger.info("\n" + "=" * 80)
        logger.info("✓ TAXONOMY FIX COMPLETE!")
        logger.info("=" * 80)
        logger.info(f"Fixed: {stats['fixed']} genomes")
        logger.info(f"Still unmatched: {stats['still_unmatched']} genomes")

    finally:
        fixer.close()


if __name__ == '__main__':
    main()
