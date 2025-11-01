#!/usr/bin/env python3
"""
Parse ICTV Virus Metadata Resource (VMR) and map to RefSeq genomes.

This script reads the ICTV VMR Excel file, extracts complete taxonomy
hierarchy, and maps it to our RefSeq viral genomes for database population.

Usage:
    python scripts/parse_ictv_taxonomy.py --vmr data/ictv/VMR_current.xlsx --output data/ictv
    python scripts/parse_ictv_taxonomy.py --vmr data/ictv/VMR_current.xlsx --database viroforge/data/viral_genomes.db

Author: ViroForge Development Team
Date: November 1, 2025
"""

import argparse
import openpyxl
import json
import csv
import sqlite3
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ICTVTaxonomyParser:
    """Parse ICTV VMR and map taxonomy to viral genomes."""

    # Taxonomy rank columns in VMR
    TAXONOMY_COLUMNS = {
        'realm': 3,         # Column D (index 3)
        'subrealm': 4,
        'kingdom': 5,
        'subkingdom': 6,
        'phylum': 7,
        'subphylum': 8,
        'class': 9,
        'subclass': 10,
        'order': 11,
        'suborder': 12,
        'family': 13,
        'subfamily': 14,
        'genus': 15,
        'subgenus': 16,
        'species': 17
    }

    # Other important columns
    VIRUS_NAME_COL = 20     # Column U (Virus name(s))
    GENBANK_ACC_COL = 23    # Column X (Virus GENBANK accession)
    HOST_SOURCE_COL = 26    # Column AA (Host source)

    def __init__(self, vmr_path: str, output_dir: str):
        """
        Initialize ICTV taxonomy parser.

        Parameters
        ----------
        vmr_path : str
            Path to ICTV VMR Excel file
        output_dir : str
            Directory to save parsed taxonomy data
        """
        self.vmr_path = Path(vmr_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if not self.vmr_path.exists():
            raise FileNotFoundError(f"VMR file not found: {vmr_path}")

    def extract_accession_from_hyperlink(self, cell_value: str) -> Optional[str]:
        """
        Extract VMR ID from hyperlink formula.

        Parameters
        ----------
        cell_value : str
            Cell value (may be hyperlink formula)

        Returns
        -------
        Optional[str]
            Extracted VMR ID or None
        """
        if not cell_value:
            return None

        # Check if it's a HYPERLINK formula
        if isinstance(cell_value, str) and cell_value.startswith('=HYPERLINK'):
            # Extract VMR ID from: =HYPERLINK("https://ictv.global/id/VMR1028127","VMR1028127")
            match = re.search(r'VMR\d+', cell_value)
            if match:
                return match.group(0)

        return str(cell_value) if cell_value else None

    def clean_taxonomy_value(self, value: any) -> Optional[str]:
        """
        Clean taxonomy value (remove None, empty strings).

        Parameters
        ----------
        value : any
            Raw taxonomy value

        Returns
        -------
        Optional[str]
            Cleaned value or None
        """
        if value is None or value == '' or str(value).lower() == 'none':
            return None
        return str(value).strip()

    def parse_vmr(self) -> List[Dict[str, any]]:
        """
        Parse ICTV VMR Excel file.

        Returns
        -------
        List[Dict[str, any]]
            List of virus taxonomy records
        """
        logger.info(f"Parsing ICTV VMR: {self.vmr_path}")

        wb = openpyxl.load_workbook(self.vmr_path, read_only=True, data_only=True)
        sheet = wb["VMR MSL40"]

        records = []
        headers_row = None

        for i, row in enumerate(sheet.iter_rows(values_only=True), 1):
            # Skip header row
            if i == 1:
                headers_row = row
                continue

            # Extract taxonomy hierarchy
            taxonomy = {}
            for rank, col_idx in self.TAXONOMY_COLUMNS.items():
                value = row[col_idx] if col_idx < len(row) else None
                taxonomy[rank] = self.clean_taxonomy_value(value)

            # Extract other metadata
            vmr_id = self.extract_accession_from_hyperlink(row[0]) if len(row) > 0 else None
            virus_name = row[self.VIRUS_NAME_COL] if len(row) > self.VIRUS_NAME_COL else None
            genbank_acc = row[self.GENBANK_ACC_COL] if len(row) > self.GENBANK_ACC_COL else None
            host_source = row[self.HOST_SOURCE_COL] if len(row) > self.HOST_SOURCE_COL else None

            # Only include if we have at least a species name
            if taxonomy['species']:
                record = {
                    'vmr_id': vmr_id,
                    'virus_name': virus_name,
                    'genbank_accession': genbank_acc,
                    'host_source': host_source,
                    **taxonomy  # Unpack all taxonomy ranks
                }
                records.append(record)

            # Progress update
            if i % 1000 == 0:
                logger.info(f"Parsed {i} rows...")

        wb.close()

        logger.info(f"✓ Parsed {len(records)} virus taxonomy records")
        return records

    def save_taxonomy_data(self, records: List[Dict[str, any]]) -> None:
        """
        Save parsed taxonomy data to JSON and TSV files.

        Parameters
        ----------
        records : List[Dict[str, any]]
            Parsed taxonomy records
        """
        # Save as JSON
        json_file = self.output_dir / "ictv_taxonomy.json"
        with open(json_file, 'w') as f:
            json.dump(records, f, indent=2)
        logger.info(f"✓ Saved taxonomy JSON: {json_file}")

        # Save as TSV
        tsv_file = self.output_dir / "ictv_taxonomy.tsv"
        if records:
            fieldnames = records[0].keys()
            with open(tsv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(records)
            logger.info(f"✓ Saved taxonomy TSV: {tsv_file}")

    def build_taxonomy_lookup(self, records: List[Dict[str, any]]) -> Dict[str, Dict]:
        """
        Build taxonomy lookup by species name and virus name.

        Parameters
        ----------
        records : List[Dict[str, any]]
            Parsed taxonomy records

        Returns
        -------
        Dict[str, Dict]
            Taxonomy lookup by species name and virus common name
        """
        lookup = {}
        for record in records:
            species = record['species']
            virus_name = record['virus_name']

            if species:
                # Add by species name
                if species not in lookup:
                    lookup[species] = record

                # Also add by virus common name if available
                if virus_name:
                    # virus_name format: "cowpox virus"
                    virus_name_clean = str(virus_name).lower().strip()
                    if virus_name_clean not in lookup:
                        lookup[virus_name_clean] = record

                    # Also try without "virus" suffix
                    virus_name_noVirus = virus_name_clean.replace(' virus', '').strip()
                    if virus_name_noVirus and virus_name_noVirus not in lookup:
                        lookup[virus_name_noVirus] = record

        logger.info(f"✓ Built taxonomy lookup for {len(lookup)} entries")
        return lookup

    def get_taxonomy_statistics(self, records: List[Dict[str, any]]) -> Dict[str, any]:
        """
        Calculate taxonomy statistics.

        Parameters
        ----------
        records : List[Dict[str, any]]
            Parsed taxonomy records

        Returns
        -------
        Dict[str, any]
            Statistics
        """
        stats = {
            'total_records': len(records),
            'unique_species': len(set(r['species'] for r in records if r['species'])),
            'unique_genera': len(set(r['genus'] for r in records if r['genus'])),
            'unique_families': len(set(r['family'] for r in records if r['family'])),
            'unique_orders': len(set(r['order'] for r in records if r['order'])),
            'unique_realms': len(set(r['realm'] for r in records if r['realm']))
        }

        # Count records per realm
        realm_counts = defaultdict(int)
        for record in records:
            if record['realm']:
                realm_counts[record['realm']] += 1

        stats['realm_distribution'] = dict(realm_counts)

        # Count records per family (top 10)
        family_counts = defaultdict(int)
        for record in records:
            if record['family']:
                family_counts[record['family']] += 1

        top_families = sorted(family_counts.items(), key=lambda x: -x[1])[:10]
        stats['top_10_families'] = dict(top_families)

        return stats

    def map_to_database_genomes(
        self,
        records: List[Dict[str, any]],
        db_path: str
    ) -> Dict[str, int]:
        """
        Map ICTV taxonomy to genomes in database.

        Parameters
        ----------
        records : List[Dict[str, any]]
            ICTV taxonomy records
        db_path : str
            Path to genome database

        Returns
        -------
        Dict[str, int]
            Mapping statistics
        """
        logger.info(f"Mapping taxonomy to database: {db_path}")

        # Build lookup by species name
        taxonomy_lookup = self.build_taxonomy_lookup(records)

        # Connect to database
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Get all genomes with their current taxonomy
        cursor.execute("SELECT genome_id, species, ncbi_taxid FROM taxonomy")
        db_genomes = cursor.fetchall()

        logger.info(f"Found {len(db_genomes)} genomes in database")

        stats = {
            'total_genomes': len(db_genomes),
            'matched': 0,
            'unmatched': 0,
            'updated': 0,
            'failed': 0
        }

        matched_species = []
        unmatched_genomes = []

        for genome_id, species_name, ncbi_taxid in db_genomes:
            # Try multiple matching strategies
            matched = False
            matched_record = None

            species_lower = species_name.lower().strip()

            # Strategy 1: Direct match (case-insensitive)
            if species_lower in taxonomy_lookup:
                matched_record = taxonomy_lookup[species_lower]
                matched = True

            # Strategy 2: Without "virus" suffix
            if not matched:
                species_noVirus = species_lower.replace(' virus', '').strip()
                if species_noVirus in taxonomy_lookup:
                    matched_record = taxonomy_lookup[species_noVirus]
                    matched = True

            # Strategy 3: Try capitalized version (for binomial names)
            if not matched:
                if species_name in taxonomy_lookup:
                    matched_record = taxonomy_lookup[species_name]
                    matched = True

            # Record result
            if matched:
                matched_species.append((genome_id, matched_record['species']))
                stats['matched'] += 1
            else:
                unmatched_genomes.append((genome_id, species_name))
                stats['unmatched'] += 1

        conn.close()

        logger.info("=" * 60)
        logger.info("Mapping Statistics:")
        logger.info(f"  Total genomes: {stats['total_genomes']}")
        logger.info(f"  Matched: {stats['matched']} ({stats['matched']/stats['total_genomes']*100:.1f}%)")
        logger.info(f"  Unmatched: {stats['unmatched']} ({stats['unmatched']/stats['total_genomes']*100:.1f}%)")
        logger.info("=" * 60)

        # Save mapping results
        mapping_file = self.output_dir / "taxonomy_mapping.tsv"
        with open(mapping_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['genome_id', 'matched_species'])
            writer.writerows(matched_species)
        logger.info(f"✓ Saved mapping: {mapping_file}")

        # Save unmatched genomes
        if unmatched_genomes:
            unmatched_file = self.output_dir / "unmatched_genomes.tsv"
            with open(unmatched_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['genome_id', 'species_name'])
                writer.writerows(unmatched_genomes)
            logger.info(f"✓ Saved unmatched: {unmatched_file}")

        return stats

    def update_database_taxonomy(
        self,
        records: List[Dict[str, any]],
        db_path: str
    ) -> Dict[str, int]:
        """
        Update database taxonomy table with ICTV data.

        Parameters
        ----------
        records : List[Dict[str, any]]
            ICTV taxonomy records
        db_path : str
            Path to genome database

        Returns
        -------
        Dict[str, int]
            Update statistics
        """
        logger.info(f"Updating database taxonomy: {db_path}")

        # Build lookup
        taxonomy_lookup = self.build_taxonomy_lookup(records)

        # Connect to database
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Get all genomes
        cursor.execute("SELECT genome_id, species FROM taxonomy")
        db_genomes = cursor.fetchall()

        stats = {
            'total': len(db_genomes),
            'updated': 0,
            'failed': 0,
            'unmatched': 0
        }

        for genome_id, species_name in db_genomes:
            # Try multiple matching strategies
            ictv_record = None
            species_lower = species_name.lower().strip()

            # Strategy 1: Direct match (case-insensitive)
            if species_lower in taxonomy_lookup:
                ictv_record = taxonomy_lookup[species_lower]

            # Strategy 2: Without "virus" suffix
            if not ictv_record:
                species_noVirus = species_lower.replace(' virus', '').strip()
                if species_noVirus in taxonomy_lookup:
                    ictv_record = taxonomy_lookup[species_noVirus]

            # Strategy 3: Try capitalized version (for binomial names)
            if not ictv_record:
                if species_name in taxonomy_lookup:
                    ictv_record = taxonomy_lookup[species_name]

            if ictv_record:
                # Update taxonomy
                try:
                    cursor.execute("""
                        UPDATE taxonomy SET
                            realm = ?,
                            kingdom = ?,
                            phylum = ?,
                            class = ?,
                            order_name = ?,
                            family = ?,
                            subfamily = ?,
                            genus = ?,
                            species = ?
                        WHERE genome_id = ?
                    """, (
                        ictv_record['realm'],
                        ictv_record['kingdom'],
                        ictv_record['phylum'],
                        ictv_record['class'],
                        ictv_record['order'],
                        ictv_record['family'],
                        ictv_record['subfamily'],
                        ictv_record['genus'],
                        ictv_record['species'],
                        genome_id
                    ))
                    stats['updated'] += 1
                except Exception as e:
                    logger.error(f"Failed to update {genome_id}: {e}")
                    stats['failed'] += 1
            else:
                stats['unmatched'] += 1

            # Progress
            if (stats['updated'] + stats['failed'] + stats['unmatched']) % 10 == 0:
                total_processed = stats['updated'] + stats['failed'] + stats['unmatched']
                logger.info(f"Progress: {total_processed}/{stats['total']} | "
                           f"Updated: {stats['updated']} | Unmatched: {stats['unmatched']}")

        # Commit changes
        conn.commit()
        conn.close()

        logger.info("=" * 60)
        logger.info("Database Update Complete!")
        logger.info(f"  Total genomes: {stats['total']}")
        logger.info(f"  Updated: {stats['updated']} ({stats['updated']/stats['total']*100:.1f}%)")
        logger.info(f"  Unmatched: {stats['unmatched']} ({stats['unmatched']/stats['total']*100:.1f}%)")
        logger.info(f"  Failed: {stats['failed']}")
        logger.info("=" * 60)

        return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Parse ICTV VMR and map taxonomy to genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse VMR and save taxonomy data
  python scripts/parse_ictv_taxonomy.py --vmr data/ictv/VMR_current.xlsx --output data/ictv

  # Parse and map to database (without updating)
  python scripts/parse_ictv_taxonomy.py --vmr data/ictv/VMR_current.xlsx --output data/ictv \\
      --database viroforge/data/viral_genomes.db --map-only

  # Parse and update database taxonomy
  python scripts/parse_ictv_taxonomy.py --vmr data/ictv/VMR_current.xlsx --output data/ictv \\
      --database viroforge/data/viral_genomes.db --update-db
        """
    )

    parser.add_argument(
        '--vmr',
        type=str,
        required=True,
        help='Path to ICTV VMR Excel file'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='data/ictv',
        help='Output directory for parsed data (default: data/ictv)'
    )

    parser.add_argument(
        '--database',
        type=str,
        default=None,
        help='Path to genome database for mapping/updating'
    )

    parser.add_argument(
        '--map-only',
        action='store_true',
        help='Only map taxonomy to database (do not update)'
    )

    parser.add_argument(
        '--update-db',
        action='store_true',
        help='Update database taxonomy table with ICTV data'
    )

    args = parser.parse_args()

    # Initialize parser
    ictv_parser = ICTVTaxonomyParser(args.vmr, args.output)

    # Parse VMR
    logger.info("Step 1: Parsing ICTV VMR...")
    records = ictv_parser.parse_vmr()

    # Save taxonomy data
    logger.info("Step 2: Saving taxonomy data...")
    ictv_parser.save_taxonomy_data(records)

    # Get statistics
    logger.info("Step 3: Calculating statistics...")
    stats = ictv_parser.get_taxonomy_statistics(records)

    logger.info("\n" + "=" * 60)
    logger.info("ICTV VMR STATISTICS")
    logger.info("=" * 60)
    logger.info(f"Total records: {stats['total_records']}")
    logger.info(f"Unique species: {stats['unique_species']}")
    logger.info(f"Unique genera: {stats['unique_genera']}")
    logger.info(f"Unique families: {stats['unique_families']}")
    logger.info(f"Unique orders: {stats['unique_orders']}")
    logger.info(f"Unique realms: {stats['unique_realms']}")
    logger.info("\nRealm distribution:")
    for realm, count in sorted(stats['realm_distribution'].items(), key=lambda x: -x[1]):
        logger.info(f"  {realm}: {count}")
    logger.info("\nTop 10 families:")
    for family, count in stats['top_10_families'].items():
        logger.info(f"  {family}: {count}")
    logger.info("=" * 60)

    # Map to database if requested
    if args.database:
        logger.info("\nStep 4: Mapping taxonomy to database...")
        mapping_stats = ictv_parser.map_to_database_genomes(records, args.database)

        if args.update_db and not args.map_only:
            logger.info("\nStep 5: Updating database taxonomy...")
            update_stats = ictv_parser.update_database_taxonomy(records, args.database)

    logger.info("\n✓ ICTV taxonomy parsing complete!")


if __name__ == "__main__":
    main()
