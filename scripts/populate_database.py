#!/usr/bin/env python3
"""
Populate ViroForge genome database with RefSeq viral genomes.

This script reads parsed genome data, applies quality filters, and inserts
genomes into the SQLite database.

Usage:
    python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db
    python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db --limit 100

Author: ViroForge Development Team
Date: November 1, 2025
"""

import argparse
import json
import sqlite3
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class QualityFilter:
    """Quality filter for viral genomes."""

    def __init__(
        self,
        min_length: int = 1000,
        max_length: int = 500000,
        min_gc: float = 0.15,
        max_gc: float = 0.75,
        allowed_genome_types: Optional[List[str]] = None
    ):
        """
        Initialize quality filter.

        Parameters
        ----------
        min_length : int
            Minimum genome length (bp)
        max_length : int
            Maximum genome length (bp)
        min_gc : float
            Minimum GC content (0.0 to 1.0)
        max_gc : float
            Maximum GC content (0.0 to 1.0)
        allowed_genome_types : Optional[List[str]]
            Allowed genome types (None = all types allowed)
        """
        self.min_length = min_length
        self.max_length = max_length
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.allowed_genome_types = allowed_genome_types

    def passes_filter(self, genome: Dict) -> Tuple[bool, Optional[str]]:
        """
        Check if genome passes quality filters.

        Parameters
        ----------
        genome : Dict
            Genome data dictionary

        Returns
        -------
        Tuple[bool, Optional[str]]
            (passes filter?, reason for failure)
        """
        # Length filter
        length = genome.get('length', 0)
        if length < self.min_length:
            return False, f"Length too short: {length} < {self.min_length}"
        if length > self.max_length:
            return False, f"Length too long: {length} > {self.max_length}"

        # GC content filter
        gc = genome.get('gc_content', 0.5)
        if gc < self.min_gc:
            return False, f"GC content too low: {gc:.3f} < {self.min_gc}"
        if gc > self.max_gc:
            return False, f"GC content too high: {gc:.3f} > {self.max_gc}"

        # Genome type filter
        if self.allowed_genome_types:
            genome_type = genome.get('genome_type', 'unknown')
            if genome_type not in self.allowed_genome_types:
                return False, f"Genome type not allowed: {genome_type}"

        # Sequence quality checks
        sequence = genome.get('sequence', '')
        if not sequence:
            return False, "Empty sequence"

        # Check for excessive N's (ambiguous bases)
        n_count = sequence.upper().count('N')
        n_fraction = n_count / len(sequence) if sequence else 1.0
        if n_fraction > 0.05:  # More than 5% N's
            return False, f"Too many ambiguous bases: {n_fraction:.1%}"

        return True, None


class DatabasePopulator:
    """Populate ViroForge genome database."""

    def __init__(
        self,
        input_dir: str,
        db_path: str,
        quality_filter: Optional[QualityFilter] = None
    ):
        """
        Initialize database populator.

        Parameters
        ----------
        input_dir : str
            Directory containing parsed genome data
        db_path : str
            Path to SQLite database
        quality_filter : Optional[QualityFilter]
            Quality filter to apply (None = use defaults)
        """
        self.input_dir = Path(input_dir)
        self.db_path = Path(db_path)

        self.metadata_file = self.input_dir / "parsed_genomes.json"
        self.genomes_dir = self.input_dir / "genomes"

        self.quality_filter = quality_filter or QualityFilter()

    def load_genome_metadata(self) -> List[Dict]:
        """
        Load parsed genome metadata.

        Returns
        -------
        List[Dict]
            List of genome metadata dictionaries
        """
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {self.metadata_file}")

        with open(self.metadata_file, 'r') as f:
            genomes = json.load(f)

        logger.info(f"Loaded metadata for {len(genomes)} genomes")
        return genomes

    def load_genome_sequence(self, genome_id: str) -> str:
        """
        Load genome sequence from FASTA file.

        Parameters
        ----------
        genome_id : str
            Genome ID

        Returns
        -------
        str
            Genome sequence
        """
        fasta_file = self.genomes_dir / f"{genome_id}.fasta"

        if not fasta_file.exists():
            raise FileNotFoundError(f"Sequence file not found: {fasta_file}")

        with open(fasta_file, 'r') as f:
            lines = f.readlines()

        # Skip header, join sequence lines
        sequence = ''.join(line.strip() for line in lines[1:])
        return sequence

    def apply_quality_filters(
        self,
        genomes: List[Dict]
    ) -> Tuple[List[Dict], List[Tuple[str, str]]]:
        """
        Apply quality filters to genomes.

        Parameters
        ----------
        genomes : List[Dict]
            List of genome metadata

        Returns
        -------
        Tuple[List[Dict], List[Tuple[str, str]]]
            (passed genomes, [(failed genome_id, reason)])
        """
        logger.info("Applying quality filters...")

        passed = []
        failed = []

        for genome in genomes:
            # Load sequence for quality checks
            try:
                sequence = self.load_genome_sequence(genome['genome_id'])
                genome['sequence'] = sequence
            except FileNotFoundError:
                failed.append((genome['genome_id'], "Sequence file not found"))
                continue

            # Check filters
            passes, reason = self.quality_filter.passes_filter(genome)

            if passes:
                passed.append(genome)
            else:
                failed.append((genome['genome_id'], reason))

        logger.info(f"✓ Quality filtering complete:")
        logger.info(f"  Passed: {len(passed)}")
        logger.info(f"  Failed: {len(failed)}")

        if failed:
            logger.info("  Failure reasons:")
            failure_counts = {}
            for _, reason in failed:
                # Extract main reason (before colon)
                main_reason = reason.split(':')[0]
                failure_counts[main_reason] = failure_counts.get(main_reason, 0) + 1

            for reason, count in sorted(failure_counts.items(), key=lambda x: -x[1]):
                logger.info(f"    {reason}: {count}")

        return passed, failed

    def insert_genome(self, conn: sqlite3.Connection, genome: Dict) -> bool:
        """
        Insert a single genome into database.

        Parameters
        ----------
        conn : sqlite3.Connection
            Database connection
        genome : Dict
            Genome data dictionary

        Returns
        -------
        bool
            True if successful, False otherwise
        """
        try:
            cursor = conn.cursor()

            # Insert into genomes table
            cursor.execute("""
                INSERT INTO genomes (
                    genome_id, genome_name, sequence, length, gc_content,
                    genome_type, genome_structure, n_segments, assembly_level,
                    quality_score, source_database, refseq_category,
                    genbank_accession, date_added, date_modified, version
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                genome['genome_id'],
                genome.get('genome_name', ''),
                genome['sequence'],
                genome['length'],
                genome['gc_content'],
                genome['genome_type'],
                genome.get('genome_structure', 'linear'),
                genome.get('n_segments', 1),
                genome.get('assembly_level', 'Complete Genome'),
                None,  # quality_score - to be calculated later
                genome.get('source_database', 'RefSeq'),
                genome.get('refseq_category', ''),
                genome.get('genbank_accession', ''),
                genome.get('date_added', datetime.now().isoformat()),
                datetime.now().isoformat(),
                genome.get('version', 1)
            ))

            # Insert into taxonomy table if we have organism info
            if genome.get('organism_name'):
                # For now, just insert basic taxonomy
                # We'll populate full ICTV taxonomy later
                cursor.execute("""
                    INSERT INTO taxonomy (
                        genome_id, family, species, ncbi_taxid
                    ) VALUES (?, ?, ?, ?)
                """, (
                    genome['genome_id'],
                    'Unknown',  # Will be updated with ICTV data
                    genome.get('organism_name', ''),
                    genome.get('taxid')
                ))

            return True

        except sqlite3.IntegrityError as e:
            logger.warning(f"Integrity error for {genome['genome_id']}: {e}")
            return False
        except Exception as e:
            logger.error(f"Error inserting {genome['genome_id']}: {e}")
            return False

    def populate_database(
        self,
        genomes: List[Dict],
        batch_size: int = 100
    ) -> Dict[str, int]:
        """
        Populate database with genomes.

        Parameters
        ----------
        genomes : List[Dict]
            List of genome data dictionaries
        batch_size : int
            Number of genomes to insert per transaction

        Returns
        -------
        Dict[str, int]
            Statistics: inserted, failed
        """
        logger.info(f"Populating database: {self.db_path}")
        logger.info(f"Inserting {len(genomes)} genomes...")

        conn = sqlite3.connect(self.db_path)

        stats = {
            'inserted': 0,
            'failed': 0
        }

        try:
            for i in range(0, len(genomes), batch_size):
                batch = genomes[i:i+batch_size]

                for genome in batch:
                    success = self.insert_genome(conn, genome)
                    if success:
                        stats['inserted'] += 1
                    else:
                        stats['failed'] += 1

                # Commit batch
                conn.commit()

                # Progress update
                progress = min(i + batch_size, len(genomes))
                logger.info(
                    f"Progress: {progress}/{len(genomes)} ({progress/len(genomes)*100:.1f}%) | "
                    f"Inserted: {stats['inserted']} | Failed: {stats['failed']}"
                )

            # Update database metadata
            cursor = conn.cursor()
            cursor.execute(
                "UPDATE database_metadata SET value = ? WHERE key = 'last_updated'",
                (datetime.now().isoformat(),)
            )
            cursor.execute(
                "UPDATE database_metadata SET value = ? WHERE key = 'total_genomes'",
                (str(stats['inserted']),)
            )
            conn.commit()

        except Exception as e:
            conn.rollback()
            logger.error(f"Error populating database: {e}")
            raise
        finally:
            conn.close()

        logger.info("=" * 60)
        logger.info("Database population complete!")
        logger.info(f"  Inserted: {stats['inserted']}")
        logger.info(f"  Failed: {stats['failed']}")
        logger.info("=" * 60)

        return stats

    def get_database_stats(self) -> Dict[str, any]:
        """
        Get database statistics after population.

        Returns
        -------
        Dict[str, any]
            Database statistics
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        stats = {}

        # Genome count
        cursor.execute("SELECT COUNT(*) FROM genomes")
        stats['total_genomes'] = cursor.fetchone()[0]

        # Genome type distribution
        cursor.execute("""
            SELECT genome_type, COUNT(*)
            FROM genomes
            GROUP BY genome_type
        """)
        stats['genome_types'] = dict(cursor.fetchall())

        # Length statistics
        cursor.execute("""
            SELECT
                AVG(length) as mean_length,
                MIN(length) as min_length,
                MAX(length) as max_length,
                AVG(gc_content) as mean_gc
            FROM genomes
        """)
        row = cursor.fetchone()
        stats['mean_length'] = row[0]
        stats['min_length'] = row[1]
        stats['max_length'] = row[2]
        stats['mean_gc'] = row[3]

        conn.close()
        return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Populate ViroForge genome database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Populate database with all parsed genomes
  python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db

  # Populate with custom quality filters
  python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db \\
      --min-length 5000 --max-length 300000 --min-gc 0.2 --max-gc 0.7

  # Populate with limit (testing)
  python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db --limit 100
        """
    )

    parser.add_argument(
        '--input',
        type=str,
        default='data/parsed',
        help='Input directory containing parsed genomes (default: data/parsed)'
    )

    parser.add_argument(
        '--database',
        type=str,
        default='viroforge/data/viral_genomes.db',
        help='Path to SQLite database (default: viroforge/data/viral_genomes.db)'
    )

    parser.add_argument(
        '--min-length',
        type=int,
        default=1000,
        help='Minimum genome length (bp) (default: 1000)'
    )

    parser.add_argument(
        '--max-length',
        type=int,
        default=500000,
        help='Maximum genome length (bp) (default: 500000)'
    )

    parser.add_argument(
        '--min-gc',
        type=float,
        default=0.15,
        help='Minimum GC content (default: 0.15)'
    )

    parser.add_argument(
        '--max-gc',
        type=float,
        default=0.75,
        help='Maximum GC content (default: 0.75)'
    )

    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='Maximum number of genomes to insert (default: None = all)'
    )

    parser.add_argument(
        '--batch-size',
        type=int,
        default=100,
        help='Batch size for database inserts (default: 100)'
    )

    parser.add_argument(
        '--create-db',
        action='store_true',
        help='Create database if it does not exist'
    )

    args = parser.parse_args()

    # Create database if requested
    if args.create_db:
        if not Path(args.database).exists():
            logger.info(f"Creating database: {args.database}")
            from viroforge.data.database_schema import create_database
            create_database(args.database)

    # Check database exists
    if not Path(args.database).exists():
        logger.error(f"Database not found: {args.database}")
        logger.error("Use --create-db to create it automatically")
        return

    # Initialize quality filter
    quality_filter = QualityFilter(
        min_length=args.min_length,
        max_length=args.max_length,
        min_gc=args.min_gc,
        max_gc=args.max_gc
    )

    # Initialize populator
    populator = DatabasePopulator(
        input_dir=args.input,
        db_path=args.database,
        quality_filter=quality_filter
    )

    # Load genomes
    logger.info("Step 1: Loading genome metadata...")
    genomes = populator.load_genome_metadata()

    if args.limit:
        genomes = genomes[:args.limit]
        logger.info(f"Limited to {args.limit} genomes")

    # Apply quality filters
    logger.info("Step 2: Applying quality filters...")
    passed_genomes, failed_genomes = populator.apply_quality_filters(genomes)

    if not passed_genomes:
        logger.error("No genomes passed quality filters!")
        return

    # Populate database
    logger.info("Step 3: Populating database...")
    stats = populator.populate_database(passed_genomes, batch_size=args.batch_size)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("DATABASE SUMMARY")
    logger.info("=" * 60)

    db_stats = populator.get_database_stats()
    logger.info(f"Total genomes in database: {db_stats['total_genomes']}")
    logger.info(f"Mean genome length: {db_stats['mean_length']:.0f} bp")
    logger.info(f"Length range: {db_stats['min_length']:,} - {db_stats['max_length']:,} bp")
    logger.info(f"Mean GC content: {db_stats['mean_gc']:.3f}")
    logger.info("\nGenome type distribution:")
    for gtype, count in sorted(db_stats['genome_types'].items()):
        logger.info(f"  {gtype}: {count}")

    logger.info("=" * 60)
    logger.info(f"✓ Database populated: {args.database}")


if __name__ == "__main__":
    main()
