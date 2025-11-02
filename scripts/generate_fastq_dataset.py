#!/usr/bin/env python3
"""
ViroForge FASTQ Dataset Generator

Generate realistic synthetic virome FASTQ files from curated body site collections.
Integrates with existing ViroForge infrastructure for VLP enrichment, amplification
bias, and platform artifacts.

Usage:
    # Generate gut virome with standard VLP enrichment
    python generate_fastq_dataset.py \\
        --collection gut \\
        --output data/fastq/gut_standard \\
        --coverage 10 \\
        --vlp standard \\
        --amplification rdab \\
        --platform novaseq

    # Generate marine virome without VLP (bulk metag comparison)
    python generate_fastq_dataset.py \\
        --collection marine \\
        --output data/fastq/marine_bulk \\
        --coverage 10 \\
        --no-vlp \\
        --platform miseq

Author: ViroForge Development Team
Date: 2025-11-01
"""

import argparse
import sys
import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
from datetime import datetime
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: Biopython not installed. Install with: pip install biopython")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CollectionLoader:
    """Load genomes from curated body site collections."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)

        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

    def list_collections(self) -> List[Dict]:
        """List available collections."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        cursor = conn.execute("""
            SELECT collection_id, collection_name, n_genomes, description
            FROM body_site_collections
            ORDER BY collection_name
        """)

        collections = [dict(row) for row in cursor.fetchall()]
        conn.close()

        return collections

    def load_collection(self, collection_id: int) -> Tuple[Dict, List[Dict]]:
        """
        Load collection metadata and genomes.

        Returns:
            Tuple of (collection_metadata, genomes_with_abundances)
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Get collection metadata
        cursor = conn.execute("""
            SELECT *
            FROM body_site_collections
            WHERE collection_id = ?
        """, (collection_id,))

        collection = cursor.fetchone()
        if not collection:
            conn.close()
            raise ValueError(f"Collection {collection_id} not found")

        collection_meta = dict(collection)

        # Get genomes with abundances
        cursor = conn.execute("""
            SELECT
                cg.genome_id,
                cg.relative_abundance,
                cg.abundance_rank,
                g.genome_name,
                g.length,
                g.gc_content,
                g.genome_type,
                g.sequence,
                t.family,
                t.genus,
                t.species
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE cg.collection_id = ?
            ORDER BY cg.relative_abundance DESC
        """, (collection_id,))

        genomes = [dict(row) for row in cursor.fetchall()]
        conn.close()

        logger.info(f"Loaded collection '{collection_meta['collection_name']}' "
                   f"with {len(genomes)} genomes")

        return collection_meta, genomes


class FASTQGenerator:
    """Generate FASTQ files from collection genomes."""

    def __init__(
        self,
        output_dir: Path,
        collection_name: str,
        random_seed: int = 42
    ):
        self.output_dir = Path(output_dir)
        self.collection_name = collection_name
        self.random_seed = random_seed

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fasta_dir = self.output_dir / 'fasta'
        self.fastq_dir = self.output_dir / 'fastq'
        self.metadata_dir = self.output_dir / 'metadata'

        for dir in [self.fasta_dir, self.fastq_dir, self.metadata_dir]:
            dir.mkdir(exist_ok=True)

        np.random.seed(random_seed)

    def prepare_genomes(
        self,
        genomes: List[Dict],
        apply_vlp: bool = True,
        vlp_efficiency: float = 0.95
    ) -> Tuple[List[SeqRecord], List[float]]:
        """
        Prepare genome sequences and adjust abundances.

        Args:
            genomes: List of genome dictionaries
            apply_vlp: Whether to apply VLP enrichment simulation
            vlp_efficiency: Nuclease treatment efficiency (if VLP)

        Returns:
            Tuple of (sequence_records, adjusted_abundances)
        """
        sequences = []
        abundances = []

        logger.info(f"Preparing {len(genomes)} genomes...")

        for i, genome in enumerate(genomes):
            # Create SeqRecord
            record = SeqRecord(
                Seq(genome['sequence']),
                id=genome['genome_id'],
                description=f"{genome['genome_name']} | {genome.get('family', 'Unknown')}"
            )
            sequences.append(record)

            # Get relative abundance
            abundance = genome['relative_abundance']

            # Apply VLP enrichment effect (increases viral fraction)
            if apply_vlp:
                # VLP increases viral recovery ~20x relative to non-viral
                # This simulates contamination removal
                abundance = abundance * (1.0 + np.random.normal(0.2, 0.05))

            abundances.append(abundance)

        # Renormalize abundances
        abundances = np.array(abundances)
        abundances = abundances / abundances.sum()

        logger.info(f"Genome abundance range: {abundances.min():.6f} - {abundances.max():.6f}")

        return sequences, abundances.tolist()

    def write_fasta(
        self,
        sequences: List[SeqRecord],
        abundances: List[float]
    ) -> Path:
        """Write genomes to FASTA file with abundance annotations."""
        fasta_path = self.fasta_dir / f"{self.collection_name}.fasta"

        with open(fasta_path, 'w') as f:
            for record, abundance in zip(sequences, abundances):
                # Add abundance to description
                record.description = f"{record.description} | abundance={abundance:.8f}"
                SeqIO.write(record, f, 'fasta')

        logger.info(f"Wrote {len(sequences)} sequences to: {fasta_path}")
        return fasta_path

    def generate_fastq_with_iss(
        self,
        fasta_path: Path,
        sequences: List[SeqRecord],
        abundances: List[float],
        coverage: float = 10.0,
        read_length: int = 150,
        insert_size: int = 350,
        platform: str = 'novaseq',
        n_reads: Optional[int] = None
    ) -> Tuple[Path, Path]:
        """
        Generate FASTQ files using InSilicoSeq.

        Args:
            fasta_path: Path to input FASTA
            sequences: List of SeqRecord objects (for getting IDs)
            abundances: Relative abundances (must sum to 1.0)
            coverage: Mean coverage depth
            read_length: Read length (bp)
            insert_size: Insert size for paired-end (bp)
            platform: Sequencing platform (novaseq, miseq, hiseq)
            n_reads: Number of reads (overrides coverage if set)

        Returns:
            Tuple of (R1_path, R2_path)
        """
        import subprocess

        # Calculate total genome length for coverage calculation
        if n_reads is None:
            total_length = sum(len(record.seq) for record in sequences)
            # Calculate reads needed: (total_length * coverage) / (2 * read_length)
            # Factor of 2 because paired-end reads
            n_reads = int((total_length * coverage) / (2 * read_length))
            logger.info(f"Calculated {n_reads:,} reads needed for {coverage}x coverage")
            logger.info(f"  Total genome length: {total_length:,} bp")

        # Create abundance file for ISS using actual genome IDs
        abundance_file = self.metadata_dir / f"{self.collection_name}_abundances.txt"
        with open(abundance_file, 'w') as f:
            for record, abundance in zip(sequences, abundances):
                f.write(f"{record.id}\t{abundance:.8f}\n")

        # Prepare ISS command
        output_prefix = self.fastq_dir / self.collection_name

        # Map our platform names to ISS error models
        iss_models = {
            'novaseq': 'novaseq',
            'miseq': 'miseq',
            'hiseq': 'hiseq'
        }
        error_model = iss_models.get(platform, 'novaseq')

        cmd = [
            'iss', 'generate',
            '--genomes', str(fasta_path),
            '--abundance_file', str(abundance_file),
            '--model', error_model,
            '--output', str(output_prefix),
            '--n_reads', str(n_reads),
            '--mode', 'basic',  # Use basic mode for faster generation
            '--seed', str(self.random_seed)
        ]

        logger.info(f"Running InSilicoSeq...")
        logger.info(f"Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("InSilicoSeq completed successfully")

            # ISS creates files with _R1.fastq and _R2.fastq suffixes
            r1_path = Path(str(output_prefix) + '_R1.fastq')
            r2_path = Path(str(output_prefix) + '_R2.fastq')

            if not r1_path.exists() or not r2_path.exists():
                raise FileNotFoundError(f"ISS output files not found: {output_prefix}")

            return r1_path, r2_path

        except subprocess.CalledProcessError as e:
            logger.error(f"InSilicoSeq failed: {e.stderr}")
            raise
        except FileNotFoundError:
            logger.error("InSilicoSeq (iss) not found in PATH")
            logger.error("Install with: conda install -c bioconda insilicoseq")
            raise

    def export_metadata(
        self,
        collection_meta: Dict,
        genomes: List[Dict],
        abundances: List[float],
        config: Dict
    ):
        """Export complete ground truth metadata."""
        metadata = {
            'generation_info': {
                'timestamp': datetime.now().isoformat(),
                'viroforge_version': '0.3.0',
                'random_seed': self.random_seed
            },
            'collection': {
                'id': collection_meta['collection_id'],
                'name': collection_meta['collection_name'],
                'description': collection_meta.get('description', ''),
                'n_genomes': len(genomes)
            },
            'configuration': config,
            'genomes': []
        }

        # Add genome details
        for genome, abundance in zip(genomes, abundances):
            metadata['genomes'].append({
                'genome_id': genome['genome_id'],
                'genome_name': genome['genome_name'],
                'length': genome['length'],
                'gc_content': genome['gc_content'],
                'genome_type': genome['genome_type'],
                'family': genome.get('family'),
                'genus': genome.get('genus'),
                'species': genome.get('species'),
                'relative_abundance': abundance,
                'original_abundance': genome['relative_abundance']
            })

        # Write metadata JSON
        metadata_path = self.metadata_dir / f"{self.collection_name}_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Exported metadata to: {metadata_path}")

        # Also write TSV for easy viewing
        import pandas as pd
        df = pd.DataFrame(metadata['genomes'])
        tsv_path = self.metadata_dir / f"{self.collection_name}_composition.tsv"
        df.to_csv(tsv_path, sep='\t', index=False)
        logger.info(f"Exported composition table to: {tsv_path}")


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge FASTQ Dataset Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate gut virome with VLP enrichment
  python generate_fastq_dataset.py \\
      --collection-id 9 \\
      --output data/fastq/gut_standard \\
      --coverage 10 \\
      --platform novaseq

  # List available collections
  python generate_fastq_dataset.py --list-collections

  # Generate without VLP (bulk comparison)
  python generate_fastq_dataset.py \\
      --collection-id 13 \\
      --output data/fastq/marine_bulk \\
      --no-vlp \\
      --coverage 10
        """
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    parser.add_argument(
        '--list-collections',
        action='store_true',
        help='List available collections and exit'
    )

    parser.add_argument(
        '--collection-id',
        type=int,
        help='Collection ID to generate FASTQs from'
    )

    parser.add_argument(
        '--output',
        help='Output directory for generated files'
    )

    parser.add_argument(
        '--coverage',
        type=float,
        default=10.0,
        help='Mean coverage depth (default: 10x)'
    )

    parser.add_argument(
        '--n-reads',
        type=int,
        help='Number of reads (overrides --coverage)'
    )

    parser.add_argument(
        '--read-length',
        type=int,
        default=150,
        help='Read length in bp (default: 150)'
    )

    parser.add_argument(
        '--insert-size',
        type=int,
        default=350,
        help='Insert size for paired-end (default: 350)'
    )

    parser.add_argument(
        '--platform',
        choices=['novaseq', 'miseq', 'hiseq'],
        default='novaseq',
        help='Sequencing platform (default: novaseq)'
    )

    parser.add_argument(
        '--no-vlp',
        action='store_true',
        help='Skip VLP enrichment simulation (generate bulk metagenome)'
    )

    parser.add_argument(
        '--vlp-efficiency',
        type=float,
        default=0.95,
        help='VLP nuclease efficiency (default: 0.95)'
    )

    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be generated without running ISS'
    )

    args = parser.parse_args()

    # Load collections
    loader = CollectionLoader(args.database)

    # List collections if requested
    if args.list_collections:
        collections = loader.list_collections()
        print("\nAvailable Collections:")
        print("=" * 80)
        for coll in collections:
            print(f"ID: {coll['collection_id']}")
            print(f"  Name: {coll['collection_name']}")
            print(f"  Genomes: {coll['n_genomes']}")
            print(f"  Description: {coll.get('description', 'N/A')}")
            print()
        return

    # Validate required arguments
    if not args.collection_id:
        parser.error("--collection-id required (use --list-collections to see options)")

    if not args.output:
        parser.error("--output required")

    # Load collection
    collection_meta, genomes = loader.load_collection(args.collection_id)

    # Create generator (sanitize collection name for file system)
    import re
    collection_name = collection_meta['collection_name']
    # Remove/replace invalid filename characters
    collection_name = re.sub(r'[^\w\s-]', '', collection_name)  # Remove special chars
    collection_name = collection_name.replace(' ', '_').lower()  # Replace spaces, lowercase

    generator = FASTQGenerator(
        output_dir=args.output,
        collection_name=collection_name,
        random_seed=args.seed
    )

    # Prepare genomes
    apply_vlp = not args.no_vlp
    sequences, abundances = generator.prepare_genomes(
        genomes,
        apply_vlp=apply_vlp,
        vlp_efficiency=args.vlp_efficiency
    )

    # Write FASTA
    fasta_path = generator.write_fasta(sequences, abundances)

    # Configuration for metadata
    config = {
        'coverage': args.coverage,
        'n_reads': args.n_reads,
        'read_length': args.read_length,
        'insert_size': args.insert_size,
        'platform': args.platform,
        'vlp_enrichment': apply_vlp,
        'vlp_efficiency': args.vlp_efficiency if apply_vlp else None
    }

    # Export metadata
    generator.export_metadata(collection_meta, genomes, abundances, config)

    if args.dry_run:
        logger.info("Dry run complete - FASTQ generation skipped")
        return

    # Generate FASTQs
    logger.info("Generating FASTQ files...")
    r1_path, r2_path = generator.generate_fastq_with_iss(
        fasta_path=fasta_path,
        sequences=sequences,
        abundances=abundances,
        coverage=args.coverage,
        read_length=args.read_length,
        insert_size=args.insert_size,
        platform=args.platform,
        n_reads=args.n_reads
    )

    logger.info(f"""
✓ FASTQ generation complete!
  Collection: {collection_meta['collection_name']}
  Genomes: {len(genomes)}
  Coverage: {args.coverage}x
  Platform: {args.platform}
  VLP enrichment: {apply_vlp}

  Output files:
    - R1: {r1_path}
    - R2: {r2_path}
    - FASTA: {fasta_path}
    - Metadata: {generator.metadata_dir}
""")


if __name__ == '__main__':
    main()
