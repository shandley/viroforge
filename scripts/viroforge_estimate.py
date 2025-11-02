#!/usr/bin/env python3
"""
ViroForge FASTQ Size and Runtime Estimator

Estimate output file sizes, read counts, runtime, and resource requirements
before generating mock virome datasets.

Usage:
    # Basic estimation
    viroforge-estimate --genomes 500 --coverage 10x

    # Detailed estimation with all parameters
    viroforge-estimate \
      --genomes 500 \
      --coverage 10x \
      --read-length 150 \
      --fragment-length 350 \
      --compression gzip

    # From collection
    viroforge-estimate --collection gut --coverage 10x

    # Multiple coverage levels
    viroforge-estimate --genomes 500 --coverage 5x,10x,20x,50x

Author: ViroForge Development Team
Date: 2025-11-01
"""

import argparse
import sys
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional


class FastqEstimator:
    """Estimate FASTQ file sizes and resource requirements."""

    # Empirical compression ratios
    COMPRESSION_RATIOS = {
        'none': 1.0,
        'gzip': 0.25,    # ~4x compression
        'bzip2': 0.20,   # ~5x compression
        'xz': 0.18       # ~5.5x compression
    }

    # Bytes per read (overhead + sequence + quality)
    # @read_id\nSEQ\n+\nQUAL\n = ~4 + read_length + 4 + read_length + 4
    BYTES_PER_READ_OVERHEAD = 12  # @, +, 3x newlines, avg read_id length

    # Runtime estimates (seconds per million reads)
    # Based on benchmarking
    TIME_PER_MILLION_READS = {
        'generation': 60,      # Read generation
        'vlp_enrichment': 1,   # VLP enrichment (fast)
        'amplification': 2,    # Amplification bias
        'artifacts': 5,        # Platform artifacts
        'writing': 30          # FASTQ writing + compression
    }

    # Memory usage (MB per genome)
    MEMORY_PER_GENOME = 2  # Average 2 MB per genome in memory

    # Disk usage multiplier for temporary files
    TEMP_FILE_MULTIPLIER = 1.5

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = Path(db_path)

    def get_collection_info(self, collection_id: str) -> Optional[Dict]:
        """Get genome count and mean length from collection."""
        if not self.db_path.exists():
            return None

        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Check schema
        cursor = conn.execute("PRAGMA table_info(body_site_collections)")
        columns = [row[1] for row in cursor.fetchall()]

        # Find collection
        if 'collection_name' in columns:
            cursor = conn.execute("""
                SELECT collection_id, n_genomes as genome_count
                FROM body_site_collections
                WHERE collection_id LIKE ? OR collection_name LIKE ?
            """, (f"%{collection_id}%", f"%{collection_id}%"))
        else:
            cursor = conn.execute("""
                SELECT collection_id, genome_count
                FROM body_site_collections
                WHERE collection_id LIKE ? OR name LIKE ?
            """, (f"%{collection_id}%", f"%{collection_id}%"))

        result = cursor.fetchone()
        if not result:
            conn.close()
            return None

        result = dict(result)
        coll_id = result['collection_id']

        # Get mean genome length
        cursor = conn.execute("""
            SELECT AVG(g.length) as mean_length
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            WHERE cg.collection_id = ?
        """, (coll_id,))

        row = cursor.fetchone()
        mean_length = int(row['mean_length']) if row and row['mean_length'] else 50000

        conn.close()

        return {
            'genome_count': result['genome_count'],
            'mean_length': mean_length
        }

    def estimate(
        self,
        n_genomes: int,
        coverage: float,
        read_length: int = 150,
        fragment_length: int = 350,
        compression: str = 'gzip',
        mean_genome_length: int = 50000
    ) -> Dict:
        """
        Estimate FASTQ file sizes and resource requirements.

        Args:
            n_genomes: Number of genomes in community
            coverage: Sequencing coverage (e.g., 10.0 for 10x)
            read_length: Read length in bp (default: 150)
            fragment_length: Mean fragment length (default: 350)
            compression: Compression type (none, gzip, bzip2, xz)
            mean_genome_length: Mean genome length (default: 50000)

        Returns:
            Dictionary with estimates
        """
        # Calculate total community size
        total_bases = n_genomes * mean_genome_length

        # Calculate required reads for coverage
        # Coverage = (reads × read_length) / genome_length
        # For paired-end: each read pair contributes 2 × read_length
        bases_needed = total_bases * coverage
        read_pairs_needed = int(bases_needed / (2 * read_length))
        total_reads = read_pairs_needed * 2

        # Calculate uncompressed FASTQ size
        bytes_per_read = self.BYTES_PER_READ_OVERHEAD + (2 * read_length)
        uncompressed_bytes = total_reads * bytes_per_read

        # Apply compression
        compression_ratio = self.COMPRESSION_RATIOS.get(compression, 0.25)
        compressed_bytes = int(uncompressed_bytes * compression_ratio)

        # Calculate runtime
        million_reads = total_reads / 1_000_000
        runtime_seconds = sum(
            million_reads * time
            for time in self.TIME_PER_MILLION_READS.values()
        )

        # Calculate memory
        memory_mb = n_genomes * self.MEMORY_PER_GENOME
        memory_mb += 500  # Base overhead
        if total_reads > 10_000_000:
            memory_mb += 500  # Extra for large datasets

        # Calculate disk space
        disk_space_gb = (compressed_bytes * 2) / (1024**3)  # R1 + R2
        temp_space_gb = disk_space_gb * self.TEMP_FILE_MULTIPLIER
        total_disk_gb = disk_space_gb + temp_space_gb

        return {
            'input': {
                'n_genomes': n_genomes,
                'mean_genome_length': mean_genome_length,
                'total_community_size': total_bases,
                'coverage': coverage,
                'read_length': read_length,
                'fragment_length': fragment_length,
                'compression': compression
            },
            'output': {
                'read_pairs': read_pairs_needed,
                'total_reads': total_reads,
                'uncompressed_size_gb': uncompressed_bytes / (1024**3),
                'compressed_size_gb': compressed_bytes / (1024**3),
                'fastq_r1_size_gb': (compressed_bytes / 2) / (1024**3),
                'fastq_r2_size_gb': (compressed_bytes / 2) / (1024**3),
                'total_fastq_size_gb': disk_space_gb,
                'compression_ratio': compression_ratio
            },
            'resources': {
                'estimated_runtime_min': runtime_seconds / 60,
                'estimated_runtime_hr': runtime_seconds / 3600,
                'memory_required_mb': memory_mb,
                'memory_required_gb': memory_mb / 1024,
                'disk_required_gb': total_disk_gb,
                'disk_final_gb': disk_space_gb,
                'disk_temp_gb': temp_space_gb
            }
        }

    def print_estimate(self, estimate: Dict):
        """Print formatted estimate."""
        inp = estimate['input']
        out = estimate['output']
        res = estimate['resources']

        print("=" * 70)
        print("ViroForge FASTQ Size & Runtime Estimation")
        print("=" * 70)
        print()

        # Input parameters
        print("Input Parameters:")
        print("-" * 70)
        print(f"  Genomes in community:      {inp['n_genomes']:>10,}")
        print(f"  Mean genome length:        {inp['mean_genome_length']:>10,} bp")
        print(f"  Total community size:      {inp['total_community_size']:>10,} bp ({inp['total_community_size']/1e6:.1f} Mbp)")
        print(f"  Target coverage:           {inp['coverage']:>10.1f}x")
        print(f"  Read length:               {inp['read_length']:>10} bp (paired-end)")
        print(f"  Fragment length:           {inp['fragment_length']:>10} bp")
        print(f"  Compression:               {inp['compression']:>10}")
        print()

        # Output files
        print("Output Files:")
        print("-" * 70)
        print(f"  Read pairs:                {out['read_pairs']:>10,}")
        print(f"  Total reads:               {out['total_reads']:>10,}")
        print()
        print(f"  Uncompressed size:         {out['uncompressed_size_gb']:>10.2f} GB")
        print(f"  Compressed size (total):   {out['total_fastq_size_gb']:>10.2f} GB")
        print(f"    - R1 FASTQ:              {out['fastq_r1_size_gb']:>10.2f} GB")
        print(f"    - R2 FASTQ:              {out['fastq_r2_size_gb']:>10.2f} GB")
        print(f"  Compression ratio:         {out['compression_ratio']:>10.1%}")
        print()

        # Resource requirements
        print("Resource Requirements:")
        print("-" * 70)

        # Runtime
        if res['estimated_runtime_hr'] >= 1.0:
            print(f"  Estimated runtime:         {res['estimated_runtime_hr']:>10.1f} hours")
        else:
            print(f"  Estimated runtime:         {res['estimated_runtime_min']:>10.0f} minutes")

        # Memory
        if res['memory_required_gb'] >= 1.0:
            print(f"  Memory required:           {res['memory_required_gb']:>10.1f} GB")
        else:
            print(f"  Memory required:           {res['memory_required_mb']:>10.0f} MB")

        # Disk
        print(f"  Disk space required:       {res['disk_required_gb']:>10.1f} GB")
        print(f"    - Final output:          {res['disk_final_gb']:>10.1f} GB")
        print(f"    - Temp files:            {res['disk_temp_gb']:>10.1f} GB")
        print()

        # Recommendations
        print("Recommendations:")
        print("-" * 70)

        # Memory warning
        if res['memory_required_gb'] > 8:
            print("  ⚠️  High memory usage - ensure adequate RAM available")
        else:
            print("  ✓ Memory usage within normal range")

        # Disk warning
        if res['disk_required_gb'] > 100:
            print("  ⚠️  Large disk space required - verify free space")
        elif res['disk_required_gb'] > 50:
            print("  ℹ️  Moderate disk space required")
        else:
            print("  ✓ Disk space within normal range")

        # Runtime warning
        if res['estimated_runtime_hr'] > 2:
            print("  ⚠️  Long runtime - consider running in background")
        elif res['estimated_runtime_min'] > 30:
            print("  ℹ️  Moderate runtime expected")
        else:
            print("  ✓ Runtime within normal range")

        # Compression recommendation
        if inp['compression'] == 'none' and out['uncompressed_size_gb'] > 10:
            print("  💡 Consider using gzip compression to save disk space (4x smaller)")

        print()


def parse_coverage(coverage_str: str) -> List[float]:
    """Parse coverage string (e.g., '10x' or '5x,10x,20x')."""
    coverages = []
    for part in coverage_str.split(','):
        part = part.strip().lower().rstrip('x')
        try:
            coverages.append(float(part))
        except ValueError:
            print(f"⚠️  Invalid coverage: {part}", file=sys.stderr)
    return coverages


def format_size(bytes_val: float) -> str:
    """Format bytes as human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.1f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.1f} PB"


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge FASTQ Size & Runtime Estimator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic estimation
  viroforge-estimate --genomes 500 --coverage 10x

  # Detailed estimation
  viroforge-estimate \\
    --genomes 500 \\
    --coverage 10x \\
    --read-length 150 \\
    --fragment-length 350 \\
    --compression gzip

  # From collection
  viroforge-estimate --collection gut --coverage 10x

  # Multiple coverage levels
  viroforge-estimate --genomes 500 --coverage 5x,10x,20x,50x

  # Different compression
  viroforge-estimate --genomes 500 --coverage 10x --compression none
  viroforge-estimate --genomes 500 --coverage 10x --compression xz

  # Quick comparison table
  viroforge-estimate --genomes 500 --coverage 5x,10x,20x --format table
        """
    )

    # Input source (mutually exclusive)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        '--genomes',
        type=int,
        help='Number of genomes in community'
    )
    source_group.add_argument(
        '--collection',
        help='Collection ID (e.g., gut, oral)'
    )

    # Sequencing parameters
    parser.add_argument(
        '--coverage',
        required=True,
        help='Target coverage (e.g., 10x or 5x,10x,20x for multiple)'
    )
    parser.add_argument(
        '--read-length',
        type=int,
        default=150,
        help='Read length in bp (default: 150)'
    )
    parser.add_argument(
        '--fragment-length',
        type=int,
        default=350,
        help='Mean fragment length in bp (default: 350)'
    )
    parser.add_argument(
        '--compression',
        choices=['none', 'gzip', 'bzip2', 'xz'],
        default='gzip',
        help='Compression type (default: gzip)'
    )
    parser.add_argument(
        '--mean-length',
        type=int,
        default=50000,
        help='Mean genome length in bp (default: 50000, ignored if --collection)'
    )

    # Output format
    parser.add_argument(
        '--format',
        choices=['detailed', 'table', 'json'],
        default='detailed',
        help='Output format (default: detailed)'
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    args = parser.parse_args()

    # Create estimator
    estimator = FastqEstimator(args.database)

    # Get genome count and mean length
    if args.collection:
        info = estimator.get_collection_info(args.collection)
        if not info:
            print(f"❌ Error: Collection '{args.collection}' not found", file=sys.stderr)
            sys.exit(1)
        n_genomes = info['genome_count']
        mean_length = info['mean_length']
        print(f"Using collection '{args.collection}': {n_genomes} genomes, mean length {mean_length:,} bp")
        print()
    else:
        n_genomes = args.genomes
        mean_length = args.mean_length

    # Parse coverage levels
    coverages = parse_coverage(args.coverage)

    if not coverages:
        print("❌ Error: Invalid coverage specification", file=sys.stderr)
        sys.exit(1)

    # Generate estimates
    estimates = []
    for coverage in coverages:
        estimate = estimator.estimate(
            n_genomes=n_genomes,
            coverage=coverage,
            read_length=args.read_length,
            fragment_length=args.fragment_length,
            compression=args.compression,
            mean_genome_length=mean_length
        )
        estimates.append(estimate)

    # Output
    if args.format == 'detailed':
        for estimate in estimates:
            estimator.print_estimate(estimate)
            if len(estimates) > 1:
                print("\n" + "=" * 70 + "\n")

    elif args.format == 'table':
        print("Coverage Comparison Table:")
        print("=" * 100)
        print(f"{'Coverage':<10} {'Read Pairs':>12} {'Total Reads':>12} {'FASTQ Size':>12} {'Runtime':>12} {'Memory':>10}")
        print("-" * 100)

        for estimate in estimates:
            inp = estimate['input']
            out = estimate['output']
            res = estimate['resources']

            runtime = f"{res['estimated_runtime_hr']:.1f}h" if res['estimated_runtime_hr'] >= 1 else f"{res['estimated_runtime_min']:.0f}m"
            memory = f"{res['memory_required_gb']:.1f}GB" if res['memory_required_gb'] >= 1 else f"{res['memory_required_mb']:.0f}MB"

            print(f"{inp['coverage']:>8.1f}x {out['read_pairs']:>12,} {out['total_reads']:>12,} "
                  f"{out['total_fastq_size_gb']:>11.1f}GB {runtime:>12} {memory:>10}")

        print()

    elif args.format == 'json':
        import json
        print(json.dumps(estimates, indent=2))


if __name__ == '__main__':
    main()
