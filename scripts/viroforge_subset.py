#!/usr/bin/env python3
"""
ViroForge Genome Subset Selector

Create small, representative genome subsets for testing and development.
Useful for quick pipeline testing before running full datasets.

Usage:
    # Random subset
    viroforge-subset --count 50 --output test_subset.txt

    # Taxonomically diverse subset
    viroforge-subset --count 50 --strategy diverse --output diverse_subset.txt

    # Representative of each family
    viroforge-subset --strategy family_representatives --output family_reps.txt

    # Subset for specific body site
    viroforge-subset --count 20 --family Siphoviridae --host Streptococcus

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import sys
import random
from pathlib import Path
from typing import Dict, List, Optional
from collections import Counter, defaultdict


class SubsetSelector:
    """Select representative genome subsets."""

    def __init__(self, db_path: str, random_seed: int = 42):
        self.db_path = Path(db_path)
        self.random_seed = random_seed
        random.seed(random_seed)

        if not self.db_path.exists():
            print(f"❌ Error: Database not found at {db_path}", file=sys.stderr)
            sys.exit(1)

    def get_connection(self):
        """Get database connection."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def random_subset(self, count: int, **filters) -> List[Dict]:
        """Select random subset with optional filters."""
        conn = self.get_connection()

        # Build query
        query = """
            SELECT
                g.genome_id,
                g.genome_name,
                g.length,
                g.gc_content,
                g.genome_type,
                t.family,
                t.genus
            FROM genomes g
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE 1=1
        """

        params = []

        # Apply filters
        if filters.get('family'):
            query += " AND t.family LIKE ?"
            params.append(f"%{filters['family']}%")

        if filters.get('host'):
            query += """ AND EXISTS (
                SELECT 1 FROM host_associations h
                WHERE h.genome_id = g.genome_id
                AND h.host_name LIKE ?
            )"""
            params.append(f"%{filters['host']}%")

        if filters.get('genome_type'):
            query += " AND g.genome_type = ?"
            params.append(filters['genome_type'])

        if filters.get('length_min'):
            query += " AND g.length >= ?"
            params.append(filters['length_min'])

        if filters.get('length_max'):
            query += " AND g.length <= ?"
            params.append(filters['length_max'])

        query += " ORDER BY RANDOM() LIMIT ?"
        params.append(count)

        cursor = conn.execute(query, params)
        results = [dict(row) for row in cursor.fetchall()]

        conn.close()
        return results

    def diverse_subset(self, count: int) -> List[Dict]:
        """
        Select taxonomically diverse subset.

        Strategy: Sample evenly across taxonomic families to maximize diversity.
        """
        conn = self.get_connection()

        # Get all families
        cursor = conn.execute("""
            SELECT DISTINCT t.family
            FROM taxonomy t
            WHERE t.family IS NOT NULL
            ORDER BY t.family
        """)
        families = [row['family'] for row in cursor.fetchall()]

        if not families:
            conn.close()
            return self.random_subset(count)

        # Calculate genomes per family
        genomes_per_family = max(1, count // len(families))
        remainder = count % len(families)

        selected = []

        for i, family in enumerate(families):
            # Get count for this family (add 1 for first N families to use remainder)
            family_count = genomes_per_family + (1 if i < remainder else 0)

            cursor = conn.execute("""
                SELECT
                    g.genome_id,
                    g.genome_name,
                    g.length,
                    g.gc_content,
                    g.genome_type,
                    t.family,
                    t.genus
                FROM genomes g
                JOIN taxonomy t ON g.genome_id = t.genome_id
                WHERE t.family = ?
                ORDER BY RANDOM()
                LIMIT ?
            """, (family, family_count))

            selected.extend([dict(row) for row in cursor.fetchall()])

            if len(selected) >= count:
                break

        conn.close()

        # Shuffle final selection
        random.shuffle(selected)
        return selected[:count]

    def family_representatives(self, per_family: int = 1) -> List[Dict]:
        """
        Select representative genomes from each family.

        Strategy: Pick N genomes from each family (default: 1).
        """
        conn = self.get_connection()

        # Get all families with counts
        cursor = conn.execute("""
            SELECT t.family, COUNT(*) as count
            FROM taxonomy t
            WHERE t.family IS NOT NULL
            GROUP BY t.family
            ORDER BY count DESC
        """)
        families = cursor.fetchall()

        selected = []

        for row in families:
            family = row['family']

            cursor = conn.execute("""
                SELECT
                    g.genome_id,
                    g.genome_name,
                    g.length,
                    g.gc_content,
                    g.genome_type,
                    t.family,
                    t.genus
                FROM genomes g
                JOIN taxonomy t ON g.genome_id = t.genome_id
                WHERE t.family = ?
                ORDER BY RANDOM()
                LIMIT ?
            """, (family, per_family))

            selected.extend([dict(row) for row in cursor.fetchall()])

        conn.close()
        return selected

    def size_stratified_subset(self, count: int) -> List[Dict]:
        """
        Select genomes stratified by size.

        Strategy: Sample evenly across size bins (small, medium, large).
        """
        conn = self.get_connection()

        # Define size bins (based on typical viral genome sizes)
        bins = [
            ('small', 0, 10000),
            ('medium', 10000, 50000),
            ('large', 50000, 1000000)
        ]

        genomes_per_bin = count // len(bins)
        remainder = count % len(bins)

        selected = []

        for i, (bin_name, min_len, max_len) in enumerate(bins):
            bin_count = genomes_per_bin + (1 if i < remainder else 0)

            cursor = conn.execute("""
                SELECT
                    g.genome_id,
                    g.genome_name,
                    g.length,
                    g.gc_content,
                    g.genome_type,
                    t.family,
                    t.genus
                FROM genomes g
                LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
                WHERE g.length >= ? AND g.length < ?
                ORDER BY RANDOM()
                LIMIT ?
            """, (min_len, max_len, bin_count))

            selected.extend([dict(row) for row in cursor.fetchall()])

        conn.close()
        random.shuffle(selected)
        return selected

    def balanced_genometype_subset(self, count: int) -> List[Dict]:
        """
        Select genomes balanced by genome type (dsDNA, ssRNA, etc).

        Strategy: Sample proportionally from each genome type.
        """
        conn = self.get_connection()

        # Get genome type distribution
        cursor = conn.execute("""
            SELECT genome_type, COUNT(*) as count
            FROM genomes
            GROUP BY genome_type
            ORDER BY count DESC
        """)
        types = cursor.fetchall()

        total = sum(row['count'] for row in types)
        selected = []

        for row in types:
            gtype = row['genome_type']
            proportion = row['count'] / total
            type_count = max(1, int(count * proportion))

            cursor = conn.execute("""
                SELECT
                    g.genome_id,
                    g.genome_name,
                    g.length,
                    g.gc_content,
                    g.genome_type,
                    t.family,
                    t.genus
                FROM genomes g
                LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
                WHERE g.genome_type = ?
                ORDER BY RANDOM()
                LIMIT ?
            """, (gtype, type_count))

            selected.extend([dict(row) for row in cursor.fetchall()])

        conn.close()
        random.shuffle(selected)
        return selected[:count]

    def export_subset(self, genomes: List[Dict], output_file: Path, format: str = 'list'):
        """Export subset to file."""
        with open(output_file, 'w') as f:
            if format == 'list':
                # Simple list of genome IDs
                for genome in genomes:
                    f.write(f"{genome['genome_id']}\n")

            elif format == 'tsv':
                # Tab-separated with metadata
                f.write("genome_id\tgenome_name\tlength\tgc_content\tgenome_type\tfamily\tgenus\n")
                for genome in genomes:
                    f.write(f"{genome['genome_id']}\t")
                    f.write(f"{genome['genome_name']}\t")
                    f.write(f"{genome['length']}\t")
                    f.write(f"{genome.get('gc_content', '')}\t")
                    f.write(f"{genome['genome_type']}\t")
                    f.write(f"{genome.get('family', '')}\t")
                    f.write(f"{genome.get('genus', '')}\n")

            elif format == 'json':
                import json
                json.dump(genomes, f, indent=2)

    def print_summary(self, genomes: List[Dict]):
        """Print summary of selected genomes."""
        print(f"\nSelected {len(genomes)} genomes")
        print()

        # Length statistics
        lengths = [g['length'] for g in genomes]
        print("Length distribution:")
        print(f"  Min:    {min(lengths):>10,} bp")
        print(f"  Mean:   {sum(lengths)//len(lengths):>10,} bp")
        print(f"  Max:    {max(lengths):>10,} bp")
        print()

        # Genome types
        types = Counter(g['genome_type'] for g in genomes)
        print("Genome types:")
        for gtype, count in types.most_common():
            pct = count / len(genomes)
            print(f"  {gtype:15s} {count:>4} ({pct:>5.1%})")
        print()

        # Families
        families = Counter(g.get('family', 'Unknown') for g in genomes)
        print(f"Families (top 10 of {len(families)}):")
        for family, count in families.most_common(10):
            pct = count / len(genomes)
            print(f"  {family[:35]:35s} {count:>4} ({pct:>5.1%})")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge Genome Subset Selector',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Subset Strategies:
  random                Random selection (default)
  diverse               Taxonomically diverse (even sampling across families)
  family_reps           One or more representatives per family
  size_stratified       Even sampling across size bins (small/medium/large)
  balanced_type         Proportional to genome type distribution

Examples:
  # Random 50 genomes
  viroforge-subset --count 50 --output test_subset.txt

  # Diverse 100 genomes
  viroforge-subset --count 100 --strategy diverse --output diverse_100.txt

  # One representative per family
  viroforge-subset --strategy family_reps --output family_reps.txt

  # Size-stratified 60 genomes
  viroforge-subset --count 60 --strategy size_stratified --output size_subset.txt

  # Subset with filters
  viroforge-subset --count 20 --family Siphoviridae --output sipho_subset.txt

  # Export as TSV with metadata
  viroforge-subset --count 50 --format tsv --output subset.tsv
        """
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    parser.add_argument(
        '--count',
        type=int,
        default=50,
        help='Number of genomes to select (default: 50)'
    )

    parser.add_argument(
        '--strategy',
        choices=['random', 'diverse', 'family_reps', 'size_stratified', 'balanced_type'],
        default='random',
        help='Subset selection strategy (default: random)'
    )

    parser.add_argument(
        '--per-family',
        type=int,
        default=1,
        help='Genomes per family (for family_reps strategy, default: 1)'
    )

    # Filters (only apply to random strategy)
    parser.add_argument('--family', help='Filter by family')
    parser.add_argument('--host', help='Filter by host')
    parser.add_argument('--genome-type', help='Filter by genome type')
    parser.add_argument('--length-min', type=int, help='Minimum length')
    parser.add_argument('--length-max', type=int, help='Maximum length')

    parser.add_argument(
        '--output',
        required=True,
        help='Output file path'
    )

    parser.add_argument(
        '--format',
        choices=['list', 'tsv', 'json'],
        default='list',
        help='Output format (default: list)'
    )

    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )

    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress summary output'
    )

    args = parser.parse_args()

    # Create selector
    selector = SubsetSelector(args.database, args.seed)

    # Select genomes based on strategy
    print(f"Selecting {args.count} genomes using '{args.strategy}' strategy...")

    if args.strategy == 'random':
        filters = {}
        if args.family:
            filters['family'] = args.family
        if args.host:
            filters['host'] = args.host
        if args.genome_type:
            filters['genome_type'] = args.genome_type
        if args.length_min:
            filters['length_min'] = args.length_min
        if args.length_max:
            filters['length_max'] = args.length_max

        genomes = selector.random_subset(args.count, **filters)

    elif args.strategy == 'diverse':
        genomes = selector.diverse_subset(args.count)

    elif args.strategy == 'family_reps':
        genomes = selector.family_representatives(args.per_family)

    elif args.strategy == 'size_stratified':
        genomes = selector.size_stratified_subset(args.count)

    elif args.strategy == 'balanced_type':
        genomes = selector.balanced_genometype_subset(args.count)

    if not genomes:
        print("❌ Error: No genomes match criteria", file=sys.stderr)
        sys.exit(1)

    # Export
    selector.export_subset(genomes, Path(args.output), args.format)

    # Print summary
    if not args.quiet:
        selector.print_summary(genomes)

    print(f"✓ Subset saved to: {args.output}")
    print(f"  Format: {args.format}")
    print(f"  Count: {len(genomes)} genomes")


if __name__ == '__main__':
    main()
