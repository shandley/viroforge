#!/usr/bin/env python3
"""
ViroForge Database Quality Control & Quick Stats

Rapid database health checks, quality metrics, and statistics generation.
Designed for quick validation and monitoring.

Usage:
    # Quick stats
    viroforge-qc

    # Detailed quality report
    viroforge-qc --detailed

    # Export to file
    viroforge-qc --output qc_report.txt

    # Check specific aspects
    viroforge-qc --check taxonomy
    viroforge-qc --check collections
    viroforge-qc --check quality

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional
from collections import Counter
import time


class DatabaseQC:
    """Quick quality control checks for ViroForge database."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)

        if not self.db_path.exists():
            print(f"❌ Error: Database not found at {db_path}", file=sys.stderr)
            sys.exit(1)

    def get_connection(self):
        """Get database connection."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def quick_stats(self) -> Dict:
        """Get quick statistics (fast queries only)."""
        conn = self.get_connection()
        stats = {}

        # Total genomes
        cursor = conn.execute("SELECT COUNT(*) as count FROM genomes")
        stats['total_genomes'] = cursor.fetchone()['count']

        # Database file size
        stats['db_size_mb'] = self.db_path.stat().st_size / (1024**2)

        # Tables present
        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table'
            ORDER BY name
        """)
        stats['tables'] = [row['name'] for row in cursor.fetchall()]

        # Taxonomy coverage
        cursor = conn.execute("""
            SELECT COUNT(*) as count FROM taxonomy WHERE realm IS NOT NULL
        """)
        stats['with_realm'] = cursor.fetchone()['count']

        # Collections
        if 'body_site_collections' in stats['tables']:
            cursor = conn.execute("SELECT COUNT(*) FROM body_site_collections")
            stats['n_collections'] = cursor.fetchone()[0]
        else:
            stats['n_collections'] = 0

        conn.close()
        return stats

    def taxonomy_check(self) -> Dict:
        """Check taxonomy completeness and quality."""
        conn = self.get_connection()

        # Check taxonomy table exists
        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name='taxonomy'
        """)
        if not cursor.fetchone():
            conn.close()
            return {'status': 'MISSING', 'message': 'Taxonomy table not found'}

        total_genomes = conn.execute("SELECT COUNT(*) FROM genomes").fetchone()[0]

        # Count by rank
        ranks = {}
        rank_columns = {
            'realm': 'realm',
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'class': 'class',
            'order': 'order_name',
            'family': 'family',
            'genus': 'genus',
            'species': 'species'
        }

        for rank_name, column_name in rank_columns.items():
            cursor = conn.execute(f"""
                SELECT COUNT(*) FROM taxonomy
                WHERE {column_name} IS NOT NULL AND {column_name} != ''
            """)
            ranks[rank_name] = cursor.fetchone()[0]

        # Calculate coverage percentages
        coverage = {
            rank: (count / total_genomes) if total_genomes > 0 else 0
            for rank, count in ranks.items()
        }

        # Top realms
        cursor = conn.execute("""
            SELECT realm, COUNT(*) as count
            FROM taxonomy
            WHERE realm IS NOT NULL
            GROUP BY realm
            ORDER BY count DESC
            LIMIT 5
        """)
        top_realms = dict(cursor.fetchall())

        # Quality assessment
        realm_coverage = coverage.get('realm', 0)
        if realm_coverage >= 0.70:
            status = 'EXCELLENT'
            message = f'Taxonomy coverage: {realm_coverage:.1%} (≥70%)'
        elif realm_coverage >= 0.50:
            status = 'GOOD'
            message = f'Taxonomy coverage: {realm_coverage:.1%} (50-70%)'
        elif realm_coverage >= 0.30:
            status = 'FAIR'
            message = f'Taxonomy coverage: {realm_coverage:.1%} (30-50%)'
        else:
            status = 'POOR'
            message = f'Taxonomy coverage: {realm_coverage:.1%} (<30%)'

        conn.close()

        return {
            'status': status,
            'message': message,
            'total_genomes': total_genomes,
            'ranks': ranks,
            'coverage': coverage,
            'top_realms': top_realms
        }

    def quality_check(self) -> Dict:
        """Check genome quality metrics."""
        conn = self.get_connection()

        # Length distribution
        cursor = conn.execute("""
            SELECT
                MIN(length) as min_length,
                AVG(length) as mean_length,
                MAX(length) as max_length
            FROM genomes
        """)
        row = cursor.fetchone()
        length_stats = dict(row)

        # GC content distribution
        cursor = conn.execute("""
            SELECT
                MIN(gc_content) as min_gc,
                AVG(gc_content) as mean_gc,
                MAX(gc_content) as max_gc
            FROM genomes
            WHERE gc_content IS NOT NULL
        """)
        row = cursor.fetchone()
        gc_stats = dict(row) if row['mean_gc'] else {}

        # Count potential issues
        issues = []

        # Very short genomes
        cursor = conn.execute("SELECT COUNT(*) FROM genomes WHERE length < 1000")
        very_short = cursor.fetchone()[0]
        if very_short > 0:
            issues.append(f"{very_short} genomes < 1kb")

        # Very long genomes (potential assembly artifacts)
        cursor = conn.execute("SELECT COUNT(*) FROM genomes WHERE length > 500000")
        very_long = cursor.fetchone()[0]
        if very_long > 0:
            issues.append(f"{very_long} genomes > 500kb")

        # Extreme GC content
        cursor = conn.execute("""
            SELECT COUNT(*) FROM genomes
            WHERE gc_content IS NOT NULL
            AND (gc_content < 0.15 OR gc_content > 0.75)
        """)
        extreme_gc = cursor.fetchone()[0]
        if extreme_gc > 0:
            issues.append(f"{extreme_gc} genomes with extreme GC (<15% or >75%)")

        # Missing sequences
        cursor = conn.execute("""
            SELECT COUNT(*) FROM genomes
            WHERE sequence IS NULL OR sequence = ''
        """)
        missing_seq = cursor.fetchone()[0]
        if missing_seq > 0:
            issues.append(f"{missing_seq} genomes missing sequences")

        # Overall status
        if not issues:
            status = 'EXCELLENT'
            message = 'All quality checks passed'
        elif len(issues) <= 2:
            status = 'GOOD'
            message = f'{len(issues)} minor issues found'
        else:
            status = 'NEEDS_ATTENTION'
            message = f'{len(issues)} issues found'

        conn.close()

        return {
            'status': status,
            'message': message,
            'length_stats': length_stats,
            'gc_stats': gc_stats,
            'issues': issues
        }

    def collection_check(self) -> Dict:
        """Check collections if present."""
        conn = self.get_connection()

        # Check table exists
        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name='body_site_collections'
        """)
        if not cursor.fetchone():
            conn.close()
            return {
                'status': 'NOT_FOUND',
                'message': 'No collections in database'
            }

        # Check schema
        cursor = conn.execute("PRAGMA table_info(body_site_collections)")
        columns = [row[1] for row in cursor.fetchall()]

        # Count collections
        cursor = conn.execute("SELECT COUNT(*) FROM body_site_collections")
        n_collections = cursor.fetchone()[0]

        if n_collections == 0:
            conn.close()
            return {
                'status': 'EMPTY',
                'message': 'Collection table exists but empty'
            }

        # Get collection info
        if 'collection_name' in columns:
            cursor = conn.execute("""
                SELECT
                    collection_id,
                    collection_name as name,
                    n_genomes as genome_count
                FROM body_site_collections
                ORDER BY n_genomes DESC
            """)
        else:
            cursor = conn.execute("""
                SELECT
                    collection_id,
                    name,
                    genome_count
                FROM body_site_collections
                ORDER BY genome_count DESC
            """)

        collections = [dict(row) for row in cursor.fetchall()]

        # Check genome associations
        cursor = conn.execute("SELECT COUNT(*) FROM collection_genomes")
        n_associations = cursor.fetchone()[0]

        conn.close()

        status = 'EXCELLENT' if n_collections >= 5 else 'GOOD'
        message = f'{n_collections} collections with {n_associations} genome associations'

        return {
            'status': status,
            'message': message,
            'n_collections': n_collections,
            'n_associations': n_associations,
            'collections': collections
        }

    def run_all_checks(self) -> Dict:
        """Run all quality checks."""
        print("Running ViroForge database quality checks...")
        print()

        start_time = time.time()

        checks = {
            'quick_stats': self.quick_stats(),
            'taxonomy': self.taxonomy_check(),
            'quality': self.quality_check(),
            'collections': self.collection_check()
        }

        checks['runtime_seconds'] = time.time() - start_time

        return checks

    def print_report(self, checks: Dict, detailed: bool = False):
        """Print formatted QC report."""
        stats = checks['quick_stats']
        tax = checks['taxonomy']
        qual = checks['quality']
        coll = checks['collections']

        print("=" * 70)
        print("ViroForge Database Quality Control Report")
        print("=" * 70)
        print()

        # Database overview
        print("Database Overview:")
        print("-" * 70)
        print(f"  Database file:         {self.db_path}")
        print(f"  File size:             {stats['db_size_mb']:.1f} MB")
        print(f"  Total genomes:         {stats['total_genomes']:,}")
        print(f"  Tables present:        {len(stats['tables'])}")
        print()

        # Taxonomy check
        print("Taxonomy Check:")
        print("-" * 70)
        status_symbol = self._get_status_symbol(tax['status'])
        print(f"  Status:                {status_symbol} {tax['status']}")
        print(f"  {tax['message']}")
        if detailed and 'coverage' in tax:
            print()
            print("  Coverage by rank:")
            for rank, pct in tax['coverage'].items():
                print(f"    {rank:12s} {pct:>6.1%}")
        print()

        # Quality check
        print("Quality Check:")
        print("-" * 70)
        status_symbol = self._get_status_symbol(qual['status'])
        print(f"  Status:                {status_symbol} {qual['status']}")
        print(f"  {qual['message']}")

        if qual['length_stats']:
            print()
            print(f"  Length statistics:")
            print(f"    Min:    {qual['length_stats']['min_length']:>10,} bp")
            print(f"    Mean:   {qual['length_stats']['mean_length']:>10,.0f} bp")
            print(f"    Max:    {qual['length_stats']['max_length']:>10,} bp")

        if qual['gc_stats']:
            print()
            print(f"  GC content statistics:")
            print(f"    Min:    {qual['gc_stats']['min_gc']:>10.1%}")
            print(f"    Mean:   {qual['gc_stats']['mean_gc']:>10.1%}")
            print(f"    Max:    {qual['gc_stats']['max_gc']:>10.1%}")

        if qual['issues']:
            print()
            print("  Issues found:")
            for issue in qual['issues']:
                print(f"    ⚠️  {issue}")

        print()

        # Collection check
        print("Collection Check:")
        print("-" * 70)
        status_symbol = self._get_status_symbol(coll['status'])
        print(f"  Status:                {status_symbol} {coll['status']}")
        print(f"  {coll['message']}")

        if coll.get('collections'):
            print()
            print("  Collections:")
            for c in coll['collections']:
                name = c['name'][:45]
                count = c['genome_count'] if 'genome_count' in c else c.get('n_genomes', 0)
                print(f"    - {name:45s} {count:>5,} genomes")

        print()

        # Summary
        print("Summary:")
        print("-" * 70)

        all_statuses = [
            tax['status'],
            qual['status'],
            coll['status'] if coll['status'] not in ['NOT_FOUND', 'EMPTY'] else 'OK'
        ]

        if all(s in ['EXCELLENT', 'OK'] for s in all_statuses):
            print("  ✓ Database is in excellent condition")
        elif all(s in ['EXCELLENT', 'GOOD', 'OK'] for s in all_statuses):
            print("  ✓ Database is in good condition")
        elif any(s == 'NEEDS_ATTENTION' for s in all_statuses):
            print("  ⚠️  Database needs attention")
        else:
            print("  ⚠️  Database may have issues")

        print()
        print(f"Quality check completed in {checks['runtime_seconds']:.2f} seconds")
        print()

    def _get_status_symbol(self, status: str) -> str:
        """Get status symbol."""
        symbols = {
            'EXCELLENT': '✓',
            'GOOD': '✓',
            'FAIR': '⚠️',
            'POOR': '❌',
            'NEEDS_ATTENTION': '⚠️',
            'OK': '✓',
            'NOT_FOUND': 'ℹ️',
            'EMPTY': 'ℹ️',
            'MISSING': '❌'
        }
        return symbols.get(status, '?')


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge Database Quality Control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick stats
  viroforge-qc

  # Detailed report
  viroforge-qc --detailed

  # Save to file
  viroforge-qc --output qc_report.txt

  # Check specific aspect
  viroforge-qc --check taxonomy
  viroforge-qc --check quality
  viroforge-qc --check collections
        """
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    parser.add_argument(
        '--detailed',
        action='store_true',
        help='Show detailed report'
    )

    parser.add_argument(
        '--check',
        choices=['taxonomy', 'quality', 'collections', 'all'],
        help='Check specific aspect only'
    )

    parser.add_argument(
        '--output',
        help='Save report to file'
    )

    parser.add_argument(
        '--format',
        choices=['text', 'json'],
        default='text',
        help='Output format (default: text)'
    )

    args = parser.parse_args()

    # Create QC checker
    qc = DatabaseQC(args.database)

    # Run checks
    if args.check and args.check != 'all':
        # Run specific check
        if args.check == 'taxonomy':
            result = qc.taxonomy_check()
        elif args.check == 'quality':
            result = qc.quality_check()
        elif args.check == 'collections':
            result = qc.collection_check()

        if args.format == 'json':
            import json
            print(json.dumps(result, indent=2))
        else:
            print(f"\n{args.check.capitalize()} Check:")
            print("-" * 70)
            print(f"Status: {qc._get_status_symbol(result['status'])} {result['status']}")
            print(f"{result['message']}")
            print()

    else:
        # Run all checks
        checks = qc.run_all_checks()

        if args.format == 'json':
            import json
            print(json.dumps(checks, indent=2, default=str))
        else:
            if args.output:
                # Save to file
                original_stdout = sys.stdout
                with open(args.output, 'w') as f:
                    sys.stdout = f
                    qc.print_report(checks, args.detailed)
                sys.stdout = original_stdout
                print(f"✓ QC report saved to: {args.output}")
            else:
                qc.print_report(checks, args.detailed)


if __name__ == '__main__':
    main()
