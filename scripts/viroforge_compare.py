#!/usr/bin/env python3
"""
ViroForge Collection Comparison Tool

Compare taxonomic composition, abundance distributions, and diversity metrics
between body site collections.

Usage:
    # Compare two collections
    viroforge-compare gut oral

    # Compare multiple collections
    viroforge-compare gut oral skin respiratory

    # Generate detailed report
    viroforge-compare gut oral --report comparison_report.txt

    # Export comparison data
    viroforge-compare gut oral skin --export comparison.tsv

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
import numpy as np


class CollectionComparator:
    """Compare ViroForge body site collections."""

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

    def find_collection(self, collection_query: str) -> Optional[str]:
        """Find collection ID by partial match."""
        conn = self.get_connection()

        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name='body_site_collections'
        """)
        if not cursor.fetchone():
            print(f"⚠️  No collections found in database", file=sys.stderr)
            print(f"   Run curation scripts first", file=sys.stderr)
            conn.close()
            return None

        # Check schema
        cursor = conn.execute("PRAGMA table_info(body_site_collections)")
        columns = [row[1] for row in cursor.fetchall()]

        if 'collection_name' in columns:
            # Old schema
            cursor = conn.execute("""
                SELECT collection_id, collection_name as name
                FROM body_site_collections
                WHERE collection_id LIKE ? OR collection_name LIKE ?
            """, (f"%{collection_query}%", f"%{collection_query}%"))
        else:
            # New schema
            cursor = conn.execute("""
                SELECT collection_id, name
                FROM body_site_collections
                WHERE collection_id LIKE ? OR name LIKE ?
            """, (f"%{collection_query}%", f"%{collection_query}%"))

        results = cursor.fetchall()
        conn.close()

        if not results:
            return None
        if len(results) > 1:
            print(f"⚠️  Multiple collections match '{collection_query}':", file=sys.stderr)
            for row in results:
                print(f"   - {row['collection_id']}: {row['name']}", file=sys.stderr)
            return None

        return results[0]['collection_id']

    def get_collection_info(self, collection_id: str) -> Dict:
        """Get collection metadata and genomes."""
        conn = self.get_connection()

        # Check schema
        cursor = conn.execute("PRAGMA table_info(body_site_collections)")
        columns = [row[1] for row in cursor.fetchall()]

        if 'collection_name' in columns:
            # Old schema
            cursor = conn.execute("""
                SELECT
                    collection_id,
                    collection_name as name,
                    '' as body_site,
                    n_genomes as genome_count,
                    NULL as phage_fraction
                FROM body_site_collections
                WHERE collection_id = ?
            """, (collection_id,))
        else:
            # New schema
            cursor = conn.execute("""
                SELECT *
                FROM body_site_collections
                WHERE collection_id = ?
            """, (collection_id,))

        meta = cursor.fetchone()
        if not meta:
            conn.close()
            return None

        meta = dict(meta)

        # Get genomes with abundances
        cursor = conn.execute("""
            SELECT
                cg.genome_id,
                cg.abundance,
                g.genome_name as species_name,
                g.length,
                g.gc_content,
                g.genome_type,
                t.realm,
                t.family,
                t.genus
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE cg.collection_id = ?
        """, (collection_id,))

        genomes = [dict(row) for row in cursor.fetchall()]
        conn.close()

        return {
            'metadata': meta,
            'genomes': genomes
        }

    def calculate_diversity(self, abundances: np.ndarray) -> Dict:
        """Calculate diversity metrics."""
        # Remove zeros
        abundances = abundances[abundances > 0]

        # Shannon diversity
        shannon = -np.sum(abundances * np.log(abundances))

        # Simpson diversity
        simpson = 1.0 - np.sum(abundances ** 2)

        # Evenness
        n_species = len(abundances)
        evenness = shannon / np.log(n_species) if n_species > 1 else 1.0

        # Dominance (max abundance)
        dominance = abundances.max()

        return {
            'shannon': shannon,
            'simpson': simpson,
            'evenness': evenness,
            'dominance': dominance,
            'n_species': n_species
        }

    def get_taxonomic_composition(self, genomes: List[Dict], rank: str = 'family') -> Counter:
        """Get taxonomic composition at specified rank."""
        composition = Counter()

        for genome in genomes:
            taxon = genome.get(rank)
            if taxon:
                composition[taxon] += 1
            else:
                composition['Unknown'] += 1

        return composition

    def calculate_overlap(self, genomes1: List[Dict], genomes2: List[Dict]) -> Dict:
        """Calculate genome overlap between collections."""
        ids1 = set(g['genome_id'] for g in genomes1)
        ids2 = set(g['genome_id'] for g in genomes2)

        shared = ids1 & ids2
        unique1 = ids1 - ids2
        unique2 = ids2 - ids1

        jaccard = len(shared) / len(ids1 | ids2) if (ids1 | ids2) else 0

        return {
            'shared': len(shared),
            'unique_to_first': len(unique1),
            'unique_to_second': len(unique2),
            'jaccard_index': jaccard
        }

    def compare_collections(self, collection_ids: List[str]) -> Dict:
        """Compare multiple collections."""
        # Get collection info
        collections = {}
        for coll_id in collection_ids:
            info = self.get_collection_info(coll_id)
            if not info:
                print(f"❌ Error: Collection '{coll_id}' not found", file=sys.stderr)
                sys.exit(1)
            collections[coll_id] = info

        comparison = {
            'collections': collections,
            'diversity': {},
            'taxonomy': {},
            'genome_stats': {},
            'pairwise_overlap': {}
        }

        # Calculate diversity metrics for each
        for coll_id, info in collections.items():
            abundances = np.array([g['abundance'] for g in info['genomes']])
            comparison['diversity'][coll_id] = self.calculate_diversity(abundances)

            # Taxonomic composition
            comparison['taxonomy'][coll_id] = {
                'realm': self.get_taxonomic_composition(info['genomes'], 'realm'),
                'family': self.get_taxonomic_composition(info['genomes'], 'family'),
                'genus': self.get_taxonomic_composition(info['genomes'], 'genus')
            }

            # Genome statistics
            lengths = [g['length'] for g in info['genomes']]
            gc_contents = [g['gc_content'] for g in info['genomes'] if g['gc_content'] is not None]
            genome_types = Counter(g['genome_type'] for g in info['genomes'])

            comparison['genome_stats'][coll_id] = {
                'n_genomes': len(info['genomes']),
                'mean_length': np.mean(lengths),
                'median_length': np.median(lengths),
                'mean_gc': np.mean(gc_contents) if gc_contents else None,
                'genome_types': genome_types
            }

        # Pairwise overlaps
        if len(collection_ids) == 2:
            coll1, coll2 = collection_ids
            overlap = self.calculate_overlap(
                collections[coll1]['genomes'],
                collections[coll2]['genomes']
            )
            comparison['pairwise_overlap'] = {
                f"{coll1}_vs_{coll2}": overlap
            }

        return comparison

    def print_comparison(self, comparison: Dict):
        """Print formatted comparison report."""
        collections = comparison['collections']
        coll_ids = list(collections.keys())

        print("=" * 80)
        print("ViroForge Collection Comparison")
        print("=" * 80)
        print()

        # Overview table
        print("Collection Overview:")
        print("-" * 80)
        print(f"{'Collection':<40} {'Genomes':>10} {'Body Site':<20}")
        print("-" * 80)

        for coll_id in coll_ids:
            meta = collections[coll_id]['metadata']
            name = meta.get('name', meta.get('collection_name', coll_id))[:38]
            body_site = meta.get('body_site', '-')
            genome_count = meta.get('genome_count', meta.get('n_genomes', len(collections[coll_id]['genomes'])))
            print(f"{name:<40} {genome_count:>10} {body_site:<20}")

        print()

        # Diversity comparison
        print("Diversity Metrics:")
        print("-" * 80)
        print(f"{'Collection':<40} {'Shannon':>10} {'Simpson':>10} {'Evenness':>10}")
        print("-" * 80)

        for coll_id in coll_ids:
            meta = collections[coll_id]['metadata']
            name = meta.get('name', meta.get('collection_name', coll_id))[:38]
            div = comparison['diversity'][coll_id]
            print(f"{name:<40} {div['shannon']:>10.2f} {div['simpson']:>10.4f} {div['evenness']:>10.3f}")

        print()

        # Genome statistics
        print("Genome Characteristics:")
        print("-" * 80)
        print(f"{'Collection':<40} {'Mean Length':>12} {'Mean GC':>10}")
        print("-" * 80)

        for coll_id in coll_ids:
            meta = collections[coll_id]['metadata']
            name = meta.get('name', meta.get('collection_name', coll_id))[:38]
            stats = comparison['genome_stats'][coll_id]
            gc_str = f"{stats['mean_gc']:.1%}" if stats['mean_gc'] else "N/A"
            print(f"{name:<40} {stats['mean_length']:>12,.0f} bp {gc_str:>10}")

        print()

        # Taxonomic composition comparison (top 10 families)
        print("Taxonomic Composition (Top 10 Families per Collection):")
        print("-" * 80)

        for coll_id in coll_ids:
            meta = collections[coll_id]['metadata']
            name = meta.get('name', meta.get('collection_name', coll_id))
            print(f"\n{name}:")

            family_comp = comparison['taxonomy'][coll_id]['family']
            total = sum(family_comp.values())

            for family, count in family_comp.most_common(10):
                pct = count / total
                print(f"  {family[:45]:45s} {count:>4} ({pct:>5.1%})")

        print()

        # Pairwise overlap (if comparing 2 collections)
        if comparison['pairwise_overlap']:
            print("Genome Overlap:")
            print("-" * 80)

            for pair, overlap in comparison['pairwise_overlap'].items():
                coll1, coll2 = coll_ids
                meta1 = collections[coll1]['metadata']
                meta2 = collections[coll2]['metadata']
                name1 = meta1.get('name', meta1.get('collection_name', coll1))
                name2 = meta2.get('name', meta2.get('collection_name', coll2))

                print(f"\n{name1} vs {name2}:")
                print(f"  Shared genomes:        {overlap['shared']:>6,}")
                print(f"  Unique to first:       {overlap['unique_to_first']:>6,}")
                print(f"  Unique to second:      {overlap['unique_to_second']:>6,}")
                print(f"  Jaccard index:         {overlap['jaccard_index']:>6.3f}")

            print()


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge Collection Comparison Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two collections
  viroforge-compare gut oral

  # Compare multiple collections
  viroforge-compare gut oral skin respiratory

  # Generate detailed report
  viroforge-compare gut oral --report comparison_report.txt
        """
    )

    parser.add_argument(
        'collections',
        nargs='+',
        help='Collections to compare (partial names supported)'
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    parser.add_argument(
        '--report',
        help='Save comparison report to file'
    )

    args = parser.parse_args()

    # Create comparator
    comparator = CollectionComparator(args.database)

    # Find collections
    collection_ids = []
    for query in args.collections:
        coll_id = comparator.find_collection(query)
        if not coll_id:
            print(f"❌ Error: Collection '{query}' not found", file=sys.stderr)
            sys.exit(1)
        collection_ids.append(coll_id)

    print(f"Comparing {len(collection_ids)} collections...")
    print()

    # Compare collections
    comparison = comparator.compare_collections(collection_ids)

    # Print comparison
    if args.report:
        # Redirect to file
        original_stdout = sys.stdout
        with open(args.report, 'w') as f:
            sys.stdout = f
            comparator.print_comparison(comparison)
        sys.stdout = original_stdout
        print(f"✓ Comparison report saved to: {args.report}")
    else:
        comparator.print_comparison(comparison)


if __name__ == '__main__':
    main()
