#!/usr/bin/env python3
"""
Body Site Collection Validation Script

Validates curated body site collections against specifications and generates
comprehensive quality metrics.

Usage:
    python scripts/validate_body_site_collections.py --collection gut
    python scripts/validate_body_site_collections.py --all

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
from collections import Counter, defaultdict


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CollectionValidator:
    """Validates body site virome collections."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)

        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

    def validate_collection(self, collection_id: str) -> Dict:
        """
        Validate a collection and generate comprehensive metrics.

        Returns:
            Dictionary with validation results and statistics
        """
        logger.info(f"Validating collection: {collection_id}")

        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Get collection metadata
        collection_meta = self._get_collection_metadata(conn, collection_id)

        if not collection_meta:
            raise ValueError(f"Collection not found: {collection_id}")

        # Get genomes in collection
        genomes = self._get_collection_genomes(conn, collection_id)

        # Get taxonomy for genomes
        taxonomy = self._get_genome_taxonomy(conn, [g['genome_id'] for g in genomes])

        # Get genome metadata
        genome_meta = self._get_genome_metadata(conn, [g['genome_id'] for g in genomes])

        conn.close()

        # Validate
        validation = {
            'collection_id': collection_id,
            'collection_meta': collection_meta,
            'genome_count': len(genomes),
            'expected_count': collection_meta['n_genomes'],
            'count_match': len(genomes) == collection_meta['n_genomes']
        }

        # Taxonomic composition
        validation['taxonomy'] = self._analyze_taxonomy(genomes, taxonomy)

        # Abundance distribution
        validation['abundances'] = self._analyze_abundances(genomes)

        # Genome characteristics
        validation['genome_stats'] = self._analyze_genomes(genome_meta)

        # Quality checks
        validation['quality_checks'] = self._run_quality_checks(validation)

        return validation

    def _get_collection_metadata(self, conn: sqlite3.Connection, collection_id: str) -> Dict:
        """Get collection metadata."""
        cursor = conn.execute("""
            SELECT * FROM body_site_collections
            WHERE collection_id = ?
        """, (collection_id,))

        row = cursor.fetchone()
        if not row:
            return None

        return dict(row)

    def _get_collection_genomes(self, conn: sqlite3.Connection, collection_id: str) -> List[Dict]:
        """Get genomes in collection with abundances."""
        cursor = conn.execute("""
            SELECT genome_id, relative_abundance as abundance, abundance_rank
            FROM collection_genomes
            WHERE collection_id = ?
            ORDER BY relative_abundance DESC
        """, (collection_id,))

        return [dict(row) for row in cursor.fetchall()]

    def _get_genome_taxonomy(self, conn: sqlite3.Connection, genome_ids: List[str]) -> Dict:
        """Get taxonomy for genomes."""
        placeholders = ','.join('?' * len(genome_ids))
        cursor = conn.execute(f"""
            SELECT genome_id, realm, kingdom, phylum, class, order_name,
                   family, subfamily, genus, species
            FROM taxonomy
            WHERE genome_id IN ({placeholders})
        """, genome_ids)

        taxonomy = {}
        for row in cursor.fetchall():
            taxonomy[row['genome_id']] = dict(row)

        return taxonomy

    def _get_genome_metadata(self, conn: sqlite3.Connection, genome_ids: List[str]) -> Dict:
        """Get genome metadata."""
        placeholders = ','.join('?' * len(genome_ids))
        cursor = conn.execute(f"""
            SELECT genome_id, genome_name, genome_type, length, gc_content
            FROM genomes
            WHERE genome_id IN ({placeholders})
        """, genome_ids)

        metadata = {}
        for row in cursor.fetchall():
            metadata[row['genome_id']] = dict(row)

        return metadata

    def _analyze_taxonomy(self, genomes: List[Dict], taxonomy: Dict) -> Dict:
        """Analyze taxonomic composition."""
        # Count by rank
        realm_counts = Counter()
        class_counts = Counter()
        family_counts = Counter()
        genus_counts = Counter()

        for genome in genomes:
            genome_id = genome['genome_id']
            if genome_id not in taxonomy:
                continue

            tax = taxonomy[genome_id]
            if tax['realm']:
                realm_counts[tax['realm']] += 1
            if tax['class']:
                class_counts[tax['class']] += 1
            if tax['family']:
                family_counts[tax['family']] += 1
            if tax['genus']:
                genus_counts[tax['genus']] += 1

        return {
            'realm_distribution': dict(realm_counts.most_common(10)),
            'class_distribution': dict(class_counts.most_common(10)),
            'family_distribution': dict(family_counts.most_common(10)),
            'genus_distribution': dict(genus_counts.most_common(10)),
            'taxonomy_coverage': len([g for g in genomes if g['genome_id'] in taxonomy]) / len(genomes)
        }

    def _analyze_abundances(self, genomes: List[Dict]) -> Dict:
        """Analyze abundance distribution."""
        abundances = np.array([g['abundance'] for g in genomes])

        # Identify tiers based on abundance
        dominant = abundances >= 0.10  # 10%+
        common = (abundances >= 0.01) & (abundances < 0.10)  # 1-10%
        moderate = (abundances >= 0.001) & (abundances < 0.01)  # 0.1-1%
        rare = (abundances >= 0.0001) & (abundances < 0.001)  # 0.01-0.1%
        very_rare = abundances < 0.0001  # <0.01%

        return {
            'total_genomes': len(abundances),
            'sum_abundances': float(abundances.sum()),
            'max_abundance': float(abundances.max()),
            'mean_abundance': float(abundances.mean()),
            'median_abundance': float(np.median(abundances)),
            'min_abundance': float(abundances.min()),
            'shannon_diversity': float(self._calculate_shannon(abundances)),
            'simpson_diversity': float(self._calculate_simpson(abundances)),
            'tier_counts': {
                'dominant_10_30': int(dominant.sum()),
                'common_1_10': int(common.sum()),
                'moderate_01_1': int(moderate.sum()),
                'rare_001_01': int(rare.sum()),
                'very_rare_lt001': int(very_rare.sum())
            },
            'tier_abundances': {
                'dominant_10_30': float(abundances[dominant].sum()),
                'common_1_10': float(abundances[common].sum()),
                'moderate_01_1': float(abundances[moderate].sum()),
                'rare_001_01': float(abundances[rare].sum()),
                'very_rare_lt001': float(abundances[very_rare].sum())
            }
        }

    def _analyze_genomes(self, genome_meta: Dict) -> Dict:
        """Analyze genome characteristics."""
        lengths = [meta['length'] for meta in genome_meta.values()]
        gc_contents = [meta['gc_content'] for meta in genome_meta.values() if meta['gc_content'] is not None]
        genome_types = Counter([meta['genome_type'] for meta in genome_meta.values()])

        return {
            'genome_count': len(genome_meta),
            'length_stats': {
                'mean': float(np.mean(lengths)),
                'median': float(np.median(lengths)),
                'min': int(np.min(lengths)),
                'max': int(np.max(lengths)),
                'std': float(np.std(lengths))
            },
            'gc_stats': {
                'mean': float(np.mean(gc_contents)) if gc_contents else None,
                'median': float(np.median(gc_contents)) if gc_contents else None,
                'min': float(np.min(gc_contents)) if gc_contents else None,
                'max': float(np.max(gc_contents)) if gc_contents else None,
                'std': float(np.std(gc_contents)) if gc_contents else None
            },
            'genome_type_distribution': dict(genome_types)
        }

    def _run_quality_checks(self, validation: Dict) -> Dict:
        """Run quality checks on collection."""
        checks = {}

        # Check genome count matches expected
        checks['genome_count'] = {
            'status': 'PASS' if validation['count_match'] else 'FAIL',
            'expected': validation['expected_count'],
            'actual': validation['genome_count']
        }

        # Check abundance sum is approximately 1.0
        abundance_sum = validation['abundances']['sum_abundances']
        checks['abundance_sum'] = {
            'status': 'PASS' if 0.99 <= abundance_sum <= 1.01 else 'FAIL',
            'expected': 1.0,
            'actual': abundance_sum
        }

        # Check taxonomy coverage
        tax_coverage = validation['taxonomy']['taxonomy_coverage']
        checks['taxonomy_coverage'] = {
            'status': 'PASS' if tax_coverage >= 0.70 else 'WARN' if tax_coverage >= 0.50 else 'FAIL',
            'threshold': 0.70,
            'actual': tax_coverage
        }

        # Check diversity metrics
        shannon = validation['abundances']['shannon_diversity']
        checks['shannon_diversity'] = {
            'status': 'PASS' if shannon >= 3.0 else 'WARN' if shannon >= 2.0 else 'FAIL',
            'threshold': 3.0,
            'actual': shannon
        }

        # Check tier distribution
        tier_counts = validation['abundances']['tier_counts']
        total_genomes = validation['genome_count']

        # At least 1% in dominant tier
        dominant_pct = tier_counts['dominant_10_30'] / total_genomes
        checks['dominant_tier'] = {
            'status': 'PASS' if dominant_pct >= 0.01 else 'WARN',
            'threshold': 0.01,
            'actual': dominant_pct
        }

        # Most genomes in moderate/rare tiers
        mid_tier_pct = (tier_counts['moderate_01_1'] + tier_counts['rare_001_01']) / total_genomes
        checks['mid_tier_distribution'] = {
            'status': 'PASS' if mid_tier_pct >= 0.50 else 'WARN',
            'threshold': 0.50,
            'actual': mid_tier_pct
        }

        return checks

    def _calculate_shannon(self, abundances: np.ndarray) -> float:
        """Calculate Shannon diversity index."""
        abundances = abundances[abundances > 0]
        return -np.sum(abundances * np.log(abundances))

    def _calculate_simpson(self, abundances: np.ndarray) -> float:
        """Calculate Simpson diversity index."""
        return 1.0 - np.sum(abundances ** 2)

    def generate_report(self, validation: Dict, output_path: Path):
        """Generate comprehensive validation report."""
        with open(output_path, 'w') as f:
            meta = validation['collection_meta']

            f.write("=" * 80 + "\n")
            f.write("BODY SITE VIROME COLLECTION VALIDATION REPORT\n")
            f.write("=" * 80 + "\n\n")

            # Collection metadata
            f.write("Collection Information:\n")
            f.write(f"  ID: {validation['collection_id']}\n")
            f.write(f"  Name: {meta['collection_name']}\n")
            f.write(f"  Description: {meta.get('description', 'N/A')}\n")
            f.write(f"  Selection Criteria: {meta.get('selection_criteria', 'N/A')}\n")
            f.write(f"  Curation Date: {meta['curation_date']}\n")
            f.write(f"  Expected Genome Count: {meta['n_genomes']}\n\n")

            # Quality checks
            f.write("Quality Checks:\n")
            f.write("-" * 80 + "\n")
            checks = validation['quality_checks']

            for check_name, check_data in checks.items():
                status_symbol = "✓" if check_data['status'] == 'PASS' else \
                               "⚠" if check_data['status'] == 'WARN' else "✗"

                f.write(f"  {status_symbol} {check_name}: {check_data['status']}\n")

                if 'expected' in check_data:
                    f.write(f"      Expected: {check_data['expected']}\n")
                if 'threshold' in check_data:
                    f.write(f"      Threshold: {check_data['threshold']}\n")
                if 'actual' in check_data:
                    if isinstance(check_data['actual'], float):
                        f.write(f"      Actual: {check_data['actual']:.4f}\n")
                    else:
                        f.write(f"      Actual: {check_data['actual']}\n")
                f.write("\n")

            # Abundance statistics
            f.write("\nAbundance Distribution:\n")
            f.write("-" * 80 + "\n")
            abund = validation['abundances']

            f.write(f"  Total Genomes: {abund['total_genomes']}\n")
            f.write(f"  Abundance Sum: {abund['sum_abundances']:.6f}\n")
            f.write(f"  Max Abundance: {abund['max_abundance']:.4%}\n")
            f.write(f"  Mean Abundance: {abund['mean_abundance']:.4%}\n")
            f.write(f"  Median Abundance: {abund['median_abundance']:.4%}\n")
            f.write(f"  Min Abundance: {abund['min_abundance']:.6%}\n\n")

            f.write(f"  Shannon Diversity: {abund['shannon_diversity']:.2f}\n")
            f.write(f"  Simpson Diversity: {abund['simpson_diversity']:.4f}\n\n")

            f.write("  Tier Distribution (genome counts):\n")
            tier_counts = abund['tier_counts']
            tier_abunds = abund['tier_abundances']

            f.write(f"    Dominant (10-30%):   {tier_counts['dominant_10_30']:4d} genomes ")
            f.write(f"({tier_abunds['dominant_10_30']:.1%} of total abundance)\n")

            f.write(f"    Common (1-10%):      {tier_counts['common_1_10']:4d} genomes ")
            f.write(f"({tier_abunds['common_1_10']:.1%} of total abundance)\n")

            f.write(f"    Moderate (0.1-1%):   {tier_counts['moderate_01_1']:4d} genomes ")
            f.write(f"({tier_abunds['moderate_01_1']:.1%} of total abundance)\n")

            f.write(f"    Rare (0.01-0.1%):    {tier_counts['rare_001_01']:4d} genomes ")
            f.write(f"({tier_abunds['rare_001_01']:.1%} of total abundance)\n")

            f.write(f"    Very Rare (<0.01%):  {tier_counts['very_rare_lt001']:4d} genomes ")
            f.write(f"({tier_abunds['very_rare_lt001']:.1%} of total abundance)\n\n")

            # Taxonomic composition
            f.write("\nTaxonomic Composition:\n")
            f.write("-" * 80 + "\n")
            tax = validation['taxonomy']

            f.write(f"  Taxonomy Coverage: {tax['taxonomy_coverage']:.1%}\n\n")

            f.write("  Top 10 Realms:\n")
            for realm, count in tax['realm_distribution'].items():
                pct = count / abund['total_genomes']
                f.write(f"    {realm:40s} {count:4d} ({pct:.1%})\n")

            f.write("\n  Top 10 Classes:\n")
            for cls, count in tax['class_distribution'].items():
                pct = count / abund['total_genomes']
                f.write(f"    {cls:40s} {count:4d} ({pct:.1%})\n")

            f.write("\n  Top 10 Families:\n")
            for family, count in tax['family_distribution'].items():
                pct = count / abund['total_genomes']
                f.write(f"    {family:40s} {count:4d} ({pct:.1%})\n")

            # Genome characteristics
            f.write("\n\nGenome Characteristics:\n")
            f.write("-" * 80 + "\n")
            stats = validation['genome_stats']

            f.write("  Genome Lengths:\n")
            f.write(f"    Mean:   {stats['length_stats']['mean']:>10,.0f} bp\n")
            f.write(f"    Median: {stats['length_stats']['median']:>10,.0f} bp\n")
            f.write(f"    Range:  {stats['length_stats']['min']:>10,d} - {stats['length_stats']['max']:>10,d} bp\n")
            f.write(f"    Std Dev: {stats['length_stats']['std']:>10,.0f} bp\n\n")

            if stats['gc_stats']['mean'] is not None:
                f.write("  GC Content:\n")
                f.write(f"    Mean:   {stats['gc_stats']['mean']:>6.2%}\n")
                f.write(f"    Median: {stats['gc_stats']['median']:>6.2%}\n")
                f.write(f"    Range:  {stats['gc_stats']['min']:>6.2%} - {stats['gc_stats']['max']:>6.2%}\n")
                f.write(f"    Std Dev: {stats['gc_stats']['std']:>6.2%}\n\n")

            f.write("  Genome Type Distribution:\n")
            for gtype, count in stats['genome_type_distribution'].items():
                pct = count / stats['genome_count']
                f.write(f"    {gtype:15s} {count:4d} ({pct:.1%})\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("End of Validation Report\n")
            f.write("=" * 80 + "\n")

        logger.info(f"✓ Validation report written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Validate body site virome collections'
    )
    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to viral genomes database'
    )
    parser.add_argument(
        '--collection',
        help='Collection ID to validate'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Validate all collections in database'
    )
    parser.add_argument(
        '--output',
        default='data/body_site_collections/validation',
        help='Output directory for validation reports'
    )

    args = parser.parse_args()

    if not args.collection and not args.all:
        parser.error("Must specify either --collection or --all")

    # Create validator
    validator = CollectionValidator(args.database)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get collection IDs
    if args.all:
        conn = sqlite3.connect(args.database)
        cursor = conn.execute("SELECT collection_id FROM body_site_collections")
        collection_ids = [row[0] for row in cursor.fetchall()]
        conn.close()

        if not collection_ids:
            logger.error("No collections found in database")
            return
    else:
        collection_ids = [args.collection]

    logger.info(f"Validating {len(collection_ids)} collection(s)")

    # Validate each collection
    for collection_id in collection_ids:
        try:
            logger.info(f"\nValidating: {collection_id}")
            validation = validator.validate_collection(collection_id)

            # Generate report
            report_path = output_dir / f"{collection_id}_validation.txt"
            validator.generate_report(validation, report_path)

            # Print summary
            checks = validation['quality_checks']
            passed = sum(1 for c in checks.values() if c['status'] == 'PASS')
            warned = sum(1 for c in checks.values() if c['status'] == 'WARN')
            failed = sum(1 for c in checks.values() if c['status'] == 'FAIL')

            logger.info(f"  Quality Checks: {passed} passed, {warned} warnings, {failed} failed")

        except Exception as e:
            logger.error(f"✗ Failed to validate {collection_id}: {e}")
            continue

    logger.info(f"\nValidation complete!")


if __name__ == '__main__':
    main()
