#!/usr/bin/env python3
"""
ViroForge Hybrid Composition Validator

Validate that short-read and long-read datasets have matching genome compositions
for hybrid assembly. Checks that both datasets were generated from the same
collection with the same random seed.

Usage:
    python validate_hybrid_composition.py \\
        --short data/gut_hybrid/short_reads/metadata/gut_virome_composition.tsv \\
        --long data/gut_hybrid/long_reads/metadata/gut_virome_composition.tsv

    python validate_hybrid_composition.py \\
        --hybrid-dir data/gut_hybrid

Author: ViroForge Development Team
Date: 2025-11-10
"""

import argparse
import sys
import json
from pathlib import Path
import pandas as pd


def load_composition(tsv_path):
    """Load composition TSV file."""
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading {tsv_path}: {e}", file=sys.stderr)
        return None


def load_metadata(json_path):
    """Load metadata JSON file."""
    try:
        with open(json_path) as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {json_path}: {e}", file=sys.stderr)
        return None


def validate_compositions(short_df, long_df, tolerance=1e-6):
    """
    Validate that two composition dataframes match.

    Args:
        short_df: Short-read composition DataFrame
        long_df: Long-read composition DataFrame
        tolerance: Tolerance for abundance comparison

    Returns:
        Tuple of (is_valid, issues)
    """
    issues = []

    # Check genome IDs match
    short_ids = set(short_df['genome_id'])
    long_ids = set(long_df['genome_id'])

    if short_ids != long_ids:
        missing_in_long = short_ids - long_ids
        missing_in_short = long_ids - short_ids

        if missing_in_long:
            issues.append(f"Genomes in short but not long: {missing_in_long}")
        if missing_in_short:
            issues.append(f"Genomes in long but not short: {missing_in_short}")

        return False, issues

    # Check abundances match (allowing for small numerical differences)
    merged = short_df.merge(long_df, on='genome_id', suffixes=('_short', '_long'))

    abundance_diff = abs(merged['relative_abundance_short'] - merged['relative_abundance_long'])
    max_diff = abundance_diff.max()

    if max_diff > tolerance:
        issues.append(f"Abundance mismatch: max difference = {max_diff:.10f} (tolerance = {tolerance})")

        # Report top mismatches
        top_mismatches = merged.nlargest(5, abundance_diff.index)[['genome_id', 'relative_abundance_short', 'relative_abundance_long']]
        issues.append("Top 5 mismatches:")
        for _, row in top_mismatches.iterrows():
            diff = abs(row['relative_abundance_short'] - row['relative_abundance_long'])
            issues.append(f"  {row['genome_id']}: {row['relative_abundance_short']:.8f} vs {row['relative_abundance_long']:.8f} (diff={diff:.10f})")

        return False, issues

    return True, issues


def main():
    parser = argparse.ArgumentParser(
        description='Validate hybrid dataset composition matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate using composition TSV files
  python validate_hybrid_composition.py \\
      --short data/gut_hybrid/short_reads/metadata/*_composition.tsv \\
      --long data/gut_hybrid/long_reads/metadata/*_composition.tsv

  # Validate using hybrid directory
  python validate_hybrid_composition.py \\
      --hybrid-dir data/gut_hybrid
        """
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--hybrid-dir',
        help='Hybrid dataset directory (will auto-find composition files)'
    )
    group.add_argument(
        '--short',
        help='Short-read composition TSV file'
    )

    parser.add_argument(
        '--long',
        help='Long-read composition TSV file (required if using --short)'
    )

    parser.add_argument(
        '--tolerance',
        type=float,
        default=1e-6,
        help='Tolerance for abundance comparison (default: 1e-6)'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show detailed comparison'
    )

    args = parser.parse_args()

    # Handle hybrid directory mode
    if args.hybrid_dir:
        hybrid_dir = Path(args.hybrid_dir)

        if not hybrid_dir.exists():
            print(f"Error: Hybrid directory not found: {hybrid_dir}", file=sys.stderr)
            sys.exit(1)

        # Find composition files
        short_comps = list((hybrid_dir / "short_reads" / "metadata").glob("*_composition.tsv"))
        long_comps = list((hybrid_dir / "long_reads" / "metadata").glob("*_composition.tsv"))

        if not short_comps:
            print(f"Error: No composition file found in {hybrid_dir}/short_reads/metadata/", file=sys.stderr)
            sys.exit(1)

        if not long_comps:
            print(f"Error: No composition file found in {hybrid_dir}/long_reads/metadata/", file=sys.stderr)
            sys.exit(1)

        short_path = short_comps[0]
        long_path = long_comps[0]

        print(f"Found short-read composition: {short_path}")
        print(f"Found long-read composition: {long_path}")

    else:
        # Handle explicit file paths
        if not args.long:
            parser.error("--long required when using --short")

        short_path = Path(args.short)
        long_path = Path(args.long)

        if not short_path.exists():
            print(f"Error: File not found: {short_path}", file=sys.stderr)
            sys.exit(1)

        if not long_path.exists():
            print(f"Error: File not found: {long_path}", file=sys.stderr)
            sys.exit(1)

    print(f"\n{'='*80}")
    print("ViroForge Hybrid Composition Validator")
    print(f"{'='*80}\n")

    # Load compositions
    print("Loading compositions...")
    short_df = load_composition(short_path)
    long_df = load_composition(long_path)

    if short_df is None or long_df is None:
        print("\n✗ Failed to load composition files", file=sys.stderr)
        sys.exit(1)

    print(f"  Short-read: {len(short_df)} genomes")
    print(f"  Long-read:  {len(long_df)} genomes")

    # Validate
    print("\nValidating compositions...")
    is_valid, issues = validate_compositions(short_df, long_df, tolerance=args.tolerance)

    if is_valid:
        print(f"\n{'='*80}")
        print("✓ VALIDATION PASSED")
        print(f"{'='*80}")
        print(f"\nGenome compositions match!")
        print(f"  • {len(short_df)} genomes in both datasets")
        print(f"  • Abundances match within tolerance ({args.tolerance})")
        print(f"\nThese datasets are suitable for hybrid assembly.")

        if args.verbose:
            print(f"\nGenome composition:")
            print(short_df[['genome_id', 'genome_type', 'relative_abundance']].to_string(index=False))

    else:
        print(f"\n{'='*80}")
        print("✗ VALIDATION FAILED")
        print(f"{'='*80}")
        print(f"\nIssues found:")
        for issue in issues:
            print(f"  • {issue}")

        print(f"\nThese datasets are NOT suitable for hybrid assembly.")
        print(f"Ensure both were generated with:")
        print(f"  • Same --collection-id")
        print(f"  • Same --seed")
        print(f"  • Same --vlp-protocol (or both --no-vlp)")

        sys.exit(1)

    # Try to load and compare metadata if available
    print(f"\n{'-'*80}")
    print("Metadata check:")
    print(f"{'-'*80}")

    short_meta_path = short_path.parent / short_path.name.replace('_composition.tsv', '_metadata.json')
    long_meta_path = long_path.parent / long_path.name.replace('_composition.tsv', '_metadata.json')

    if short_meta_path.exists() and long_meta_path.exists():
        short_meta = load_metadata(short_meta_path)
        long_meta = load_metadata(long_meta_path)

        if short_meta and long_meta:
            # Check random seed
            short_seed = short_meta.get('generation_info', {}).get('random_seed')
            long_seed = long_meta.get('generation_info', {}).get('random_seed')

            if short_seed == long_seed:
                print(f"✓ Random seed matches: {short_seed}")
            else:
                print(f"⚠ Random seed mismatch: short={short_seed}, long={long_seed}")

            # Check collection
            short_coll = short_meta.get('collection', {}).get('id')
            long_coll = long_meta.get('collection', {}).get('id')

            if short_coll == long_coll:
                print(f"✓ Collection ID matches: {short_coll}")
            else:
                print(f"⚠ Collection ID mismatch: short={short_coll}, long={long_coll}")

            # Check VLP protocol
            short_vlp = short_meta.get('configuration', {}).get('vlp_protocol', 'none')
            long_vlp = long_meta.get('configuration', {}).get('vlp_protocol', 'none')

            if short_vlp == long_vlp:
                print(f"✓ VLP protocol matches: {short_vlp}")
            else:
                print(f"⚠ VLP protocol mismatch: short={short_vlp}, long={long_vlp}")

    print()


if __name__ == '__main__':
    main()
