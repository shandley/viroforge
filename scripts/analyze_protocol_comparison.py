#!/usr/bin/env python3
"""
Analyze VLP protocol comparison results

Extracts statistics from generated datasets and compares against literature.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List


def extract_protocol_stats(dataset_dir: Path) -> Dict:
    """Extract stats from a single dataset"""
    metadata_dir = dataset_dir / 'metadata'
    metadata_files = list(metadata_dir.glob('*_metadata.json'))

    if not metadata_files:
        return None

    with open(metadata_files[0]) as f:
        metadata = json.load(f)

    enrichment_stats = metadata['enrichment_stats']

    # Extract key metrics
    stats = {
        'protocol': enrichment_stats['vlp_protocol'],
        'viral_fraction': enrichment_stats['viral_fraction'] * 100,  # As percentage
        'contamination_fraction': enrichment_stats['contamination_fraction'] * 100,
        'n_viral_genomes': enrichment_stats['n_viral_genomes'],
        'n_contaminants': enrichment_stats['n_contaminants']
    }

    # Add viral enrichment stats if present
    if 'viral_enrichment' in enrichment_stats and enrichment_stats['viral_enrichment']:
        viral_enrich = enrichment_stats['viral_enrichment']
        stats['mean_enrichment_factor'] = viral_enrich.get('mean_enrichment_factor')
        stats['recovery_rate'] = viral_enrich.get('viral_recovery_rate')

        if 'size_bias' in viral_enrich:
            stats['size_correlation'] = viral_enrich['size_bias'].get('size_enrichment_correlation')

    # Add contamination reduction stats
    if 'contamination_reduction' in enrichment_stats and enrichment_stats['contamination_reduction']:
        contam_red = enrichment_stats['contamination_reduction']
        stats['original_contamination'] = contam_red.get('original_total_contamination', 0) * 100
        stats['reduction_factor'] = contam_red.get('overall_reduction_factor', 0) * 100

        # Extract by-type reductions
        if 'reduction_by_type' in contam_red:
            by_type = contam_red['reduction_by_type']
            stats['host_dna_reduction'] = by_type.get('host_dna', {}).get('reduction_pct', 0)
            stats['rrna_reduction'] = by_type.get('rrna', {}).get('reduction_pct', 0)
            stats['bacteria_reduction'] = by_type.get('reagent_bacteria', {}).get('reduction_pct', 0)

    return stats


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_protocol_comparison.py <comparison_dir>")
        sys.exit(1)

    comparison_dir = Path(sys.argv[1])

    if not comparison_dir.exists():
        print(f"Error: Directory not found: {comparison_dir}")
        sys.exit(1)

    # Find all dataset directories
    datasets = sorted([d for d in comparison_dir.iterdir() if d.is_dir() and d.name.startswith('collection_')])

    if not datasets:
        print(f"No datasets found in {comparison_dir}")
        sys.exit(1)

    print("=" * 100)
    print("VLP Protocol Comparison Analysis")
    print("=" * 100)
    print()

    all_stats = []
    for dataset_dir in datasets:
        stats = extract_protocol_stats(dataset_dir)
        if stats:
            all_stats.append(stats)

    if not all_stats:
        print("No valid metadata found")
        sys.exit(1)

    # Print summary table
    print("SUMMARY TABLE")
    print("-" * 100)
    print(f"{'Protocol':<25} {'Viral %':<12} {'Contam %':<12} {'Reduction %':<15} {'Recovery':<12}")
    print("-" * 100)

    for stats in all_stats:
        protocol = stats['protocol'] if stats['protocol'] != 'none' else 'Bulk (no VLP)'
        viral = stats['viral_fraction']
        contam = stats['contamination_fraction']
        reduction = stats.get('reduction_factor', 0)
        recovery = stats.get('recovery_rate', 1.0)

        print(f"{protocol:<25} {viral:>10.2f}% {contam:>10.2f}% {reduction:>13.1f}% {recovery:>10.1f}%")

    print("-" * 100)
    print()

    # Print detailed contamination reduction
    print("CONTAMINATION REDUCTION BY TYPE")
    print("-" * 100)
    print(f"{'Protocol':<25} {'Host DNA %':<15} {'rRNA %':<15} {'Bacteria %':<15}")
    print("-" * 100)

    for stats in all_stats:
        if stats['protocol'] == 'none':
            continue  # Skip bulk

        protocol = stats['protocol']
        host = stats.get('host_dna_reduction', 0)
        rrna = stats.get('rrna_reduction', 0)
        bacteria = stats.get('bacteria_reduction', 0)

        print(f"{protocol:<25} {host:>13.1f}% {rrna:>13.1f}% {bacteria:>13.1f}%")

    print("-" * 100)
    print()

    # Literature comparison
    print("LITERATURE VALIDATION")
    print("-" * 100)
    print()

    print("Expected Ranges (from literature):")
    print("  VLP-enriched viral fraction: >90% (Roux et al. 2016)")
    print("  Bulk metagenome viral fraction: 10-50% (Roux et al. 2016)")
    print("  Nuclease efficiency (host DNA): >85% (Thurber et al. 2009)")
    print("  Bacterial filtration (0.2 μm): >90% (Shkoporov et al. 2018)")
    print()

    print("ViroForge Results:")
    for stats in all_stats:
        protocol = stats['protocol'] if stats['protocol'] != 'none' else 'Bulk'
        viral = stats['viral_fraction']

        if stats['protocol'] == 'none':
            status = "✓" if 10 <= viral <= 50 else "✗"
            print(f"  {protocol:<20}: {viral:>6.2f}% viral  {status}")
        else:
            status = "✓" if viral > 90 else "✗"
            print(f"  {protocol:<20}: {viral:>6.2f}% viral  {status}")

            if 'host_dna_reduction' in stats:
                host_red = stats['host_dna_reduction']
                host_status = "✓" if host_red > 85 else "✗"
                print(f"    - Host DNA reduction: {host_red:>5.1f}%  {host_status}")

            if 'bacteria_reduction' in stats:
                bact_red = stats['bacteria_reduction']
                bact_status = "✓" if bact_red > 90 else "✗"
                print(f"    - Bacteria reduction:  {bact_red:>5.1f}%  {bact_status}")

    print()
    print("=" * 100)


if __name__ == '__main__':
    main()
