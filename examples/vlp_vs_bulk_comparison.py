#!/usr/bin/env python3
"""
VLP vs Bulk Metagenome Comparison
==================================

This example demonstrates the difference between VLP-enriched virome
sequencing and bulk metagenomic sequencing by creating two parallel
samples from the same starting material.

Key Insight: VLP enrichment dramatically changes the composition,
while bulk metagenomics preserves the original community structure.
"""

from copy import deepcopy
from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp, no_enrichment


def print_composition_stats(name, composition):
    """Print detailed composition statistics."""
    stats = composition.get_summary_stats()

    print(f"\n{name}:")
    print(f"  Viral genomes:      {stats['n_viral_genomes']}")
    print(f"  Contaminants:       {stats['n_contaminants']}")
    print(f"  Viral fraction:     {stats['viral_fraction']:.1%}")
    print(f"  Contam. fraction:   {stats['contamination_fraction']:.1%}")

    if 'contamination_by_type' in stats:
        print(f"  Contamination breakdown:")
        for ctype, count in stats['contamination_by_type'].items():
            print(f"    - {ctype}: {count}")


def main():
    print("="*70)
    print("ViroForge: VLP vs Bulk Metagenome Comparison")
    print("="*70)
    print()
    print("Simulating a paired sequencing experiment:")
    print("  Sample A: VLP-enriched virome")
    print("  Sample B: Bulk metagenome (no enrichment)")
    print()

    # Create starting material
    print("Creating starting material...")
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=25,
        random_seed=42
    )

    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )

    # Create two identical compositions (one for VLP, one for bulk)
    vlp_composition = MockViromeComposition(
        name='sample_A_vlp',
        viral_community=deepcopy(viral_community),
        contamination_profile=deepcopy(contamination),
        viral_fraction=0.50  # 50% viral starting material
    )

    bulk_composition = MockViromeComposition(
        name='sample_B_bulk',
        viral_community=deepcopy(viral_community),
        contamination_profile=deepcopy(contamination),
        viral_fraction=0.50  # Identical starting material
    )

    print(f"  ✓ Created paired samples from identical starting material")
    print(f"  Starting viral fraction: {vlp_composition.viral_fraction:.1%}")
    print()

    # Process Sample A: VLP enrichment
    print("Processing Sample A: VLP Enrichment")
    print("-" * 70)
    vlp_protocol = standard_vlp()
    vlp_protocol.apply(vlp_composition)
    print(f"  ✓ Applied standard VLP enrichment protocol")

    # Process Sample B: No enrichment
    print("\nProcessing Sample B: Bulk Metagenome (No Enrichment)")
    print("-" * 70)
    bulk_protocol = no_enrichment()
    bulk_protocol.apply(bulk_composition)
    print(f"  ✓ No enrichment applied (preserves original composition)")

    # Compare results
    print("\n" + "="*70)
    print("RESULTS: Side-by-Side Comparison")
    print("="*70)

    print_composition_stats("Sample A: VLP-Enriched Virome", vlp_composition)
    print_composition_stats("Sample B: Bulk Metagenome", bulk_composition)

    # Calculate differences
    print("\n" + "="*70)
    print("ENRICHMENT IMPACT")
    print("="*70)

    vlp_viral = vlp_composition.viral_fraction
    bulk_viral = bulk_composition.viral_fraction
    enrichment = vlp_viral / bulk_viral

    print(f"\nViral Enrichment: {enrichment:.2f}x")
    print(f"  VLP sample:  {vlp_viral:.1%} viral")
    print(f"  Bulk sample: {bulk_viral:.1%} viral")
    print(f"  Difference:  {(vlp_viral - bulk_viral)*100:.1f} percentage points")

    # Sequencing implications
    print("\n" + "="*70)
    print("SEQUENCING IMPLICATIONS")
    print("="*70)

    total_reads = 10_000_000
    vlp_viral_reads = total_reads * vlp_viral
    bulk_viral_reads = total_reads * bulk_viral

    print(f"\nFor {total_reads:,} sequencing reads:")
    print(f"  VLP enriched:   {vlp_viral_reads:,.0f} viral reads")
    print(f"  Bulk metag:     {bulk_viral_reads:,.0f} viral reads")
    print(f"  Gain:           {vlp_viral_reads - bulk_viral_reads:,.0f} additional viral reads")

    # Cost efficiency
    print("\n" + "="*70)
    print("COST EFFICIENCY")
    print("="*70)

    cost_per_million_reads = 50  # Example cost

    # To get same number of viral reads
    bulk_reads_needed = vlp_viral_reads / bulk_viral
    bulk_cost = (bulk_reads_needed / 1_000_000) * cost_per_million_reads
    vlp_cost = (total_reads / 1_000_000) * cost_per_million_reads

    print(f"\nTo obtain {vlp_viral_reads:,.0f} viral reads:")
    print(f"  VLP approach:   {total_reads:,} total reads (${vlp_cost:.2f})")
    print(f"  Bulk approach:  {bulk_reads_needed:,.0f} total reads (${bulk_cost:.2f})")
    print(f"  Cost savings:   ${bulk_cost - vlp_cost:.2f} ({((bulk_cost - vlp_cost)/bulk_cost)*100:.0f}% cheaper)")

    # Research questions
    print("\n" + "="*70)
    print("RESEARCH APPLICATIONS")
    print("="*70)

    print("\nBest use cases:")
    print("  VLP-enriched virome:")
    print("    ✓ Viral diversity studies")
    print("    ✓ Viral genome assembly")
    print("    ✓ Rare virus detection")
    print("    ✓ Viral evolution analysis")
    print()
    print("  Bulk metagenome:")
    print("    ✓ Virus-host interaction ratios")
    print("    ✓ Community structure (bacteria + viruses)")
    print("    ✓ Metabolic pathway analysis")
    print("    ✓ Multi-kingdom studies")

    print("\n" + "="*70)
    print("✓ Comparison complete!")
    print("="*70)


if __name__ == '__main__':
    main()
