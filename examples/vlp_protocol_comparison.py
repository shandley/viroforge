#!/usr/bin/env python3
"""
VLP Enrichment Protocol Comparison
===================================

This example compares different VLP enrichment protocols to show how
methodology affects the final viral composition. Different labs use
different protocols, and this impacts the resulting data.

Protocols compared:
1. Standard VLP (most common)
2. Iron chloride VLP (Conceição-Neto et al.)
3. Ultracentrifugation-based
4. Syringe filter (sharp cutoff)
"""

from copy import deepcopy
from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import (
    standard_vlp,
    iron_chloride_vlp,
    ultracentrifuge_vlp,
    syringe_filter_vlp,
    no_enrichment
)


def create_test_composition():
    """Create a consistent test composition."""
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=20,
        random_seed=42
    )

    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )

    return MockViromeComposition(
        name='test_sample',
        viral_community=viral_community,
        contamination_profile=contamination,
        viral_fraction=0.50
    )


def main():
    print("="*70)
    print("ViroForge: VLP Protocol Comparison")
    print("="*70)
    print()

    # Create base composition for each protocol
    print("Creating identical starting samples...")
    base_composition = create_test_composition()
    initial_viral = base_composition.viral_fraction

    protocols = {
        'Bulk Metagenome (No VLP)': no_enrichment(),
        'Standard VLP': standard_vlp(),
        'Iron Chloride VLP': iron_chloride_vlp(),
        'Ultracentrifugation VLP': ultracentrifuge_vlp(),
        'Syringe Filter VLP': syringe_filter_vlp(),
    }

    results = {}

    print(f"  ✓ Created base composition: {initial_viral:.1%} viral\n")

    # Apply each protocol
    print("Applying VLP enrichment protocols...\n")
    print("-" * 70)

    for protocol_name, protocol in protocols.items():
        # Create fresh composition for each protocol
        composition = create_test_composition()

        # Apply protocol
        protocol.apply(composition)

        # Store results
        results[protocol_name] = {
            'viral_fraction': composition.viral_fraction,
            'contam_fraction': composition.contamination_fraction,
            'enrichment': composition.viral_fraction / initial_viral,
            'composition': composition
        }

        print(f"{protocol_name}:")
        print(f"  Final viral fraction: {composition.viral_fraction:.1%}")
        print(f"  Enrichment:           {results[protocol_name]['enrichment']:.2f}x")
        print()

    # Detailed comparison
    print("="*70)
    print("DETAILED PROTOCOL COMPARISON")
    print("="*70)
    print()

    # Table header
    print(f"{'Protocol':<30} {'Viral %':>10} {'Enrichment':>12} {'Contam %':>10}")
    print("-" * 70)

    for protocol_name, result in results.items():
        print(f"{protocol_name:<30} "
              f"{result['viral_fraction']:>9.1%} "
              f"{result['enrichment']:>11.2f}x "
              f"{result['contam_fraction']:>9.1%}")

    # Protocol descriptions
    print("\n" + "="*70)
    print("PROTOCOL DETAILS")
    print("="*70)

    print("\n1. Bulk Metagenome (No VLP)")
    print("   - No filtration")
    print("   - No nuclease treatment")
    print("   - Preserves original composition")
    print("   - Use case: Virus-host ratio studies")

    print("\n2. Standard VLP (Most Common)")
    print("   - 0.45 μm pre-filtration")
    print("   - 0.2 μm tangential flow filtration")
    print("   - 95% nuclease efficiency")
    print("   - Use case: General virome studies")

    print("\n3. Iron Chloride VLP (Conceição-Neto et al.)")
    print("   - FeCl3 precipitation for virus aggregation")
    print("   - 0.2 μm filtration")
    print("   - 98% nuclease efficiency (enhanced)")
    print("   - Use case: Maximum contamination removal")

    print("\n4. Ultracentrifugation VLP")
    print("   - Density gradient ultracentrifugation")
    print("   - 0.2 μm filtration")
    print("   - 90% nuclease efficiency")
    print("   - Use case: High-quality viral DNA")

    print("\n5. Syringe Filter VLP")
    print("   - 0.2 μm syringe filter (sharp cutoff)")
    print("   - 93% nuclease efficiency")
    print("   - Use case: Field studies, low sample volume")

    # Recommendations
    print("\n" + "="*70)
    print("PROTOCOL SELECTION GUIDE")
    print("="*70)

    print("\nChoose based on your research goal:")
    print()
    print("  Maximum viral purity:")
    print("    → Iron Chloride VLP or Standard VLP")
    print()
    print("  Standardized/reproducible:")
    print("    → Standard VLP (most common in literature)")
    print()
    print("  Low sample volume:")
    print("    → Syringe Filter VLP")
    print()
    print("  Virus-host interactions:")
    print("    → Bulk Metagenome (no enrichment)")
    print()
    print("  Highest quality DNA:")
    print("    → Ultracentrifugation VLP")

    # Show which protocol achieved best enrichment
    best_protocol = max(results.items(),
                       key=lambda x: x[1]['viral_fraction'] if x[0] != 'Bulk Metagenome (No VLP)' else 0)

    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    print(f"\nBest viral enrichment: {best_protocol[0]}")
    print(f"  Achieved: {best_protocol[1]['viral_fraction']:.1%} viral")
    print(f"  Enrichment: {best_protocol[1]['enrichment']:.2f}x over bulk")

    print("\n" + "="*70)
    print("✓ Protocol comparison complete!")
    print("="*70)


if __name__ == '__main__':
    main()
