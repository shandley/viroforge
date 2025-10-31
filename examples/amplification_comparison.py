#!/usr/bin/env python3
"""
Amplification Bias Comparison Example
======================================

This example demonstrates how different library preparation amplification methods
introduce different biases into virome sequencing data.

Comparison of:
1. RdAB (Random RT + dsDNA + PCR) - Most common, length + GC bias
2. MDA (Multiple Displacement Amplification) - Extreme GC bias, high stochasticity
3. Linker (Adapter ligation + PCR) - Minimal bias, modern protocols
4. No Amplification - Control (high-biomass samples)

Author: ViroForge Development Team
"""

from viroforge.core import create_body_site_profile
from viroforge.amplification import (
    rdab_40_cycles,
    mda_standard,
    linker_standard,
    no_amplification
)


def main():
    print("="*70)
    print("ViroForge: Amplification Bias Comparison")
    print("="*70)
    print()
    print("This example shows how different amplification methods introduce")
    print("different biases into virome composition.")
    print()

    # Create a test viral community
    print("Creating gut virome community with diverse genomes...")
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=20,
        random_seed=42
    )

    print(f"  ✓ Created community with {len(viral_community.genomes)} genomes")
    print()

    # Show initial composition
    print("Initial Composition (before amplification):")
    print("-"*70)
    print(f"{'Genome ID':<25} {'Length':>10} {'GC%':>8} {'Abundance':>12}")
    print("-"*70)

    # Sort by abundance for display
    sorted_genomes = sorted(viral_community.genomes,
                          key=lambda g: g.abundance,
                          reverse=True)[:5]

    for genome in sorted_genomes:
        print(f"{genome.genome_id:<25} {genome.length:>10,} "
              f"{genome.gc_content*100:>7.1f}% {genome.abundance:>12.6f}")
    print(f"{'... (showing top 5)':<25}")
    print()

    # Test each amplification method
    methods = {
        'No Amplification (Control)': no_amplification(),
        'RdAB (40 cycles)': rdab_40_cycles(),
        'Linker (20 cycles)': linker_standard(),
        'MDA (4 hours)': mda_standard()
    }

    print("="*70)
    print("AMPLIFICATION METHOD COMPARISON")
    print("="*70)
    print()

    results = {}

    for method_name, method in methods.items():
        print(f"{method_name}:")
        print("-"*70)

        # Create fresh community for each method
        test_community = create_body_site_profile(
            body_site='gut',
            n_genomes=20,
            random_seed=42
        )

        # Import MockViromeComposition
        from viroforge.utils.composition import MockViromeComposition

        # Create composition
        composition = MockViromeComposition(
            name=method_name,
            viral_community=test_community,
            contamination_profile=None,
            viral_fraction=1.0
        )

        # Apply amplification
        method.apply(composition)

        # Show top genomes after amplification
        sorted_genomes = sorted(composition.viral_community.genomes,
                              key=lambda g: g.abundance,
                              reverse=True)[:5]

        print(f"{'Genome ID':<25} {'Length':>10} {'GC%':>8} {'Final Abund':>13}")
        for genome in sorted_genomes:
            print(f"{genome.genome_id:<25} {genome.length:>10,} "
                  f"{genome.gc_content*100:>7.1f}% {genome.abundance:>13.6f}")

        # Calculate statistics
        import numpy as np
        abundances = [g.abundance for g in composition.viral_community.genomes]

        stats = {
            'mean': np.mean(abundances),
            'std': np.std(abundances),
            'cv': np.std(abundances) / np.mean(abundances) if np.mean(abundances) > 0 else 0,
            'max/min': max(abundances) / min(abundances) if min(abundances) > 0 else float('inf')
        }

        results[method_name] = stats

        print(f"\nStatistics:")
        print(f"  Mean abundance:        {stats['mean']:.6f}")
        print(f"  Std deviation:         {stats['std']:.6f}")
        print(f"  Coefficient variation: {stats['cv']:.2f}")
        print(f"  Max/Min ratio:         {stats['max/min']:.1f}x")
        print()

    # Summary comparison
    print("="*70)
    print("SUMMARY: BIAS COMPARISON")
    print("="*70)
    print()

    print("Coefficient of Variation (higher = more bias):")
    print("-"*70)
    for method_name, stats in results.items():
        bar_length = int(stats['cv'] * 20)  # Scale for visualization
        bar = '█' * bar_length
        print(f"{method_name:<30} {stats['cv']:>6.2f} {bar}")
    print()

    print("Max/Min Abundance Ratio (higher = stronger enrichment/depletion):")
    print("-"*70)
    for method_name, stats in results.items():
        if stats['max/min'] == float('inf'):
            print(f"{method_name:<30}   ∞")
        else:
            bar_length = int(min(stats['max/min'] / 10, 50))
            bar = '█' * bar_length
            print(f"{method_name:<30} {stats['max/min']:>6.1f}x {bar}")
    print()

    # Guidance
    print("="*70)
    print("INTERPRETATION & GUIDANCE")
    print("="*70)
    print()

    print("No Amplification (Control):")
    print("  • Minimal bias, preserves true composition")
    print("  • Only possible with high-biomass samples (rare)")
    print("  • Use when: >10 ng viral DNA available")
    print()

    print("RdAB Amplification:")
    print("  • Most common method in virome studies")
    print("  • Moderate bias: favors short genomes and optimal GC")
    print("  • Use when: Standard virome study, 0.1-10 ng DNA")
    print()

    print("Linker Amplification:")
    print("  • Modern library prep (Nextera, TruSeq)")
    print("  • Reduced bias compared to RdAB")
    print("  • Use when: Want less bias, have modern kit")
    print()

    print("MDA Amplification:")
    print("  • For ultra-low biomass samples")
    print("  • Strong GC bias and high stochasticity")
    print("  • Use when: <0.1 ng DNA, single-cell viruses")
    print()

    print("="*70)
    print("✓ Comparison complete!")
    print("="*70)


if __name__ == '__main__':
    main()
