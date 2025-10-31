#!/usr/bin/env python3
"""
Basic VLP Enrichment Example
=============================

This example demonstrates how to apply VLP (Virus-Like Particle) enrichment
to a mock virome composition, simulating the common laboratory process of
enriching viral DNA while removing host and bacterial contamination.

VLP enrichment is the defining feature of viromics that differentiates it
from bulk metagenomics.
"""

from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp


def main():
    print("="*70)
    print("ViroForge: Basic VLP Enrichment Example")
    print("="*70)
    print()

    # Step 1: Create a viral community
    print("Step 1: Creating gut virome community...")
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=20,
        random_seed=42
    )
    print(f"  ✓ Created community with {len(viral_community.genomes)} viral genomes")
    print()

    # Step 2: Create contamination profile
    print("Step 2: Creating contamination profile (realistic level)...")
    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )
    print(f"  ✓ Created contamination with {len(contamination.contaminants)} contaminants")
    print()

    # Step 3: Create bulk metagenome composition (50% viral, 50% contamination)
    print("Step 3: Creating bulk metagenome composition...")
    composition = MockViromeComposition(
        name='gut_sample',
        viral_community=viral_community,
        contamination_profile=contamination,
        viral_fraction=0.50  # 50% viral (typical bulk metagenome)
    )

    initial_viral_fraction = composition.viral_fraction
    initial_contam_fraction = composition.contamination_fraction

    print(f"  Composition before enrichment:")
    print(f"    - Viral fraction:         {initial_viral_fraction:.1%}")
    print(f"    - Contamination fraction: {initial_contam_fraction:.1%}")
    print()

    # Step 4: Apply VLP enrichment
    print("Step 4: Applying VLP enrichment...")
    print("  Standard VLP protocol:")
    print("    - 0.45 μm pre-filtration")
    print("    - 0.2 μm tangential flow filtration")
    print("    - DNase/RNase treatment (95% efficiency)")
    print()

    vlp = standard_vlp()
    vlp.apply(composition)

    final_viral_fraction = composition.viral_fraction
    final_contam_fraction = composition.contamination_fraction

    print(f"  Composition after enrichment:")
    print(f"    - Viral fraction:         {final_viral_fraction:.1%}")
    print(f"    - Contamination fraction: {final_contam_fraction:.1%}")
    print()

    # Step 5: Calculate enrichment metrics
    print("Step 5: Enrichment Results")
    print("="*70)
    enrichment_fold = final_viral_fraction / initial_viral_fraction
    contam_reduction = initial_contam_fraction / final_contam_fraction

    print(f"  Viral enrichment:       {enrichment_fold:.2f}x increase")
    print(f"  Contamination removal:  {contam_reduction:.1f}x reduction")
    print(f"  Final viral purity:     {final_viral_fraction:.1%}")
    print()

    # Step 6: Show what this means for sequencing
    print("Step 6: Impact on Sequencing")
    print("="*70)
    print("  If you sequence 10 million reads:")

    bulk_viral_reads = 10_000_000 * initial_viral_fraction
    vlp_viral_reads = 10_000_000 * final_viral_fraction

    print(f"    Bulk metagenome:  {bulk_viral_reads:,.0f} viral reads ({initial_viral_fraction:.1%})")
    print(f"    VLP enriched:     {vlp_viral_reads:,.0f} viral reads ({final_viral_fraction:.1%})")
    print(f"    Gain:             {vlp_viral_reads - bulk_viral_reads:,.0f} additional viral reads")
    print()

    print("  This means:")
    print(f"    - {enrichment_fold:.1f}x more viral sequence data")
    print(f"    - {enrichment_fold:.1f}x better viral genome coverage")
    print(f"    - {enrichment_fold:.1f}x more power to detect rare viruses")
    print()

    print("="*70)
    print("✓ VLP enrichment complete!")
    print("="*70)


if __name__ == '__main__':
    main()
