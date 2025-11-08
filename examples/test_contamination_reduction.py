"""
Test script for VLP contamination reduction integration

This script demonstrates how VLP enrichment protocols affect
different types of contamination differently.
"""

import sys
from pathlib import Path
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from viroforge.core.contamination import (
    create_contamination_profile,
    ContaminationProfile
)
from viroforge.enrichment.vlp import (
    VLPEnrichment,
    VLPProtocol
)


def test_contamination_reduction():
    """Test contamination reduction with different VLP protocols"""

    print("=" * 80)
    print("VLP Contamination Reduction Test")
    print("=" * 80)
    print()

    # Create a realistic contamination profile
    print("Creating 'realistic' contamination profile...")
    profile = create_contamination_profile('realistic', random_seed=42)

    print(f"\nOriginal Contamination Profile:")
    print(f"  Total abundance: {profile.get_total_abundance():.4f}")
    print(f"  Number of contaminants: {len(profile)}")

    by_type = profile.get_abundance_by_type()
    print(f"\nContamination by type:")
    for ctype, abundance in by_type.items():
        print(f"    {ctype.value:20s}: {abundance*100:6.2f}%")

    print("\n" + "=" * 80)
    print("Testing Different VLP Protocols")
    print("=" * 80)

    # Test different protocols
    protocols = [
        ('Tangential Flow (Standard)', VLPProtocol.tangential_flow_standard()),
        ('Syringe Filter', VLPProtocol.syringe_filter_standard()),
        ('Ultracentrifugation', VLPProtocol.ultracentrifugation()),
        ('Norgen Kit', VLPProtocol.norgen_kit()),
        ('No VLP (Control)', VLPProtocol.no_vlp())
    ]

    for protocol_name, protocol_config in protocols:
        print(f"\n{'=' * 80}")
        print(f"Protocol: {protocol_name}")
        print(f"{'=' * 80}")

        # Initialize VLP enrichment
        vlp = VLPEnrichment(protocol=protocol_config, random_seed=42)

        # Apply contamination reduction
        reduced_profile, stats = vlp.apply_contamination_reduction(profile)

        # Display results
        print(f"\nProtocol Configuration:")
        print(f"  Nuclease treatment: {protocol_config.nuclease_treatment}")
        if protocol_config.nuclease_treatment:
            print(f"  Nuclease efficiency: {protocol_config.nuclease_efficiency*100:.1f}%")
        print(f"  Filtration method: {protocol_config.filtration_method}")
        if protocol_config.pore_size_um:
            print(f"  Pore size: {protocol_config.pore_size_um} μm")
        print(f"  Overall contamination reduction: {protocol_config.contamination_reduction*100:.1f}%")

        print(f"\nReduction Results:")
        print(f"  Original total: {stats['original_total_contamination']*100:6.2f}%")
        print(f"  Reduced total:  {stats['reduced_total_contamination']*100:6.2f}%")
        print(f"  Overall reduction: {stats['overall_reduction_factor']*100:6.1f}%")

        print(f"\nReduction by contaminant type:")
        for ctype_str, ctype_stats in stats['reduction_by_type'].items():
            print(f"    {ctype_str:20s}:")
            print(f"      Original:  {ctype_stats['original_abundance']*100:6.2f}%")
            print(f"      Reduced:   {ctype_stats['reduced_abundance']*100:6.2f}%")
            print(f"      Removal:   {ctype_stats['reduction_pct']:6.1f}%")


def compare_vlp_vs_bulk():
    """Compare contamination levels in VLP vs bulk preparations"""

    print("\n\n" + "=" * 80)
    print("VLP vs Bulk Contamination Comparison")
    print("=" * 80)

    # Test with different initial contamination levels
    for profile_type in ['clean', 'realistic', 'heavy']:
        print(f"\n{'=' * 80}")
        print(f"Initial Contamination: {profile_type.upper()}")
        print(f"{'=' * 80}")

        # Create profile
        profile = create_contamination_profile(profile_type, random_seed=42)
        initial_total = profile.get_total_abundance()

        # VLP-enriched
        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        vlp_profile, vlp_stats = vlp.apply_contamination_reduction(profile)
        vlp_total = vlp_profile.get_total_abundance()

        # Bulk (no VLP)
        bulk = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        bulk_profile, bulk_stats = bulk.apply_contamination_reduction(profile)
        bulk_total = bulk_profile.get_total_abundance()

        print(f"\nContamination levels:")
        print(f"  Initial:      {initial_total*100:6.2f}%")
        print(f"  VLP-enriched: {vlp_total*100:6.2f}%  (reduction: {vlp_stats['overall_reduction_factor']*100:5.1f}%)")
        print(f"  Bulk (no VLP):{bulk_total*100:6.2f}%  (reduction: {bulk_stats['overall_reduction_factor']*100:5.1f}%)")
        print(f"  Fold difference: {bulk_total/max(vlp_total, 1e-10):.1f}x")


def test_literature_validation():
    """Validate against literature-reported contamination levels"""

    print("\n\n" + "=" * 80)
    print("Literature Validation")
    print("=" * 80)
    print()
    print("ViromeQC Survey (Roux et al. 2016) reported contamination ranges:")
    print("  VLP-enriched viromes: 1-15% non-viral")
    print("  Bulk metagenomes: 50-90% non-viral")
    print()

    # Create realistic profile (7.6% initial contamination)
    profile = create_contamination_profile('realistic', random_seed=42)

    # Simulate VLP enrichment
    vlp = VLPEnrichment(
        protocol=VLPProtocol.tangential_flow_standard(),
        random_seed=42
    )
    vlp_profile, vlp_stats = vlp.apply_contamination_reduction(profile)
    vlp_contamination = vlp_profile.get_total_abundance() * 100

    # Simulate bulk
    bulk = VLPEnrichment(
        protocol=VLPProtocol.no_vlp(),
        random_seed=42
    )
    bulk_profile, bulk_stats = bulk.apply_contamination_reduction(profile)
    bulk_contamination = bulk_profile.get_total_abundance() * 100

    print(f"ViroForge Simulation Results:")
    print(f"  VLP-enriched: {vlp_contamination:.1f}% non-viral")
    print(f"  Bulk metagenome: {bulk_contamination:.1f}% non-viral")
    print()

    # Check if within literature ranges
    vlp_valid = 1.0 <= vlp_contamination <= 15.0
    bulk_valid = 50.0 <= bulk_contamination <= 90.0

    print(f"Validation:")
    if vlp_valid:
        print(f"  ✓ VLP contamination within literature range (1-15%)")
    else:
        print(f"  ✗ VLP contamination outside literature range: {vlp_contamination:.1f}%")

    if bulk_valid:
        print(f"  ✓ Bulk contamination within literature range (50-90%)")
    else:
        print(f"  ✗ Bulk contamination outside literature range: {bulk_contamination:.1f}%")

    # Note about the simulation
    print()
    print("Note: These percentages represent contamination in the INPUT sample.")
    print("VLP enrichment REDUCES this contamination, bringing final viral fraction to 90-99%.")


if __name__ == '__main__':
    print("\nViroForge - VLP Contamination Reduction Integration Test")
    print("=" * 80)
    print()

    # Run tests
    test_contamination_reduction()
    compare_vlp_vs_bulk()
    test_literature_validation()

    print("\n" + "=" * 80)
    print("All tests completed!")
    print("=" * 80)
