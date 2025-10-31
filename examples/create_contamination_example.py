#!/usr/bin/env python3
"""
Example: Creating Contamination Profiles with ViroForge

This script demonstrates how to use the core.contamination module to create
realistic contamination profiles and combine them with viral communities to
create complete mock virome datasets.
"""

from viroforge.core import (
    create_body_site_profile,
    create_contamination_profile,
    ContaminationProfile,
    add_host_contamination,
    add_rrna_contamination,
    add_reagent_contamination,
    add_phix_control,
)
from viroforge.utils import create_mock_virome
from pathlib import Path


def example_1_predefined_profiles():
    """Example 1: Create pre-defined contamination profiles."""
    print("=" * 70)
    print("Example 1: Pre-defined Contamination Profiles")
    print("=" * 70)

    profile_types = ['clean', 'realistic', 'heavy', 'failed']

    for profile_type in profile_types:
        print(f"\n{profile_type.upper()} CONTAMINATION PROFILE:")
        print("-" * 40)

        profile = create_contamination_profile(profile_type, random_seed=42)

        stats = profile.get_summary_stats()
        print(f"  Total contaminants: {stats['n_contaminants']}")
        print(f"  Total contamination: {stats['total_abundance']*100:.2f}%")
        print(f"\n  Breakdown by type:")
        for ctype, abundance in stats['by_type'].items():
            print(f"    {ctype}: {abundance*100:.2f}%")


def example_2_custom_profile():
    """Example 2: Build a custom contamination profile."""
    print("\n" + "=" * 70)
    print("Example 2: Custom Contamination Profile")
    print("=" * 70)

    # Create empty profile
    profile = ContaminationProfile(name="my_custom_profile")

    # Add specific contaminants with custom levels
    print("\nAdding contaminants:")
    print("  - 8% human host DNA")
    print("  - 10% rRNA")
    print("  - 1% reagent bacteria")
    print("  - 0.5% PhiX")

    add_host_contamination(
        profile,
        host_organism="human",
        abundance_pct=8.0,
        n_fragments=100,
        random_seed=42
    )

    add_rrna_contamination(
        profile,
        abundance_pct=10.0,
        n_sequences=50,
        random_seed=42
    )

    add_reagent_contamination(
        profile,
        abundance_pct=1.0,
        n_genomes=5,
        random_seed=42
    )

    add_phix_control(
        profile,
        abundance_pct=0.5
    )

    print(f"\nCreated: {profile}")
    print(f"Summary: {profile.get_summary_stats()}")


def example_3_host_organisms():
    """Example 3: Different host organisms."""
    print("\n" + "=" * 70)
    print("Example 3: Host DNA from Different Organisms")
    print("=" * 70)

    organisms = ['human', 'mouse', 'rat']

    for organism in organisms:
        profile = ContaminationProfile(name=f"{organism}_host")

        add_host_contamination(
            profile,
            host_organism=organism,
            abundance_pct=5.0,
            n_fragments=50,
            random_seed=42
        )

        print(f"\n{organism.upper()} host DNA:")
        table = profile.get_contamination_table()
        print(f"  Sequences: {len(table)}")
        print(f"  Total abundance: {profile.get_total_abundance()*100:.2f}%")
        print(f"  Mean GC content: {table['gc_content'].mean():.1f}%")
        print(f"  Organism: {table['organism'].iloc[0]}")


def example_4_complete_mock_virome():
    """Example 4: Create a complete mock virome with viral + contamination."""
    print("\n" + "=" * 70)
    print("Example 4: Complete Mock Virome Composition")
    print("=" * 70)

    # Create a gut virome with realistic contamination
    mock_virome = create_mock_virome(
        name="gut_virome_realistic",
        body_site="gut",
        contamination_level="realistic",
        n_viral_genomes=50,
        viral_fraction=0.90,  # 90% viral, 10% contamination
        random_seed=42
    )

    print(f"\n{mock_virome}")

    print("\nComposition Summary:")
    stats = mock_virome.get_summary_stats()
    print(f"  Viral genomes: {stats['n_viral_genomes']}")
    print(f"  Viral abundance: {stats['viral_abundance']*100:.2f}%")
    print(f"  Contaminants: {stats['n_contaminants']}")
    print(f"  Contamination abundance: {stats['contamination_abundance']*100:.2f}%")

    print("\n  Contamination breakdown:")
    for ctype, abundance in stats['contamination_by_type'].items():
        print(f"    {ctype}: {abundance*100:.2f}%")

    # Get composition table
    comp_table = mock_virome.get_composition_table()

    print(f"\n  Total sequences: {len(comp_table)}")
    print("\n  Top 10 most abundant sequences:")
    top_10 = comp_table.nlargest(10, 'abundance')[
        ['genome_id', 'type', 'organism', 'abundance']
    ]
    for idx, row in top_10.iterrows():
        print(f"    {row['genome_id'][:30]:30s} | {row['type']:15s} | "
              f"{row['abundance']*100:6.2f}%")


def example_5_vlp_comparison():
    """Example 5: Compare VLP success vs failure."""
    print("\n" + "=" * 70)
    print("Example 5: VLP Enrichment Success vs Failure")
    print("=" * 70)

    scenarios = {
        'Successful VLP': 'clean',
        'Typical VLP': 'realistic',
        'Poor VLP': 'heavy',
        'Failed VLP': 'failed'
    }

    print("\nComparing different VLP enrichment outcomes:\n")
    print(f"{'Scenario':<20} {'Viral %':<10} {'Host %':<10} {'rRNA %':<10} "
          f"{'Reagent %':<10} {'PhiX %':<10}")
    print("-" * 70)

    for scenario_name, profile_type in scenarios.items():
        mock_virome = create_mock_virome(
            name=scenario_name,
            body_site="gut",
            contamination_level=profile_type,
            n_viral_genomes=30,
            viral_fraction=0.95 if profile_type == 'clean' else
                          0.90 if profile_type == 'realistic' else
                          0.75 if profile_type == 'heavy' else 0.60,
            random_seed=42
        )

        stats = mock_virome.get_summary_stats()
        by_type = stats['contamination_by_type']

        print(f"{scenario_name:<20} "
              f"{stats['viral_abundance']*100:>9.1f} "
              f"{by_type.get('host_dna', 0)*100:>9.2f} "
              f"{by_type.get('rrna', 0)*100:>9.2f} "
              f"{by_type.get('reagent_bacteria', 0)*100:>9.2f} "
              f"{by_type.get('phix', 0)*100:>9.2f}")


def example_6_export_data():
    """Example 6: Export contamination data to files."""
    print("\n" + "=" * 70)
    print("Example 6: Exporting Contamination Data")
    print("=" * 70)

    # Create a contamination profile
    profile = create_contamination_profile('realistic', random_seed=42)

    # Create output directory
    output_dir = Path("example_output")
    output_dir.mkdir(exist_ok=True)

    # Export contamination table
    contam_table_file = output_dir / "contamination_table.tsv"
    profile.export_contamination_table(contam_table_file)
    print(f"\n  Exported contamination table to: {contam_table_file}")

    # Export contaminant sequences (FASTA)
    contam_fasta_file = output_dir / "contaminants.fasta"
    profile.export_contaminants_fasta(contam_fasta_file)
    print(f"  Exported contaminant sequences to: {contam_fasta_file}")

    # Create complete mock virome and export
    mock_virome = create_mock_virome(
        name="export_example",
        body_site="gut",
        contamination_level="realistic",
        n_viral_genomes=30,
        random_seed=42
    )

    # Export complete composition
    composition_file = output_dir / "complete_composition.tsv"
    mock_virome.export_composition_table(composition_file)
    print(f"  Exported complete composition to: {composition_file}")

    print(f"\n  All files created in: {output_dir.absolute()}")


if __name__ == "__main__":
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 13 + "ViroForge Contamination Module Examples" + " " * 15 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    # Run examples
    example_1_predefined_profiles()
    example_2_custom_profile()
    example_3_host_organisms()
    example_4_complete_mock_virome()
    example_5_vlp_comparison()
    example_6_export_data()

    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
    print()
