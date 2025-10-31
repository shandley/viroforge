#!/usr/bin/env python3
"""
Example: Creating Viral Communities with ViroForge

This script demonstrates how to use the core.community module to create
realistic viral community compositions with different abundance distributions
and body-site specific profiles.
"""

from viroforge.core import (
    ViralGenome,
    ViralCommunity,
    create_abundance_profile,
    create_body_site_profile,
)
import numpy as np


def example_1_create_custom_community():
    """Example 1: Create a custom viral community from scratch."""
    print("=" * 70)
    print("Example 1: Creating a Custom Viral Community")
    print("=" * 70)

    # Create a new community
    community = ViralCommunity(name="my_custom_virome")

    # Add individual genomes
    for i in range(20):
        genome = ViralGenome(
            genome_id=f"viral_genome_{i:03d}",
            sequence="ATCGATCGATCG" * 500,  # 6kb genome
            taxonomy=f"Siphoviridae;Genus_{i % 5};Species_{i}",
            family="Siphoviridae",
            genus=f"Genus_{i % 5}",
            species=f"Species_{i}",
            description=f"Example viral genome {i}"
        )
        community.add_genome(genome)

    # Apply log-normal abundance distribution
    print(f"\nCreated community with {len(community)} genomes")
    community.apply_abundance_distribution('lognormal', sigma=2.0)

    # Show summary
    print(f"\n{community}")
    print("\nSummary Statistics:")
    for key, value in community.get_summary_stats().items():
        print(f"  {key}: {value}")

    # Show top 5 most abundant genomes
    print("\nTop 5 Most Abundant Genomes:")
    abundance_df = community.get_abundance_table()
    top_5 = abundance_df.nlargest(5, 'abundance')[['genome_id', 'taxonomy', 'abundance']]
    print(top_5.to_string(index=False))

    return community


def example_2_body_site_profiles():
    """Example 2: Create body-site specific communities."""
    print("\n" + "=" * 70)
    print("Example 2: Creating Body-Site Specific Communities")
    print("=" * 70)

    body_sites = ['gut', 'oral', 'skin']

    for body_site in body_sites:
        print(f"\n{body_site.upper()} VIROME:")
        print("-" * 40)

        # Create body-site specific community
        community = create_body_site_profile(
            body_site=body_site,
            n_genomes=30,
            random_seed=42
        )

        # Show family distribution
        abundance_df = community.get_abundance_table()
        family_abundance = abundance_df.groupby('family')['abundance'].sum().sort_values(ascending=False)

        print(f"  Total genomes: {len(community)}")
        print(f"  Total abundance: {community.get_total_abundance():.4f}")
        print(f"\n  Family Distribution:")
        for family, abundance in family_abundance.head(5).items():
            print(f"    {family}: {abundance:.2%}")


def example_3_abundance_distributions():
    """Example 3: Compare different abundance distributions."""
    print("\n" + "=" * 70)
    print("Example 3: Comparing Different Abundance Distributions")
    print("=" * 70)

    distributions = {
        'lognormal': {'sigma': 2.0},
        'powerlaw': {'alpha': 1.5},
        'even': {}
    }

    for dist_name, params in distributions.items():
        print(f"\n{dist_name.upper()} DISTRIBUTION:")
        print("-" * 40)

        # Create community
        community = ViralCommunity(name=f"{dist_name}_community")

        # Add 50 genomes
        for i in range(50):
            genome = ViralGenome(
                genome_id=f"genome_{i:03d}",
                sequence="ATCG" * 1000,
                taxonomy=f"Family_{i % 10};Genus_{i % 20};Species_{i}",
                family=f"Family_{i % 10}"
            )
            community.add_genome(genome)

        # Apply distribution
        community.apply_abundance_distribution(dist_name, **params)

        # Calculate statistics
        abundances = [g.abundance for g in community.genomes]
        print(f"  Mean abundance: {np.mean(abundances):.4f}")
        print(f"  Median abundance: {np.median(abundances):.4f}")
        print(f"  Max abundance: {np.max(abundances):.4f}")
        print(f"  Min abundance: {np.min(abundances):.6f}")
        print(f"  Std dev: {np.std(abundances):.4f}")


def example_4_export_data():
    """Example 4: Export community data to files."""
    print("\n" + "=" * 70)
    print("Example 4: Exporting Community Data")
    print("=" * 70)

    # Create a gut virome community
    community = create_body_site_profile('gut', n_genomes=50, random_seed=42)

    # Export abundance table
    output_dir = Path("example_output")
    output_dir.mkdir(exist_ok=True)

    abundance_file = output_dir / "gut_virome_abundance.tsv"
    community.export_abundance_table(abundance_file)
    print(f"\n  Exported abundance table to: {abundance_file}")

    # Export genomes FASTA
    fasta_file = output_dir / "gut_virome_genomes.fasta"
    community.export_genomes_fasta(fasta_file)
    print(f"  Exported genomes FASTA to: {fasta_file}")

    print(f"\n  Files created in: {output_dir.absolute()}")


if __name__ == "__main__":
    from pathlib import Path

    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "ViroForge Community Module Examples" + " " * 17 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    # Run examples
    example_1_create_custom_community()
    example_2_body_site_profiles()
    example_3_abundance_distributions()
    example_4_export_data()

    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
    print()
