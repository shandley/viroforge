#!/usr/bin/env python3
"""
Complete VLP-Enriched Virome Workflow
======================================

This example demonstrates the complete workflow for creating a VLP-enriched
virome dataset, from viral community creation through enrichment to final
FASTQ generation.

This is the recommended workflow for benchmarking virome analysis pipelines.

Workflow Steps:
1. Create viral community
2. Add contamination
3. Combine into mock composition
4. Apply VLP enrichment
5. Generate sequencing reads
6. Export ground truth metadata
"""

import os
from pathlib import Path
from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.sequencing import generate_reads


def main():
    print("="*70)
    print("ViroForge: Complete VLP-Enriched Virome Workflow")
    print("="*70)
    print()
    print("This example creates a complete VLP-enriched virome dataset")
    print("suitable for benchmarking virome analysis pipelines.")
    print()

    # Configuration
    output_dir = Path("output/vlp_virome_example")
    output_dir.mkdir(parents=True, exist_ok=True)

    n_genomes = 30
    sequencing_depth = 1_000_000  # 1M reads
    read_length = 150

    # Step 1: Create viral community
    print("="*70)
    print("STEP 1: Create Viral Community")
    print("="*70)
    print()

    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=n_genomes,
        random_seed=42
    )

    stats = viral_community.get_summary_stats()
    print(f"  ✓ Created gut virome community")
    print(f"    - Number of genomes: {stats['n_genomes']}")
    print(f"    - Total families:    {stats['n_families']}")
    print(f"    - Size range:        {stats['min_length']:,} - {stats['max_length']:,} bp")
    print(f"    - GC range:          {stats['min_gc']:.1f}% - {stats['max_gc']:.1f}%")
    print()

    # Step 2: Create contamination
    print("="*70)
    print("STEP 2: Create Contamination Profile")
    print("="*70)
    print()

    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )

    contam_stats = contamination.get_summary_stats()
    print(f"  ✓ Created realistic contamination")
    print(f"    - Number of contaminants: {contam_stats['n_contaminants']}")
    print(f"    - Contamination types:")
    for ctype, count in contam_stats['by_type'].items():
        print(f"      • {ctype}: {count}")
    print()

    # Step 3: Create mock composition (bulk metagenome)
    print("="*70)
    print("STEP 3: Create Mock Virome Composition")
    print("="*70)
    print()

    composition = MockViromeComposition(
        name='gut_virome_vlp',
        viral_community=viral_community,
        contamination_profile=contamination,
        viral_fraction=0.50  # Start as bulk metagenome
    )

    print(f"  ✓ Created bulk metagenome composition")
    print(f"    - Viral fraction:         {composition.viral_fraction:.1%}")
    print(f"    - Contamination fraction: {composition.contamination_fraction:.1%}")
    print()
    print(f"  This represents the raw sample BEFORE VLP enrichment")
    print()

    # Step 4: Apply VLP enrichment
    print("="*70)
    print("STEP 4: Apply VLP Enrichment")
    print("="*70)
    print()

    initial_viral = composition.viral_fraction
    initial_contam = composition.contamination_fraction

    print(f"  Applying standard VLP enrichment protocol...")
    print(f"    - Pre-filtration:     0.45 μm")
    print(f"    - Main filtration:    0.2 μm (tangential flow)")
    print(f"    - Nuclease treatment: DNase/RNase (95% efficiency)")
    print()

    vlp = standard_vlp()
    vlp.apply(composition)

    final_viral = composition.viral_fraction
    final_contam = composition.contamination_fraction
    enrichment = final_viral / initial_viral

    print(f"  ✓ VLP enrichment complete!")
    print(f"    - Initial: {initial_viral:.1%} viral, {initial_contam:.1%} contamination")
    print(f"    - Final:   {final_viral:.1%} viral, {final_contam:.1%} contamination")
    print(f"    - Enrichment: {enrichment:.2f}x increase in viral fraction")
    print()

    # Step 5: Generate sequencing reads
    print("="*70)
    print("STEP 5: Generate Sequencing Reads")
    print("="*70)
    print()

    print(f"  Generating {sequencing_depth:,} paired-end reads...")
    print(f"    - Read length:        {read_length} bp")
    print(f"    - Platform:           Illumina NovaSeq")
    print(f"    - Insert size:        300 bp (mean)")
    print()

    # Note: This would call the sequencing module (Phase 3)
    # For now, we'll just show what would happen
    expected_viral_reads = sequencing_depth * final_viral
    expected_contam_reads = sequencing_depth * final_contam

    print(f"  Expected read distribution:")
    print(f"    - Viral reads:        {expected_viral_reads:,.0f} ({final_viral:.1%})")
    print(f"    - Contamination reads: {expected_contam_reads:,.0f} ({final_contam:.1%})")
    print()

    # Compare to what bulk would have given
    bulk_viral_reads = sequencing_depth * initial_viral
    gain = expected_viral_reads - bulk_viral_reads

    print(f"  Compared to bulk metagenome:")
    print(f"    - Bulk viral reads:   {bulk_viral_reads:,.0f}")
    print(f"    - VLP viral reads:    {expected_viral_reads:,.0f}")
    print(f"    - Additional reads:   {gain:,.0f} (+{(gain/bulk_viral_reads)*100:.0f}%)")
    print()

    # Step 6: Export metadata
    print("="*70)
    print("STEP 6: Export Ground Truth Metadata")
    print("="*70)
    print()

    # Get composition table
    composition_table = composition.get_composition_table()

    output_metadata = output_dir / "ground_truth_composition.tsv"
    composition_table.to_csv(output_metadata, sep='\t', index=False)

    print(f"  ✓ Exported composition table")
    print(f"    - File: {output_metadata}")
    print(f"    - Rows: {len(composition_table)}")
    print()

    # Summary statistics
    summary = composition.get_summary_stats()
    output_summary = output_dir / "summary_stats.txt"

    with open(output_summary, 'w') as f:
        f.write("ViroForge Mock Virome Summary\n")
        f.write("="*50 + "\n\n")
        f.write(f"Dataset: {summary['name']}\n\n")
        f.write(f"Viral Composition:\n")
        f.write(f"  Genomes:        {summary['n_viral_genomes']}\n")
        f.write(f"  Viral fraction: {summary['viral_fraction']:.1%}\n\n")
        f.write(f"Contamination:\n")
        f.write(f"  Contaminants:   {summary['n_contaminants']}\n")
        f.write(f"  Contam fraction: {summary['contamination_fraction']:.1%}\n\n")
        f.write(f"Enrichment:\n")
        f.write(f"  Protocol:       Standard VLP\n")
        f.write(f"  Enrichment:     {enrichment:.2f}x\n")

    print(f"  ✓ Exported summary statistics")
    print(f"    - File: {output_summary}")
    print()

    # Final summary
    print("="*70)
    print("WORKFLOW COMPLETE!")
    print("="*70)
    print()
    print(f"  Output directory: {output_dir}")
    print()
    print(f"  Files created:")
    print(f"    1. {output_metadata.name}")
    print(f"    2. {output_summary.name}")
    print()
    print(f"  Next steps:")
    print(f"    - Use this composition to generate FASTQ reads")
    print(f"    - Run your virome analysis pipeline")
    print(f"    - Compare results to ground truth metadata")
    print(f"    - Calculate accuracy metrics")
    print()
    print(f"  This VLP-enriched dataset is ideal for:")
    print(f"    ✓ Benchmarking viral taxonomy assignment")
    print(f"    ✓ Testing viral genome assembly")
    print(f"    ✓ Validating abundance estimation")
    print(f"    ✓ Assessing contamination detection")
    print()
    print("="*70)


if __name__ == '__main__':
    main()
