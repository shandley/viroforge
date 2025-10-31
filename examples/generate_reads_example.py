#!/usr/bin/env python3
"""
Example: Generating Illumina Reads from Mock Virome Compositions

This script demonstrates how to generate realistic Illumina sequencing reads
from mock virome compositions using ViroForge.

Requirements:
    - InSilicoSeq must be installed:
      conda install -c bioconda insilicoseq
      OR
      pip install InSilicoSeq

Examples demonstrate:
    1. Quick generation (one-liner)
    2. Custom composition with full control
    3. Multiple contamination levels comparison
    4. Different sequencing platforms (NovaSeq, MiSeq, HiSeq)
    5. File size estimation
"""

import logging
from pathlib import Path
from viroforge.utils import create_mock_virome
from viroforge.simulators import (
    generate_reads,
    quick_generate,
    estimate_file_size,
    check_insilicoseq_installed
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def example_1_quick_generate():
    """
    Example 1: Quick generation - one-liner for complete dataset.

    This is the simplest way to generate a complete mock virome.
    """
    print("\n" + "="*70)
    print("Example 1: Quick Generate - One-Liner")
    print("="*70)

    # Check if InSilicoSeq is installed
    if not check_insilicoseq_installed():
        print("⚠️  InSilicoSeq is not installed!")
        print("   Install with: conda install -c bioconda insilicoseq")
        print("   Skipping this example...")
        return

    print("\nGenerating gut virome with realistic contamination...")
    print("This will create:")
    print("  - 100,000 paired-end reads")
    print("  - NovaSeq error model (150bp reads)")
    print("  - Gut body-site profile (50 viral genomes)")
    print("  - Realistic contamination (~7% total)")
    print()

    # Estimate file size
    print(f"Estimated output size: {estimate_file_size(100_000)}")
    print()

    try:
        output = quick_generate(
            body_site='gut',
            contamination_level='realistic',
            n_viral_genomes=50,
            n_reads=100_000,  # Small dataset for quick demo
            output_prefix='example_output/quick_gut_virome',
            random_seed=42
        )

        print("\n✓ Generation complete!")
        print(f"  R1 FASTQ: {output['r1']}")
        print(f"  R2 FASTQ: {output['r2']}")
        print(f"  Ground truth: {output['ground_truth']}")

    except Exception as e:
        print(f"\n✗ Error: {e}")


def example_2_custom_composition():
    """
    Example 2: Custom composition with full control.

    Demonstrates creating a custom composition and generating reads
    with specific parameters.
    """
    print("\n" + "="*70)
    print("Example 2: Custom Composition with Full Control")
    print("="*70)

    if not check_insilicoseq_installed():
        print("⚠️  InSilicoSeq not installed, skipping...")
        return

    print("\nCreating custom oral virome composition...")

    # Create custom composition
    composition = create_mock_virome(
        name='oral_virome_custom',
        body_site='oral',
        contamination_level='clean',
        n_viral_genomes=30,
        viral_fraction=0.95,  # 95% viral, 5% contamination
        random_seed=123
    )

    print(f"\nComposition summary:")
    print(composition.summary())

    print("\nGenerating reads with NovaSeq error model...")

    try:
        output = generate_reads(
            composition=composition,
            output_prefix='example_output/oral_virome_custom',
            n_reads=100_000,
            model='NovaSeq',
            cpus=2,
            compress=False,
            gc_bias=True,  # Enable GC bias simulation
            validate_output=True,  # Validate FASTQ files
            random_seed=123,
            keep_temp_files=True  # Keep intermediate files for inspection
        )

        print("\n✓ Generation complete!")
        print(f"\nOutput files:")
        for key, path in output.items():
            print(f"  {key}: {path}")

    except Exception as e:
        print(f"\n✗ Error: {e}")


def example_3_contamination_comparison():
    """
    Example 3: Generate datasets with different contamination levels.

    Useful for testing QC pipeline performance with different
    contamination scenarios.
    """
    print("\n" + "="*70)
    print("Example 3: Contamination Level Comparison")
    print("="*70)

    if not check_insilicoseq_installed():
        print("⚠️  InSilicoSeq not installed, skipping...")
        return

    print("\nGenerating 4 datasets with different contamination levels:")
    print("  1. Clean (0.7% contamination)")
    print("  2. Realistic (7.4% contamination)")
    print("  3. Heavy (26.5% contamination)")
    print("  4. Failed (39.3% contamination)")
    print()

    contamination_levels = ['clean', 'realistic', 'heavy', 'failed']

    for level in contamination_levels:
        print(f"\nGenerating {level} profile...")

        try:
            output = quick_generate(
                body_site='gut',
                contamination_level=level,
                n_viral_genomes=50,
                n_reads=50_000,  # Small datasets for comparison
                output_prefix=f'example_output/gut_{level}',
                random_seed=42
            )

            print(f"  ✓ {level}: {output['r1'].name}, {output['r2'].name}")

        except Exception as e:
            print(f"  ✗ {level} failed: {e}")

    print("\n✓ All contamination levels generated!")
    print("\nUse these datasets to test your QC pipeline:")
    print("  - 'clean' should PASS all QC checks")
    print("  - 'realistic' should PASS most QC checks")
    print("  - 'heavy' should WARN or FAIL enrichment checks")
    print("  - 'failed' should FAIL enrichment checks")


def example_4_platform_comparison():
    """
    Example 4: Different Illumina platforms (NovaSeq vs MiSeq).

    Demonstrates generating reads with different sequencing platforms.
    Each platform has different read lengths and error profiles.
    """
    print("\n" + "="*70)
    print("Example 4: Sequencing Platform Comparison")
    print("="*70)

    if not check_insilicoseq_installed():
        print("⚠️  InSilicoSeq not installed, skipping...")
        return

    print("\nGenerating reads with different platforms:")
    print("  - NovaSeq: 150bp PE, ~350bp insert")
    print("  - MiSeq:   250bp PE, ~450bp insert")
    print("  - HiSeq:   125bp PE, ~350bp insert")
    print()

    # Create one composition to use for all platforms
    composition = create_mock_virome(
        name='platform_comparison',
        body_site='gut',
        contamination_level='realistic',
        n_viral_genomes=50,
        random_seed=42
    )

    platforms = ['NovaSeq', 'MiSeq', 'HiSeq']

    for platform in platforms:
        print(f"\nGenerating {platform} reads...")

        try:
            output = generate_reads(
                composition=composition,
                output_prefix=f'example_output/{platform.lower()}_reads',
                n_reads=50_000,
                model=platform,
                cpus=2,
                random_seed=42
            )

            print(f"  ✓ {platform}: {output['r1'].name}")

        except Exception as e:
            print(f"  ✗ {platform} failed: {e}")

    print("\n✓ All platforms generated!")
    print("\nCompare these datasets to see platform-specific differences:")
    print("  - Read length distributions")
    print("  - Quality score profiles")
    print("  - Error rates and patterns")


def example_5_file_size_estimation():
    """
    Example 5: Estimate file sizes before generation.

    Useful for planning storage and ensuring you have enough disk space.
    """
    print("\n" + "="*70)
    print("Example 5: File Size Estimation")
    print("="*70)

    print("\nEstimated file sizes for different read counts:")
    print()

    read_counts = [100_000, 1_000_000, 10_000_000, 50_000_000]

    for n_reads in read_counts:
        size_uncompressed = estimate_file_size(n_reads, compress=False)
        size_compressed = estimate_file_size(n_reads, compress=True)

        print(f"{n_reads:>12,} reads:")
        print(f"  Uncompressed: {size_uncompressed}")
        print(f"  Compressed:   {size_compressed}")
        print()

    print("💡 Tips:")
    print("  - Use compress=True for large datasets to save ~75% space")
    print("  - Plan for ~2x the final size during generation (temp files)")
    print("  - Validation requires decompressing (can't validate .gz files)")


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("ViroForge: Illumina Read Generation Examples")
    print("="*70)

    # Create output directory
    output_dir = Path('example_output')
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir.absolute()}")

    # Check if InSilicoSeq is installed
    if check_insilicoseq_installed():
        print("✓ InSilicoSeq is installed")
    else:
        print("⚠️  InSilicoSeq is NOT installed")
        print("\nTo run these examples, install InSilicoSeq:")
        print("  conda install -c bioconda insilicoseq")
        print("  OR")
        print("  pip install InSilicoSeq")
        print("\nSome examples will be skipped...\n")

    # Run examples
    example_1_quick_generate()
    example_2_custom_composition()
    example_3_contamination_comparison()
    example_4_platform_comparison()
    example_5_file_size_estimation()

    print("\n" + "="*70)
    print("Examples complete!")
    print("="*70)
    print("\nNext steps:")
    print("  1. Examine the generated FASTQ files")
    print("  2. Check the ground truth TSV files")
    print("  3. Run your virome analysis pipeline on these datasets")
    print("  4. Validate pipeline performance against ground truth")
    print()


if __name__ == '__main__':
    main()
