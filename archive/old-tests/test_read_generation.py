#!/usr/bin/env python3
"""
End-to-end test of FASTQ read generation.

This script tests the complete ViroForge workflow:
1. Create mock virome composition
2. Generate FASTQ reads
3. Validate output files
4. Verify ground truth accuracy
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_1_quick_generate():
    """Test 1: Quick generate with small dataset."""
    logger.info("="*70)
    logger.info("Test 1: Quick Generate (Small Dataset)")
    logger.info("="*70)

    from viroforge.simulators import quick_generate, check_insilicoseq_installed

    # Verify InSilicoSeq is installed
    if not check_insilicoseq_installed():
        logger.error("InSilicoSeq is not installed!")
        return False

    logger.info("✓ InSilicoSeq 2.0.1 installed")

    try:
        logger.info("\nGenerating gut virome with realistic contamination...")
        logger.info("Parameters:")
        logger.info("  - Body site: gut")
        logger.info("  - Contamination: realistic (~7% total)")
        logger.info("  - Viral genomes: 50")
        logger.info("  - Reads: 10,000 (small test dataset)")
        logger.info("  - Model: NovaSeq")
        logger.info("  - Random seed: 42")

        start_time = datetime.now()

        output = quick_generate(
            body_site='gut',
            contamination_level='realistic',
            n_viral_genomes=50,
            n_reads=10_000,  # Small dataset for quick testing
            output_prefix='test_output/test1_gut_realistic',
            random_seed=42
        )

        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        logger.info(f"\n✓ Generation completed in {duration:.1f} seconds")
        logger.info("\nOutput files:")
        for key, path in output.items():
            size_mb = path.stat().st_size / 1024 / 1024
            logger.info(f"  {key}: {path} ({size_mb:.2f} MB)")

        return True, output

    except Exception as e:
        logger.error(f"\n✗ Test 1 failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_2_validate_fastq(output):
    """Test 2: Validate FASTQ files."""
    logger.info("\n" + "="*70)
    logger.info("Test 2: Validate FASTQ Files")
    logger.info("="*70)

    from viroforge.utils.validation import verify_fastq_file
    from Bio import SeqIO

    try:
        r1_path = output['r1']
        r2_path = output['r2']

        logger.info(f"\nValidating {r1_path}...")
        n_reads_r1 = verify_fastq_file(r1_path)
        logger.info(f"✓ R1 validated: {n_reads_r1:,} reads")

        logger.info(f"\nValidating {r2_path}...")
        n_reads_r2 = verify_fastq_file(r2_path)
        logger.info(f"✓ R2 validated: {n_reads_r2:,} reads")

        if n_reads_r1 != n_reads_r2:
            logger.error(f"✗ Read count mismatch: R1={n_reads_r1}, R2={n_reads_r2}")
            return False

        logger.info(f"\n✓ Read counts match: {n_reads_r1:,} reads in both files")

        # Sample some reads to check quality
        logger.info("\nSampling reads to check format...")
        r1_records = list(SeqIO.parse(r1_path, 'fastq'))
        r2_records = list(SeqIO.parse(r2_path, 'fastq'))

        logger.info(f"Sample R1 record:")
        logger.info(f"  ID: {r1_records[0].id}")
        logger.info(f"  Seq: {str(r1_records[0].seq)[:50]}...")
        logger.info(f"  Qual: {len(r1_records[0].letter_annotations['phred_quality'])} scores")
        logger.info(f"  Seq length: {len(r1_records[0].seq)} bp")

        # Check first 10 records for seq/qual length match
        logger.info("\nChecking first 10 records for seq/qual length match...")
        for i in range(min(10, len(r1_records))):
            r1_seq_len = len(r1_records[i].seq)
            r1_qual_len = len(r1_records[i].letter_annotations['phred_quality'])
            r2_seq_len = len(r2_records[i].seq)
            r2_qual_len = len(r2_records[i].letter_annotations['phred_quality'])

            if r1_seq_len != r1_qual_len:
                logger.error(f"✗ R1 record {i}: seq={r1_seq_len}, qual={r1_qual_len}")
                return False
            if r2_seq_len != r2_qual_len:
                logger.error(f"✗ R2 record {i}: seq={r2_seq_len}, qual={r2_qual_len}")
                return False

        logger.info("✓ All sampled records have matching seq/qual lengths")

        return True

    except Exception as e:
        logger.error(f"\n✗ Test 2 failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_3_verify_ground_truth(output):
    """Test 3: Verify ground truth metadata."""
    logger.info("\n" + "="*70)
    logger.info("Test 3: Verify Ground Truth Metadata")
    logger.info("="*70)

    try:
        gt_path = output['ground_truth']

        logger.info(f"\nReading ground truth: {gt_path}")
        df = pd.read_csv(gt_path, sep='\t')

        logger.info(f"✓ Ground truth loaded: {len(df)} genomes")

        # Check required columns
        required_cols = ['genome_id', 'genome_type', 'taxonomy', 'length',
                        'gc_content', 'abundance', 'source']
        missing_cols = set(required_cols) - set(df.columns)
        if missing_cols:
            logger.error(f"✗ Missing columns: {missing_cols}")
            return False

        logger.info("✓ All required columns present")

        # Check genome types
        genome_types = df['genome_type'].value_counts()
        logger.info("\nGenome type distribution:")
        for gtype, count in genome_types.items():
            logger.info(f"  {gtype}: {count} genomes")

        # Check abundances sum to 1.0
        total_abundance = df['abundance'].sum()
        logger.info(f"\nTotal abundance: {total_abundance:.6f}")
        if abs(total_abundance - 1.0) > 1e-6:
            logger.warning(f"⚠ Abundance sum is {total_abundance}, not 1.0")
        else:
            logger.info("✓ Abundances sum to 1.0")

        # Check viral vs contamination fractions
        viral_abundance = df[df['source'] == 'viral_community']['abundance'].sum()
        contam_abundance = df[df['source'] == 'contamination']['abundance'].sum()

        logger.info(f"\nViral abundance: {viral_abundance:.4f} ({viral_abundance*100:.1f}%)")
        logger.info(f"Contamination: {contam_abundance:.4f} ({contam_abundance*100:.1f}%)")

        # For realistic profile, expect ~90-95% viral
        if 0.85 < viral_abundance < 0.95:
            logger.info("✓ Viral fraction in expected range for 'realistic' profile")
        else:
            logger.warning(f"⚠ Viral fraction {viral_abundance:.2f} outside expected range")

        # Show sample of ground truth
        logger.info("\nSample ground truth entries:")
        logger.info(df.head(5).to_string())

        return True

    except Exception as e:
        logger.error(f"\n✗ Test 3 failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_4_contamination_levels():
    """Test 4: Generate datasets with different contamination levels."""
    logger.info("\n" + "="*70)
    logger.info("Test 4: Different Contamination Levels")
    logger.info("="*70)

    from viroforge.simulators import quick_generate

    contamination_levels = ['clean', 'realistic', 'heavy']
    results = []

    for level in contamination_levels:
        logger.info(f"\nGenerating {level} profile...")

        try:
            start_time = datetime.now()

            output = quick_generate(
                body_site='gut',
                contamination_level=level,
                n_viral_genomes=50,
                n_reads=5_000,  # Very small for quick testing
                output_prefix=f'test_output/test4_{level}',
                random_seed=42
            )

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            # Read ground truth
            df = pd.read_csv(output['ground_truth'], sep='\t')
            viral_frac = df[df['source'] == 'viral_community']['abundance'].sum()
            contam_frac = df[df['source'] == 'contamination']['abundance'].sum()

            results.append({
                'level': level,
                'duration': duration,
                'viral_fraction': viral_frac,
                'contamination': contam_frac
            })

            logger.info(f"  ✓ {level}: viral={viral_frac:.1%}, contamination={contam_frac:.1%}, time={duration:.1f}s")

        except Exception as e:
            logger.error(f"  ✗ {level} failed: {e}")
            return False

    # Summary
    logger.info("\n" + "-"*70)
    logger.info("Contamination Level Summary:")
    logger.info("-"*70)
    for r in results:
        logger.info(f"{r['level']:10s}: {r['contamination']:6.1%} contamination, {r['duration']:5.1f}s")

    logger.info("\n✓ All contamination levels tested successfully")
    return True


def test_5_platform_comparison():
    """Test 5: Different Illumina platforms."""
    logger.info("\n" + "="*70)
    logger.info("Test 5: Platform Comparison")
    logger.info("="*70)

    from viroforge.utils import create_mock_virome
    from viroforge.simulators import generate_reads
    from Bio import SeqIO

    # Create one composition to use for all platforms
    logger.info("\nCreating composition...")
    composition = create_mock_virome(
        name='platform_test',
        body_site='gut',
        contamination_level='clean',
        n_viral_genomes=30,
        random_seed=42
    )

    platforms = ['NovaSeq', 'MiSeq']  # Skip HiSeq for speed
    results = []

    for platform in platforms:
        logger.info(f"\nTesting {platform}...")

        try:
            start_time = datetime.now()

            output = generate_reads(
                composition=composition,
                output_prefix=f'test_output/test5_{platform.lower()}',
                n_reads=5_000,
                model=platform,
                cpus=2,
                random_seed=42
            )

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            # Check read lengths
            records = list(SeqIO.parse(output['r1'], 'fastq'))
            read_length = len(records[0].seq)

            results.append({
                'platform': platform,
                'read_length': read_length,
                'duration': duration
            })

            logger.info(f"  ✓ {platform}: {read_length}bp reads, {duration:.1f}s")

        except Exception as e:
            logger.error(f"  ✗ {platform} failed: {e}")
            return False

    # Summary
    logger.info("\n" + "-"*70)
    logger.info("Platform Comparison Summary:")
    logger.info("-"*70)
    for r in results:
        logger.info(f"{r['platform']:10s}: {r['read_length']}bp reads, {r['duration']:5.1f}s generation time")

    logger.info("\n✓ All platforms tested successfully")
    return True


def main():
    """Run all tests."""
    logger.info("\n" + "="*70)
    logger.info("ViroForge End-to-End Testing")
    logger.info("="*70)

    # Create output directory
    output_dir = Path('test_output')
    output_dir.mkdir(exist_ok=True)
    logger.info(f"Output directory: {output_dir.absolute()}")

    # Track results
    test_results = []

    # Test 1: Quick generate
    logger.info("\n")
    success, output = test_1_quick_generate()
    test_results.append(('Quick Generate', success))

    if not success:
        logger.error("\n✗ Test 1 failed, stopping here")
        return 1

    # Test 2: Validate FASTQ
    logger.info("\n")
    success = test_2_validate_fastq(output)
    test_results.append(('FASTQ Validation', success))

    if not success:
        logger.error("\n✗ Test 2 failed, stopping here")
        return 1

    # Test 3: Ground truth
    logger.info("\n")
    success = test_3_verify_ground_truth(output)
    test_results.append(('Ground Truth', success))

    # Test 4: Contamination levels
    logger.info("\n")
    success = test_4_contamination_levels()
    test_results.append(('Contamination Levels', success))

    # Test 5: Platforms
    logger.info("\n")
    success = test_5_platform_comparison()
    test_results.append(('Platform Comparison', success))

    # Final summary
    logger.info("\n" + "="*70)
    logger.info("TEST SUMMARY")
    logger.info("="*70)

    for test_name, success in test_results:
        status = "✓ PASS" if success else "✗ FAIL"
        logger.info(f"{test_name:25s}: {status}")

    passed = sum(1 for _, success in test_results if success)
    total = len(test_results)

    logger.info("-"*70)
    logger.info(f"Results: {passed}/{total} tests passed")

    if passed == total:
        logger.info("\n🎉 ALL TESTS PASSED!")
        return 0
    else:
        logger.info(f"\n⚠️  {total - passed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
