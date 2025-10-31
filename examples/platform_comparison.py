#!/usr/bin/env python3
"""
Platform Artifact Comparison Example
=====================================

This example demonstrates how different Illumina sequencing platforms introduce
different artifacts that can impact virome analysis.

Platforms compared:
1. NovaSeq 6000 (patterned flow cell, high throughput)
2. NextSeq 2000 (patterned flow cell, mid throughput)
3. MiSeq (cluster-based flow cell, low throughput)
4. HiSeq 2500 (cluster-based flow cell, legacy platform)
5. Ideal (no artifacts - control)

Key artifacts:
- PolyG tails: Patterned flow cells only (NovaSeq, NextSeq)
- Optical duplicates: All platforms (rate varies)
- Index hopping: All platforms (higher in patterned flow cells)

Author: ViroForge Development Team
"""

from viroforge.artifacts import (
    ReadPair,
    novaseq_6000,
    nextseq_2000,
    miseq,
    hiseq_2500,
    no_artifacts
)


def create_test_reads(n=10000, seed=42):
    """
    Create test reads for demonstration.

    Args:
        n: Number of reads to create
        seed: Random seed for reproducibility

    Returns:
        List of ReadPair objects
    """
    import random
    random.seed(seed)

    reads = []
    for i in range(n):
        # Simulate variety of read lengths and GC content
        read_length = random.choice([75, 100, 125, 150])

        # Generate random sequences
        bases = ['A', 'C', 'G', 'T']
        forward_seq = ''.join(random.choices(bases, k=read_length))
        reverse_seq = ''.join(random.choices(bases, k=read_length))

        read = ReadPair(
            read_id=f"read_{i:06d}",
            forward_seq=forward_seq,
            reverse_seq=reverse_seq,
            sample_index="SAMPLE01",
            genome_id=f"genome_{i % 100}"  # 100 different source genomes
        )
        reads.append(read)

    return reads


def analyze_reads(reads, platform_name):
    """
    Analyze reads for artifacts.

    Args:
        reads: List of ReadPair objects
        platform_name: Name of platform

    Returns:
        Dictionary of statistics
    """
    stats = {}

    # Total reads
    stats['total_reads'] = len(reads)

    # PolyG detection (reads ending in >10 G's)
    polyg_r1 = sum(1 for r in reads if r.forward_seq.endswith('G' * 10))
    polyg_r2 = sum(1 for r in reads if r.reverse_seq and r.reverse_seq.endswith('G' * 10))
    stats['polyg_r1_count'] = polyg_r1
    stats['polyg_r2_count'] = polyg_r2
    stats['polyg_rate'] = (polyg_r1 + polyg_r2) / (len(reads) * 2) * 100

    # Optical duplicates (reads with _optdup in ID)
    optical_dups = sum(1 for r in reads if '_optdup' in r.read_id)
    stats['optical_dup_count'] = optical_dups
    stats['optical_dup_rate'] = optical_dups / (len(reads) - optical_dups) * 100 if len(reads) > optical_dups else 0

    # Index hopping (reads with different index from SAMPLE01)
    hopped = sum(1 for r in reads if r.sample_index != "SAMPLE01")
    stats['index_hop_count'] = hopped
    stats['index_hop_rate'] = hopped / len(reads) * 100

    # Average read length
    all_lengths = [len(r.forward_seq) for r in reads]
    if reads[0].reverse_seq:
        all_lengths.extend([len(r.reverse_seq) for r in reads if r.reverse_seq])

    stats['avg_read_length'] = sum(all_lengths) / len(all_lengths)
    stats['max_read_length'] = max(all_lengths)

    return stats


def main():
    print("="*70)
    print("ViroForge: Platform Artifact Comparison")
    print("="*70)
    print()
    print("This example shows how different Illumina platforms introduce")
    print("different artifacts into sequencing data.")
    print()

    # Create baseline reads
    print("Creating 10,000 baseline read pairs...")
    baseline_reads = create_test_reads(n=10000, seed=42)
    print(f"  ✓ Created {len(baseline_reads)} read pairs")
    print(f"  Average read length: {sum(len(r.forward_seq) for r in baseline_reads)/len(baseline_reads):.1f} bp")
    print()

    # Test each platform
    platforms = {
        'NovaSeq 6000': novaseq_6000(),
        'NextSeq 2000': nextseq_2000(),
        'MiSeq': miseq(),
        'HiSeq 2500': hiseq_2500(),
        'Ideal (No Artifacts)': no_artifacts()
    }

    results = {}

    print("="*70)
    print("APPLYING PLATFORM ARTIFACTS")
    print("="*70)
    print()

    for platform_name, platform in platforms.items():
        print(f"{platform_name}:")
        print("-"*70)

        # Apply platform artifacts
        modified_reads = platform.apply(baseline_reads.copy(), random_seed=42)

        # Analyze
        stats = analyze_reads(modified_reads, platform_name)
        results[platform_name] = stats

        # Print summary
        print(f"  Total reads:          {stats['total_reads']:>10,}")
        print(f"  Optical duplicates:   {stats['optical_dup_count']:>10,} ({stats['optical_dup_rate']:>5.2f}%)")
        print(f"  PolyG tails (R1):     {stats['polyg_r1_count']:>10,}")
        print(f"  PolyG tails (R2):     {stats['polyg_r2_count']:>10,}")
        print(f"  PolyG rate:           {stats['polyg_rate']:>10.2f}%")
        print(f"  Index hopping:        {stats['index_hop_count']:>10,} ({stats['index_hop_rate']:>5.2f}%)")
        print(f"  Max read length:      {stats['max_read_length']:>10.0f} bp")
        print()

    # Comparison table
    print("="*70)
    print("PLATFORM COMPARISON TABLE")
    print("="*70)
    print()

    print(f"{'Platform':<25} {'Flow Cell':<12} {'Dup Rate':<10} {'PolyG%':<10} {'Hop%':<10}")
    print("-"*70)

    flow_cell_types = {
        'NovaSeq 6000': 'Patterned',
        'NextSeq 2000': 'Patterned',
        'MiSeq': 'Cluster',
        'HiSeq 2500': 'Cluster',
        'Ideal (No Artifacts)': 'N/A'
    }

    for platform_name, stats in results.items():
        flow_cell = flow_cell_types[platform_name]
        dup_rate = f"{stats['optical_dup_rate']:.2f}%"
        polyg_rate = f"{stats['polyg_rate']:.2f}%" if stats['polyg_rate'] > 0 else "None"
        hop_rate = f"{stats['index_hop_rate']:.2f}%"

        print(f"{platform_name:<25} {flow_cell:<12} {dup_rate:<10} {polyg_rate:<10} {hop_rate:<10}")

    print()

    # Key insights
    print("="*70)
    print("KEY INSIGHTS")
    print("="*70)
    print()

    print("PolyG Tails:")
    print("  • Only present in PATTERNED flow cells (NovaSeq, NextSeq)")
    print("  • NOT present in cluster-based flow cells (MiSeq, HiSeq)")
    print("  • Caused by incomplete fluorophore quenching")
    print("  • Can affect adapter trimming and alignment")
    print()

    print("Optical Duplicates:")
    print("  • Present in ALL platforms (rate varies)")
    print("  • Higher in patterned flow cells (8-10%)")
    print("  • Lower in cluster-based flow cells (2-5%)")
    print("  • Can inflate coverage and skew abundance estimates")
    print()

    print("Index Hopping:")
    print("  • More common in patterned flow cells (1-2%)")
    print("  • Minimal in cluster-based flow cells (~0.1-0.2%)")
    print("  • Critical for multiplexed samples")
    print("  • Can cause cross-contamination between samples")
    print()

    # Recommendations
    print("="*70)
    print("PLATFORM SELECTION GUIDE")
    print("="*70)
    print()

    print("NovaSeq 6000:")
    print("  ✓ Best for: Large cohort studies, high throughput needs")
    print("  ✓ Advantages: Lowest cost per base, fastest turnaround")
    print("  ⚠ Considerations: PolyG tails, higher optical dups, index hopping")
    print("  → Solution: Trim polyG, mark duplicates, unique dual indexes")
    print()

    print("NextSeq 2000:")
    print("  ✓ Best for: Medium-scale studies, clinical samples")
    print("  ✓ Advantages: Balanced throughput and cost")
    print("  ⚠ Considerations: PolyG tails (moderate), moderate artifacts")
    print()

    print("MiSeq:")
    print("  ✓ Best for: Pilot studies, method development, small cohorts")
    print("  ✓ Advantages: Long reads (2x300), NO polyG, lowest artifacts")
    print("  ✓ Ideal for: High-accuracy virome studies")
    print("  ⚠ Considerations: Lower throughput, higher cost per base")
    print()

    print("HiSeq 2500:")
    print("  ✓ Best for: Legacy data comparison")
    print("  ✓ Advantages: Well-characterized, NO polyG")
    print("  ⚠ Considerations: Older platform, being phased out")
    print()

    # Impact on virome analysis
    print("="*70)
    print("IMPACT ON VIROME ANALYSIS")
    print("="*70)
    print()

    print("Artifact Mitigation Strategies:")
    print()

    print("1. PolyG Tails (NovaSeq/NextSeq):")
    print("   • Use fastp or cutadapt with --polyg-trim")
    print("   • Set minimum length after trimming")
    print("   • Monitor trim rates in QC reports")
    print()

    print("2. Optical Duplicates (All Platforms):")
    print("   • Use Picard MarkDuplicates")
    print("   • Consider tile/coordinate-based detection")
    print("   • Monitor duplicate rates (>15% indicates issues)")
    print()

    print("3. Index Hopping (Patterned Flow Cells):")
    print("   • Use unique dual indexes (not combinatorial)")
    print("   • Filter reads with low-frequency barcodes (<0.1%)")
    print("   • Be cautious with low-abundance samples")
    print()

    # Benchmarking recommendations
    print("="*70)
    print("BENCHMARKING RECOMMENDATIONS")
    print("="*70)
    print()

    print("When benchmarking virome analysis pipelines:")
    print()
    print("  1. Test with MULTIPLE platforms (especially NovaSeq vs MiSeq)")
    print("  2. Include Ideal (no artifacts) as baseline")
    print("  3. Verify artifact removal steps work correctly")
    print("  4. Check if results are robust to platform choice")
    print("  5. Document which platform was used for real data")
    print()

    print("ViroForge allows you to:")
    print("  ✓ Generate matched datasets across platforms")
    print("  ✓ Test artifact removal pipelines")
    print("  ✓ Validate cross-platform reproducibility")
    print("  ✓ Optimize for specific platform artifacts")
    print()

    print("="*70)
    print("✓ Platform comparison complete!")
    print("="*70)


if __name__ == '__main__':
    main()
