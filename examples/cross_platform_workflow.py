#!/usr/bin/env python3
"""
Cross-Platform Comparison Workflow
===================================

This example demonstrates how the SAME viral community produces different
sequencing data on different platforms (NovaSeq vs MiSeq).

This is critical for:
- Cross-platform reproducibility studies
- Platform selection for virome projects
- Validating that analysis pipelines are robust to platform choice
- Testing artifact removal strategies

Platforms compared:
1. NovaSeq 6000: Patterned flow cell, high throughput, more artifacts
2. MiSeq: Cluster flow cell, low throughput, minimal artifacts

Author: ViroForge Development Team
"""

from pathlib import Path
from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000, miseq, ReadPair


def create_mock_reads(composition: MockViromeComposition, n_reads: int, seed: int = 42) -> list:
    """Create mock reads based on composition abundances."""
    import random
    random.seed(seed)

    all_sequences = composition.get_all_sequences()
    genome_ids = [seq[0] for seq in all_sequences]
    abundances = [seq[2] for seq in all_sequences]

    total_abundance = sum(abundances)
    probabilities = [a / total_abundance for a in abundances]

    reads = []
    for i in range(n_reads):
        genome_id = random.choices(genome_ids, weights=probabilities, k=1)[0]

        read_length = 150
        bases = ['A', 'C', 'G', 'T']
        forward_seq = ''.join(random.choices(bases, k=read_length))
        reverse_seq = ''.join(random.choices(bases, k=read_length))

        read = ReadPair(
            read_id=f"read_{i:07d}",
            forward_seq=forward_seq,
            reverse_seq=reverse_seq,
            sample_index="SAMPLE01",
            genome_id=genome_id
        )
        reads.append(read)

    return reads


def analyze_reads(reads: list, platform_name: str) -> dict:
    """Analyze reads for artifacts and composition."""
    from collections import Counter

    stats = {
        'platform': platform_name,
        'total_reads': len(reads),
    }

    # Artifact detection
    stats['polyg_r1'] = sum(1 for r in reads if len(r.forward_seq) > 150)
    stats['polyg_r2'] = sum(1 for r in reads if r.reverse_seq and len(r.reverse_seq) > 150)
    stats['optical_dups'] = sum(1 for r in reads if '_optdup' in r.read_id)
    stats['index_hopped'] = sum(1 for r in reads if r.sample_index != "SAMPLE01")

    # Genome distribution
    genome_counts = Counter(r.genome_id for r in reads)
    stats['n_unique_genomes'] = len(genome_counts)
    stats['top_genome'] = genome_counts.most_common(1)[0] if genome_counts else (None, 0)

    return stats


def main():
    print("="*70)
    print("ViroForge: Cross-Platform Comparison Workflow")
    print("="*70)
    print()
    print("This example generates the SAME viral community on two different")
    print("sequencing platforms to demonstrate platform-specific artifacts")
    print("and their impact on virome analysis.")
    print()

    # Configuration
    output_dir = Path("output/cross_platform")
    output_dir.mkdir(parents=True, exist_ok=True)

    n_viral_genomes = 30
    n_reads = 50_000  # 50k reads per platform

    # ========================================================================
    # STEP 1: Create Shared Viral Community
    # ========================================================================
    print("="*70)
    print("STEP 1: Create Shared Viral Community")
    print("="*70)
    print()

    print(f"Creating gut virome community (shared across platforms)...")
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=n_viral_genomes,
        random_seed=42
    )

    print(f"  ✓ Created viral community")
    print(f"    - Genomes:     {len(viral_community.genomes)}")
    print(f"    - Mean length: {sum(g.length for g in viral_community.genomes)/len(viral_community.genomes):,.0f} bp")
    print()
    print(f"  This community will be sequenced on BOTH platforms.")
    print()

    # ========================================================================
    # STEP 2: Create Shared Contamination
    # ========================================================================
    print("="*70)
    print("STEP 2: Create Shared Contamination")
    print("="*70)
    print()

    print("Creating contamination profile (shared across platforms)...")
    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )

    print(f"  ✓ Created contamination profile")
    print(f"    - Contaminants: {len(contamination.contaminants)}")
    print()

    # ========================================================================
    # PLATFORM A: NovaSeq 6000
    # ========================================================================
    print("="*70)
    print("PLATFORM A: NovaSeq 6000 Workflow")
    print("="*70)
    print()

    # Create composition (separate for each platform to avoid shared state)
    print("Creating NovaSeq composition...")
    composition_novaseq = MockViromeComposition(
        name='gut_virome_novaseq',
        viral_community=create_body_site_profile('gut', n_genomes=n_viral_genomes, random_seed=42),
        contamination_profile=create_contamination_profile('realistic', random_seed=42),
        viral_fraction=0.50
    )

    print(f"  ✓ Initial composition: {composition_novaseq.viral_fraction:.1%} viral")
    print()

    # Apply VLP enrichment
    print("Applying VLP enrichment...")
    vlp = standard_vlp()
    vlp.apply(composition_novaseq)
    print(f"  ✓ After VLP: {composition_novaseq.viral_fraction:.1%} viral")
    print()

    # Apply amplification
    print("Applying RdAB amplification...")
    amplification = rdab_40_cycles()
    amplification.apply(composition_novaseq)
    print(f"  ✓ After amplification")
    print()

    # Generate reads
    print(f"Generating {n_reads:,} reads...")
    reads_novaseq = create_mock_reads(composition_novaseq, n_reads=n_reads, seed=42)
    print(f"  ✓ Generated {len(reads_novaseq):,} reads")
    print()

    # Apply NovaSeq artifacts
    print("Applying NovaSeq 6000 artifacts...")
    print("  - PolyG tails, optical duplicates, index hopping")
    platform_novaseq = novaseq_6000()
    reads_novaseq_final = platform_novaseq.apply(reads_novaseq, random_seed=42)

    novaseq_stats = analyze_reads(reads_novaseq_final, "NovaSeq 6000")
    print(f"  ✓ NovaSeq dataset complete")
    print(f"    - Total reads:        {novaseq_stats['total_reads']:,}")
    print(f"    - PolyG tails (R1):   {novaseq_stats['polyg_r1']:,}")
    print(f"    - PolyG tails (R2):   {novaseq_stats['polyg_r2']:,}")
    print(f"    - Optical duplicates: {novaseq_stats['optical_dups']:,}")
    print(f"    - Index hopping:      {novaseq_stats['index_hopped']:,}")
    print()

    # ========================================================================
    # PLATFORM B: MiSeq
    # ========================================================================
    print("="*70)
    print("PLATFORM B: MiSeq Workflow")
    print("="*70)
    print()

    # Create composition (separate instance)
    print("Creating MiSeq composition...")
    composition_miseq = MockViromeComposition(
        name='gut_virome_miseq',
        viral_community=create_body_site_profile('gut', n_genomes=n_viral_genomes, random_seed=42),
        contamination_profile=create_contamination_profile('realistic', random_seed=42),
        viral_fraction=0.50
    )

    print(f"  ✓ Initial composition: {composition_miseq.viral_fraction:.1%} viral")
    print()

    # Apply VLP enrichment (same protocol)
    print("Applying VLP enrichment...")
    vlp = standard_vlp()
    vlp.apply(composition_miseq)
    print(f"  ✓ After VLP: {composition_miseq.viral_fraction:.1%} viral")
    print()

    # Apply amplification (same protocol)
    print("Applying RdAB amplification...")
    amplification = rdab_40_cycles()
    amplification.apply(composition_miseq)
    print(f"  ✓ After amplification")
    print()

    # Generate reads (same seed for fair comparison)
    print(f"Generating {n_reads:,} reads...")
    reads_miseq = create_mock_reads(composition_miseq, n_reads=n_reads, seed=42)
    print(f"  ✓ Generated {len(reads_miseq):,} reads")
    print()

    # Apply MiSeq artifacts
    print("Applying MiSeq artifacts...")
    print("  - NO polyG tails (cluster flow cell)")
    print("  - Optical duplicates (lower rate)")
    print("  - Index hopping (minimal)")
    platform_miseq = miseq()
    reads_miseq_final = platform_miseq.apply(reads_miseq, random_seed=42)

    miseq_stats = analyze_reads(reads_miseq_final, "MiSeq")
    print(f"  ✓ MiSeq dataset complete")
    print(f"    - Total reads:        {miseq_stats['total_reads']:,}")
    print(f"    - PolyG tails (R1):   {miseq_stats['polyg_r1']:,}")
    print(f"    - PolyG tails (R2):   {miseq_stats['polyg_r2']:,}")
    print(f"    - Optical duplicates: {miseq_stats['optical_dups']:,}")
    print(f"    - Index hopping:      {miseq_stats['index_hopped']:,}")
    print()

    # ========================================================================
    # COMPARISON
    # ========================================================================
    print("="*70)
    print("CROSS-PLATFORM COMPARISON")
    print("="*70)
    print()

    # Comparison table
    print(f"{'Metric':<25} {'NovaSeq 6000':>15} {'MiSeq':>15} {'Difference':>15}")
    print("-"*70)

    metrics = [
        ('Total reads', novaseq_stats['total_reads'], miseq_stats['total_reads']),
        ('PolyG (R1)', novaseq_stats['polyg_r1'], miseq_stats['polyg_r1']),
        ('PolyG (R2)', novaseq_stats['polyg_r2'], miseq_stats['polyg_r2']),
        ('Optical dups', novaseq_stats['optical_dups'], miseq_stats['optical_dups']),
        ('Index hopping', novaseq_stats['index_hopped'], miseq_stats['index_hopped']),
    ]

    for metric_name, novaseq_val, miseq_val in metrics:
        diff = novaseq_val - miseq_val
        diff_str = f"+{diff:,}" if diff > 0 else f"{diff:,}"
        print(f"{metric_name:<25} {novaseq_val:>15,} {miseq_val:>15,} {diff_str:>15}")

    print()

    # Key insights
    print("="*70)
    print("KEY INSIGHTS")
    print("="*70)
    print()

    polyg_diff = novaseq_stats['polyg_r1'] + novaseq_stats['polyg_r2'] - miseq_stats['polyg_r1'] - miseq_stats['polyg_r2']
    dup_diff = novaseq_stats['optical_dups'] - miseq_stats['optical_dups']
    hop_diff = novaseq_stats['index_hopped'] - miseq_stats['index_hopped']

    print(f"1. PolyG Tails:")
    print(f"   NovaSeq has {polyg_diff:,} MORE reads with polyG tails")
    print(f"   This is expected for patterned flow cells")
    print(f"   → Solution: Use fastp --polyg-trim")
    print()

    print(f"2. Optical Duplicates:")
    print(f"   NovaSeq has {dup_diff:,} MORE optical duplicates")
    print(f"   This is due to higher cluster density")
    print(f"   → Solution: Use Picard MarkDuplicates")
    print()

    print(f"3. Index Hopping:")
    print(f"   NovaSeq has {hop_diff:,} MORE index hopping events")
    print(f"   This is due to ExAmp chemistry in patterned flow cells")
    print(f"   → Solution: Use unique dual indexes, filter low-frequency barcodes")
    print()

    print(f"4. Data Quality:")
    print(f"   MiSeq produces cleaner data with fewer artifacts")
    print(f"   NovaSeq produces more data but requires more QC")
    print(f"   → Choose based on throughput needs vs artifact tolerance")
    print()

    # ========================================================================
    # EXPORT RESULTS
    # ========================================================================
    print("="*70)
    print("EXPORTING RESULTS")
    print("="*70)
    print()

    # Export comparison summary
    output_comparison = output_dir / "platform_comparison_summary.txt"
    with open(output_comparison, 'w') as f:
        f.write("ViroForge Cross-Platform Comparison Summary\n")
        f.write("="*70 + "\n\n")

        f.write("Study Design:\n")
        f.write(f"  - Same viral community sequenced on both platforms\n")
        f.write(f"  - Same VLP enrichment protocol\n")
        f.write(f"  - Same amplification protocol (RdAB, 40 cycles)\n")
        f.write(f"  - {n_reads:,} reads per platform\n\n")

        f.write("Platform A: NovaSeq 6000\n")
        f.write(f"  Total reads:        {novaseq_stats['total_reads']:,}\n")
        f.write(f"  PolyG tails (R1):   {novaseq_stats['polyg_r1']:,} ({novaseq_stats['polyg_r1']/novaseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  PolyG tails (R2):   {novaseq_stats['polyg_r2']:,} ({novaseq_stats['polyg_r2']/novaseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  Optical duplicates: {novaseq_stats['optical_dups']:,} ({novaseq_stats['optical_dups']/novaseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  Index hopping:      {novaseq_stats['index_hopped']:,} ({novaseq_stats['index_hopped']/novaseq_stats['total_reads']*100:.2f}%)\n\n")

        f.write("Platform B: MiSeq\n")
        f.write(f"  Total reads:        {miseq_stats['total_reads']:,}\n")
        f.write(f"  PolyG tails (R1):   {miseq_stats['polyg_r1']:,} ({miseq_stats['polyg_r1']/miseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  PolyG tails (R2):   {miseq_stats['polyg_r2']:,} ({miseq_stats['polyg_r2']/miseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  Optical duplicates: {miseq_stats['optical_dups']:,} ({miseq_stats['optical_dups']/miseq_stats['total_reads']*100:.2f}%)\n")
        f.write(f"  Index hopping:      {miseq_stats['index_hopped']:,} ({miseq_stats['index_hopped']/miseq_stats['total_reads']*100:.2f}%)\n\n")

        f.write("Conclusion:\n")
        f.write(f"  - NovaSeq has {polyg_diff:,} more polyG tails\n")
        f.write(f"  - NovaSeq has {dup_diff:,} more optical duplicates\n")
        f.write(f"  - NovaSeq has {hop_diff:,} more index hopping events\n")
        f.write(f"  - MiSeq produces cleaner data with fewer artifacts\n")
        f.write(f"  - Both platforms are suitable for virome studies with proper QC\n")

    print(f"  ✓ Exported comparison summary")
    print(f"    - File: {output_comparison}")
    print()

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    print("="*70)
    print("WORKFLOW COMPLETE!")
    print("="*70)
    print()
    print(f"  Output directory: {output_dir}")
    print()
    print(f"  This cross-platform comparison demonstrates:")
    print(f"    ✓ Same viral community on different platforms")
    print(f"    ✓ Platform-specific artifact profiles")
    print(f"    ✓ Impact of flow cell technology (patterned vs cluster)")
    print(f"    ✓ Importance of platform-aware QC strategies")
    print()
    print(f"  Use this for:")
    print(f"    • Testing cross-platform reproducibility")
    print(f"    • Validating artifact removal pipelines")
    print(f"    • Platform selection for new studies")
    print(f"    • Training bioinformaticians on platform differences")
    print()
    print("="*70)


if __name__ == '__main__':
    main()
