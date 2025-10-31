#!/usr/bin/env python3
"""
Complete Integrated Virome Workflow
====================================

This example demonstrates the complete ViroForge pipeline, integrating all
Phase 2 features into a realistic virome sequencing simulation.

Pipeline Steps:
1. Create viral community (body-site specific)
2. Add realistic contamination (host DNA, bacteria, etc.)
3. Apply VLP enrichment (virus purification)
4. Apply amplification bias (library preparation)
5. Generate reads with platform artifacts (sequencing)
6. Export ground truth metadata

This workflow simulates a realistic gut virome study using NovaSeq sequencing.

Author: ViroForge Development Team
"""

from pathlib import Path
from viroforge.core import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000, ReadPair


def create_mock_reads(composition: MockViromeComposition, n_reads: int, seed: int = 42) -> list:
    """
    Create mock reads based on composition abundances.

    In a real workflow, this would call InSilicoSeq or similar.
    For demonstration, we create ReadPair objects directly.

    Args:
        composition: MockViromeComposition with abundances
        n_reads: Number of read pairs to generate
        seed: Random seed

    Returns:
        List of ReadPair objects
    """
    import random
    random.seed(seed)

    # Get all genomes with their abundances
    all_sequences = composition.get_all_sequences()

    # Create cumulative abundance for sampling
    genome_ids = [seq[0] for seq in all_sequences]
    abundances = [seq[2] for seq in all_sequences]

    # Normalize abundances
    total_abundance = sum(abundances)
    probabilities = [a / total_abundance for a in abundances]

    reads = []
    for i in range(n_reads):
        # Sample genome based on abundance
        genome_id = random.choices(genome_ids, weights=probabilities, k=1)[0]

        # Generate random sequences (in real workflow, these come from genomes)
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


def analyze_composition(composition: MockViromeComposition, stage: str) -> dict:
    """
    Analyze composition at a specific pipeline stage.

    Args:
        composition: MockViromeComposition to analyze
        stage: Pipeline stage name

    Returns:
        Dictionary of statistics
    """
    stats = {
        'stage': stage,
        'n_viral_genomes': len(composition.viral_community.genomes),
        'viral_fraction': composition.viral_fraction,
        'contamination_fraction': composition.contamination_fraction
    }

    if composition.contamination_profile:
        stats['n_contaminants'] = len(composition.contamination_profile.contaminants)

        # Count by contamination type
        from collections import Counter
        contam_types = Counter(c.contaminant_type.value for c in composition.contamination_profile.contaminants)
        stats['contamination_breakdown'] = dict(contam_types)

    return stats


def main():
    print("="*70)
    print("ViroForge: Complete Integrated Virome Workflow")
    print("="*70)
    print()
    print("This example demonstrates the complete ViroForge pipeline:")
    print("  Community → Contamination → VLP → Amplification → Sequencing")
    print()

    # Configuration
    output_dir = Path("output/complete_workflow")
    output_dir.mkdir(parents=True, exist_ok=True)

    n_viral_genomes = 50
    n_reads = 100_000  # 100k reads for demonstration

    # ========================================================================
    # STEP 1: Create Viral Community
    # ========================================================================
    print("="*70)
    print("STEP 1: Create Viral Community")
    print("="*70)
    print()

    print(f"Creating gut virome community with {n_viral_genomes} genomes...")
    viral_community = create_body_site_profile(
        body_site='gut',
        n_genomes=n_viral_genomes,
        random_seed=42
    )

    stats = viral_community.get_summary_stats()
    print(f"  ✓ Created gut virome community")
    print(f"    - Genomes:      {stats['n_genomes']}")
    print(f"    - Families:     {stats['families']}")
    print(f"    - Mean length:  {stats['mean_length']:,.0f} bp")
    print(f"    - Mean GC:      {stats['mean_gc']*100:.1f}%")
    print()

    # ========================================================================
    # STEP 2: Add Contamination
    # ========================================================================
    print("="*70)
    print("STEP 2: Add Realistic Contamination")
    print("="*70)
    print()

    print("Creating realistic contamination profile...")
    contamination = create_contamination_profile(
        contamination_level='realistic',
        random_seed=42
    )

    contam_stats = contamination.get_summary_stats()
    print(f"  ✓ Created contamination profile")
    print(f"    - Contaminants: {contam_stats['n_contaminants']}")
    print(f"    - Types:")
    for ctype, count in contam_stats['by_type'].items():
        print(f"      • {ctype}: {count}")
    print()

    # ========================================================================
    # STEP 3: Create Initial Composition (Bulk Metagenome)
    # ========================================================================
    print("="*70)
    print("STEP 3: Create Initial Composition")
    print("="*70)
    print()

    print("Combining viral community and contamination...")
    composition = MockViromeComposition(
        name='gut_virome_novaseq',
        viral_community=viral_community,
        contamination_profile=contamination,
        viral_fraction=0.50  # Start with 50% viral (typical bulk metagenome)
    )

    initial_stats = analyze_composition(composition, 'Initial (Bulk Metagenome)')
    print(f"  ✓ Created bulk metagenome composition")
    print(f"    - Viral fraction:         {initial_stats['viral_fraction']:.1%}")
    print(f"    - Contamination fraction: {initial_stats['contamination_fraction']:.1%}")
    print()
    print(f"  This represents the raw sample BEFORE any processing.")
    print()

    # ========================================================================
    # STEP 4: Apply VLP Enrichment
    # ========================================================================
    print("="*70)
    print("STEP 4: Apply VLP Enrichment")
    print("="*70)
    print()

    print("Applying standard VLP enrichment protocol...")
    print("  - Pre-filtration:     0.45 μm")
    print("  - Main filtration:    0.2 μm (tangential flow)")
    print("  - Nuclease treatment: DNase/RNase (95% efficiency)")
    print()

    vlp = standard_vlp()
    vlp.apply(composition)

    vlp_stats = analyze_composition(composition, 'After VLP Enrichment')
    enrichment_factor = vlp_stats['viral_fraction'] / initial_stats['viral_fraction']

    print(f"  ✓ VLP enrichment complete")
    print(f"    - Viral fraction:         {vlp_stats['viral_fraction']:.1%}")
    print(f"    - Contamination fraction: {vlp_stats['contamination_fraction']:.1%}")
    print(f"    - Enrichment factor:      {enrichment_factor:.2f}x")
    print()

    # ========================================================================
    # STEP 5: Apply Amplification Bias
    # ========================================================================
    print("="*70)
    print("STEP 5: Apply Amplification Bias")
    print("="*70)
    print()

    print("Applying RdAB amplification (40 cycles)...")
    print("  - Method:        Random RT + dsDNA synthesis + PCR")
    print("  - Cycles:        40")
    print("  - Length bias:   Favors short genomes")
    print("  - GC bias:       Optimal ~50% GC")
    print()

    amplification = rdab_40_cycles()
    amplification.apply(composition)

    amp_stats = analyze_composition(composition, 'After Amplification')
    print(f"  ✓ Amplification complete")
    print(f"    - Viral fraction:         {amp_stats['viral_fraction']:.1%}")
    print(f"    - Contamination fraction: {amp_stats['contamination_fraction']:.1%}")
    print()
    print(f"  Note: Amplification may slightly change viral/contam ratio")
    print(f"        due to differential bias on different genome sizes.")
    print()

    # ========================================================================
    # STEP 6: Generate Reads
    # ========================================================================
    print("="*70)
    print("STEP 6: Generate Sequencing Reads")
    print("="*70)
    print()

    print(f"Generating {n_reads:,} paired-end reads...")
    print(f"  - Read length:  150 bp")
    print(f"  - Insert size:  ~300 bp (simulated)")
    print()

    reads = create_mock_reads(composition, n_reads=n_reads, seed=42)

    print(f"  ✓ Generated {len(reads):,} read pairs")
    print()

    # Count reads by source
    from collections import Counter
    genome_counts = Counter(r.genome_id for r in reads)

    # Separate viral from contamination
    viral_genomes = {g.genome_id for g in composition.viral_community.genomes}
    viral_reads = sum(count for genome_id, count in genome_counts.items() if genome_id in viral_genomes)
    contam_reads = len(reads) - viral_reads

    print(f"  Read distribution:")
    print(f"    - Viral reads:        {viral_reads:>8,} ({viral_reads/len(reads)*100:>5.1f}%)")
    print(f"    - Contamination reads: {contam_reads:>8,} ({contam_reads/len(reads)*100:>5.1f}%)")
    print()

    # ========================================================================
    # STEP 7: Apply Platform Artifacts
    # ========================================================================
    print("="*70)
    print("STEP 7: Apply Platform Artifacts (NovaSeq 6000)")
    print("="*70)
    print()

    print("Applying NovaSeq 6000 platform artifacts...")
    print("  - PolyG tails:       2.5% of reads")
    print("  - Optical duplicates: 9% duplication rate")
    print("  - Index hopping:     1.5% misassignment")
    print()

    platform = novaseq_6000()
    reads_with_artifacts = platform.apply(reads, random_seed=42)

    print(f"  ✓ Platform artifacts applied")
    print(f"    - Input reads:   {len(reads):>8,}")
    print(f"    - Output reads:  {len(reads_with_artifacts):>8,}")
    print(f"    - Increase:      {len(reads_with_artifacts) - len(reads):>8,} "
          f"(optical duplicates)")
    print()

    # Analyze artifacts
    polyg_r1 = sum(1 for r in reads_with_artifacts if len(r.forward_seq) > 150)
    polyg_r2 = sum(1 for r in reads_with_artifacts if r.reverse_seq and len(r.reverse_seq) > 150)
    optical_dups = sum(1 for r in reads_with_artifacts if '_optdup' in r.read_id)
    hopped = sum(1 for r in reads_with_artifacts if r.sample_index != "SAMPLE01")

    print(f"  Artifact statistics:")
    print(f"    - PolyG tails (R1):   {polyg_r1:>6} ({polyg_r1/len(reads_with_artifacts)*100:.2f}%)")
    print(f"    - PolyG tails (R2):   {polyg_r2:>6} ({polyg_r2/len(reads_with_artifacts)*100:.2f}%)")
    print(f"    - Optical duplicates: {optical_dups:>6} ({optical_dups/len(reads_with_artifacts)*100:.2f}%)")
    print(f"    - Index hopping:      {hopped:>6} ({hopped/len(reads_with_artifacts)*100:.2f}%)")
    print()

    # ========================================================================
    # STEP 8: Export Results and Ground Truth
    # ========================================================================
    print("="*70)
    print("STEP 8: Export Results and Ground Truth")
    print("="*70)
    print()

    # Export composition table
    composition_table = composition.get_composition_table()
    output_composition = output_dir / "ground_truth_composition.tsv"
    composition_table.to_csv(output_composition, sep='\t', index=False)

    print(f"  ✓ Exported composition table")
    print(f"    - File: {output_composition}")
    print(f"    - Entries: {len(composition_table)}")
    print()

    # Export read-to-genome mapping (ground truth)
    output_read_mapping = output_dir / "ground_truth_read_mapping.tsv"
    with open(output_read_mapping, 'w') as f:
        f.write("read_id\tgenome_id\tsample_index\thas_polyg_r1\thas_polyg_r2\tis_optical_dup\n")
        for read in reads_with_artifacts[:1000]:  # First 1000 for demonstration
            has_polyg_r1 = "TRUE" if len(read.forward_seq) > 150 else "FALSE"
            has_polyg_r2 = "TRUE" if (read.reverse_seq and len(read.reverse_seq) > 150) else "FALSE"
            is_dup = "TRUE" if '_optdup' in read.read_id else "FALSE"
            f.write(f"{read.read_id}\t{read.genome_id}\t{read.sample_index}\t"
                   f"{has_polyg_r1}\t{has_polyg_r2}\t{is_dup}\n")

    print(f"  ✓ Exported read mapping (first 1000 reads)")
    print(f"    - File: {output_read_mapping}")
    print()

    # Export pipeline summary
    output_summary = output_dir / "pipeline_summary.txt"
    with open(output_summary, 'w') as f:
        f.write("ViroForge Complete Workflow Summary\n")
        f.write("="*70 + "\n\n")

        f.write("Pipeline Configuration:\n")
        f.write(f"  Body site:         Gut\n")
        f.write(f"  Viral genomes:     {n_viral_genomes}\n")
        f.write(f"  Contamination:     Realistic\n")
        f.write(f"  VLP enrichment:    Standard (0.2 μm, 95% nuclease)\n")
        f.write(f"  Amplification:     RdAB (40 cycles)\n")
        f.write(f"  Platform:          NovaSeq 6000\n")
        f.write(f"  Total reads:       {len(reads_with_artifacts):,}\n\n")

        f.write("Pipeline Stages:\n\n")

        f.write(f"1. Initial Composition (Bulk Metagenome)\n")
        f.write(f"   Viral fraction:    {initial_stats['viral_fraction']:.1%}\n")
        f.write(f"   Contam fraction:   {initial_stats['contamination_fraction']:.1%}\n\n")

        f.write(f"2. After VLP Enrichment\n")
        f.write(f"   Viral fraction:    {vlp_stats['viral_fraction']:.1%}\n")
        f.write(f"   Contam fraction:   {vlp_stats['contamination_fraction']:.1%}\n")
        f.write(f"   Enrichment:        {enrichment_factor:.2f}x\n\n")

        f.write(f"3. After Amplification (RdAB)\n")
        f.write(f"   Viral fraction:    {amp_stats['viral_fraction']:.1%}\n")
        f.write(f"   Contam fraction:   {amp_stats['contamination_fraction']:.1%}\n\n")

        f.write(f"4. Final Read Statistics\n")
        f.write(f"   Total reads:       {len(reads_with_artifacts):,}\n")
        f.write(f"   Viral reads:       {viral_reads:,} ({viral_reads/len(reads)*100:.1f}%)\n")
        f.write(f"   Contam reads:      {contam_reads:,} ({contam_reads/len(reads)*100:.1f}%)\n\n")

        f.write(f"5. Platform Artifacts (NovaSeq 6000)\n")
        f.write(f"   PolyG tails (R1):  {polyg_r1} ({polyg_r1/len(reads_with_artifacts)*100:.2f}%)\n")
        f.write(f"   PolyG tails (R2):  {polyg_r2} ({polyg_r2/len(reads_with_artifacts)*100:.2f}%)\n")
        f.write(f"   Optical dups:      {optical_dups} ({optical_dups/len(reads_with_artifacts)*100:.2f}%)\n")
        f.write(f"   Index hopping:     {hopped} ({hopped/len(reads_with_artifacts)*100:.2f}%)\n\n")

        f.write("Ground Truth Files:\n")
        f.write(f"  - {output_composition.name}\n")
        f.write(f"  - {output_read_mapping.name}\n")
        f.write(f"  - {output_summary.name}\n")

    print(f"  ✓ Exported pipeline summary")
    print(f"    - File: {output_summary}")
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
    print(f"  Key Results:")
    print(f"    - Initial viral fraction:  {initial_stats['viral_fraction']:.1%}")
    print(f"    - Final viral fraction:    {amp_stats['viral_fraction']:.1%}")
    print(f"    - VLP enrichment:          {enrichment_factor:.2f}x")
    print(f"    - Total reads generated:   {len(reads_with_artifacts):,}")
    print(f"    - Optical duplicates:      {optical_dups:,}")
    print()
    print(f"  This dataset is ready for:")
    print(f"    ✓ Virome analysis pipeline benchmarking")
    print(f"    ✓ Taxonomy assignment validation")
    print(f"    ✓ Abundance estimation accuracy testing")
    print(f"    ✓ Artifact removal pipeline validation")
    print(f"    ✓ Cross-platform comparison studies")
    print()
    print(f"  Next steps:")
    print(f"    1. Run your virome analysis pipeline on the generated reads")
    print(f"    2. Compare results to ground truth (ground_truth_composition.tsv)")
    print(f"    3. Calculate accuracy metrics (precision, recall, abundance correlation)")
    print(f"    4. Test artifact removal (fastp for polyG, Picard for duplicates)")
    print()
    print("="*70)


if __name__ == '__main__':
    main()
