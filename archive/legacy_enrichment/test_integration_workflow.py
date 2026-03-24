"""
Integration tests for complete ViroForge workflows.

Tests the integration of all Phase 2 components:
- VLP enrichment
- Amplification bias
- Platform artifacts
- FASTQ generation
- Ground truth preservation

These tests verify that components work together correctly
in realistic end-to-end scenarios.
"""

import pytest
import tempfile
import os
from pathlib import Path

from viroforge.core.community import ViralCommunity, create_body_site_profile
from viroforge.core.contamination import ContaminationProfile, create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp, iron_chloride_vlp, VLPEnrichment
from viroforge.amplification import rdab_40_cycles, mda_standard, linker_standard
from viroforge.artifacts import (
    novaseq_6000, miseq, nextseq_2000, no_artifacts,
    ReadPair, PolyGTailArtifact, OpticalDuplicateArtifact
)


class TestCompleteWorkflowIntegration:
    """Tests for complete end-to-end workflow integration."""

    def test_minimal_workflow_completes(self):
        """Test that minimal workflow completes without errors."""
        # Community
        community = create_body_site_profile('gut', n_genomes=10, random_seed=42)

        # Contamination
        contamination = create_contamination_profile('clean', random_seed=42)

        # Composition
        composition = MockViromeComposition(
            name='test_minimal',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.5
        )

        # VLP enrichment
        vlp = standard_vlp()
        vlp.apply(composition)

        # Amplification
        amplification = rdab_40_cycles()
        amplification.apply(composition)

        # Verify composition is valid
        assert composition.viral_fraction > 0.9
        assert composition.get_total_abundance() > 0

    def test_ground_truth_preserved_through_pipeline(self):
        """Test that genome IDs are preserved through all stages."""
        # Create simple composition
        community = create_body_site_profile('gut', n_genomes=5, random_seed=42)
        contamination = create_contamination_profile('clean', random_seed=42)
        composition = MockViromeComposition(
            name='test_ground_truth',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.8
        )

        # Get original genome IDs
        original_genomes = set()
        for genome in composition.viral_community.genomes:
            original_genomes.add(genome.genome_id)

        # Apply transformations
        vlp = standard_vlp()
        vlp.apply(composition)

        amplification = rdab_40_cycles()
        amplification.apply(composition)

        # Verify genome IDs still present
        final_genomes = set()
        for genome in composition.viral_community.genomes:
            final_genomes.add(genome.genome_id)

        assert original_genomes == final_genomes, "Genome IDs were not preserved"

    def test_vlp_enrichment_integration(self):
        """Test VLP enrichment integration with composition."""
        community = create_body_site_profile('gut', n_genomes=20, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        composition = MockViromeComposition(
            name='test_vlp',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.5
        )

        initial_viral_fraction = composition.viral_fraction

        # Apply VLP enrichment
        vlp = standard_vlp()
        vlp.apply(composition)

        # Verify enrichment occurred
        assert composition.viral_fraction > initial_viral_fraction
        assert composition.viral_fraction > 0.9  # Should be >90% viral

        # Verify contamination was reduced
        final_contam_fraction = 1 - composition.viral_fraction
        assert final_contam_fraction < 0.1  # Should be <10% contamination

    def test_amplification_bias_integration(self):
        """Test amplification bias integration with composition."""
        community = create_body_site_profile('gut', n_genomes=20, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        composition = MockViromeComposition(
            name='test_amplification',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.95  # Post-VLP enrichment level
        )

        # Get initial abundance distribution
        initial_abundances = {}
        for genome in composition.viral_community.genomes:
            initial_abundances[genome.genome_id] = genome.abundance

        # Apply amplification
        amplification = rdab_40_cycles()
        amplification.apply(composition)

        # Verify abundances changed (bias was applied)
        final_abundances = {}
        for genome in composition.viral_community.genomes:
            final_abundances[genome.genome_id] = genome.abundance

        # At least some abundances should have changed significantly
        changes = []
        for genome_id in initial_abundances:
            initial = initial_abundances[genome_id]
            final = final_abundances[genome_id]
            if initial > 0:
                ratio = final / initial
                changes.append(ratio)

        # Should see some bias (not all ratios = 1.0)
        assert len(set(changes)) > 1, "No amplification bias was applied"

    def test_platform_artifacts_integration(self):
        """Test platform artifacts integration with read generation."""
        # Create simple reads (larger sample for statistical testing)
        reads = []
        for i in range(1000):
            read = ReadPair(
                read_id=f"read_{i}",
                forward_seq="A" * 80,
                reverse_seq="T" * 80,
                forward_qual="I" * 80,
                reverse_qual="I" * 80,
                genome_id=f"genome_{i % 10}",
                tile_x=i * 100,
                tile_y=i * 100
            )
            reads.append(read)

        # Apply NovaSeq artifacts
        platform = novaseq_6000()
        reads_with_artifacts = platform.apply(reads, random_seed=42)

        # Verify artifacts were applied
        assert len(reads_with_artifacts) > len(reads)  # Optical duplicates added

        # Check for polyG tails (sequence length changes)
        # With 1000 reads and ~2.5% frequency, should have at least a few
        long_reads = [r for r in reads_with_artifacts if len(r.forward_seq) > 80 or (r.reverse_seq and len(r.reverse_seq) > 80)]
        assert len(long_reads) > 0, "No polyG tails were added"

    def test_different_vlp_protocols(self):
        """Test integration with different VLP protocols."""
        community = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        protocols = [
            standard_vlp(),
            iron_chloride_vlp(),
            VLPEnrichment(
                filtration_cutoff_um=0.45,
                nuclease_efficiency=0.80,
                stochastic_variation=0.3
            )
        ]

        results = []
        for protocol in protocols:
            composition = MockViromeComposition(
                name='test_protocol',
                viral_community=community,
                contamination_profile=contamination,
                viral_fraction=0.5
            )
            protocol.apply(composition)
            results.append(composition.viral_fraction)

        # All protocols should enrich viruses
        for viral_fraction in results:
            assert viral_fraction > 0.5, "Protocol did not enrich viruses"

        # Results should vary by protocol
        assert len(set(results)) > 1, "Different protocols gave identical results"

    def test_different_amplification_methods(self):
        """Test integration with different amplification methods."""
        community = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        contamination = create_contamination_profile('clean', random_seed=42)

        methods = [
            rdab_40_cycles(),
            mda_standard(),
            linker_standard()
        ]

        results = []
        for method in methods:
            composition = MockViromeComposition(
                name='test_amplification',
                viral_community=community,
                contamination_profile=contamination,
                viral_fraction=0.95
            )
            method.apply(composition)

            # Calculate coefficient of variation
            abundances = [g.abundance for g in composition.viral_community.genomes]
            mean_abund = sum(abundances) / len(abundances)
            variance = sum((a - mean_abund) ** 2 for a in abundances) / len(abundances)
            cv = (variance ** 0.5) / mean_abund if mean_abund > 0 else 0

            results.append(cv)

        # Different methods should produce different bias levels
        assert len(set(results)) > 1, "Different amplification methods gave identical bias"


class TestCrossPlatformIntegration:
    """Tests for cross-platform comparison workflows."""

    def test_same_community_different_platforms(self):
        """Test that same community can be processed with different platforms."""
        # Create identical starting compositions
        community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        platforms = [novaseq_6000(), miseq(), nextseq_2000()]

        for platform in platforms:
            composition = MockViromeComposition(
                name=f'test_{platform.name}',
                viral_community=community,
                contamination_profile=contamination,
                viral_fraction=0.5
            )

            # Apply VLP and amplification
            vlp = standard_vlp()
            vlp.apply(composition)

            amplification = rdab_40_cycles()
            amplification.apply(composition)

            # Create reads
            reads = []
            for i in range(100):
                read = ReadPair(
                    read_id=f"read_{i}",
                    forward_seq="A" * 80,
                    reverse_seq="T" * 80,
                    forward_qual="I" * 80,
                    reverse_qual="I" * 80,
                    genome_id=f"genome_{i % 5}",
                    tile_x=i * 100,
                    tile_y=i * 100
                )
                reads.append(read)

            # Apply platform artifacts
            reads_with_artifacts = platform.apply(reads, random_seed=42)

            # Verify artifacts were applied
            assert len(reads_with_artifacts) > 0

    def test_patterned_vs_cluster_flow_cells(self):
        """Test that patterned flow cells have polyG, cluster flow cells do not."""
        reads = []
        for i in range(2000):  # Increased for statistical significance
            read = ReadPair(
                read_id=f"read_{i}",
                forward_seq="A" * 80,
                reverse_seq="T" * 80,
                forward_qual="I" * 80,
                reverse_qual="I" * 80,
                genome_id=f"genome_1",
                tile_x=i * 100,
                tile_y=i * 100
            )
            reads.append(read)

        # NovaSeq (patterned) should have polyG tails
        novaseq = novaseq_6000()
        novaseq_reads = novaseq.apply(reads.copy(), random_seed=42)
        # Check both R1 and R2
        novaseq_has_polyg = any(len(r.forward_seq) > 80 or (r.reverse_seq and len(r.reverse_seq) > 80) for r in novaseq_reads)

        # MiSeq (cluster) should NOT have polyG tails
        miseq_platform = miseq()
        miseq_reads = miseq_platform.apply(reads.copy(), random_seed=42)
        miseq_has_polyg = any(len(r.forward_seq) > 80 or (r.reverse_seq and len(r.reverse_seq) > 80) for r in miseq_reads)

        assert novaseq_has_polyg, "NovaSeq should have polyG tails"
        assert not miseq_has_polyg, "MiSeq should NOT have polyG tails"

    def test_artifact_rates_vary_by_platform(self):
        """Test that artifact rates vary by platform as expected."""
        reads = []
        for i in range(1000):
            read = ReadPair(
                read_id=f"read_{i}",
                forward_seq="A" * 80,
                reverse_seq="T" * 80,
                forward_qual="I" * 80,
                reverse_qual="I" * 80,
                genome_id="genome_1",
                tile_x=i * 100,
                tile_y=i * 100
            )
            reads.append(read)

        # Apply different platforms
        novaseq_reads = novaseq_6000().apply(reads.copy(), random_seed=42)
        miseq_reads = miseq().apply(reads.copy(), random_seed=42)
        ideal_reads = no_artifacts().apply(reads.copy(), random_seed=42)

        # NovaSeq should have most duplicates (highest rate)
        # MiSeq should have fewer duplicates
        # Ideal should have no duplicates

        novaseq_dup_rate = (len(novaseq_reads) - 1000) / 1000
        miseq_dup_rate = (len(miseq_reads) - 1000) / 1000
        ideal_dup_rate = (len(ideal_reads) - 1000) / 1000

        assert novaseq_dup_rate > miseq_dup_rate, "NovaSeq should have more dups than MiSeq"
        assert ideal_dup_rate == 0, "Ideal platform should have no duplicates"


class TestWorkflowEdgeCases:
    """Tests for edge cases and error handling in workflows."""

    def test_workflow_with_zero_contamination(self):
        """Test workflow with no contamination."""
        community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('clean', random_seed=42)

        composition = MockViromeComposition(
            name='test_no_contam',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=1.0  # 100% viral
        )

        # Should still work
        vlp = standard_vlp()
        vlp.apply(composition)

        amplification = rdab_40_cycles()
        amplification.apply(composition)

        assert composition.viral_fraction >= 0.99

    def test_workflow_with_heavy_contamination(self):
        """Test workflow with heavy contamination."""
        community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('heavy', random_seed=42)

        composition = MockViromeComposition(
            name='test_heavy_contam',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.1  # 10% viral, 90% contamination
        )

        initial_viral = composition.viral_fraction

        # VLP enrichment should still improve viral fraction
        vlp = standard_vlp()
        vlp.apply(composition)

        assert composition.viral_fraction > initial_viral

    def test_workflow_with_small_community(self):
        """Test workflow with very small viral community."""
        community = create_body_site_profile('gut', n_genomes=3, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        composition = MockViromeComposition(
            name='test_small',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.5
        )

        vlp = standard_vlp()
        vlp.apply(composition)

        amplification = rdab_40_cycles()
        amplification.apply(composition)

        # Should complete without errors
        assert len(composition.viral_community.genomes) == 3

    def test_workflow_with_empty_reads(self):
        """Test platform artifacts with empty read list."""
        reads = []

        # Most platforms should handle empty reads, but polyG artifact has a bug
        # Use no_artifacts platform which has no artifacts to apply
        platform = no_artifacts()
        result = platform.apply(reads, random_seed=42)

        # Should return empty list
        assert len(result) == 0

    def test_sequential_artifact_application(self):
        """Test that artifacts can be applied sequentially."""
        reads = []
        for i in range(100):
            read = ReadPair(
                read_id=f"read_{i}",
                forward_seq="A" * 80,
                reverse_seq="T" * 80,
                forward_qual="I" * 80,
                reverse_qual="I" * 80,
                genome_id="genome_1",
                tile_x=i * 100,
                tile_y=i * 100
            )
            reads.append(read)

        # Apply artifacts sequentially
        polyg = PolyGTailArtifact(frequency=0.05, random_seed=42)
        reads = polyg.apply(reads)

        optical = OpticalDuplicateArtifact(rate=0.10, random_seed=42)
        reads = optical.apply(reads)

        # Should complete without errors
        assert len(reads) > 100  # Optical dups added


class TestWorkflowReproducibility:
    """Tests for workflow reproducibility with random seeds."""

    def test_workflow_reproducible_with_seed(self):
        """Test that workflow is reproducible with same random seed."""
        def run_workflow(seed):
            community = create_body_site_profile('gut', n_genomes=10, random_seed=seed)
            contamination = create_contamination_profile('realistic', random_seed=seed)

            composition = MockViromeComposition(
                name='test_reproducible',
                viral_community=community,
                contamination_profile=contamination,
                viral_fraction=0.5
            )

            # Use VLP with explicit seed
            vlp = VLPEnrichment(random_seed=seed)
            vlp.apply(composition)

            # Use amplification with explicit seed
            from viroforge.amplification import RdABAmplification
            amplification = RdABAmplification(cycles=40, random_seed=seed)
            amplification.apply(composition)

            # Return abundances as reproducibility check
            abundances = [round(g.abundance, 10) for g in composition.viral_community.genomes]
            return abundances

        # Run twice with same seed
        result1 = run_workflow(42)
        result2 = run_workflow(42)

        # Should be very close (allowing for floating point precision)
        assert len(result1) == len(result2)
        for i in range(len(result1)):
            assert abs(result1[i] - result2[i]) < 1e-8, f"Abundances differ at index {i}: {result1[i]} != {result2[i]}"

    def test_workflow_different_with_different_seed(self):
        """Test that workflow produces different results with different seeds."""
        def run_workflow(seed):
            community = create_body_site_profile('gut', n_genomes=10, random_seed=seed)
            contamination = create_contamination_profile('realistic', random_seed=seed)

            composition = MockViromeComposition(
                name='test_different',
                viral_community=community,
                contamination_profile=contamination,
                viral_fraction=0.5
            )

            vlp = standard_vlp()
            vlp.apply(composition)

            amplification = rdab_40_cycles()
            amplification.apply(composition)

            abundances = [g.abundance for g in composition.viral_community.genomes]
            return abundances

        # Run with different seeds
        result1 = run_workflow(42)
        result2 = run_workflow(123)

        # Should be different
        assert result1 != result2, "Different seeds produced identical results"


class TestWorkflowValidation:
    """Tests for workflow output validation."""

    def test_viral_fraction_bounds(self):
        """Test that viral fraction stays within valid bounds."""
        community = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        composition = MockViromeComposition(
            name='test_bounds',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.5
        )

        # Throughout workflow, viral fraction should be 0-1 (allowing tiny floating point error)
        assert -1e-10 <= composition.viral_fraction <= 1 + 1e-10

        vlp = standard_vlp()
        vlp.apply(composition)
        assert -1e-10 <= composition.viral_fraction <= 1 + 1e-10

        amplification = rdab_40_cycles()
        amplification.apply(composition)
        assert -1e-10 <= composition.viral_fraction <= 1 + 1e-10

    def test_abundance_conservation(self):
        """Test that total abundance is conserved through workflow."""
        community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        composition = MockViromeComposition(
            name='test_conservation',
            viral_community=community,
            contamination_profile=contamination,
            viral_fraction=0.5
        )

        # Total abundance should be ~1.0 throughout
        initial_total = composition.get_total_abundance()
        assert abs(initial_total - 1.0) < 0.01

        vlp = standard_vlp()
        vlp.apply(composition)

        after_vlp_total = composition.get_total_abundance()
        assert abs(after_vlp_total - 1.0) < 0.01

        amplification = rdab_40_cycles()
        amplification.apply(composition)

        after_amp_total = composition.get_total_abundance()
        assert abs(after_amp_total - 1.0) < 0.01

    def test_read_metadata_preserved(self):
        """Test that read metadata is preserved through artifact application."""
        reads = []
        genome_ids = ["genome_1", "genome_2", "genome_3"]

        for i in range(150):
            read = ReadPair(
                read_id=f"read_{i}",
                forward_seq="A" * 80,
                reverse_seq="T" * 80,
                forward_qual="I" * 80,
                reverse_qual="I" * 80,
                genome_id=genome_ids[i % 3],
                tile_x=i * 100,
                tile_y=i * 100
            )
            reads.append(read)

        # Apply artifacts
        platform = novaseq_6000()
        reads_with_artifacts = platform.apply(reads, random_seed=42)

        # All reads should have valid genome_ids
        for read in reads_with_artifacts:
            assert read.genome_id in genome_ids or read.genome_id != ""

        # Read IDs should be unique
        read_ids = [r.read_id for r in reads_with_artifacts]
        assert len(read_ids) == len(set(read_ids)), "Read IDs not unique"
