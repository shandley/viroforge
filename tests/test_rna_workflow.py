"""
Comprehensive tests for RNA virome workflow.

Tests for Phase 8.2 RNA workflow components:
- Reverse transcription with virus-type specific efficiency
- rRNA depletion (Ribo-Zero/RiboMinus)
- RNA degradation and fragmentation
- Complete RNA virome workflow integration

Author: ViroForge Development Team
Date: 2025-11-09
"""

import unittest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viroforge.workflows.rna_virome import (
    RNAViromeWorkflow,
    ReverseTranscription,
    RiboDepletion,
    RNADegradation,
    RNAVirusType,
    PrimerType,
    RiboDepleteMethod,
    infer_virus_type_from_taxonomy
)


class TestReverseTranscription(unittest.TestCase):
    """Test reverse transcription modeling."""

    def setUp(self):
        """Set up test fixtures."""
        self.rt = ReverseTranscription(random_seed=42)

    def test_rt_initialization(self):
        """Test RT object initialization."""
        self.assertEqual(self.rt.primer_type, PrimerType.RANDOM_HEXAMER)
        self.assertEqual(self.rt.base_efficiency, 0.80)
        self.assertIsNotNone(self.rt.rng)

    def test_rt_efficiency_by_virus_type(self):
        """Test RT efficiency varies by virus type."""
        # ssRNA+ should have higher efficiency than ssRNA- and dsRNA
        eff_positive = self.rt.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
        eff_negative = self.rt.get_efficiency(RNAVirusType.SSRNA_NEGATIVE)
        eff_dsrna = self.rt.get_efficiency(RNAVirusType.DSRNA)

        # Check efficiency ranges
        self.assertGreaterEqual(eff_positive, 0.70)
        self.assertLessEqual(eff_positive, 0.90)

        self.assertGreaterEqual(eff_negative, 0.50)
        self.assertLessEqual(eff_negative, 0.70)

        self.assertGreaterEqual(eff_dsrna, 0.40)
        self.assertLessEqual(eff_dsrna, 0.80)

        # ssRNA+ should generally be more efficient than ssRNA-
        # (not guaranteed due to randomness, but test average)
        efficiencies_pos = [
            self.rt.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
            for _ in range(100)
        ]
        efficiencies_neg = [
            self.rt.get_efficiency(RNAVirusType.SSRNA_NEGATIVE)
            for _ in range(100)
        ]

        self.assertGreater(np.mean(efficiencies_pos), np.mean(efficiencies_neg))

    def test_rt_length_penalty(self):
        """Test RT efficiency decreases for very long genomes."""
        short_eff = self.rt.get_efficiency(
            RNAVirusType.SSRNA_POSITIVE,
            genome_length=5000
        )
        long_eff = self.rt.get_efficiency(
            RNAVirusType.SSRNA_POSITIVE,
            genome_length=20000
        )

        # Long genomes should have lower efficiency
        self.assertLess(long_eff, short_eff)

    def test_rt_template_switching(self):
        """Test template switching creates chimeric sequences."""
        sequences = [
            SeqRecord(Seq("A" * 1000), id=f"seq_{i}")
            for i in range(100)
        ]

        modified, n_chimeric = self.rt.apply_template_switching(sequences)

        # Should create some chimeric sequences (around 2%)
        self.assertGreater(n_chimeric, 0)
        self.assertLess(n_chimeric, 10)  # With 100 seqs, expect ~2
        self.assertEqual(len(modified), len(sequences))

    def test_rt_truncation(self):
        """Test RT truncation creates shorter sequences."""
        sequences = [
            SeqRecord(Seq("A" * 1000), id=f"seq_{i}")
            for i in range(100)
        ]

        modified, truncation_stats = self.rt.apply_truncation(sequences)

        # Should truncate some sequences (around 15%)
        total_truncated = (
            truncation_stats['5_prime_truncated'] +
            truncation_stats['3_prime_truncated'] -
            truncation_stats['both_truncated']
        )
        self.assertGreater(total_truncated, 0)
        self.assertLess(total_truncated, 30)  # Expect ~15/100

        # Check that some sequences are shorter
        original_lengths = [len(seq.seq) for seq in sequences]
        modified_lengths = [len(seq.seq) for seq in modified]

        self.assertLess(np.mean(modified_lengths), np.mean(original_lengths))

    def test_rt_primer_type_affects_efficiency(self):
        """Test different primer types affect RT efficiency."""
        rt_hexamer = ReverseTranscription(
            primer_type=PrimerType.RANDOM_HEXAMER,
            random_seed=42
        )
        rt_octamer = ReverseTranscription(
            primer_type=PrimerType.RANDOM_OCTAMER,
            random_seed=42
        )
        rt_specific = ReverseTranscription(
            primer_type=PrimerType.SPECIFIC,
            random_seed=42
        )

        # Get average efficiencies
        hexamer_effs = [
            rt_hexamer.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
            for _ in range(50)
        ]
        octamer_effs = [
            rt_octamer.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
            for _ in range(50)
        ]
        specific_effs = [
            rt_specific.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
            for _ in range(50)
        ]

        # Specific should be best, octamer better than hexamer
        self.assertGreater(np.mean(specific_effs), np.mean(hexamer_effs))
        self.assertGreater(np.mean(octamer_effs), np.mean(hexamer_effs))


class TestRiboDepletion(unittest.TestCase):
    """Test rRNA depletion modeling."""

    def setUp(self):
        """Set up test fixtures."""
        self.ribo = RiboDepletion(random_seed=42)

    def test_ribo_initialization(self):
        """Test RiboDepletion initialization."""
        self.assertEqual(self.ribo.method, RiboDepleteMethod.RIBO_ZERO)
        self.assertIsNotNone(self.ribo.efficiency)
        self.assertGreaterEqual(self.ribo.efficiency, 0.90)
        self.assertLessEqual(self.ribo.efficiency, 0.95)

    def test_ribo_depletion_efficiency(self):
        """Test rRNA depletion reduces rRNA abundance."""
        rrna_before = 0.90  # 90% rRNA

        rrna_after, stats = self.ribo.apply_depletion(rrna_before)

        # Should dramatically reduce rRNA
        self.assertLess(rrna_after, 0.15)  # Should be <15%
        self.assertGreater(rrna_after, 0.0)  # But not 0

        # Check stats
        self.assertIn('rrna_removed', stats)
        self.assertIn('viral_enrichment', stats)
        self.assertGreater(stats['viral_enrichment'], 10.0)  # Should enrich >10x

    def test_ribo_viral_enrichment_calculation(self):
        """Test viral enrichment is calculated correctly."""
        rrna_before = 0.90  # 90% rRNA, 10% viral

        rrna_after, stats = self.ribo.apply_depletion(rrna_before)

        # Before: 10% viral out of 100%
        # After: 10% viral out of ~15% total (assuming 92.5% rRNA removed)
        # Enrichment should be ~(10/15) / (10/100) = 6.67x minimum

        viral_enrichment = stats['viral_enrichment']
        self.assertGreater(viral_enrichment, 6.0)
        self.assertLess(viral_enrichment, 25.0)  # Maximum reasonable enrichment

    def test_ribo_method_comparison(self):
        """Test different rRNA depletion methods have different efficiencies."""
        ribo_zero = RiboDepletion(method=RiboDepleteMethod.RIBO_ZERO, random_seed=42)
        ribominus = RiboDepletion(method=RiboDepleteMethod.RIBOMINUS, random_seed=43)
        none = RiboDepletion(method=RiboDepleteMethod.NONE, random_seed=44)

        rrna_before = 0.90

        rrna_after_zero, _ = ribo_zero.apply_depletion(rrna_before)
        rrna_after_minus, _ = ribominus.apply_depletion(rrna_before)
        rrna_after_none, _ = none.apply_depletion(rrna_before)

        # Ribo-Zero should be most efficient
        self.assertLess(rrna_after_zero, rrna_after_minus)
        self.assertLess(rrna_after_minus, rrna_after_none)

        # None should not remove any rRNA
        self.assertEqual(rrna_after_none, rrna_before)

    def test_ribo_extreme_contamination(self):
        """Test rRNA depletion works with extreme contamination."""
        rrna_before = 0.95  # 95% rRNA (heavy contamination)

        rrna_after, stats = self.ribo.apply_depletion(rrna_before)

        # Should still reduce substantially
        self.assertLess(rrna_after, 0.20)

        # Viral enrichment should be dramatic
        self.assertGreater(stats['viral_enrichment'], 15.0)


class TestRNADegradation(unittest.TestCase):
    """Test RNA degradation modeling."""

    def setUp(self):
        """Set up test fixtures."""
        self.rna_deg = RNADegradation(degradation_rate=0.30, random_seed=42)

    def test_rna_degradation_initialization(self):
        """Test RNADegradation initialization."""
        self.assertEqual(self.rna_deg.degradation_rate, 0.30)
        self.assertEqual(self.rna_deg.rnase_contamination, 0.10)
        self.assertIsNotNone(self.rna_deg.rng)

    def test_rna_degradation_fragments_sequences(self):
        """Test RNA degradation fragments some sequences."""
        sequences = [
            SeqRecord(Seq("A" * 2000), id=f"seq_{i}")
            for i in range(100)
        ]

        modified, degradation_stats = self.rna_deg.apply_degradation(sequences)

        # Should fragment some sequences (~30%)
        self.assertGreater(degradation_stats['n_fragmented'], 0)
        self.assertLess(degradation_stats['n_fragmented'], 50)

        # Total sequences should increase (fragments)
        self.assertGreater(len(modified), len(sequences))

    def test_rna_degradation_creates_smaller_fragments(self):
        """Test degradation creates smaller fragments."""
        sequences = [
            SeqRecord(Seq("A" * 3000), id=f"seq_{i}")
            for i in range(50)
        ]

        modified, degradation_stats = self.rna_deg.apply_degradation(sequences)

        # Average fragment size should be less than original
        if degradation_stats['n_fragmented'] > 0:
            mean_fragment_size = degradation_stats['mean_fragment_size']
            self.assertLess(mean_fragment_size, 3000)
            self.assertGreater(mean_fragment_size, 100)  # Minimum kept

    def test_rna_degradation_rate_affects_fragmentation(self):
        """Test higher degradation rate causes more fragmentation."""
        sequences = [
            SeqRecord(Seq("A" * 2000), id=f"seq_{i}")
            for i in range(100)
        ]

        low_deg = RNADegradation(degradation_rate=0.10, random_seed=42)
        high_deg = RNADegradation(degradation_rate=0.50, random_seed=43)

        _, low_stats = low_deg.apply_degradation(sequences.copy())
        _, high_stats = high_deg.apply_degradation(sequences.copy())

        # Higher rate should fragment more sequences
        self.assertGreater(
            high_stats['n_fragmented'],
            low_stats['n_fragmented']
        )


class TestRNAViromeWorkflow(unittest.TestCase):
    """Test complete RNA virome workflow."""

    def setUp(self):
        """Set up test fixtures."""
        self.workflow = RNAViromeWorkflow(random_seed=42)

        # Create test sequences
        self.sequences = [
            SeqRecord(Seq("A" * 1000), id=f"genome_{i:03d}")
            for i in range(50)
        ]

        # Create virus type mapping
        self.virus_types = {
            seq.id: RNAVirusType.SSRNA_POSITIVE
            for seq in self.sequences
        }

    def test_workflow_initialization(self):
        """Test RNA workflow initializes correctly."""
        self.assertIsNotNone(self.workflow.reverse_transcription)
        self.assertIsNotNone(self.workflow.ribo_depletion)
        self.assertIsNotNone(self.workflow.rna_degradation)

    def test_workflow_applies_all_steps(self):
        """Test workflow applies all steps in correct order."""
        processed, stats = self.workflow.apply(
            sequences=self.sequences,
            virus_types=self.virus_types,
            rrna_abundance_before=0.90
        )

        # Check all stats are present
        self.assertIn('degradation_stats', stats)
        self.assertIn('rt_stats', stats)
        self.assertIn('ribo_depletion_stats', stats)
        self.assertIn('overall_recovery', stats)

        # Some sequences should survive workflow
        self.assertGreater(len(processed), 0)
        self.assertLessEqual(len(processed), len(self.sequences))

    def test_workflow_recovery_rate(self):
        """Test overall recovery rate is reasonable."""
        processed, stats = self.workflow.apply(
            sequences=self.sequences,
            virus_types=self.virus_types,
            rrna_abundance_before=0.90
        )

        # Recovery should be 40-90% (degradation + RT efficiency)
        recovery = stats['overall_recovery']
        self.assertGreater(recovery, 0.30)
        self.assertLess(recovery, 1.0)

    def test_workflow_viral_enrichment(self):
        """Test workflow produces strong viral enrichment."""
        processed, stats = self.workflow.apply(
            sequences=self.sequences,
            virus_types=self.virus_types,
            rrna_abundance_before=0.90
        )

        # Viral enrichment should be >10x
        viral_enrichment = stats['ribo_depletion_stats']['viral_enrichment']
        self.assertGreater(viral_enrichment, 10.0)

    def test_workflow_with_different_virus_types(self):
        """Test workflow handles mixed virus types."""
        # Create mixed virus types
        mixed_types = {}
        for i, seq in enumerate(self.sequences):
            if i % 3 == 0:
                mixed_types[seq.id] = RNAVirusType.SSRNA_POSITIVE
            elif i % 3 == 1:
                mixed_types[seq.id] = RNAVirusType.SSRNA_NEGATIVE
            else:
                mixed_types[seq.id] = RNAVirusType.DSRNA

        processed, stats = self.workflow.apply(
            sequences=self.sequences,
            virus_types=mixed_types,
            rrna_abundance_before=0.90
        )

        # Should handle mixed types without error
        self.assertGreater(len(processed), 0)
        self.assertIn('mean_efficiency', stats['rt_stats'])

    def test_workflow_reproducibility(self):
        """Test workflow is reproducible with same seed."""
        workflow1 = RNAViromeWorkflow(random_seed=12345)
        workflow2 = RNAViromeWorkflow(random_seed=12345)

        sequences1 = [
            SeqRecord(Seq("A" * 1000), id=f"genome_{i:03d}")
            for i in range(30)
        ]
        sequences2 = [
            SeqRecord(Seq("A" * 1000), id=f"genome_{i:03d}")
            for i in range(30)
        ]

        virus_types = {seq.id: RNAVirusType.SSRNA_POSITIVE for seq in sequences1}

        processed1, stats1 = workflow1.apply(sequences1, virus_types, 0.90)
        processed2, stats2 = workflow2.apply(sequences2, virus_types, 0.90)

        # Should get same number of sequences
        self.assertEqual(len(processed1), len(processed2))

        # Stats should be very similar
        self.assertAlmostEqual(
            stats1['overall_recovery'],
            stats2['overall_recovery'],
            places=2
        )


class TestInferVirusType(unittest.TestCase):
    """Test virus type inference from taxonomy."""

    def test_infer_ssrna_positive(self):
        """Test inference of ssRNA+ viruses."""
        # Picornaviridae (ssRNA+)
        taxonomy = {
            'family': 'Picornaviridae',
            'genus': 'Enterovirus',
            'species': 'Enterovirus A',
            'genome_name': 'Enterovirus A strain'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.SSRNA_POSITIVE)

        # Caliciviridae (norovirus, ssRNA+)
        taxonomy = {
            'family': 'Caliciviridae',
            'genus': 'Norovirus',
            'species': 'Norwalk virus'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.SSRNA_POSITIVE)

    def test_infer_ssrna_negative(self):
        """Test inference of ssRNA- viruses."""
        # Orthomyxoviridae (influenza, ssRNA-)
        taxonomy = {
            'family': 'Orthomyxoviridae',
            'genus': 'Alphainfluenzavirus',
            'species': 'Influenza A virus'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.SSRNA_NEGATIVE)

        # Pneumoviridae (RSV, ssRNA-)
        taxonomy = {
            'family': 'Pneumoviridae',
            'genus': 'Orthopneumovirus',
            'species': 'Human orthopneumovirus'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.SSRNA_NEGATIVE)

    def test_infer_dsrna(self):
        """Test inference of dsRNA viruses."""
        # Reoviridae (dsRNA)
        taxonomy = {
            'family': 'Reoviridae',
            'genus': 'Orthoreovirus',
            'species': 'Mammalian orthoreovirus'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.DSRNA)

        # Sedoreoviridae (rotavirus, dsRNA)
        taxonomy = {
            'family': 'Sedoreoviridae',
            'genus': 'Rotavirus',
            'species': 'Rotavirus A'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        self.assertEqual(virus_type, RNAVirusType.DSRNA)

    def test_infer_unknown_family_defaults_positive(self):
        """Test unknown families default to ssRNA+."""
        taxonomy = {
            'family': 'UnknownViridae',
            'genus': 'Unknown',
            'species': 'Unknown virus'
        }

        virus_type = infer_virus_type_from_taxonomy(taxonomy)
        # Should default to ssRNA+ (most common)
        self.assertEqual(virus_type, RNAVirusType.SSRNA_POSITIVE)


if __name__ == '__main__':
    unittest.main()
