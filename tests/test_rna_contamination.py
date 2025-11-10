"""
Comprehensive tests for RNA virome contamination profiles.

Tests for Phase 8.2 RNA contamination functions:
- Host RNA contamination (post-Ribo-Zero)
- Bacterial RNA contamination
- RNA-specific contamination profiles

Author: ViroForge Development Team
Date: 2025-11-09
"""

import unittest
import numpy as np

from viroforge.core.contamination import (
    ContaminationProfile,
    ContaminantType,
    add_host_rna_contamination,
    add_bacterial_rna_contamination,
    create_rna_contamination_profile
)


class TestHostRNAContamination(unittest.TestCase):
    """Test host RNA contamination functions."""

    def setUp(self):
        """Set up test fixtures."""
        self.profile = ContaminationProfile(name="test_rna_profile")

    def test_host_rna_basic_addition(self):
        """Test adding host RNA contamination."""
        add_host_rna_contamination(
            self.profile,
            host_organism="human",
            abundance_pct_before_depletion=90.0,
            abundance_pct_after_depletion=10.0,
            random_seed=42
        )

        # Should have contaminants
        self.assertGreater(len(self.profile), 0)

        # Total abundance should be ~10% (after depletion)
        total_abundance = self.profile.get_total_abundance()
        self.assertAlmostEqual(total_abundance, 0.10, places=1)

    def test_host_rna_rrna_and_mrna(self):
        """Test host RNA includes both rRNA and mRNA."""
        add_host_rna_contamination(
            self.profile,
            host_organism="human",
            abundance_pct_after_depletion=10.0,
            rrna_fraction=0.95,
            random_seed=42
        )

        # Check for rRNA contaminants
        rrna_contaminants = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.RRNA
        ]
        self.assertGreater(len(rrna_contaminants), 0)

        # Check for mRNA/transcript contaminants (stored as HOST_DNA type)
        mrna_contaminants = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.HOST_DNA
            and 'mrna' in c.genome_id
        ]
        self.assertGreater(len(mrna_contaminants), 0)

    def test_host_rna_rrna_dominates(self):
        """Test rRNA is dominant fraction of host RNA."""
        add_host_rna_contamination(
            self.profile,
            host_organism="human",
            abundance_pct_after_depletion=10.0,
            rrna_fraction=0.95,  # 95% rRNA
            random_seed=42
        )

        # Calculate rRNA vs mRNA abundances
        rrna_abundance = sum(
            c.abundance for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.RRNA
        )
        mrna_abundance = sum(
            c.abundance for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.HOST_DNA
        )

        # rRNA should be dominant
        self.assertGreater(rrna_abundance, mrna_abundance * 10)

    def test_host_rna_different_organisms(self):
        """Test different host organisms."""
        profile_human = ContaminationProfile()
        add_host_rna_contamination(
            profile_human,
            host_organism="human",
            random_seed=42
        )

        profile_mouse = ContaminationProfile()
        add_host_rna_contamination(
            profile_mouse,
            host_organism="mouse",
            random_seed=42
        )

        # Both should work
        self.assertGreater(len(profile_human), 0)
        self.assertGreater(len(profile_mouse), 0)

        # Check organism names
        human_organisms = set(c.organism for c in profile_human.contaminants)
        mouse_organisms = set(c.organism for c in profile_mouse.contaminants)

        self.assertIn('Homo sapiens', human_organisms)
        self.assertIn('Mus musculus', mouse_organisms)

    def test_host_rna_rrna_lengths(self):
        """Test rRNA sequences have correct lengths."""
        add_host_rna_contamination(
            self.profile,
            host_organism="human",
            random_seed=42
        )

        rrna_contaminants = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.RRNA
        ]

        # Check lengths are reasonable for rRNA
        for rrna in rrna_contaminants:
            # rRNA should be 100-5000 bp range
            self.assertGreater(rrna.length, 100)
            self.assertLess(rrna.length, 5000)

            # Check for expected rRNA types in description
            self.assertTrue(
                any(rna_type in rrna.description
                    for rna_type in ['18S', '28S', '5.8S', '5S'])
            )


class TestBacterialRNAContamination(unittest.TestCase):
    """Test bacterial RNA contamination functions."""

    def setUp(self):
        """Set up test fixtures."""
        self.profile = ContaminationProfile(name="test_bacterial_rna")

    def test_bacterial_rna_basic_addition(self):
        """Test adding bacterial RNA contamination."""
        add_bacterial_rna_contamination(
            self.profile,
            abundance_pct=5.0,
            microbiome_type="gut",
            random_seed=42
        )

        # Should have contaminants
        self.assertGreater(len(self.profile), 0)

        # Total abundance should be ~5%
        total_abundance = self.profile.get_total_abundance()
        self.assertAlmostEqual(total_abundance, 0.05, places=2)

    def test_bacterial_rna_rrna_and_mrna(self):
        """Test bacterial RNA includes both rRNA and mRNA."""
        add_bacterial_rna_contamination(
            self.profile,
            abundance_pct=5.0,
            rrna_fraction=0.80,
            random_seed=42
        )

        # Check for bacterial rRNA
        bacterial_rrna = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.RRNA
        ]
        self.assertGreater(len(bacterial_rrna), 0)

        # Check for bacterial mRNA (stored as REAGENT_BACTERIA type)
        bacterial_mrna = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.REAGENT_BACTERIA
        ]
        self.assertGreater(len(bacterial_mrna), 0)

    def test_bacterial_rna_microbiome_types(self):
        """Test different microbiome types have different taxa."""
        profile_gut = ContaminationProfile()
        add_bacterial_rna_contamination(
            profile_gut,
            microbiome_type="gut",
            random_seed=42
        )

        profile_oral = ContaminationProfile()
        add_bacterial_rna_contamination(
            profile_oral,
            microbiome_type="oral",
            random_seed=42
        )

        profile_skin = ContaminationProfile()
        add_bacterial_rna_contamination(
            profile_skin,
            microbiome_type="skin",
            random_seed=42
        )

        # All should work
        self.assertGreater(len(profile_gut), 0)
        self.assertGreater(len(profile_oral), 0)
        self.assertGreater(len(profile_skin), 0)

        # Check for expected taxa
        gut_organisms = set(c.organism for c in profile_gut.contaminants)
        oral_organisms = set(c.organism for c in profile_oral.contaminants)
        skin_organisms = set(c.organism for c in profile_skin.contaminants)

        # Gut should have Bacteroides
        self.assertTrue(any('Bacteroides' in org for org in gut_organisms))

        # Oral should have Streptococcus
        self.assertTrue(any('Streptococcus' in org for org in oral_organisms))

        # Skin should have Staphylococcus
        self.assertTrue(any('Staphylococcus' in org for org in skin_organisms))

    def test_bacterial_rna_16s_23s_rrna(self):
        """Test bacterial rRNA includes 16S and 23S."""
        add_bacterial_rna_contamination(
            self.profile,
            abundance_pct=5.0,
            random_seed=42
        )

        rrna_contaminants = [
            c for c in self.profile.contaminants
            if c.contaminant_type == ContaminantType.RRNA
        ]

        # Check for 16S and 23S in descriptions
        has_16s = any('16S' in c.description for c in rrna_contaminants)
        has_23s = any('23S' in c.description for c in rrna_contaminants)

        self.assertTrue(has_16s)
        self.assertTrue(has_23s)

        # Check lengths are appropriate
        for rrna in rrna_contaminants:
            # Bacterial rRNA should be in expected range
            self.assertGreater(rrna.length, 1000)
            self.assertLess(rrna.length, 4000)


class TestRNAContaminationProfiles(unittest.TestCase):
    """Test pre-defined RNA contamination profiles."""

    def test_create_clean_rna_profile(self):
        """Test creating clean RNA contamination profile."""
        profile = create_rna_contamination_profile(
            profile_type='clean',
            random_seed=42
        )

        # Should have contaminants
        self.assertGreater(len(profile), 0)

        # Total contamination should be low (~7-8%)
        total = profile.get_total_abundance()
        self.assertLess(total, 0.10)  # <10%
        self.assertGreater(total, 0.05)  # >5%

    def test_create_realistic_rna_profile(self):
        """Test creating realistic RNA contamination profile."""
        profile = create_rna_contamination_profile(
            profile_type='realistic',
            random_seed=42
        )

        # Should have contaminants
        self.assertGreater(len(profile), 0)

        # Total contamination should be moderate (~15-16%)
        total = profile.get_total_abundance()
        self.assertLess(total, 0.20)  # <20%
        self.assertGreater(total, 0.10)  # >10%

    def test_create_heavy_rna_profile(self):
        """Test creating heavy RNA contamination profile."""
        profile = create_rna_contamination_profile(
            profile_type='heavy',
            random_seed=42
        )

        # Should have contaminants
        self.assertGreater(len(profile), 0)

        # Total contamination should be high (~30-35%)
        total = profile.get_total_abundance()
        self.assertGreater(total, 0.25)  # >25%

    def test_create_failed_rna_profile(self):
        """Test creating failed RNA profile (Ribo-Zero failed)."""
        profile = create_rna_contamination_profile(
            profile_type='failed',
            random_seed=42
        )

        # Should have VERY high contamination (~95-96%)
        total = profile.get_total_abundance()
        self.assertGreater(total, 0.90)  # >90%

    def test_rna_profile_with_ribo_depletion(self):
        """Test RNA profile with vs without Ribo-Zero."""
        profile_with = create_rna_contamination_profile(
            profile_type='realistic',
            ribo_depletion=True,
            random_seed=42
        )

        profile_without = create_rna_contamination_profile(
            profile_type='realistic',
            ribo_depletion=False,
            random_seed=42
        )

        # Profile without Ribo-Zero should have MUCH more contamination
        total_with = profile_with.get_total_abundance()
        total_without = profile_without.get_total_abundance()

        self.assertGreater(total_without, total_with * 5)  # >5x more

    def test_rna_profile_with_without_microbiome(self):
        """Test RNA profile with vs without microbiome contamination."""
        profile_with_microbiome = create_rna_contamination_profile(
            profile_type='realistic',
            microbiome_rich=True,
            random_seed=42
        )

        profile_without_microbiome = create_rna_contamination_profile(
            profile_type='realistic',
            microbiome_rich=False,
            random_seed=42
        )

        # Profile with microbiome should have more contaminants
        self.assertGreater(
            len(profile_with_microbiome),
            len(profile_without_microbiome)
        )

    def test_rna_profile_composition(self):
        """Test RNA profile has expected composition."""
        profile = create_rna_contamination_profile(
            profile_type='realistic',
            ribo_depletion=True,
            microbiome_rich=True,
            random_seed=42
        )

        # Check for expected contamination types
        type_abundances = profile.get_abundance_by_type()

        # Should have rRNA
        self.assertIn(ContaminantType.RRNA, type_abundances)

        # Should have host DNA/RNA (HOST_DNA type reused)
        self.assertIn(ContaminantType.HOST_DNA, type_abundances)

        # Should have reagent bacteria
        self.assertIn(ContaminantType.REAGENT_BACTERIA, type_abundances)

        # Should have PhiX
        self.assertIn(ContaminantType.PHIX, type_abundances)

    def test_rna_profile_override_parameters(self):
        """Test overriding default RNA profile parameters."""
        profile = create_rna_contamination_profile(
            profile_type='realistic',
            host_rna_after=20.0,  # Override: 20% instead of 10%
            bacterial_rna=10.0,   # Override: 10% instead of 5%
            random_seed=42
        )

        # Total should be higher due to overrides
        total = profile.get_total_abundance()
        self.assertGreater(total, 0.25)  # Should be >25%

    def test_rna_vs_dna_contamination_comparison(self):
        """Test RNA contamination is dramatically different from DNA."""
        from viroforge.core.contamination import create_contamination_profile

        rna_profile = create_rna_contamination_profile(
            profile_type='realistic',
            ribo_depletion=False,  # No Ribo-Zero for fair comparison
            random_seed=42
        )

        dna_profile = create_contamination_profile(
            profile_type='realistic',
            random_seed=42
        )

        rna_total = rna_profile.get_total_abundance()
        dna_total = dna_profile.get_total_abundance()

        # RNA should have MUCH more contamination (without Ribo-Zero)
        self.assertGreater(rna_total, dna_total * 10)  # >10x more

        # After Ribo-Zero, RNA should be more comparable
        rna_profile_depleted = create_rna_contamination_profile(
            profile_type='realistic',
            ribo_depletion=True,
            random_seed=42
        )

        rna_depleted_total = rna_profile_depleted.get_total_abundance()

        # After depletion, RNA contamination should be similar order of magnitude
        self.assertLess(rna_depleted_total, dna_total * 5)


class TestRNAContaminationIntegration(unittest.TestCase):
    """Test RNA contamination integration with profiles."""

    def test_combined_host_and_bacterial_rna(self):
        """Test combining host and bacterial RNA contamination."""
        profile = ContaminationProfile(name="combined_rna")

        add_host_rna_contamination(
            profile,
            abundance_pct_after_depletion=10.0,
            random_seed=42
        )

        add_bacterial_rna_contamination(
            profile,
            abundance_pct=5.0,
            random_seed=42
        )

        # Total should be ~15%
        total = profile.get_total_abundance()
        self.assertAlmostEqual(total, 0.15, places=1)

    def test_rna_contamination_export(self):
        """Test exporting RNA contamination profile."""
        profile = create_rna_contamination_profile(
            profile_type='realistic',
            random_seed=42
        )

        # Get contamination table
        table = profile.get_contamination_table()

        # Should have rows
        self.assertGreater(len(table), 0)

        # Should have expected columns
        expected_columns = ['genome_id', 'type', 'organism', 'source', 'abundance']
        for col in expected_columns:
            self.assertIn(col, table.columns)

    def test_rna_contamination_summary_stats(self):
        """Test RNA contamination summary statistics."""
        profile = create_rna_contamination_profile(
            profile_type='realistic',
            random_seed=42
        )

        stats = profile.get_summary_stats()

        # Check expected keys
        self.assertIn('n_contaminants', stats)
        self.assertIn('total_abundance', stats)
        self.assertIn('by_type', stats)

        # Should have multiple contamination types
        self.assertGreater(len(stats['by_type']), 2)


if __name__ == '__main__':
    unittest.main()
