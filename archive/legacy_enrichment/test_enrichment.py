"""
Tests for VLP enrichment modeling.

Tests all components of VLP enrichment:
- Size-based filtration
- Nuclease treatment
- Family enrichment factors
- Stochastic variation
- Complete enrichment workflow
- Pre-defined protocols
"""

import pytest
import numpy as np
from unittest.mock import Mock

from viroforge.enrichment import (
    VLPEnrichment,
    standard_vlp,
    iron_chloride_vlp,
    ultracentrifuge_vlp,
    no_enrichment,
    syringe_filter_vlp,
    FAMILY_ENRICHMENT_FACTORS
)
from viroforge.core import ViralGenome, ViralCommunity, create_body_site_profile
from viroforge.core.contamination import (
    ContaminationProfile,
    ContaminantGenome,
    ContaminantType
)
from viroforge.utils import MockViromeComposition
from Bio.Seq import Seq


class TestVLPEnrichmentInit:
    """Test VLPEnrichment initialization and validation."""

    def test_default_initialization(self):
        """Test default VLP enrichment initialization."""
        vlp = VLPEnrichment()

        assert vlp.filtration_method == 'tangential_flow'
        assert vlp.filtration_cutoff_um == 0.2
        assert vlp.nuclease_treatment is True
        assert vlp.nuclease_efficiency == 0.95

    def test_custom_initialization(self):
        """Test custom VLP enrichment parameters."""
        vlp = VLPEnrichment(
            filtration_cutoff_um=0.22,
            nuclease_efficiency=0.98,
            stochastic_variation=0.1
        )

        assert vlp.filtration_cutoff_um == 0.22
        assert vlp.nuclease_efficiency == 0.98
        assert vlp.stochastic_variation == 0.1

    def test_invalid_filtration_cutoff_raises_error(self):
        """Test that invalid filtration cutoff raises error."""
        with pytest.raises(ValueError, match="filtration_cutoff_um"):
            VLPEnrichment(filtration_cutoff_um=-0.1)

        with pytest.raises(ValueError, match="filtration_cutoff_um"):
            VLPEnrichment(filtration_cutoff_um=2.0)

    def test_invalid_nuclease_efficiency_raises_error(self):
        """Test that invalid nuclease efficiency raises error."""
        with pytest.raises(ValueError, match="nuclease_efficiency"):
            VLPEnrichment(nuclease_efficiency=-0.1)

        with pytest.raises(ValueError, match="nuclease_efficiency"):
            VLPEnrichment(nuclease_efficiency=1.5)

    def test_random_seed_sets_reproducibility(self):
        """Test that random seed ensures reproducibility."""
        vlp1 = VLPEnrichment(random_seed=42)
        vlp2 = VLPEnrichment(random_seed=42)

        # Should produce same stochastic values
        val1 = vlp1.add_stochastic_variation(1.0)
        val2 = vlp2.add_stochastic_variation(1.0)

        assert val1 == val2


class TestVirusSizeEstimation:
    """Test virus size estimation from genome length."""

    def test_small_virus_size(self):
        """Test size estimation for small viruses (Microviridae)."""
        vlp = VLPEnrichment()

        # PhiX174: 5.4 kb
        size = vlp.estimate_virus_size_nm(5400)

        # Should be ~25-30 nm
        assert 20 < size < 35

    def test_medium_virus_size(self):
        """Test size estimation for medium viruses."""
        vlp = VLPEnrichment()

        # Typical phage: 50 kb
        size = vlp.estimate_virus_size_nm(50000)

        # Should be ~50-70 nm (based on power-law formula)
        assert 50 < size < 70

    def test_large_virus_size(self):
        """Test size estimation for large/jumbo viruses."""
        vlp = VLPEnrichment()

        # T4 phage: ~169 kb
        size = vlp.estimate_virus_size_nm(169000)

        # Should be ~85-100 nm (based on power-law formula)
        assert 85 < size < 100

    def test_jumbo_virus_size(self):
        """Test size estimation for jumbo phages."""
        vlp = VLPEnrichment()

        # Jumbo phage: 300 kb
        size = vlp.estimate_virus_size_nm(300000)

        # Should be ~105-115 nm (based on power-law formula)
        assert 105 < size < 115


class TestSizeBasedRetention:
    """Test size-based filtration retention."""

    def test_small_virus_high_retention(self):
        """Test that small viruses are highly retained."""
        vlp = VLPEnrichment(filtration_cutoff_um=0.2, size_retention_curve='sigmoid')

        # Small virus (5 kb, ~27 nm)
        retention = vlp.calculate_size_retention(5000)

        # Should be highly retained (>94%)
        assert retention > 0.94

    def test_medium_virus_good_retention(self):
        """Test that medium viruses are well retained."""
        vlp = VLPEnrichment(filtration_cutoff_um=0.2, size_retention_curve='sigmoid')

        # Medium virus (50 kb)
        retention = vlp.calculate_size_retention(50000)

        # Should be well retained (>80%)
        assert retention > 0.80

    def test_large_virus_partial_retention(self):
        """Test that large viruses are partially retained."""
        vlp = VLPEnrichment(filtration_cutoff_um=0.2, size_retention_curve='sigmoid')

        # Large virus (200 kb)
        retention = vlp.calculate_size_retention(200000)

        # Should be partially retained (50-90%)
        assert 0.50 < retention < 0.90

    def test_jumbo_virus_low_retention(self):
        """Test that jumbo viruses are poorly retained with small cutoff."""
        vlp = VLPEnrichment(filtration_cutoff_um=0.05, size_retention_curve='sigmoid')

        # Jumbo virus (400 kb, ~122 nm)
        retention = vlp.calculate_size_retention(400000)

        # Should be poorly retained with 50 nm cutoff (<50%)
        assert retention < 0.50

    def test_step_function_sharp_cutoff(self):
        """Test step function has sharp cutoff."""
        vlp = VLPEnrichment(filtration_cutoff_um=0.08, size_retention_curve='step')

        # Just below cutoff (100 kb → ~75 nm)
        retention_below = vlp.calculate_size_retention(100000)

        # Just above cutoff (120 kb → ~80 nm)
        retention_above = vlp.calculate_size_retention(120000)

        # Step function should have binary retention
        assert retention_below == 1.0 or retention_below > 0.9
        assert retention_above == 0.0 or retention_above < 0.1

    def test_no_filtration_full_retention(self):
        """Test that no filtration retains everything."""
        vlp = VLPEnrichment(filtration_method='none')

        retention_small = vlp.calculate_size_retention(5000)
        retention_jumbo = vlp.calculate_size_retention(500000)

        assert retention_small == 1.0
        assert retention_jumbo == 1.0


class TestFamilyEnrichment:
    """Test family-specific enrichment factors."""

    def test_microviridae_highly_enriched(self):
        """Test that Microviridae are highly enriched."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_family_enrichment('Microviridae')

        # Should be 2.5x (from ViromeQC)
        assert factor == 2.5

    def test_inoviridae_depleted(self):
        """Test that Inoviridae are depleted."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_family_enrichment('Inoviridae')

        # Should be 0.3x (filamentous, excluded by filter)
        assert factor == 0.3

    def test_siphoviridae_moderate_enrichment(self):
        """Test that Siphoviridae have moderate enrichment."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_family_enrichment('Siphoviridae')

        # Should be 1.2x
        assert factor == 1.2

    def test_unknown_family_baseline(self):
        """Test that unknown families get baseline enrichment."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_family_enrichment('UnknownViridae')

        # Should be 1.0 (baseline)
        assert factor == 1.0

    def test_all_enrichment_factors_positive(self):
        """Test that all enrichment factors are positive."""
        for family, factor in FAMILY_ENRICHMENT_FACTORS.items():
            assert factor > 0, f"{family} has non-positive factor: {factor}"


class TestStabilityFactor:
    """Test capsid stability factor calculation."""

    def test_low_gc_lower_stability(self):
        """Test that low GC gives lower stability."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_stability_factor(30.0)

        # Should be <1.0 (less stable)
        assert factor < 1.0

    def test_high_gc_higher_stability(self):
        """Test that high GC gives higher stability."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_stability_factor(70.0)

        # Should be >1.0 (more stable)
        assert factor > 1.0

    def test_medium_gc_baseline_stability(self):
        """Test that medium GC gives baseline stability."""
        vlp = VLPEnrichment()

        factor = vlp.calculate_stability_factor(50.0)

        # Should be ~1.0 (baseline)
        assert 0.99 < factor < 1.01

    def test_stability_factor_bounded(self):
        """Test that stability factor is bounded."""
        vlp = VLPEnrichment()

        # Extreme GC values should be clipped
        factor_low = vlp.calculate_stability_factor(0.0)
        factor_high = vlp.calculate_stability_factor(100.0)

        assert 0.85 <= factor_low <= 1.15
        assert 0.85 <= factor_high <= 1.15


class TestNucleaseTreatment:
    """Test nuclease treatment effects."""

    def test_viral_genome_protected(self):
        """Test that encapsidated viral genomes are protected."""
        vlp = VLPEnrichment(nuclease_efficiency=0.95)

        retention = vlp.calculate_nuclease_retention(is_viral=True, is_encapsidated=True)

        # Should be fully retained
        assert retention == 1.0

    def test_free_dna_removed(self):
        """Test that free DNA is removed by nuclease."""
        vlp = VLPEnrichment(nuclease_efficiency=0.95)

        retention = vlp.calculate_nuclease_retention(is_viral=False, is_encapsidated=False)

        # Should be 5% retained (95% removed)
        assert retention == pytest.approx(0.05)

    def test_no_nuclease_full_retention(self):
        """Test that no nuclease treatment retains everything."""
        vlp = VLPEnrichment(nuclease_treatment=False)

        retention_viral = vlp.calculate_nuclease_retention(True, True)
        retention_free = vlp.calculate_nuclease_retention(False, False)

        assert retention_viral == 1.0
        assert retention_free == 1.0

    def test_variable_nuclease_efficiency(self):
        """Test different nuclease efficiencies."""
        # 90% efficiency
        vlp_90 = VLPEnrichment(nuclease_efficiency=0.90)
        retention_90 = vlp_90.calculate_nuclease_retention(False, False)
        assert retention_90 == pytest.approx(0.10)

        # 98% efficiency
        vlp_98 = VLPEnrichment(nuclease_efficiency=0.98)
        retention_98 = vlp_98.calculate_nuclease_retention(False, False)
        assert retention_98 == pytest.approx(0.02)


class TestStochasticVariation:
    """Test stochastic variation modeling."""

    def test_no_variation_returns_input(self):
        """Test that zero variation returns input unchanged."""
        vlp = VLPEnrichment(stochastic_variation=0.0)

        retention = vlp.add_stochastic_variation(0.8)

        assert retention == 0.8

    def test_variation_changes_value(self):
        """Test that variation changes the value."""
        vlp = VLPEnrichment(stochastic_variation=0.2, random_seed=42)

        retention = vlp.add_stochastic_variation(0.8)

        # Should be different from input
        assert retention != 0.8

    def test_variation_bounded_to_valid_range(self):
        """Test that variation doesn't exceed valid range."""
        vlp = VLPEnrichment(stochastic_variation=0.5, random_seed=42)

        # Test many values
        for _ in range(100):
            retention = vlp.add_stochastic_variation(0.5)
            assert 0.0 <= retention <= 1.0

    def test_high_variation_more_spread(self):
        """Test that higher variation produces more spread."""
        vlp_low = VLPEnrichment(stochastic_variation=0.1, random_seed=42)
        vlp_high = VLPEnrichment(stochastic_variation=0.5, random_seed=43)

        values_low = [vlp_low.add_stochastic_variation(0.8) for _ in range(100)]
        values_high = [vlp_high.add_stochastic_variation(0.8) for _ in range(100)]

        # Higher variation should have larger standard deviation
        std_low = np.std(values_low)
        std_high = np.std(values_high)

        assert std_high > std_low


class TestVLPEnrichmentApplication:
    """Test complete VLP enrichment application."""

    def test_vlp_increases_viral_fraction(self):
        """Test that VLP enrichment increases viral fraction."""
        # Create composition with 50% viral (bulk metagenome)
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = ContaminationProfile(name='test_contam')

        # Add significant host DNA contamination
        for i in range(10):
            contam = ContaminantGenome(
                genome_id=f"host_{i}",
                sequence=Seq("ATCGATCG" * 1250),  # 10kb
                organism="Human",
                contaminant_type=ContaminantType.HOST_DNA,
                abundance=0.05  # 5% each, 50% total
            )
            contamination.add_contaminant(contam)

        composition = MockViromeComposition(
            name='test_vlp',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.50
        )

        initial_viral_frac = composition.viral_fraction

        # Apply VLP enrichment
        vlp = VLPEnrichment(
            filtration_cutoff_um=0.2,
            nuclease_efficiency=0.95,
            random_seed=42
        )
        vlp.apply(composition)

        final_viral_frac = composition.viral_fraction

        # Viral fraction should increase significantly
        assert final_viral_frac > initial_viral_frac
        assert final_viral_frac > 0.85  # Should be >85% viral

    def test_no_enrichment_preserves_abundances(self):
        """Test that no enrichment preserves relative abundances."""
        from viroforge.core.contamination import ContaminantGenome, ContaminantType

        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = ContaminationProfile(name='test')

        # Add some contaminants (host and bacterial DNA)
        contamination.add_contaminant(ContaminantGenome(
            genome_id='host_1',
            organism='Homo sapiens',
            sequence='ATCGATCG' * 625,  # 5000 bp
            abundance=0.3,
            contaminant_type=ContaminantType.HOST_DNA
        ))
        contamination.add_contaminant(ContaminantGenome(
            genome_id='bacteria_1',
            organism='E. coli',
            sequence='GCGCGCGC' * 625,  # 5000 bp
            abundance=0.2,
            contaminant_type=ContaminantType.REAGENT_BACTERIA
        ))

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.50
        )

        initial_viral_frac = composition.viral_fraction

        # Apply no enrichment
        bulk = no_enrichment()
        bulk.apply(composition)

        final_viral_frac = composition.viral_fraction

        # Should be approximately unchanged (within rounding error from stochastic variation)
        assert abs(final_viral_frac - initial_viral_frac) < 0.05  # Allow 5% tolerance for variation

    def test_enrichment_reduces_host_dna(self):
        """Test that enrichment dramatically reduces host DNA."""
        viral_community = create_body_site_profile('gut', n_genomes=5, random_seed=42)

        contamination = ContaminationProfile(name='test')
        # Add host DNA
        host = ContaminantGenome(
            genome_id="host",
            sequence=Seq("ATCG" * 2500),
            organism="Human",
            contaminant_type=ContaminantType.HOST_DNA,
            abundance=0.3  # 30% host DNA
        )
        contamination.add_contaminant(host)

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.70
        )

        initial_host_abundance = host.abundance

        # Apply VLP enrichment
        vlp = standard_vlp()
        vlp.apply(composition)

        final_host_abundance = host.abundance

        # Host DNA should be dramatically reduced
        assert final_host_abundance < initial_host_abundance * 0.1  # <10% of original


class TestPreDefinedProtocols:
    """Test pre-defined VLP enrichment protocols."""

    def test_standard_vlp_configuration(self):
        """Test standard VLP protocol configuration."""
        vlp = standard_vlp()

        assert vlp.filtration_method == 'tangential_flow'
        assert vlp.filtration_cutoff_um == 0.2
        assert vlp.nuclease_treatment is True
        assert vlp.nuclease_efficiency == 0.95

    def test_iron_chloride_vlp_configuration(self):
        """Test iron chloride protocol configuration."""
        vlp = iron_chloride_vlp()

        assert vlp.filtration_cutoff_um == 0.2
        assert vlp.nuclease_efficiency == 0.98  # Higher with FeCl3

    def test_ultracentrifuge_vlp_configuration(self):
        """Test ultracentrifugation protocol configuration."""
        vlp = ultracentrifuge_vlp()

        assert vlp.filtration_method == 'ultracentrifugation'
        assert vlp.nuclease_efficiency == 0.90

    def test_syringe_filter_vlp_configuration(self):
        """Test syringe filtration protocol configuration."""
        vlp = syringe_filter_vlp()

        assert vlp.filtration_method == 'syringe'
        assert vlp.size_retention_curve == 'step'  # Sharp cutoff

    def test_no_enrichment_configuration(self):
        """Test no enrichment protocol configuration."""
        bulk = no_enrichment()

        assert bulk.filtration_method == 'none'
        assert bulk.nuclease_treatment is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
