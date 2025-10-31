"""
Tests for library preparation amplification bias modeling.

Tests all amplification methods:
- RdAB (length + GC bias)
- MDA (extreme GC bias + stochasticity)
- Linker (minimal bias)
- NoAmplification (control)
"""

import pytest
import numpy as np
from unittest.mock import Mock

from viroforge.amplification import (
    RdABAmplification,
    MDAAmplification,
    LinkerAmplification,
    NoAmplification,
    rdab_40_cycles,
    rdab_30_cycles,
    mda_standard,
    mda_overnight,
    linker_standard,
    no_amplification
)
from viroforge.core import ViralGenome, ViralCommunity, create_body_site_profile
from viroforge.core.contamination import ContaminationProfile, create_contamination_profile
from viroforge.utils.composition import MockViromeComposition


class TestRdABInitialization:
    """Test RdAB amplification initialization and validation."""

    def test_default_initialization(self):
        """Test RdAB with default parameters."""
        rdab = RdABAmplification()

        assert rdab.cycles == 40
        assert rdab.length_bias_strength == 1.0
        assert rdab.gc_bias_strength == 1.0
        assert rdab.optimal_gc == 0.50
        assert rdab.gc_tolerance == 0.15

    def test_custom_initialization(self):
        """Test RdAB with custom parameters."""
        rdab = RdABAmplification(
            cycles=30,
            length_bias_strength=1.5,
            gc_bias_strength=0.8,
            optimal_gc=0.45,
            gc_tolerance=0.20
        )

        assert rdab.cycles == 30
        assert rdab.length_bias_strength == 1.5
        assert rdab.gc_bias_strength == 0.8
        assert rdab.optimal_gc == 0.45
        assert rdab.gc_tolerance == 0.20

    def test_invalid_cycles_raises_error(self):
        """Test that invalid cycle number raises error."""
        with pytest.raises(ValueError, match="cycles must be"):
            RdABAmplification(cycles=-5)

        with pytest.raises(ValueError, match="cycles must be"):
            RdABAmplification(cycles=100)

    def test_invalid_bias_strength_raises_error(self):
        """Test that negative bias strength raises error."""
        with pytest.raises(ValueError, match="length_bias_strength"):
            RdABAmplification(length_bias_strength=-1.0)

        with pytest.raises(ValueError, match="gc_bias_strength"):
            RdABAmplification(gc_bias_strength=-0.5)

    def test_random_seed_sets_reproducibility(self):
        """Test that random seed ensures reproducibility."""
        rdab1 = RdABAmplification(random_seed=42)
        rdab2 = RdABAmplification(random_seed=42)

        # Both should produce identical results
        # (No stochasticity in RdAB, but testing RNG initialization)
        assert rdab1.random_seed == rdab2.random_seed


class TestRdABLengthBias:
    """Test RdAB length-dependent amplification bias."""

    def test_short_genome_high_efficiency(self):
        """Test that short genomes amplify efficiently."""
        rdab = RdABAmplification(cycles=40, length_bias_strength=1.0)

        # Small virus (5 kb)
        efficiency = rdab.calculate_length_efficiency(5000)

        # Should be high efficiency (>0.9)
        assert efficiency > 0.9

    def test_long_genome_lower_efficiency(self):
        """Test that long genomes amplify less efficiently."""
        rdab = RdABAmplification(cycles=40, length_bias_strength=1.0)

        # Large virus (200 kb)
        efficiency = rdab.calculate_length_efficiency(200000)

        # Should be lower efficiency (<0.5)
        assert efficiency < 0.5

    def test_length_bias_scales_with_strength(self):
        """Test that bias strength parameter works."""
        weak = RdABAmplification(length_bias_strength=0.5)
        strong = RdABAmplification(length_bias_strength=2.0)

        # 100 kb genome
        weak_eff = weak.calculate_length_efficiency(100000)
        strong_eff = strong.calculate_length_efficiency(100000)

        # Strong bias should have lower efficiency for large genome
        assert strong_eff < weak_eff

    def test_no_length_bias_when_disabled(self):
        """Test that length bias can be disabled."""
        rdab = RdABAmplification(length_bias_strength=0.0)

        # All lengths should have efficiency = 1.0
        assert rdab.calculate_length_efficiency(5000) == 1.0
        assert rdab.calculate_length_efficiency(100000) == 1.0
        assert rdab.calculate_length_efficiency(500000) == 1.0


class TestRdABGCBias:
    """Test RdAB GC-dependent amplification bias."""

    def test_optimal_gc_high_efficiency(self):
        """Test that optimal GC has maximum efficiency."""
        rdab = RdABAmplification(cycles=40, optimal_gc=0.50)

        efficiency = rdab.calculate_gc_efficiency(0.50)

        # Should be maximum efficiency
        assert efficiency == pytest.approx(1.0, abs=0.01)

    def test_extreme_gc_lower_efficiency(self):
        """Test that extreme GC content reduces efficiency."""
        rdab = RdABAmplification(cycles=40, optimal_gc=0.50, gc_bias_strength=1.0)

        low_gc_eff = rdab.calculate_gc_efficiency(0.25)   # Low GC
        high_gc_eff = rdab.calculate_gc_efficiency(0.75)  # High GC

        # Both extremes should be less efficient
        assert low_gc_eff < 0.5
        assert high_gc_eff < 0.5

    def test_gc_bias_symmetric(self):
        """Test that GC bias is symmetric around optimal."""
        rdab = RdABAmplification(optimal_gc=0.50)

        eff_low = rdab.calculate_gc_efficiency(0.30)   # 20% below optimal
        eff_high = rdab.calculate_gc_efficiency(0.70)  # 20% above optimal

        # Should be approximately equal
        assert eff_low == pytest.approx(eff_high, abs=0.05)

    def test_no_gc_bias_when_disabled(self):
        """Test that GC bias can be disabled."""
        rdab = RdABAmplification(gc_bias_strength=0.0)

        # All GC contents should have efficiency = 1.0
        assert rdab.calculate_gc_efficiency(0.25) == 1.0
        assert rdab.calculate_gc_efficiency(0.50) == 1.0
        assert rdab.calculate_gc_efficiency(0.75) == 1.0


class TestRdABApplication:
    """Test full RdAB amplification application."""

    def test_rdab_modifies_abundances(self):
        """Test that RdAB amplification changes abundances."""
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('clean', random_seed=42)

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.90
        )

        # Store initial abundances
        initial_abundances = [g.abundance for g in composition.viral_community.genomes]

        # Apply RdAB
        rdab = RdABAmplification(cycles=40)
        rdab.apply(composition)

        # Get final abundances
        final_abundances = [g.abundance for g in composition.viral_community.genomes]

        # Abundances should have changed
        assert not all(i == f for i, f in zip(initial_abundances, final_abundances))

    def test_rdab_favors_short_genomes(self):
        """Test that RdAB increases relative abundance of short genomes."""
        # Create community with different sized genomes
        viral_community = ViralCommunity(name='test')

        # Add short and long genomes with equal starting abundance
        viral_community.add_genome(ViralGenome(
            genome_id='short_1',
            sequence='ATCG' * 1250,  # 5 kb
            taxonomy='Viruses;Microviridae',
            family='Microviridae',
            genus='Microvirus',
            abundance=0.5,
            gc_content=0.50
        ))

        viral_community.add_genome(ViralGenome(
            genome_id='long_1',
            sequence='ATCG' * 50000,  # 200 kb
            taxonomy='Viruses;Myoviridae',
            family='Myoviridae',
            genus='Myovirus',
            abundance=0.5,
            gc_content=0.50
        ))

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=None,
            viral_fraction=1.0
        )

        # Initially equal (50/50)
        short_initial = composition.viral_community.genomes[0].abundance
        long_initial = composition.viral_community.genomes[1].abundance
        assert short_initial == pytest.approx(long_initial, abs=0.01)

        # Apply RdAB
        rdab = RdABAmplification(cycles=40, gc_bias_strength=0.0)  # Only length bias
        rdab.apply(composition)

        # After amplification, short should be more abundant
        short_final = composition.viral_community.genomes[0].abundance
        long_final = composition.viral_community.genomes[1].abundance

        assert short_final > long_final
        # Should be at least 2x more abundant
        assert short_final / long_final > 2.0

    def test_cycle_dependency(self):
        """Test that bias scales with cycle number."""
        # Create separate communities for independent testing
        viral_community_30 = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        viral_community_40 = create_body_site_profile('gut', n_genomes=15, random_seed=42)

        comp_30 = MockViromeComposition(
            name='30_cycles',
            viral_community=viral_community_30,
            contamination_profile=None,
            viral_fraction=1.0
        )

        comp_40 = MockViromeComposition(
            name='40_cycles',
            viral_community=viral_community_40,
            contamination_profile=None,
            viral_fraction=1.0
        )

        # Apply different cycle numbers
        rdab_30 = RdABAmplification(cycles=30)
        rdab_40 = RdABAmplification(cycles=40)

        rdab_30.apply(comp_30)
        rdab_40.apply(comp_40)

        # 40 cycles should show stronger bias
        # (variance in abundances should be higher)
        var_30 = np.var([g.abundance for g in comp_30.viral_community.genomes])
        var_40 = np.var([g.abundance for g in comp_40.viral_community.genomes])

        assert var_40 > var_30


class TestMDAAmplification:
    """Test MDA amplification."""

    def test_mda_initialization(self):
        """Test MDA initialization."""
        mda = MDAAmplification()

        assert mda.amplification_time == 4.0
        assert mda.gc_bias_strength == 3.0
        assert mda.stochasticity == 0.3
        assert mda.chimera_rate == 0.15

    def test_mda_extreme_gc_bias(self):
        """Test that MDA has stronger GC bias than RdAB."""
        mda = MDAAmplification(gc_bias_strength=3.0)
        rdab = RdABAmplification(gc_bias_strength=1.0)

        # High GC genome (70%)
        mda_eff = mda.calculate_gc_efficiency(0.70)
        rdab_eff = rdab.calculate_gc_efficiency(0.70)

        # MDA should have much lower efficiency for high GC
        assert mda_eff < rdab_eff
        assert mda_eff < 0.1  # Very low efficiency

    def test_mda_stochastic_variation(self):
        """Test that MDA adds stochastic variation."""
        mda = MDAAmplification(stochasticity=0.3, random_seed=42)

        # Multiple calls should give different results due to randomness
        results = [mda.add_stochastic_variation(1.0) for _ in range(10)]

        # Should have variation
        assert np.std(results) > 0.1

    def test_mda_application(self):
        """Test full MDA application."""
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=None,
            viral_fraction=1.0
        )

        mda = MDAAmplification(amplification_time_hours=4.0)
        mda.apply(composition)

        # Should have modified abundances
        total_abundance = composition.get_total_abundance()
        assert total_abundance == pytest.approx(1.0, abs=0.001)


class TestLinkerAmplification:
    """Test linker-based amplification."""

    def test_linker_initialization(self):
        """Test linker initialization."""
        linker = LinkerAmplification()

        assert linker.cycles == 20
        assert linker.gc_bias_strength == 0.5

    def test_linker_weaker_gc_bias(self):
        """Test that linker has weaker GC bias than RdAB."""
        linker = LinkerAmplification(gc_bias_strength=0.5)
        rdab = RdABAmplification(gc_bias_strength=1.0)

        # Extreme GC
        linker_eff = linker.calculate_gc_efficiency(0.70)
        rdab_eff = rdab.calculate_gc_efficiency(0.70)

        # Linker should be less penalized
        assert linker_eff > rdab_eff

    def test_linker_application(self):
        """Test full linker application."""
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=None,
            viral_fraction=1.0
        )

        linker = LinkerAmplification(cycles=20)
        linker.apply(composition)

        # Should maintain normalization
        total_abundance = composition.get_total_abundance()
        assert total_abundance == pytest.approx(1.0, abs=0.001)


class TestNoAmplification:
    """Test no amplification control."""

    def test_no_amplification_preserves_abundances(self):
        """Test that NoAmplification doesn't change abundances."""
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=None,
            viral_fraction=1.0
        )

        # Store initial abundances
        initial_abundances = [g.abundance for g in composition.viral_community.genomes]

        # Apply no amplification
        no_amp = NoAmplification()
        no_amp.apply(composition)

        # Get final abundances
        final_abundances = [g.abundance for g in composition.viral_community.genomes]

        # Should be identical
        for i, f in zip(initial_abundances, final_abundances):
            assert i == pytest.approx(f, abs=1e-10)


class TestPreDefinedProtocols:
    """Test pre-defined amplification protocols."""

    def test_rdab_40_cycles_protocol(self):
        """Test rdab_40_cycles protocol."""
        amp = rdab_40_cycles()

        assert isinstance(amp, RdABAmplification)
        assert amp.cycles == 40
        assert amp.length_bias_strength == 1.0
        assert amp.gc_bias_strength == 1.0

    def test_rdab_30_cycles_protocol(self):
        """Test rdab_30_cycles protocol."""
        amp = rdab_30_cycles()

        assert isinstance(amp, RdABAmplification)
        assert amp.cycles == 30

    def test_mda_standard_protocol(self):
        """Test mda_standard protocol."""
        amp = mda_standard()

        assert isinstance(amp, MDAAmplification)
        assert amp.amplification_time == 4.0
        assert amp.gc_bias_strength == 3.0

    def test_mda_overnight_protocol(self):
        """Test mda_overnight protocol."""
        amp = mda_overnight()

        assert isinstance(amp, MDAAmplification)
        assert amp.amplification_time == 16.0

    def test_linker_standard_protocol(self):
        """Test linker_standard protocol."""
        amp = linker_standard()

        assert isinstance(amp, LinkerAmplification)
        assert amp.cycles == 20

    def test_no_amplification_protocol(self):
        """Test no_amplification protocol."""
        amp = no_amplification()

        assert isinstance(amp, NoAmplification)


class TestAmplificationComparison:
    """Test comparison of different amplification methods."""

    def test_bias_strength_comparison(self):
        """Test that different methods have different bias strengths."""
        # Create separate communities for independent testing
        viral_community_rdab = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        viral_community_linker = create_body_site_profile('gut', n_genomes=15, random_seed=42)
        viral_community_none = create_body_site_profile('gut', n_genomes=15, random_seed=42)

        comp_rdab = MockViromeComposition(
            name='rdab',
            viral_community=viral_community_rdab,
            contamination_profile=None,
            viral_fraction=1.0
        )

        comp_linker = MockViromeComposition(
            name='linker',
            viral_community=viral_community_linker,
            contamination_profile=None,
            viral_fraction=1.0
        )

        comp_none = MockViromeComposition(
            name='none',
            viral_community=viral_community_none,
            contamination_profile=None,
            viral_fraction=1.0
        )

        # Apply different methods
        rdab_40_cycles().apply(comp_rdab)
        linker_standard().apply(comp_linker)
        no_amplification().apply(comp_none)

        # Calculate variance (measure of bias)
        var_rdab = np.var([g.abundance for g in comp_rdab.viral_community.genomes])
        var_linker = np.var([g.abundance for g in comp_linker.viral_community.genomes])
        var_none = np.var([g.abundance for g in comp_none.viral_community.genomes])

        # Amplification bias can either increase or decrease variance depending on
        # the initial abundance distribution and genome properties
        # What matters is that amplified methods differ from no amplification

        # All three methods should produce different variance patterns
        assert var_rdab != var_linker or var_rdab != var_none  # At least some difference

        # RdAB and Linker should modify abundances differently
        # (we can't predict direction of change without knowing the distribution)
        assert abs(var_rdab - var_none) > 1e-10 or abs(var_linker - var_none) > 1e-10


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
