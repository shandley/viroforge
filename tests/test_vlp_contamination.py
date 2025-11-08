"""
Unit tests for VLP contamination reduction integration

Tests the integration between VLP enrichment protocols and
contamination reduction, validating literature-based parameters.
"""

import pytest
import numpy as np
from viroforge.core.contamination import (
    create_contamination_profile,
    ContaminationProfile,
    ContaminantType
)
from viroforge.enrichment.vlp import (
    VLPEnrichment,
    VLPProtocol
)


class TestContaminationReduction:
    """Test suite for VLP contamination reduction"""

    def test_tangential_flow_reduces_contamination(self):
        """Test that tangential flow filtration significantly reduces contamination"""
        # Create realistic contamination profile
        profile = create_contamination_profile('realistic', random_seed=42)
        initial_contamination = profile.get_total_abundance()

        # Apply VLP enrichment
        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        reduced_profile, stats = vlp.apply_contamination_reduction(profile)
        final_contamination = reduced_profile.get_total_abundance()

        # Should reduce contamination by >85%
        reduction_factor = stats['overall_reduction_factor']
        assert reduction_factor > 0.85, f"Expected >85% reduction, got {reduction_factor*100:.1f}%"

        # Final contamination should be much lower
        assert final_contamination < initial_contamination * 0.15

    def test_no_vlp_preserves_contamination(self):
        """Test that no VLP enrichment leaves contamination unchanged"""
        profile = create_contamination_profile('realistic', random_seed=42)
        initial_contamination = profile.get_total_abundance()

        # No VLP enrichment
        bulk = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        reduced_profile, stats = bulk.apply_contamination_reduction(profile)
        final_contamination = reduced_profile.get_total_abundance()

        # Should be unchanged
        assert stats['overall_reduction_factor'] == pytest.approx(0.0, abs=1e-6)
        assert final_contamination == pytest.approx(initial_contamination, rel=1e-6)

    def test_host_dna_nuclease_dependent(self):
        """Test that host DNA reduction depends on nuclease treatment"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # With nuclease (TFF)
        vlp_nuclease = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        _, stats_nuclease = vlp_nuclease.apply_contamination_reduction(profile)

        # Without nuclease (no VLP)
        vlp_no_nuclease = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        _, stats_no_nuclease = vlp_no_nuclease.apply_contamination_reduction(profile)

        # Host DNA should be dramatically reduced with nuclease
        host_dna_reduction_nuclease = stats_nuclease['reduction_by_type']['host_dna']['reduction_factor']
        host_dna_reduction_no_nuclease = stats_no_nuclease['reduction_by_type']['host_dna']['reduction_factor']

        assert host_dna_reduction_nuclease > 0.85  # >85% with nuclease
        assert host_dna_reduction_no_nuclease < 0.05  # <5% without

    def test_bacteria_filtration_dependent(self):
        """Test that bacterial contamination is primarily removed by filtration"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # TFF (0.2 μm) - should remove bacteria
        vlp_tff = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        _, stats_tff = vlp_tff.apply_contamination_reduction(profile)

        # No filtration - bacteria remain
        vlp_no_filt = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        _, stats_no_filt = vlp_no_filt.apply_contamination_reduction(profile)

        bacteria_reduction_tff = stats_tff['reduction_by_type']['reagent_bacteria']['reduction_factor']
        bacteria_reduction_no_filt = stats_no_filt['reduction_by_type']['reagent_bacteria']['reduction_factor']

        # Bacteria should be >95% removed by filtration
        assert bacteria_reduction_tff > 0.95
        # Without filtration, minimal removal
        assert bacteria_reduction_no_filt < 0.05

    def test_phix_treated_like_small_virus(self):
        """Test that PhiX (encapsidated) is retained like small viruses"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # VLP enrichment
        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        _, stats = vlp.apply_contamination_reduction(profile)

        # PhiX should have moderate loss (10-60% removal), not complete removal
        phix_reduction = stats['reduction_by_type']['phix']['reduction_factor']
        assert 0.1 < phix_reduction < 0.6, f"PhiX reduction ({phix_reduction*100:.1f}%) outside expected range"

        # PhiX should not be removed as efficiently as free DNA
        host_dna_reduction = stats['reduction_by_type']['host_dna']['reduction_factor']
        assert phix_reduction < host_dna_reduction * 0.6  # PhiX protected by capsid

    def test_different_protocols_different_efficiency(self):
        """Test that different VLP protocols have different efficiencies"""
        profile = create_contamination_profile('realistic', random_seed=42)

        protocols = [
            ('TFF', VLPProtocol.tangential_flow_standard()),
            ('Syringe', VLPProtocol.syringe_filter_standard()),
            ('UC', VLPProtocol.ultracentrifugation()),
        ]

        reduction_factors = {}
        for name, protocol_config in protocols:
            vlp = VLPEnrichment(protocol=protocol_config, random_seed=42)
            _, stats = vlp.apply_contamination_reduction(profile)
            reduction_factors[name] = stats['overall_reduction_factor']

        # All should reduce contamination
        for name, reduction in reduction_factors.items():
            assert reduction > 0.75, f"{name} reduction too low: {reduction*100:.1f}%"

        # TFF should be most efficient (highest nuclease efficiency)
        assert reduction_factors['TFF'] >= reduction_factors['Syringe']

    def test_contamination_level_affects_absolute_reduction(self):
        """Test that higher initial contamination leads to higher absolute reduction"""
        # Clean, realistic, and heavy profiles
        profiles = {
            'clean': create_contamination_profile('clean', random_seed=42),
            'realistic': create_contamination_profile('realistic', random_seed=42),
            'heavy': create_contamination_profile('heavy', random_seed=42)
        }

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )

        absolute_reductions = {}
        for name, profile in profiles.items():
            initial = profile.get_total_abundance()
            reduced_profile, _ = vlp.apply_contamination_reduction(profile)
            final = reduced_profile.get_total_abundance()
            absolute_reductions[name] = initial - final

        # Higher contamination → larger absolute reduction
        assert absolute_reductions['heavy'] > absolute_reductions['realistic']
        assert absolute_reductions['realistic'] > absolute_reductions['clean']

    def test_reduction_stats_complete(self):
        """Test that reduction statistics are complete and correct"""
        profile = create_contamination_profile('realistic', random_seed=42)

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        _, stats = vlp.apply_contamination_reduction(profile)

        # Check required keys
        required_keys = [
            'protocol',
            'original_total_contamination',
            'reduced_total_contamination',
            'overall_reduction_factor',
            'reduction_by_type',
            'mean_removal_by_type'
        ]

        for key in required_keys:
            assert key in stats, f"Missing key in stats: {key}"

        # Check contaminant types
        expected_types = ['host_dna', 'rrna', 'reagent_bacteria', 'phix']
        for ctype in expected_types:
            assert ctype in stats['reduction_by_type']
            assert ctype in stats['mean_removal_by_type']

        # Check reduction_by_type structure
        for ctype, ctype_stats in stats['reduction_by_type'].items():
            assert 'original_abundance' in ctype_stats
            assert 'reduced_abundance' in ctype_stats
            assert 'reduction_factor' in ctype_stats
            assert 'reduction_pct' in ctype_stats

    def test_stochastic_variation_present(self):
        """Test that multiple runs with different seeds produce variation"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # Run with different seeds
        reduction_factors = []
        for seed in [42, 123, 456, 789, 101112]:
            vlp = VLPEnrichment(
                protocol=VLPProtocol.tangential_flow_standard(),
                random_seed=seed
            )
            _, stats = vlp.apply_contamination_reduction(profile)
            reduction_factors.append(stats['overall_reduction_factor'])

        # Should have some variation (realistic biological/technical variation)
        # Note: Variation is modest because stochasticity is controlled (5-15% CV per contaminant)
        std_reduction = np.std(reduction_factors)
        assert 0.001 < std_reduction < 0.20, f"Variation too {'high' if std_reduction > 0.20 else 'low'}: {std_reduction:.3f}"

        # Check that not all values are identical
        assert len(set(reduction_factors)) > 1, "No variation detected across seeds"

    def test_vlp_vs_bulk_fold_difference(self):
        """Test that VLP shows substantial fold-reduction vs bulk"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # VLP
        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        vlp_profile, _ = vlp.apply_contamination_reduction(profile)
        vlp_contamination = vlp_profile.get_total_abundance()

        # Bulk
        bulk = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        bulk_profile, _ = bulk.apply_contamination_reduction(profile)
        bulk_contamination = bulk_profile.get_total_abundance()

        # Calculate fold difference
        fold_difference = bulk_contamination / max(vlp_contamination, 1e-10)

        # VLP should reduce contamination by >5-fold
        assert fold_difference > 5.0, f"Fold difference too low: {fold_difference:.1f}x"


class TestContaminationIntegration:
    """Test integration with existing contamination module"""

    def test_profile_structure_preserved(self):
        """Test that contamination profile structure is preserved"""
        profile = create_contamination_profile('realistic', random_seed=42)
        initial_n_contaminants = len(profile)

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        reduced_profile, _ = vlp.apply_contamination_reduction(profile)

        # Same number of contaminants (just reduced abundances)
        assert len(reduced_profile) == initial_n_contaminants

        # All contaminant types preserved
        initial_types = set(c.contaminant_type for c in profile.contaminants)
        reduced_types = set(c.contaminant_type for c in reduced_profile.contaminants)
        assert initial_types == reduced_types

    def test_reduced_profile_has_correct_name(self):
        """Test that reduced profile has appropriate name"""
        profile = create_contamination_profile('realistic', random_seed=42)

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        reduced_profile, _ = vlp.apply_contamination_reduction(profile)

        # Name should indicate VLP reduction
        assert 'vlp_reduced' in reduced_profile.name.lower()

    def test_contaminant_metadata_preserved(self):
        """Test that contaminant metadata (organism, source, etc.) is preserved"""
        profile = create_contamination_profile('realistic', random_seed=42)

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        reduced_profile, _ = vlp.apply_contamination_reduction(profile)

        # Check first contaminant
        original_cont = profile.contaminants[0]
        reduced_cont = reduced_profile.contaminants[0]

        assert original_cont.genome_id == reduced_cont.genome_id
        assert original_cont.organism == reduced_cont.organism
        assert original_cont.source == reduced_cont.source
        assert original_cont.contaminant_type == reduced_cont.contaminant_type


class TestLiteratureValidation:
    """Validate against literature-reported values"""

    def test_nuclease_efficiency_range(self):
        """Test that nuclease treatment removes >85% free DNA (Thurber et al. 2009)"""
        profile = create_contamination_profile('realistic', random_seed=42)

        protocols = [
            VLPProtocol.tangential_flow_standard(),
            VLPProtocol.syringe_filter_standard(),
            VLPProtocol.ultracentrifugation(),
        ]

        for protocol_config in protocols:
            if protocol_config.nuclease_treatment:
                vlp = VLPEnrichment(protocol=protocol_config, random_seed=42)
                _, stats = vlp.apply_contamination_reduction(profile)

                host_dna_removal = stats['reduction_by_type']['host_dna']['reduction_factor']
                # Literature: >85% removal
                assert host_dna_removal > 0.85, f"Nuclease efficiency too low: {host_dna_removal*100:.1f}%"

    def test_bacterial_filtration_efficiency(self):
        """Test that 0.2 μm filters remove >90% bacteria (Shkoporov et al. 2018)"""
        profile = create_contamination_profile('realistic', random_seed=42)

        # 0.2 μm filtration
        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        _, stats = vlp.apply_contamination_reduction(profile)

        bacteria_removal = stats['reduction_by_type']['reagent_bacteria']['reduction_factor']
        # Literature: >90% removal for bacteria (1-5 μm) through 0.2 μm filter
        assert bacteria_removal > 0.90, f"Bacterial filtration too low: {bacteria_removal*100:.1f}%"

    def test_vlp_enrichment_fold_change(self):
        """Test that VLP shows 5-20x contamination reduction vs bulk"""
        profile = create_contamination_profile('realistic', random_seed=42)

        vlp = VLPEnrichment(
            protocol=VLPProtocol.tangential_flow_standard(),
            random_seed=42
        )
        vlp_profile, _ = vlp.apply_contamination_reduction(profile)

        bulk = VLPEnrichment(
            protocol=VLPProtocol.no_vlp(),
            random_seed=42
        )
        bulk_profile, _ = bulk.apply_contamination_reduction(profile)

        fold_change = bulk_profile.get_total_abundance() / max(vlp_profile.get_total_abundance(), 1e-10)

        # Literature-based expectation: 5-20x reduction
        assert 5.0 < fold_change < 25.0, f"Fold change outside literature range: {fold_change:.1f}x"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
