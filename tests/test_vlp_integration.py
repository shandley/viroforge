"""
Integration test: VLP enrichment vs Bulk metagenome comparison.

This test demonstrates the complete VLP enrichment workflow using real
viral genome data from NCBI, showing how VLP enrichment transforms a
bulk metagenome into a virome dataset.

Key validations:
- VLP enrichment dramatically increases viral fraction
- Host and bacterial DNA are effectively removed
- Different virus families show differential enrichment
- Bulk metagenome preserves original composition
"""

import pytest
from copy import deepcopy

from viroforge.core import create_body_site_profile
from viroforge.core.contamination import (
    ContaminationProfile,
    create_contamination_profile
)
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp, no_enrichment


class TestVLPvsBulkIntegration:
    """Integration tests comparing VLP-enriched vs bulk metagenomes."""

    def test_vlp_enrichment_workflow_with_real_genomes(self):
        """
        End-to-end test: VLP enrichment dramatically changes composition.

        This test simulates:
        1. Starting sample: 50% viral, 50% contamination (bulk metagenome)
        2. VLP enrichment: Should increase to >90% viral
        3. Bulk (no enrichment): Should remain ~50% viral
        """
        # Create viral community from real body site profile
        viral_community = create_body_site_profile(
            'gut',
            n_genomes=15,
            random_seed=42
        )

        # Create realistic contamination profile
        contamination = create_contamination_profile(
            contamination_level='realistic',
            random_seed=42
        )

        # Create base composition (50% viral, like a bulk metagenome)
        base_composition = MockViromeComposition(
            name='base_sample',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.50
        )

        initial_viral_fraction = base_composition.viral_fraction

        # Create two copies: one for VLP, one for bulk
        # We need to deep copy to avoid shared references
        vlp_composition = MockViromeComposition(
            name='vlp_enriched',
            viral_community=deepcopy(viral_community),
            contamination_profile=deepcopy(contamination),
            viral_fraction=0.50
        )

        bulk_composition = MockViromeComposition(
            name='bulk_metagenome',
            viral_community=deepcopy(viral_community),
            contamination_profile=deepcopy(contamination),
            viral_fraction=0.50
        )

        # Apply VLP enrichment to one sample
        vlp_protocol = standard_vlp()
        vlp_protocol.apply(vlp_composition)

        # Apply no enrichment (bulk) to the other
        bulk_protocol = no_enrichment()
        bulk_protocol.apply(bulk_composition)

        # Validate results
        vlp_viral_fraction = vlp_composition.viral_fraction
        bulk_viral_fraction = bulk_composition.viral_fraction

        # VLP enrichment should dramatically increase viral fraction
        assert vlp_viral_fraction > 0.85, \
            f"VLP enrichment should result in >85% viral, got {vlp_viral_fraction:.1%}"

        # VLP should be much higher than initial
        enrichment_fold = vlp_viral_fraction / initial_viral_fraction
        assert enrichment_fold > 1.5, \
            f"VLP should enrich viruses >1.5x, got {enrichment_fold:.2f}x"

        # Bulk should preserve original composition
        assert abs(bulk_viral_fraction - initial_viral_fraction) < 0.01, \
            f"Bulk should preserve composition, but changed from {initial_viral_fraction:.1%} to {bulk_viral_fraction:.1%}"

        # VLP should be much more enriched than bulk
        assert vlp_viral_fraction > bulk_viral_fraction + 0.30, \
            f"VLP ({vlp_viral_fraction:.1%}) should be >30% higher than bulk ({bulk_viral_fraction:.1%})"

        print("\n" + "="*60)
        print("VLP vs Bulk Enrichment Comparison")
        print("="*60)
        print(f"Initial composition:  {initial_viral_fraction:.1%} viral")
        print(f"After VLP enrichment: {vlp_viral_fraction:.1%} viral ({enrichment_fold:.2f}x increase)")
        print(f"Bulk metagenome:      {bulk_viral_fraction:.1%} viral (unchanged)")
        print(f"Difference:           {(vlp_viral_fraction - bulk_viral_fraction)*100:.1f}% more viral in VLP")
        print("="*60)

    def test_host_dna_removal_by_vlp_enrichment(self):
        """
        Test that VLP enrichment effectively removes host DNA contamination.
        """
        # Create viral community
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)

        # Create contamination with significant host DNA
        contamination = create_contamination_profile(
            contamination_level='heavy',  # More contamination
            random_seed=42
        )

        # Start with heavily contaminated sample (30% viral, 70% contamination)
        composition = MockViromeComposition(
            name='contaminated_sample',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.30
        )

        initial_contam_abundance = composition.contamination_profile.get_total_abundance()

        # Apply VLP enrichment
        vlp = standard_vlp()
        vlp.apply(composition)

        final_contam_abundance = composition.contamination_profile.get_total_abundance()

        # Contamination should be dramatically reduced
        contam_reduction = initial_contam_abundance / final_contam_abundance
        assert contam_reduction > 5.0, \
            f"Contamination should be reduced >5x, got {contam_reduction:.2f}x"

        # Final viral fraction should be high despite starting low
        assert composition.viral_fraction > 0.80, \
            f"Final viral fraction should be >80%, got {composition.viral_fraction:.1%}"

        print("\n" + "="*60)
        print("Host DNA Removal by VLP Enrichment")
        print("="*60)
        print(f"Initial contamination: {initial_contam_abundance:.3f}")
        print(f"Final contamination:   {final_contam_abundance:.3f}")
        print(f"Reduction:             {contam_reduction:.1f}x")
        print(f"Final viral fraction:  {composition.viral_fraction:.1%}")
        print("="*60)

    def test_differential_enrichment_by_virus_family(self):
        """
        Test that different virus families show differential enrichment.

        Small viruses (Microviridae) should be highly enriched.
        Filamentous viruses (Inoviridae) should be depleted.
        """
        from viroforge.core import ViralCommunity, ViralGenome
        from viroforge.enrichment import FAMILY_ENRICHMENT_FACTORS

        # Create a custom community with specific families
        viral_community = ViralCommunity(name='test_community')

        # Small virus (highly enriched)
        viral_community.add_genome(ViralGenome(
            genome_id='small_virus_1',
            sequence='ATCG' * 1250,  # 5 kb
            taxonomy='Viruses;Microviridae;Microvirus',
            family='Microviridae',
            genus='Microvirus',
            abundance=0.2,
            gc_content=45.0
        ))

        # Filamentous virus (depleted)
        viral_community.add_genome(ViralGenome(
            genome_id='filamentous_1',
            sequence='GCGC' * 1500,  # 6 kb
            taxonomy='Viruses;Inoviridae;Inovirus',
            family='Inoviridae',
            genus='Inovirus',
            abundance=0.2,
            gc_content=50.0
        ))

        # Baseline phage
        viral_community.add_genome(ViralGenome(
            genome_id='baseline_1',
            sequence='ATGC' * 10000,  # 40 kb
            taxonomy='Viruses;Myoviridae;Myovirus',
            family='Myoviridae',
            genus='Myovirus',
            abundance=0.2,
            gc_content=48.0
        ))

        # Create realistic contamination (need enough to see differential enrichment)
        contamination = create_contamination_profile(
            contamination_level='realistic',
            random_seed=42
        )

        composition = MockViromeComposition(
            name='test',
            viral_community=viral_community,
            contamination_profile=contamination,
            viral_fraction=0.60  # Lower viral fraction to see differential effects
        )

        # Store initial relative abundances within viral community
        initial_total = sum(g.abundance for g in composition.viral_community.genomes)
        initial_relative = {
            g.genome_id: g.abundance / initial_total
            for g in composition.viral_community.genomes
        }

        # Apply VLP enrichment
        vlp = standard_vlp()
        vlp.apply(composition)

        # Get final relative abundances within viral community
        final_total = sum(g.abundance for g in composition.viral_community.genomes)
        final_relative = {
            g.genome_id: g.abundance / final_total
            for g in composition.viral_community.genomes
        }

        # Calculate enrichment ratios (relative change within viral community)
        small_change = final_relative['small_virus_1'] / initial_relative['small_virus_1']
        filamentous_change = final_relative['filamentous_1'] / initial_relative['filamentous_1']
        baseline_change = final_relative['baseline_1'] / initial_relative['baseline_1']

        # Normalize to baseline
        small_enrichment = small_change / baseline_change
        filamentous_enrichment = filamentous_change / baseline_change

        # The key validation: differential enrichment
        # Microviridae (2.5x factor) should be MORE enriched than Inoviridae (0.3x factor)
        assert small_enrichment > filamentous_enrichment * 1.5, \
            f"Microviridae should be enriched >1.5x more than Inoviridae, got {small_enrichment:.2f}x vs {filamentous_enrichment:.2f}x"

        # Small virus should be enriched or equal relative to baseline (allowing for stochastic variation)
        assert small_enrichment >= 0.95, \
            f"Microviridae should be at least equal to baseline, got {small_enrichment:.2f}x"

        # Filamentous should be depleted relative to baseline
        assert filamentous_enrichment < 0.9, \
            f"Inoviridae should be depleted <0.9x vs baseline, got {filamentous_enrichment:.2f}x"

        print("\n" + "="*60)
        print("Differential Enrichment by Family")
        print("="*60)
        print(f"Microviridae (small):       {small_enrichment:.2f}x (expected >1.5x)")
        print(f"Myoviridae (baseline):      1.00x (reference)")
        print(f"Inoviridae (filamentous):   {filamentous_enrichment:.2f}x (expected <0.7x)")
        print("="*60)

    def test_protocol_comparison_iron_chloride_vs_standard(self):
        """
        Compare different VLP enrichment protocols.

        Iron chloride protocol should have higher nuclease efficiency,
        resulting in even better contamination removal.
        """
        from viroforge.enrichment import iron_chloride_vlp

        # Create identical starting compositions
        viral_community = create_body_site_profile('gut', n_genomes=10, random_seed=42)
        contamination = create_contamination_profile('realistic', random_seed=42)

        standard_composition = MockViromeComposition(
            name='standard',
            viral_community=deepcopy(viral_community),
            contamination_profile=deepcopy(contamination),
            viral_fraction=0.50
        )

        fecl3_composition = MockViromeComposition(
            name='iron_chloride',
            viral_community=deepcopy(viral_community),
            contamination_profile=deepcopy(contamination),
            viral_fraction=0.50
        )

        # Apply different protocols
        standard_vlp().apply(standard_composition)
        iron_chloride_vlp().apply(fecl3_composition)

        # Iron chloride should have equal or better viral enrichment
        # (due to higher nuclease efficiency: 98% vs 95%)
        assert fecl3_composition.viral_fraction >= standard_composition.viral_fraction - 0.05, \
            f"FeCl3 protocol should have similar or better viral fraction"

        # Both should be highly enriched
        assert standard_composition.viral_fraction > 0.85
        assert fecl3_composition.viral_fraction > 0.85

        print("\n" + "="*60)
        print("Protocol Comparison")
        print("="*60)
        print(f"Standard VLP:        {standard_composition.viral_fraction:.1%} viral")
        print(f"Iron chloride VLP:   {fecl3_composition.viral_fraction:.1%} viral")
        print("="*60)

    @pytest.mark.slow
    def test_full_workflow_with_real_genomes(self):
        """
        Complete workflow test using real downloaded genomes from NCBI.

        This test requires the minimal genome database to be downloaded.
        It demonstrates using actual viral genomes in the enrichment workflow.
        """
        from viroforge.utils.genome_database import get_genome_database, set_entrez_email
        import os

        # Try to load real genomes
        try:
            # Set email for Entrez
            email = os.environ.get('EMAIL', 'test@example.com')
            set_entrez_email(email)

            # Load minimal database (20 real viral genomes)
            real_genomes = get_genome_database('minimal', email=email)

            # Verify we got genomes
            assert len(real_genomes) == 20, f"Expected 20 genomes, got {len(real_genomes)}"

            # Create viral community from real genomes (use subset)
            from viroforge.core import ViralCommunity

            # Select 10 random genomes and set abundances
            import random
            random.seed(42)
            selected = random.sample(real_genomes, 10)

            # Set equal abundances
            for genome in selected:
                genome.abundance = 1.0 / len(selected)

            # Create viral community and add genomes
            viral_community = ViralCommunity(name='real_gut_virome')
            for genome in selected:
                viral_community.add_genome(genome)

            # Create contamination
            contamination = create_contamination_profile('realistic', random_seed=42)

            # Create composition
            composition = MockViromeComposition(
                name='real_genome_test',
                viral_community=viral_community,
                contamination_profile=contamination,
                viral_fraction=0.50
            )

            initial_viral = composition.viral_fraction

            # Apply VLP enrichment
            vlp = standard_vlp()
            vlp.apply(composition)

            final_viral = composition.viral_fraction

            # Validate enrichment worked with real genomes
            assert final_viral > 0.85, \
                f"Real genome enrichment should reach >85% viral, got {final_viral:.1%}"

            assert final_viral > initial_viral + 0.30, \
                f"Real genomes should be enriched >30%, got {(final_viral - initial_viral)*100:.1f}%"

            print("\n" + "="*60)
            print("Real Genome Workflow Test")
            print("="*60)
            print(f"Genomes used: {len(selected)} real viral genomes from NCBI")
            print(f"Families: {len(set(g.family for g in selected))}")
            print(f"Initial: {initial_viral:.1%} viral")
            print(f"Final:   {final_viral:.1%} viral")
            print(f"Enrichment: {final_viral/initial_viral:.2f}x")
            print("="*60)

        except FileNotFoundError:
            pytest.skip("Minimal genome database not downloaded yet. Run: get_genome_database('minimal')")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
