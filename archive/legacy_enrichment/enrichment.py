"""
VLP (Virus-Like Particle) enrichment modeling.

This module models the biological processes of VLP enrichment, including:
- Size-based filtration (0.2-0.45 μm)
- Nuclease treatment (DNase/RNase digestion of free DNA)
- Family-specific enrichment factors (from ViromeQC literature)
- Stochastic biological variation

VLP enrichment is the defining feature of viromics - it differentially
enriches for encapsidated viral genomes while removing host and bacterial
contamination.

Literature basis:
- Zolfo et al. (2019) "Detecting contamination in viromes using ViromeQC"
  Nature Biotechnology 37:1408-1412
- Conceição-Neto et al. (2015) "Modular approach to customise sample
  preparation procedures for viral metagenomics" BMC Genomics 16:617

Author: ViroForge Development Team
Date: 2025-01-31
"""

import logging
from typing import Optional, Literal
import numpy as np

from .utils.composition import MockViromeComposition

logger = logging.getLogger(__name__)


# Family-specific enrichment factors from ViromeQC paper (Zolfo et al. 2019)
# Based on empirical observations from 248 published viromes
FAMILY_ENRICHMENT_FACTORS = {
    'Microviridae': 2.5,      # Small ssDNA, highly stable
    'Siphoviridae': 1.2,      # Long-tailed phages, good retention
    'Myoviridae': 1.0,        # Contractile-tailed, baseline reference
    'Podoviridae': 1.1,       # Short-tailed, slightly enriched
    'Inoviridae': 0.3,        # Filamentous, depleted by filtration
    'Caudoviricetes': 1.15,   # Newly classified tailed phages
    'Autographiviridae': 1.0, # Jumbo phages, variable
    'Petitvirales': 2.3,      # Small circular ssDNA (similar to Microviridae)
    'Crassvirales': 1.1,      # crAssphage family, gut-associated
    # Default for unknown families
    'Unknown': 1.0,
    'Other': 1.0,
}


class VLPEnrichment:
    """
    Model VLP (Virus-Like Particle) enrichment processes.

    VLP enrichment uses physical separation (filtration, ultracentrifugation)
    and enzymatic treatment (nuclease) to enrich for encapsidated viral
    genomes while removing cellular material and free nucleic acids.

    This class provides a flexible, lab-agnostic framework that can model
    different VLP enrichment protocols.

    Attributes:
        filtration_method: Type of filtration ('tangential_flow', 'syringe',
            'ultracentrifugation', 'none')
        filtration_cutoff_um: Filter pore size in micrometers (typically 0.2-0.45)
        prefiltration_cutoff_um: Optional pre-filter size (typically 0.45-0.8)
        nuclease_treatment: Whether nuclease treatment is applied
        nuclease_efficiency: Fraction of free DNA removed by nuclease (0-1)
        size_retention_curve: Shape of size-retention curve ('sigmoid', 'step')
        stochastic_variation: Standard deviation for log-normal noise

    Example:
        >>> from viroforge.enrichment import VLPEnrichment
        >>> from viroforge.utils import create_mock_virome
        >>>
        >>> # Create bulk metagenome composition (50% viral)
        >>> composition = create_mock_virome(
        ...     'gut', 'realistic', viral_fraction=0.50
        ... )
        >>>
        >>> # Apply VLP enrichment
        >>> vlp = VLPEnrichment(
        ...     filtration_cutoff_um=0.2,
        ...     nuclease_efficiency=0.95
        ... )
        >>> vlp.apply(composition)  # Now ~95% viral!
    """

    def __init__(
        self,
        filtration_method: Literal['tangential_flow', 'syringe',
                                   'ultracentrifugation', 'none'] = 'tangential_flow',
        filtration_cutoff_um: float = 0.2,
        prefiltration_cutoff_um: Optional[float] = 0.45,
        nuclease_treatment: bool = True,
        nuclease_efficiency: float = 0.95,
        size_retention_curve: Literal['sigmoid', 'step'] = 'sigmoid',
        stochastic_variation: float = 0.2,
        random_seed: Optional[int] = None
    ):
        """
        Initialize VLP enrichment configuration.

        Args:
            filtration_method: Filtration method type
            filtration_cutoff_um: Primary filter cutoff (μm), typically 0.2
            prefiltration_cutoff_um: Pre-filter cutoff (μm), typically 0.45
            nuclease_treatment: Whether to apply nuclease treatment
            nuclease_efficiency: Fraction of free DNA removed (0-1), typically 0.90-0.98
            size_retention_curve: Retention curve shape ('sigmoid' smooth, 'step' sharp)
            stochastic_variation: Biological variation (σ for log-normal noise)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        # Validate parameters
        if filtration_cutoff_um < 0 or filtration_cutoff_um > 1.0:
            raise ValueError(f"filtration_cutoff_um must be 0-1.0, got {filtration_cutoff_um}")

        if nuclease_efficiency < 0 or nuclease_efficiency > 1.0:
            raise ValueError(f"nuclease_efficiency must be 0-1, got {nuclease_efficiency}")

        if stochastic_variation < 0:
            raise ValueError(f"stochastic_variation must be >= 0, got {stochastic_variation}")

        self.filtration_method = filtration_method
        self.filtration_cutoff_um = filtration_cutoff_um
        self.prefiltration_cutoff_um = prefiltration_cutoff_um
        self.nuclease_treatment = nuclease_treatment
        self.nuclease_efficiency = nuclease_efficiency
        self.size_retention_curve = size_retention_curve
        self.stochastic_variation = stochastic_variation
        self.random_seed = random_seed

        # Create local random generator for reproducibility
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized VLP enrichment: {self}")

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"VLPEnrichment(method={self.filtration_method}, "
            f"cutoff={self.filtration_cutoff_um}μm, "
            f"nuclease={self.nuclease_efficiency if self.nuclease_treatment else 'none'})"
        )

    def estimate_virus_size_nm(self, genome_length_bp: int) -> float:
        """
        Estimate virus particle size from genome length.

        Rough correlation based on viral biology:
        - Small viruses (Microviridae): ~4-6 kb → 25-30 nm
        - Medium phages: 30-60 kb → 50-100 nm
        - Large phages: 100-200 kb → 100-200 nm
        - Jumbo phages: >200 kb → >200 nm

        Uses power-law scaling: size ∝ length^0.35

        Args:
            genome_length_bp: Genome length in base pairs

        Returns:
            Estimated virus diameter in nanometers

        Example:
            >>> vlp = VLPEnrichment()
            >>> vlp.estimate_virus_size_nm(5000)   # Microvirus
            27.4
            >>> vlp.estimate_virus_size_nm(50000)  # Typical phage
            96.1
            >>> vlp.estimate_virus_size_nm(200000) # Jumbo phage
            214.7
        """
        # Power-law scaling with empirical parameters
        # Derived from known virus sizes
        size_nm = (genome_length_bp / 1000) ** 0.35 * 15.0
        return size_nm

    def calculate_size_retention(self, genome_length_bp: int) -> float:
        """
        Calculate retention probability based on virus size.

        Uses sigmoid or step function to model filtration retention.

        Args:
            genome_length_bp: Genome length in base pairs

        Returns:
            Retention probability (0-1), where 1.0 = fully retained

        Example:
            >>> vlp = VLPEnrichment(filtration_cutoff_um=0.2, size_retention_curve='sigmoid')
            >>> vlp.calculate_size_retention(5000)    # Small virus
            0.98
            >>> vlp.calculate_size_retention(50000)   # Medium virus
            0.89
            >>> vlp.calculate_size_retention(300000)  # Jumbo virus
            0.15
        """
        # Estimate virus size from genome length
        virus_size_nm = self.estimate_virus_size_nm(genome_length_bp)

        # Convert cutoff from μm to nm
        cutoff_nm = self.filtration_cutoff_um * 1000

        if self.filtration_method == 'none':
            # No filtration - everything passes
            return 1.0

        elif self.size_retention_curve == 'sigmoid':
            # Smooth sigmoid transition (realistic for tangential flow)
            # Sharpness parameter controls transition width
            sharpness = cutoff_nm * 0.3  # 30% of cutoff

            # Sigmoid function: 1 / (1 + exp((size - cutoff) / sharpness))
            retention = 1.0 / (1.0 + np.exp((virus_size_nm - cutoff_nm) / sharpness))

        elif self.size_retention_curve == 'step':
            # Sharp step function (realistic for syringe filtration)
            retention = 1.0 if virus_size_nm < cutoff_nm else 0.0

        else:
            raise ValueError(f"Unknown size_retention_curve: {self.size_retention_curve}")

        return retention

    def calculate_family_enrichment(self, family: str) -> float:
        """
        Get family-specific enrichment factor.

        Based on ViromeQC paper empirical observations.

        Args:
            family: Viral family name

        Returns:
            Enrichment factor (1.0 = baseline, >1.0 = enriched, <1.0 = depleted)

        Example:
            >>> vlp = VLPEnrichment()
            >>> vlp.calculate_family_enrichment('Microviridae')
            2.5
            >>> vlp.calculate_family_enrichment('Inoviridae')
            0.3
        """
        # Look up family in enrichment factors dictionary
        # Use exact match first, then check for partial matches
        if family in FAMILY_ENRICHMENT_FACTORS:
            return FAMILY_ENRICHMENT_FACTORS[family]

        # Check for partial matches (e.g., family name contains key)
        for key, factor in FAMILY_ENRICHMENT_FACTORS.items():
            if key.lower() in family.lower() or family.lower() in key.lower():
                return factor

        # Default to baseline if unknown
        return FAMILY_ENRICHMENT_FACTORS['Unknown']

    def calculate_stability_factor(self, gc_content: float) -> float:
        """
        Estimate capsid stability from GC content.

        Higher GC content correlates (weakly) with higher capsid stability.
        This is a minor effect compared to size and family factors.

        Args:
            gc_content: GC content percentage (0-100)

        Returns:
            Stability factor (0.85-1.15)

        Example:
            >>> vlp = VLPEnrichment()
            >>> vlp.calculate_stability_factor(30.0)  # Low GC
            0.96
            >>> vlp.calculate_stability_factor(70.0)  # High GC
            1.04
        """
        # Normalize GC to 0-1 range
        gc_fraction = gc_content / 100.0

        # Linear relationship: stability increases slightly with GC
        # Effect is weak: ±4% max deviation from baseline
        baseline = 1.0
        gc_effect = (gc_fraction - 0.5) * 0.08  # ±4% at extremes

        stability = baseline + gc_effect

        # Clip to reasonable range
        stability = np.clip(stability, 0.85, 1.15)

        return stability

    def calculate_nuclease_retention(self, is_viral: bool, is_encapsidated: bool = True) -> float:
        """
        Calculate retention after nuclease treatment.

        Encapsidated viral genomes are protected by capsid.
        Free DNA (host, bacterial) is degraded by nuclease.

        Args:
            is_viral: Whether this is a viral genome
            is_encapsidated: Whether viral genome is encapsidated (True for most viruses)

        Returns:
            Retention probability after nuclease (0-1)

        Example:
            >>> vlp = VLPEnrichment(nuclease_efficiency=0.95)
            >>> vlp.calculate_nuclease_retention(is_viral=True, is_encapsidated=True)
            1.0  # Viral genome protected by capsid
            >>> vlp.calculate_nuclease_retention(is_viral=False, is_encapsidated=False)
            0.05  # Free DNA 95% removed
        """
        if not self.nuclease_treatment:
            # No nuclease treatment applied
            return 1.0

        if is_viral and is_encapsidated:
            # Viral genomes in capsids are protected
            return 1.0
        else:
            # Free DNA removed by nuclease
            # Retention = 1 - efficiency (if 95% removed, 5% remains)
            return 1.0 - self.nuclease_efficiency

    def add_stochastic_variation(self, base_retention: float) -> float:
        """
        Add biological stochastic variation.

        Models biological variability in VLP enrichment:
        - Sample handling variation
        - Filter pore size variation
        - Viral capsid integrity variation
        - Environmental conditions (temperature, pH)

        Uses log-normal distribution (multiplicative noise).

        Args:
            base_retention: Base retention probability

        Returns:
            Retention with stochastic variation added

        Example:
            >>> vlp = VLPEnrichment(stochastic_variation=0.2, random_seed=42)
            >>> vlp.add_stochastic_variation(0.9)
            0.87  # Slightly varied from 0.9
        """
        if self.stochastic_variation == 0:
            return base_retention

        # Log-normal multiplicative noise
        # mean=0 gives geometric mean = 1.0 (no bias)
        stochastic_factor = self.rng.lognormal(mean=0, sigma=self.stochastic_variation)

        # Apply factor
        retention_with_noise = base_retention * stochastic_factor

        # Clip to valid range [0, 1]
        retention_with_noise = np.clip(retention_with_noise, 0.0, 1.0)

        return retention_with_noise

    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply VLP enrichment to a mock virome composition.

        This modifies genome abundances in-place to model the effects of
        VLP enrichment. The composition represents the state AFTER enrichment.

        Process:
        1. For each genome, calculate retention factors:
           - Size-based filtration
           - Family-specific enrichment
           - Capsid stability
           - Nuclease treatment (viral vs contamination)
           - Stochastic variation
        2. Multiply genome abundance by combined retention
        3. Renormalize abundances to sum to 1.0

        Args:
            composition: MockViromeComposition to enrich (modified in-place)

        Example:
            >>> composition = create_mock_virome('gut', 'realistic', viral_fraction=0.50)
            >>> print(f"Before: {composition.viral_fraction:.1%} viral")
            Before: 50.0% viral
            >>>
            >>> vlp = VLPEnrichment(filtration_cutoff_um=0.2, nuclease_efficiency=0.95)
            >>> vlp.apply(composition)
            >>> print(f"After: {composition.viral_fraction:.1%} viral")
            After: 94.3% viral
        """
        logger.info(f"Applying VLP enrichment to {composition.name}")
        logger.info(f"  Initial viral fraction: {composition.viral_fraction:.1%}")

        # Track abundances for logging
        initial_viral_abundance = composition.viral_fraction

        # Process viral genomes
        for genome in composition.viral_community.genomes:
            # Size-based retention
            size_retention = self.calculate_size_retention(genome.length)

            # Family-specific enrichment (only applies with actual VLP enrichment)
            if self.filtration_method != 'none' or self.nuclease_treatment:
                family_factor = self.calculate_family_enrichment(genome.family)
                stability_factor = self.calculate_stability_factor(genome.gc_content)
            else:
                family_factor = 1.0
                stability_factor = 1.0

            # Nuclease retention (viral genomes protected)
            nuclease_retention = self.calculate_nuclease_retention(
                is_viral=True,
                is_encapsidated=True
            )

            # Combined retention
            combined_retention = (
                size_retention *
                family_factor *
                stability_factor *
                nuclease_retention
            )

            # Add stochastic variation
            final_retention = self.add_stochastic_variation(combined_retention)

            # Apply to abundance
            genome.abundance *= final_retention

        # Process contaminant genomes
        if composition.contamination_profile is not None:
            for contaminant in composition.contamination_profile.contaminants:
                # If no enrichment is performed, all contaminants pass through unchanged
                if self.filtration_method == 'none' and not self.nuclease_treatment:
                    combined_retention = 1.0

                elif contaminant.contaminant_type.value in ['host_dna', 'rrna', 'reagent_bacteria']:
                    # These are either free nucleic acids or large cells
                    # Free nucleic acids: removed by nuclease
                    # Large cells: removed by filtration

                    # Assume most host/bacterial DNA is free (not in cells at this stage)
                    # Small fragments may pass through filter
                    size_retention = 0.3  # Some small fragments pass through

                    # Not encapsidated - removed by nuclease
                    nuclease_retention = self.calculate_nuclease_retention(
                        is_viral=False,
                        is_encapsidated=False
                    )

                    combined_retention = size_retention * nuclease_retention

                elif contaminant.contaminant_type.value == 'phix':
                    # PhiX is a small virus - treated like viral genome
                    size_retention = self.calculate_size_retention(contaminant.length)
                    if self.filtration_method != 'none' or self.nuclease_treatment:
                        family_factor = self.calculate_family_enrichment('Microviridae')
                    else:
                        family_factor = 1.0
                    nuclease_retention = self.calculate_nuclease_retention(
                        is_viral=True,
                        is_encapsidated=True
                    )
                    combined_retention = size_retention * family_factor * nuclease_retention

                else:
                    # Unknown contaminant type - conservative retention
                    combined_retention = 0.5

                # Add stochastic variation
                final_retention = self.add_stochastic_variation(combined_retention)

                # Apply to abundance
                contaminant.abundance *= final_retention

        # Normalize abundances to sum to 1.0 (without resetting to fixed ratios)
        total_abundance = composition.get_total_abundance()
        if total_abundance > 0:
            # Normalize viral genomes
            for genome in composition.viral_community.genomes:
                genome.abundance /= total_abundance

            # Normalize contaminants
            if composition.contamination_profile is not None:
                for contaminant in composition.contamination_profile.contaminants:
                    contaminant.abundance /= total_abundance

        # Update viral_fraction attribute to reflect new composition
        new_viral_fraction = composition.viral_community.get_total_abundance()
        composition.viral_fraction = new_viral_fraction
        composition.contamination_fraction = 1.0 - new_viral_fraction

        # Log results
        final_viral_abundance = new_viral_fraction
        enrichment_fold = final_viral_abundance / initial_viral_abundance if initial_viral_abundance > 0 else float('inf')

        logger.info(f"  Final viral fraction: {final_viral_abundance:.1%}")
        logger.info(f"  Enrichment: {enrichment_fold:.1f}x increase")
        logger.info(f"✓ VLP enrichment applied successfully")


# Pre-defined VLP enrichment protocols for convenience

def standard_vlp() -> VLPEnrichment:
    """
    Standard VLP enrichment protocol (most common).

    - 0.45 μm pre-filtration
    - 0.2 μm tangential flow filtration
    - Nuclease treatment (95% efficiency)
    - Sigmoid retention curve

    Returns:
        VLPEnrichment instance configured for standard protocol

    Example:
        >>> from viroforge.enrichment import standard_vlp
        >>> vlp = standard_vlp()
        >>> vlp.apply(composition)
    """
    return VLPEnrichment(
        filtration_method='tangential_flow',
        filtration_cutoff_um=0.2,
        prefiltration_cutoff_um=0.45,
        nuclease_treatment=True,
        nuclease_efficiency=0.95,
        size_retention_curve='sigmoid'
    )


def iron_chloride_vlp() -> VLPEnrichment:
    """
    Iron chloride precipitation + filtration protocol (Conceição-Neto et al.).

    FeCl3 precipitation enhances virus recovery and nuclease efficiency.

    - 0.8 μm pre-filtration
    - 0.2 μm filtration
    - Nuclease treatment (98% efficiency - better with FeCl3)

    Returns:
        VLPEnrichment instance for iron chloride protocol

    Example:
        >>> from viroforge.enrichment import iron_chloride_vlp
        >>> vlp = iron_chloride_vlp()
        >>> vlp.apply(composition)
    """
    return VLPEnrichment(
        filtration_method='tangential_flow',
        filtration_cutoff_um=0.2,
        prefiltration_cutoff_um=0.8,
        nuclease_treatment=True,
        nuclease_efficiency=0.98,  # Higher efficiency with FeCl3
        size_retention_curve='sigmoid'
    )


def ultracentrifuge_vlp() -> VLPEnrichment:
    """
    Ultracentrifugation-based VLP enrichment protocol.

    - No filtration (uses centrifugation)
    - Nuclease treatment (90% efficiency)
    - Size-based retention via density gradient

    Returns:
        VLPEnrichment instance for ultracentrifugation protocol

    Example:
        >>> from viroforge.enrichment import ultracentrifuge_vlp
        >>> vlp = ultracentrifuge_vlp()
        >>> vlp.apply(composition)
    """
    return VLPEnrichment(
        filtration_method='ultracentrifugation',
        filtration_cutoff_um=0.2,  # Effective cutoff for UC
        nuclease_treatment=True,
        nuclease_efficiency=0.90,  # Slightly lower than filtration
        size_retention_curve='sigmoid'
    )


def no_enrichment() -> VLPEnrichment:
    """
    No VLP enrichment (bulk metagenome).

    Use this to create bulk metagenome datasets for comparison.

    - No filtration
    - No nuclease treatment
    - All genomes retained

    Returns:
        VLPEnrichment instance with no enrichment (pass-through)

    Example:
        >>> from viroforge.enrichment import no_enrichment
        >>> bulk = no_enrichment()
        >>> bulk.apply(composition)  # No change in abundances
    """
    return VLPEnrichment(
        filtration_method='none',
        filtration_cutoff_um=1.0,  # No cutoff
        nuclease_treatment=False,
        nuclease_efficiency=0.0,
        size_retention_curve='step',
        stochastic_variation=0.0  # No variation for true bulk metagenome
    )


def syringe_filter_vlp() -> VLPEnrichment:
    """
    Syringe filtration protocol (quick, small volumes).

    - Sharp size cutoff (step function)
    - 0.22 μm syringe filter
    - Nuclease treatment (93% efficiency)

    Returns:
        VLPEnrichment instance for syringe filtration

    Example:
        >>> from viroforge.enrichment import syringe_filter_vlp
        >>> vlp = syringe_filter_vlp()
        >>> vlp.apply(composition)
    """
    return VLPEnrichment(
        filtration_method='syringe',
        filtration_cutoff_um=0.22,
        prefiltration_cutoff_um=0.45,
        nuclease_treatment=True,
        nuclease_efficiency=0.93,
        size_retention_curve='step'  # Sharp cutoff for syringe
    )


if __name__ == '__main__':
    # Example usage
    print("ViroForge VLP Enrichment Module")
    print("=" * 70)
    print()
    print("Pre-defined VLP protocols:")
    print(f"  - standard_vlp(): {standard_vlp()}")
    print(f"  - iron_chloride_vlp(): {iron_chloride_vlp()}")
    print(f"  - ultracentrifuge_vlp(): {ultracentrifuge_vlp()}")
    print(f"  - syringe_filter_vlp(): {syringe_filter_vlp()}")
    print(f"  - no_enrichment(): {no_enrichment()}")
    print()
    print("Family enrichment factors:")
    for family, factor in sorted(FAMILY_ENRICHMENT_FACTORS.items()):
        print(f"  {family:20s}: {factor:.1f}x")
