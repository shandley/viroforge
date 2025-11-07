"""
Enhanced VLP (Virus-Like Particle) Enrichment Modeling

This module provides comprehensive VLP enrichment simulation based on:
- Virion size estimation from genome properties
- Size-based filtration with protocol-specific retention curves
- Nuclease treatment efficiency
- Contamination modeling

Literature references:
- Thurber et al. 2009 - VLP methodology
- Shkoporov et al. 2018 - VLP enrichment efficiency
- Roux et al. 2016 - ViromeQC survey data
- Lim et al. 2020 - VLP protocol comparison
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class VirionSizeEstimate:
    """Estimated virion physical properties"""
    genome_length: int
    genome_type: str
    estimated_diameter_nm: float
    estimated_volume_nm3: float
    confidence: str  # 'high', 'medium', 'low'


class VirionSizeEstimator:
    """
    Estimate virion size from genome properties

    Based on literature relationships between genome length,
    type, and physical virion dimensions.

    References:
    - Cui et al. 2014 - Virus capsid size database
    - Nasir et al. 2017 - Virion size scaling laws
    """

    # Empirical relationships from literature
    # diameter (nm) = a * log10(genome_length) + b
    SIZE_RELATIONSHIPS = {
        'dsDNA': {'a': 35, 'b': -70, 'variance': 15},  # Tailed phages, large dsDNA
        'ssDNA': {'a': 8, 'b': -5, 'variance': 5},     # Microviruses, circoviruses
        'dsRNA': {'a': 25, 'b': -50, 'variance': 10},  # Reoviruses
        'ssRNA': {'a': 15, 'b': -20, 'variance': 8},   # Many RNA viruses
        'ssRNA-RT': {'a': 30, 'b': -60, 'variance': 12},  # Retroviruses
        'unknown': {'a': 25, 'b': -45, 'variance': 20}  # Conservative estimate
    }

    def estimate_size(
        self,
        genome_length: int,
        genome_type: str = 'dsDNA',
        rng: Optional[np.random.Generator] = None
    ) -> VirionSizeEstimate:
        """
        Estimate virion diameter from genome length and type

        Args:
            genome_length: Genome length in base pairs
            genome_type: Genome type (dsDNA, ssDNA, dsRNA, ssRNA, etc.)
            rng: Optional random number generator for reproducibility

        Returns:
            VirionSizeEstimate with diameter, volume, and confidence
        """
        # Get relationship parameters
        if genome_type not in self.SIZE_RELATIONSHIPS:
            genome_type = 'unknown'
            confidence = 'low'
        else:
            # dsDNA/ssDNA have most empirical data (tailed phages, microviruses)
            # RNA viruses have more variable morphologies
            confidence = 'high' if genome_type in ['dsDNA', 'ssDNA'] else 'medium'

        params = self.SIZE_RELATIONSHIPS[genome_type]

        # Calculate base diameter
        log_length = np.log10(genome_length)
        diameter_nm = params['a'] * log_length + params['b']

        # Add stochastic variation
        if rng is not None:
            diameter_nm += rng.normal(0, params['variance'] * 0.3)
        else:
            diameter_nm += np.random.normal(0, params['variance'] * 0.3)

        # Ensure reasonable bounds
        diameter_nm = max(15, min(500, diameter_nm))  # 15 nm (small ssDNA) to 500 nm (giant viruses)

        # Calculate volume (assume spherical or icosahedral, close to spherical)
        volume_nm3 = (4/3) * np.pi * (diameter_nm/2)**3

        return VirionSizeEstimate(
            genome_length=genome_length,
            genome_type=genome_type,
            estimated_diameter_nm=diameter_nm,
            estimated_volume_nm3=volume_nm3,
            confidence=confidence
        )


class FiltrationCurve:
    """
    Model size-based retention during filtration

    Different filtration methods have different retention characteristics:
    - Sigmoid: Gradual transition at pore size (most filters)
    - Step: Sharp cutoff (ideal filter, doesn't exist)
    - Linear: Gradual transition (some depth filters)
    """

    @staticmethod
    def sigmoid(virion_diameter_nm: float, pore_size_um: float, steepness: float = 0.05) -> float:
        """
        Sigmoid retention curve

        Args:
            virion_diameter_nm: Virion diameter in nanometers
            pore_size_um: Filter pore size in micrometers
            steepness: Curve steepness (higher = sharper transition)

        Returns:
            Retention probability (0-1)
        """
        pore_size_nm = pore_size_um * 1000

        # Sigmoid centered at pore size
        # Small virions (<pore size) → low retention
        # Large virions (>pore size) → high retention
        retention = 1 / (1 + np.exp(-steepness * (virion_diameter_nm - pore_size_nm)))

        return retention

    @staticmethod
    def step(virion_diameter_nm: float, pore_size_um: float) -> float:
        """
        Step function retention (ideal filter)

        Returns:
            1.0 if diameter > pore size, else 0.0
        """
        pore_size_nm = pore_size_um * 1000
        return 1.0 if virion_diameter_nm > pore_size_nm else 0.0

    @staticmethod
    def linear(virion_diameter_nm: float, pore_size_um: float, transition_width_nm: float = 100) -> float:
        """
        Linear transition retention

        Args:
            virion_diameter_nm: Virion diameter in nanometers
            pore_size_um: Filter pore size in micrometers
            transition_width_nm: Width of transition region

        Returns:
            Retention probability (0-1)
        """
        pore_size_nm = pore_size_um * 1000

        if virion_diameter_nm < pore_size_nm - transition_width_nm/2:
            return 0.0
        elif virion_diameter_nm > pore_size_nm + transition_width_nm/2:
            return 1.0
        else:
            # Linear interpolation in transition region
            position = (virion_diameter_nm - (pore_size_nm - transition_width_nm/2)) / transition_width_nm
            return position


@dataclass
class VLPProtocolConfig:
    """Configuration for a specific VLP enrichment protocol"""
    name: str
    filtration_method: str  # 'tangential_flow', 'syringe', 'centrifugal', 'ultracentrifugation'
    pore_size_um: Optional[float]  # None for ultracentrifugation
    prefiltration_pore_size_um: Optional[float]
    retention_curve_type: str  # 'sigmoid', 'step', 'linear'
    retention_curve_steepness: float
    nuclease_treatment: bool
    nuclease_efficiency: float  # 0-1, fraction of free DNA/RNA degraded
    recovery_rate: float  # 0-1, fraction of virions retained
    contamination_reduction: float  # 0-1, fraction of contaminants removed


class VLPProtocol:
    """
    Pre-defined VLP enrichment protocols based on common methods

    Each protocol has specific characteristics based on literature:
    - Filtration method and pore size
    - Nuclease treatment efficiency
    - Recovery rates
    - Contamination reduction factors
    """

    @classmethod
    def tangential_flow_standard(cls) -> VLPProtocolConfig:
        """
        Standard tangential flow filtration (0.2 μm)

        High efficiency, good recovery
        Most commonly used in virome studies

        Reference: Thurber et al. 2009
        """
        return VLPProtocolConfig(
            name="Tangential Flow Filtration (0.2 μm)",
            filtration_method='tangential_flow',
            pore_size_um=0.2,
            prefiltration_pore_size_um=0.45,
            retention_curve_type='sigmoid',
            retention_curve_steepness=0.008,  # Gentle slope
            nuclease_treatment=True,
            nuclease_efficiency=0.98,  # Very high efficiency
            recovery_rate=0.85,  # Good recovery
            contamination_reduction=0.95  # Excellent decontamination
        )

    @classmethod
    def syringe_filter_standard(cls) -> VLPProtocolConfig:
        """
        Syringe filtration (0.2 μm)

        Variable efficiency, lower recovery
        Simpler but less controlled than TFF

        Reference: Lim et al. 2020
        """
        return VLPProtocolConfig(
            name="Syringe Filter (0.2 μm)",
            filtration_method='syringe',
            pore_size_um=0.2,
            prefiltration_pore_size_um=0.45,
            retention_curve_type='sigmoid',
            retention_curve_steepness=0.010,  # Sharper slope than TFF
            nuclease_treatment=True,
            nuclease_efficiency=0.90,  # Lower efficiency (less mixing)
            recovery_rate=0.60,  # Lower recovery (filter clogging)
            contamination_reduction=0.85  # Good but variable
        )

    @classmethod
    def ultracentrifugation(cls) -> VLPProtocolConfig:
        """
        Ultracentrifugation (100,000 x g)

        Size-independent (density-based)
        High recovery, moderate decontamination

        Reference: Shkoporov et al. 2018
        """
        return VLPProtocolConfig(
            name="Ultracentrifugation (100,000g)",
            filtration_method='ultracentrifugation',
            pore_size_um=None,  # Not filtration-based
            prefiltration_pore_size_um=0.45,  # Remove large debris first
            retention_curve_type='step',  # All virions retained equally
            retention_curve_steepness=1.0,
            nuclease_treatment=True,
            nuclease_efficiency=0.95,
            recovery_rate=0.90,  # High recovery
            contamination_reduction=0.75  # Moderate (density-based, not size-based)
        )

    @classmethod
    def norgen_kit(cls) -> VLPProtocolConfig:
        """
        Norgen Phage DNA Isolation Kit

        Commercial kit, column-based
        Variable performance

        Reference: Manufacturer data + Roux et al. 2016
        """
        return VLPProtocolConfig(
            name="Norgen Phage DNA Kit",
            filtration_method='column',
            pore_size_um=None,  # Column-based, not filtration
            prefiltration_pore_size_um=None,
            retention_curve_type='sigmoid',
            retention_curve_steepness=0.015,
            nuclease_treatment=True,
            nuclease_efficiency=0.92,
            recovery_rate=0.70,  # Variable
            contamination_reduction=0.88
        )

    @classmethod
    def no_vlp(cls) -> VLPProtocolConfig:
        """
        No VLP enrichment (bulk metagenome)

        For comparison - includes all cellular contamination
        """
        return VLPProtocolConfig(
            name="No VLP Enrichment (Bulk)",
            filtration_method='none',
            pore_size_um=None,
            prefiltration_pore_size_um=None,
            retention_curve_type='step',
            retention_curve_steepness=1.0,
            nuclease_treatment=False,
            nuclease_efficiency=0.0,
            recovery_rate=1.0,
            contamination_reduction=0.0  # No decontamination
        )


class VLPEnrichment:
    """
    Main VLP enrichment model

    Applies protocol-specific VLP enrichment to viral abundances,
    modeling:
    - Size-based retention (for filtration methods)
    - Nuclease treatment efficiency
    - Recovery rates
    - Contamination reduction

    Usage:
        vlp = VLPEnrichment(protocol=VLPProtocol.tangential_flow_standard())
        enriched_abundances = vlp.apply_enrichment(
            genomes, abundances, contamination_level='clean'
        )
    """

    def __init__(
        self,
        protocol: VLPProtocolConfig,
        random_seed: Optional[int] = None
    ):
        self.protocol = protocol
        self.size_estimator = VirionSizeEstimator()
        self.random_seed = random_seed

        # Use local RNG instead of global state for thread safety
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized VLP enrichment: {protocol.name}")

    def apply_enrichment(
        self,
        genomes: List[Dict],
        abundances: np.ndarray,
        add_contamination: bool = True,
        contamination_level: str = 'clean'
    ) -> Tuple[np.ndarray, Dict]:
        """
        Apply VLP enrichment to genome abundances

        Args:
            genomes: List of genome dictionaries with 'length' and 'genome_type'
            abundances: Numpy array of relative abundances (must sum to 1.0)
            add_contamination: Whether to add contamination based on protocol efficiency
            contamination_level: 'clean', 'realistic', 'heavy' (if add_contamination=True)

        Returns:
            Tuple of (enriched_abundances, enrichment_stats)
        """
        n_genomes = len(genomes)
        enriched_abundances = abundances.copy()

        enrichment_factors = np.ones(n_genomes)
        size_estimates = []

        # Apply size-based retention (if applicable)
        if self.protocol.filtration_method in ['tangential_flow', 'syringe', 'centrifugal']:
            for i, genome in enumerate(genomes):
                # Estimate virion size
                size_est = self.size_estimator.estimate_size(
                    genome['length'],
                    genome.get('genome_type', 'dsDNA'),
                    rng=self.rng
                )
                size_estimates.append(size_est)

                # Calculate retention probability
                if self.protocol.pore_size_um is not None:
                    if self.protocol.retention_curve_type == 'sigmoid':
                        retention = FiltrationCurve.sigmoid(
                            size_est.estimated_diameter_nm,
                            self.protocol.pore_size_um,
                            self.protocol.retention_curve_steepness
                        )
                    elif self.protocol.retention_curve_type == 'step':
                        retention = FiltrationCurve.step(
                            size_est.estimated_diameter_nm,
                            self.protocol.pore_size_um
                        )
                    else:  # linear
                        retention = FiltrationCurve.linear(
                            size_est.estimated_diameter_nm,
                            self.protocol.pore_size_um
                        )
                    enrichment_factors[i] = retention
                else:
                    # Column-based or other methods without pore size
                    # Use moderate size-dependent retention
                    # Smaller viruses retained less efficiently
                    size_factor = min(1.0, size_est.estimated_diameter_nm / 100.0)
                    enrichment_factors[i] = 0.5 + 0.5 * size_factor

        elif self.protocol.filtration_method == 'column':
            # Column-based methods (e.g., Norgen kit)
            # Moderate size-dependent retention
            for i, genome in enumerate(genomes):
                size_est = self.size_estimator.estimate_size(
                    genome['length'],
                    genome.get('genome_type', 'dsDNA'),
                    rng=self.rng
                )
                size_estimates.append(size_est)

                # Modest size bias (larger viruses bind better to columns)
                size_factor = min(1.0, size_est.estimated_diameter_nm / 80.0)
                enrichment_factors[i] = 0.6 + 0.4 * size_factor

        # Apply recovery rate (overall loss during processing)
        enrichment_factors *= self.protocol.recovery_rate

        # Add stochastic variation (biological and technical variability)
        variation = self.rng.normal(1.0, 0.05, n_genomes)  # 5% CV
        enrichment_factors *= variation

        # Apply enrichment
        enriched_abundances = abundances * enrichment_factors

        # Renormalize (total viral content changes, but relative viral abundances renormalized)
        viral_total = enriched_abundances.sum()
        if viral_total > 0:
            enriched_abundances /= viral_total

        # Calculate enrichment statistics
        stats = {
            'protocol': self.protocol.name,
            'mean_enrichment_factor': float(np.mean(enrichment_factors)),
            'median_enrichment_factor': float(np.median(enrichment_factors)),
            'enrichment_range': (float(np.min(enrichment_factors)), float(np.max(enrichment_factors))),
            'viral_recovery_rate': self.protocol.recovery_rate,
            'nuclease_efficiency': self.protocol.nuclease_efficiency if self.protocol.nuclease_treatment else 0.0,
            'contamination_reduction': self.protocol.contamination_reduction,
            'n_genomes_retained': int(np.sum(enriched_abundances > 0.0001)),  # >0.01% abundance
            'size_bias': self._calculate_size_bias(size_estimates, enrichment_factors)
        }

        logger.info(f"VLP enrichment applied: {stats['n_genomes_retained']}/{n_genomes} genomes retained")
        logger.info(f"Mean enrichment factor: {stats['mean_enrichment_factor']:.3f}")

        return enriched_abundances, stats

    def _calculate_size_bias(self, size_estimates: List[VirionSizeEstimate],
                            enrichment_factors: np.ndarray) -> Dict[str, float]:
        """Calculate size bias metrics"""
        if not size_estimates:
            return {}

        sizes = np.array([s.estimated_diameter_nm for s in size_estimates])

        # Correlation between size and enrichment
        if len(sizes) > 1:
            correlation = np.corrcoef(sizes, enrichment_factors)[0, 1]
        else:
            correlation = 0.0

        return {
            'size_enrichment_correlation': float(correlation),
            'mean_virion_diameter_nm': float(np.mean(sizes)),
            'size_range_nm': (float(np.min(sizes)), float(np.max(sizes)))
        }

    def apply_contamination_reduction(
        self,
        contamination_profile: 'ContaminationProfile'
    ) -> Tuple['ContaminationProfile', Dict]:
        """
        Apply protocol-specific contamination reduction

        Different contaminant types are affected differently:
        - Host DNA (free): Removed by nuclease treatment (90-98%)
        - rRNA (free/debris): Removed by nuclease treatment (85-95%)
        - Reagent bacteria (cells): Removed by filtration (75-95%)
        - PhiX (encapsidated): Treated like small virus (10-40% loss)

        Args:
            contamination_profile: Original contamination profile

        Returns:
            Tuple of (reduced_profile, reduction_stats)

        Literature basis:
        - Thurber et al. 2009: DNase reduces free DNA by >95%
        - Shkoporov et al. 2018: 0.2 μm filtration removes >90% bacteria
        - Roux et al. 2016 (ViromeQC): Typical contamination ranges
        """
        from ..core.contamination import ContaminationProfile, ContaminantType

        # Create a copy of the contamination profile
        reduced_profile = ContaminationProfile(name=f"{contamination_profile.name}_vlp_reduced")

        original_abundances = {}
        reduced_abundances = {}
        reduction_factors_by_type = {}

        for contaminant in contamination_profile.contaminants:
            ctype = contaminant.contaminant_type
            original_abundance = contaminant.abundance

            # Calculate reduction factor based on contaminant type
            if ctype == ContaminantType.HOST_DNA:
                # Free host DNA - highly susceptible to nuclease
                if self.protocol.nuclease_treatment:
                    base_removal = self.protocol.nuclease_efficiency
                    # Add some variation (biological + technical)
                    removal = base_removal * self.rng.normal(1.0, 0.05)
                    removal = np.clip(removal, 0.85, 0.99)
                else:
                    # Without nuclease, only filtration helps (minimal for free DNA)
                    removal = 0.20 * self.protocol.contamination_reduction

            elif ctype == ContaminantType.RRNA:
                # rRNA - some protection in cellular debris
                if self.protocol.nuclease_treatment:
                    base_removal = self.protocol.nuclease_efficiency * 0.92  # Slightly less efficient
                    removal = base_removal * self.rng.normal(1.0, 0.08)
                    removal = np.clip(removal, 0.80, 0.97)
                else:
                    removal = 0.15 * self.protocol.contamination_reduction

            elif ctype == ContaminantType.REAGENT_BACTERIA:
                # Bacterial cells - primarily removed by filtration
                if self.protocol.filtration_method in ['tangential_flow', 'syringe', 'centrifugal']:
                    # Bacteria are 1-5 μm, much larger than pore size
                    # Should be nearly 100% removed by 0.2 μm filter
                    base_removal = self.protocol.contamination_reduction * 0.98
                    removal = base_removal * self.rng.normal(1.0, 0.10)
                    removal = np.clip(removal, 0.70, 0.99)
                elif self.protocol.filtration_method == 'ultracentrifugation':
                    # Cells pellet differently than viruses
                    removal = 0.85 * self.rng.normal(1.0, 0.15)
                    removal = np.clip(removal, 0.60, 0.95)
                elif self.protocol.filtration_method == 'column':
                    # Column-based kits variable
                    removal = 0.88 * self.rng.normal(1.0, 0.12)
                    removal = np.clip(removal, 0.65, 0.95)
                else:  # no filtration
                    removal = 0.0

                # Nuclease also helps if cells are lysed
                if self.protocol.nuclease_treatment:
                    # Assume ~30% of bacteria are lysed
                    lysed_fraction = 0.30
                    additional_removal = lysed_fraction * self.protocol.nuclease_efficiency * 0.5
                    removal = min(0.99, removal + additional_removal)

            elif ctype == ContaminantType.PHIX:
                # PhiX is an encapsidated virus (27 nm, ssDNA, circular)
                # Treat similarly to small viruses
                # PhiX genome: 5,386 bp, ssDNA
                phix_size = self.size_estimator.estimate_size(5386, 'ssDNA', rng=self.rng)

                if self.protocol.filtration_method in ['tangential_flow', 'syringe', 'centrifugal', 'column']:
                    # PhiX is small (27 nm), some loss through filters
                    if self.protocol.pore_size_um is not None:
                        # Calculate retention like other viruses
                        if self.protocol.retention_curve_type == 'sigmoid':
                            retention = FiltrationCurve.sigmoid(
                                phix_size.estimated_diameter_nm,
                                self.protocol.pore_size_um,
                                self.protocol.retention_curve_steepness
                            )
                        elif self.protocol.retention_curve_type == 'step':
                            retention = FiltrationCurve.step(
                                phix_size.estimated_diameter_nm,
                                self.protocol.pore_size_um
                            )
                        else:  # linear
                            retention = FiltrationCurve.linear(
                                phix_size.estimated_diameter_nm,
                                self.protocol.pore_size_um
                            )

                        # Apply recovery rate
                        retention *= self.protocol.recovery_rate

                        # Small viruses have some additional loss
                        retention *= self.rng.normal(0.85, 0.10)
                        retention = np.clip(retention, 0.60, 0.95)

                        removal = 1.0 - retention
                    else:
                        # Column-based: moderate retention
                        retention = self.protocol.recovery_rate * 0.75
                        removal = 1.0 - retention
                else:
                    # No filtration: PhiX retained
                    removal = 0.0

            else:  # OTHER or unknown
                # Apply general contamination reduction
                removal = self.protocol.contamination_reduction * self.rng.normal(1.0, 0.15)
                removal = np.clip(removal, 0.50, 0.99)

            # Calculate reduced abundance
            reduced_abundance = original_abundance * (1.0 - removal)

            # Track statistics
            if ctype not in original_abundances:
                original_abundances[ctype] = 0.0
                reduced_abundances[ctype] = 0.0
                reduction_factors_by_type[ctype] = []

            original_abundances[ctype] += original_abundance
            reduced_abundances[ctype] += reduced_abundance
            reduction_factors_by_type[ctype].append(removal)

            # Create reduced contaminant (copy with new abundance)
            from ..core.contamination import ContaminantGenome
            reduced_contaminant = ContaminantGenome(
                genome_id=contaminant.genome_id,
                sequence=contaminant.sequence,
                contaminant_type=contaminant.contaminant_type,
                description=contaminant.description,
                organism=contaminant.organism,
                source=contaminant.source,
                abundance=reduced_abundance,
                gc_content=contaminant.gc_content,
                metadata=contaminant.metadata.copy()
            )

            reduced_profile.add_contaminant(reduced_contaminant)

        # Calculate reduction statistics
        stats = {
            'protocol': self.protocol.name,
            'original_total_contamination': contamination_profile.get_total_abundance(),
            'reduced_total_contamination': reduced_profile.get_total_abundance(),
            'overall_reduction_factor': 1.0 - (reduced_profile.get_total_abundance() /
                                              max(contamination_profile.get_total_abundance(), 1e-10)),
            'reduction_by_type': {},
            'mean_removal_by_type': {}
        }

        for ctype in original_abundances:
            original = original_abundances[ctype]
            reduced = reduced_abundances[ctype]

            if original > 0:
                reduction = 1.0 - (reduced / original)
            else:
                reduction = 0.0

            stats['reduction_by_type'][ctype.value] = {
                'original_abundance': float(original),
                'reduced_abundance': float(reduced),
                'reduction_factor': float(reduction),
                'reduction_pct': float(reduction * 100)
            }

            stats['mean_removal_by_type'][ctype.value] = float(np.mean(reduction_factors_by_type[ctype]))

        logger.info(f"Contamination reduction applied: {stats['overall_reduction_factor']*100:.1f}% total removal")
        logger.info(f"Original: {stats['original_total_contamination']:.4f} → "
                   f"Reduced: {stats['reduced_total_contamination']:.4f}")

        return reduced_profile, stats

    def compare_to_bulk(
        self,
        genomes: List[Dict],
        abundances: np.ndarray
    ) -> Dict:
        """
        Compare VLP-enriched to bulk (no VLP) composition

        Returns:
            Dictionary with comparison metrics
        """
        # Get VLP-enriched abundances
        vlp_abundances, vlp_stats = self.apply_enrichment(genomes, abundances)

        # Calculate fold-changes
        fold_changes = np.where(abundances > 0, vlp_abundances / abundances, 0)

        comparison = {
            'vlp_stats': vlp_stats,
            'mean_fold_change': float(np.mean(fold_changes[fold_changes > 0])),
            'median_fold_change': float(np.median(fold_changes[fold_changes > 0])),
            'max_enrichment': float(np.max(fold_changes)),
            'max_depletion': float(np.min(fold_changes[fold_changes > 0])) if np.any(fold_changes > 0) else 0.0,
            'genomes_lost': int(np.sum((abundances > 0) & (vlp_abundances < 0.0001))),  # Below detection
            'genomes_gained': 0  # VLP doesn't add genomes, only retains/depletes
        }

        return comparison
