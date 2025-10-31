"""
Library Preparation Amplification Bias Modeling
================================================

This module models amplification biases introduced during library preparation,
a critical source of bias in virome sequencing studies.

Different amplification methods introduce different biases:
- RdAB: Length and GC bias (most common for viromes)
- MDA: Extreme GC bias and stochasticity (low-biomass samples)
- Linker: Minimal bias (modern protocols)
- None: No bias (high-biomass samples)

Literature basis:
- Kim et al. (2013) "Amplification methods bias metagenomic libraries"
  Nat Methods 10:47-48
- Marine et al. (2014) "Evaluation of a transposase protocol for viral metagenomics"
  PeerJ 2:e868
- Duhaime et al. (2012) "Comparative omics and trait analyses within a genome
  of a new marine T4-like cyanophage" Environ Microbiol 14:2055-2069

Author: ViroForge Development Team
Date: 2025-01-31
"""

import logging
from abc import ABC, abstractmethod
from typing import Optional, Literal
import numpy as np

from .utils.composition import MockViromeComposition

logger = logging.getLogger(__name__)


class AmplificationMethod(ABC):
    """
    Abstract base class for library preparation amplification methods.

    All amplification methods modify genome abundances in-place to simulate
    the biases introduced during library preparation.
    """

    @abstractmethod
    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply amplification bias to a mock virome composition.

        This modifies genome abundances in-place to model the biases
        introduced by this amplification method.

        Args:
            composition: MockViromeComposition to modify

        Note:
            This method modifies the composition in-place and automatically
            renormalizes abundances to sum to 1.0.
        """
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """String representation of the amplification method."""
        pass


class RdABAmplification(AmplificationMethod):
    """
    RdAB amplification (Random RT + dsDNA synthesis + PCR).

    This is the most common amplification method for viral metagenomics,
    combining random reverse transcription, second-strand synthesis with
    Sequenase, and PCR amplification.

    Biases:
    - **Length bias**: Exponential advantage for shorter genomes
      (each PCR cycle doubles short fragments more efficiently)
    - **GC bias**: Quadratic penalty for extreme GC content
      (optimal ~50% GC, decreasing efficiency at extremes)
    - **Cycle dependency**: Bias scales with PCR cycle number
      (more cycles = stronger bias)

    Literature:
    - Standard protocol: 40 cycles typical
    - Length bias: ~1.5-2x per kb difference
    - GC bias: 2-10x difference across 30-70% GC range

    Attributes:
        cycles: Number of PCR cycles (typically 30-45)
        length_bias_strength: Intensity of length bias (0=none, 1=standard, >1=stronger)
        gc_bias_strength: Intensity of GC bias (0=none, 1=standard, >1=stronger)
        optimal_gc: GC content with maximum amplification efficiency (default 0.50)
        gc_tolerance: GC range with near-optimal efficiency (default 0.15)
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.amplification import RdABAmplification
        >>> from viroforge.utils import create_mock_virome
        >>>
        >>> composition = create_mock_virome('gut', 'realistic')
        >>> rdab = RdABAmplification(cycles=40)
        >>> rdab.apply(composition)  # Applies length and GC bias
    """

    def __init__(
        self,
        cycles: int = 40,
        length_bias_strength: float = 1.0,
        gc_bias_strength: float = 1.0,
        optimal_gc: float = 0.50,
        gc_tolerance: float = 0.15,
        random_seed: Optional[int] = None
    ):
        """
        Initialize RdAB amplification parameters.

        Args:
            cycles: PCR cycle number (30-45 typical)
            length_bias_strength: Intensity of length bias (0-2)
            gc_bias_strength: Intensity of GC bias (0-2)
            optimal_gc: GC content with peak efficiency (0-1)
            gc_tolerance: GC range for near-peak efficiency (0-0.3)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if cycles < 0 or cycles > 60:
            raise ValueError(f"cycles must be 0-60, got {cycles}")

        if length_bias_strength < 0:
            raise ValueError(f"length_bias_strength must be >= 0, got {length_bias_strength}")

        if gc_bias_strength < 0:
            raise ValueError(f"gc_bias_strength must be >= 0, got {gc_bias_strength}")

        if not 0 <= optimal_gc <= 1:
            raise ValueError(f"optimal_gc must be 0-1, got {optimal_gc}")

        if not 0 <= gc_tolerance <= 0.5:
            raise ValueError(f"gc_tolerance must be 0-0.5, got {gc_tolerance}")

        self.cycles = cycles
        self.length_bias_strength = length_bias_strength
        self.gc_bias_strength = gc_bias_strength
        self.optimal_gc = optimal_gc
        self.gc_tolerance = gc_tolerance
        self.random_seed = random_seed

        # Create local random generator for reproducibility
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized RdAB amplification: {self}")

    def calculate_length_efficiency(self, genome_length_bp: int) -> float:
        """
        Calculate length-dependent amplification efficiency.

        Shorter genomes amplify more efficiently in PCR because:
        - Faster polymerase extension
        - More complete products per cycle
        - Less secondary structure interference

        Uses exponential decay model:
        efficiency = exp(-k * length_kb)

        Args:
            genome_length_bp: Genome length in base pairs

        Returns:
            Amplification efficiency (0-1, where 1 is optimal)

        Example:
            >>> rdab = RdABAmplification(cycles=40)
            >>> rdab.calculate_length_efficiency(5000)    # Small virus
            0.95
            >>> rdab.calculate_length_efficiency(100000)  # Large virus
            0.35
        """
        if self.length_bias_strength == 0:
            return 1.0

        # Convert to kb for scaling
        length_kb = genome_length_bp / 1000.0

        # Exponential decay: shorter = better
        # Decay rate: ~0.015 per kb (empirically derived)
        decay_rate = 0.015 * self.length_bias_strength

        efficiency = np.exp(-decay_rate * length_kb)

        return efficiency

    def calculate_gc_efficiency(self, gc_content: float) -> float:
        """
        Calculate GC-dependent amplification efficiency.

        PCR efficiency is optimal around 50% GC and decreases at extremes:
        - Low GC (<30%): AT-rich regions, low melting temp, poor priming
        - High GC (>70%): Secondary structures, high melting temp, incomplete denaturation

        Uses exponential decay model centered at optimal_gc:
        efficiency = exp(-((gc - optimal) / tolerance)^2 * strength)

        Args:
            gc_content: GC content fraction (0-1)

        Returns:
            Amplification efficiency (0-1, where 1 is optimal)

        Example:
            >>> rdab = RdABAmplification(cycles=40, optimal_gc=0.50)
            >>> rdab.calculate_gc_efficiency(0.50)  # Optimal
            1.0
            >>> rdab.calculate_gc_efficiency(0.30)  # Low GC
            0.22
            >>> rdab.calculate_gc_efficiency(0.70)  # High GC
            0.22
        """
        if self.gc_bias_strength == 0:
            return 1.0

        # Calculate deviation from optimal
        gc_deviation = abs(gc_content - self.optimal_gc)

        # Exponential decay with quadratic penalty
        # Efficiency drops as we move from optimal, but never reaches exactly 0
        relative_deviation = gc_deviation / self.gc_tolerance

        # Apply exponential penalty with configurable strength
        penalty = (relative_deviation ** 2) * self.gc_bias_strength

        # Efficiency (exponential decay, approaches 0 but never reaches it)
        efficiency = np.exp(-penalty)

        return efficiency

    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply RdAB amplification bias to composition.

        Modifies genome abundances based on:
        1. Length-dependent efficiency (exponential decay)
        2. GC-dependent efficiency (quadratic around optimal)
        3. Cycle-dependent scaling (bias compounds with cycles)

        The combined bias for each genome is:
        final_abundance = initial_abundance * (length_eff * gc_eff) ^ cycles

        Args:
            composition: MockViromeComposition to modify in-place
        """
        logger.info(f"Applying RdAB amplification to {composition.name}")
        logger.info(f"  Cycles: {self.cycles}")
        logger.info(f"  Length bias strength: {self.length_bias_strength}")
        logger.info(f"  GC bias strength: {self.gc_bias_strength}")

        # Track initial state
        initial_total = composition.get_total_abundance()

        # Apply bias to viral genomes
        for genome in composition.viral_community.genomes:
            # Calculate efficiencies
            length_eff = self.calculate_length_efficiency(genome.length)
            gc_eff = self.calculate_gc_efficiency(genome.gc_content)

            # Combined efficiency per cycle
            cycle_efficiency = length_eff * gc_eff

            # Amplification factor scales exponentially with cycles
            # Each cycle multiplies by cycle_efficiency
            amplification_factor = cycle_efficiency ** self.cycles

            # Apply to abundance
            genome.abundance *= amplification_factor

        # Apply bias to contaminants (if present)
        if composition.contamination_profile is not None:
            for contaminant in composition.contamination_profile.contaminants:
                length_eff = self.calculate_length_efficiency(contaminant.length)
                gc_eff = self.calculate_gc_efficiency(contaminant.gc_content)

                cycle_efficiency = length_eff * gc_eff
                amplification_factor = cycle_efficiency ** self.cycles

                contaminant.abundance *= amplification_factor

        # Normalize abundances to sum to 1.0
        total_abundance = composition.get_total_abundance()
        if total_abundance > 0:
            # Normalize viral genomes
            for genome in composition.viral_community.genomes:
                genome.abundance /= total_abundance

            # Normalize contaminants
            if composition.contamination_profile is not None:
                for contaminant in composition.contamination_profile.contaminants:
                    contaminant.abundance /= total_abundance

        # Update viral_fraction attribute
        new_viral_fraction = composition.viral_community.get_total_abundance()
        composition.viral_fraction = new_viral_fraction
        composition.contamination_fraction = 1.0 - new_viral_fraction

        logger.info(f"  Amplification complete")
        logger.info(f"✓ RdAB amplification applied successfully")

    def __repr__(self) -> str:
        return (f"RdABAmplification(cycles={self.cycles}, "
                f"length_bias={self.length_bias_strength}, "
                f"gc_bias={self.gc_bias_strength})")


class MDAAmplification(AmplificationMethod):
    """
    MDA (Multiple Displacement Amplification) using φ29 polymerase.

    MDA is used for whole-genome amplification of low-biomass samples,
    including environmental viromes. It introduces extreme biases compared
    to PCR-based methods.

    Biases:
    - **Extreme GC bias**: 10-1000x difference across GC range
      (much stronger than PCR, φ29 struggles with high GC)
    - **High stochasticity**: Random, uneven amplification
      (some templates dominate by chance)
    - **Chimera formation**: Template switching creates artifacts (5-30%)

    Literature:
    - Lasken & Stockwell (2007) Nat Rev Microbiol 5:755-763
    - Yilmaz et al. (2010) PLoS ONE 5:e15533
    - Typical amplification time: 2-16 hours
    - GC bias stronger than any PCR method

    Attributes:
        amplification_time_hours: MDA reaction time (2-16 hours typical)
        gc_bias_strength: Intensity of GC bias (typically 2-5x stronger than PCR)
        stochasticity: Random variation (CV, typically 0.2-0.5)
        chimera_rate: Fraction of chimeric reads (0-0.3)
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.amplification import MDAAmplification
        >>>
        >>> composition = create_mock_virome('gut', 'clean')
        >>> mda = MDAAmplification(amplification_time_hours=4, gc_bias_strength=3.0)
        >>> mda.apply(composition)  # Extreme GC bias and stochasticity
    """

    def __init__(
        self,
        amplification_time_hours: float = 4.0,
        gc_bias_strength: float = 3.0,
        stochasticity: float = 0.3,
        chimera_rate: float = 0.15,
        random_seed: Optional[int] = None
    ):
        """
        Initialize MDA amplification parameters.

        Args:
            amplification_time_hours: Reaction time (2-16 hours)
            gc_bias_strength: Intensity of GC bias (2-5, stronger than RdAB)
            stochasticity: Random variation (CV, 0-1)
            chimera_rate: Fraction of chimeric products (0-0.3)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if amplification_time_hours < 0 or amplification_time_hours > 24:
            raise ValueError(f"amplification_time_hours must be 0-24, got {amplification_time_hours}")

        if gc_bias_strength < 0:
            raise ValueError(f"gc_bias_strength must be >= 0, got {gc_bias_strength}")

        if not 0 <= stochasticity <= 1:
            raise ValueError(f"stochasticity must be 0-1, got {stochasticity}")

        if not 0 <= chimera_rate <= 1:
            raise ValueError(f"chimera_rate must be 0-1, got {chimera_rate}")

        self.amplification_time = amplification_time_hours
        self.gc_bias_strength = gc_bias_strength
        self.stochasticity = stochasticity
        self.chimera_rate = chimera_rate
        self.random_seed = random_seed

        # Create local random generator
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized MDA amplification: {self}")

    def calculate_gc_efficiency(self, gc_content: float) -> float:
        """
        Calculate extreme GC-dependent efficiency for MDA.

        φ29 polymerase has much stronger GC bias than Taq polymerase:
        - Struggles with high GC content (>60%)
        - Optimal range: 30-50% GC
        - Can show 10-1000x differences

        Uses steep exponential decay for high GC:
        efficiency = exp(-k * (gc - optimal)^2)

        Args:
            gc_content: GC content fraction (0-1)

        Returns:
            Amplification efficiency (0-1)
        """
        if self.gc_bias_strength == 0:
            return 1.0

        # Optimal GC for φ29 is lower than Taq (~40%)
        optimal_gc = 0.40

        # Calculate squared deviation
        gc_deviation_sq = (gc_content - optimal_gc) ** 2

        # Steep exponential penalty
        # Stronger bias than RdAB
        penalty_rate = 15.0 * self.gc_bias_strength

        efficiency = np.exp(-penalty_rate * gc_deviation_sq)

        # Floor at very low value (not zero, some amplification always happens)
        efficiency = max(0.001, efficiency)

        return efficiency

    def add_stochastic_variation(self, base_efficiency: float) -> float:
        """
        Add random stochastic variation to efficiency.

        MDA is highly stochastic - some templates randomly amplify much more
        than others due to:
        - Random priming events
        - Template accessibility
        - Local conditions

        Uses log-normal distribution (multiplicative noise).

        Args:
            base_efficiency: Base amplification efficiency

        Returns:
            Efficiency with stochastic variation added
        """
        if self.stochasticity == 0:
            return base_efficiency

        # Log-normal multiplicative noise
        stochastic_factor = self.rng.lognormal(mean=0, sigma=self.stochasticity)

        # Apply factor
        efficiency_with_noise = base_efficiency * stochastic_factor

        # Clip to reasonable range
        efficiency_with_noise = np.clip(efficiency_with_noise, 0.0001, 10000.0)

        return efficiency_with_noise

    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply MDA amplification bias to composition.

        Modifies genome abundances based on:
        1. Extreme GC-dependent efficiency
        2. High stochastic variation
        3. Amplification time scaling

        Note: Chimera formation is tracked but chimeric reads are only
        generated during actual read simulation.

        Args:
            composition: MockViromeComposition to modify in-place
        """
        logger.info(f"Applying MDA amplification to {composition.name}")
        logger.info(f"  Amplification time: {self.amplification_time} hours")
        logger.info(f"  GC bias strength: {self.gc_bias_strength}")
        logger.info(f"  Stochasticity: {self.stochasticity}")
        logger.info(f"  Chimera rate: {self.chimera_rate}")

        # Apply bias to viral genomes
        for genome in composition.viral_community.genomes:
            # Extreme GC bias
            gc_eff = self.calculate_gc_efficiency(genome.gc_content)

            # Add stochastic variation
            final_eff = self.add_stochastic_variation(gc_eff)

            # Amplification scales with time (longer = more rounds)
            # Model as exponential growth with time-dependent exponent
            time_scaling = self.amplification_time / 4.0  # Normalize to 4h standard
            amplification_factor = final_eff ** time_scaling

            # Apply to abundance
            genome.abundance *= amplification_factor

        # Apply bias to contaminants
        if composition.contamination_profile is not None:
            for contaminant in composition.contamination_profile.contaminants:
                gc_eff = self.calculate_gc_efficiency(contaminant.gc_content)
                final_eff = self.add_stochastic_variation(gc_eff)
                time_scaling = self.amplification_time / 4.0
                amplification_factor = final_eff ** time_scaling

                contaminant.abundance *= amplification_factor

        # Normalize abundances
        total_abundance = composition.get_total_abundance()
        if total_abundance > 0:
            for genome in composition.viral_community.genomes:
                genome.abundance /= total_abundance

            if composition.contamination_profile is not None:
                for contaminant in composition.contamination_profile.contaminants:
                    contaminant.abundance /= total_abundance

        # Update viral_fraction
        new_viral_fraction = composition.viral_community.get_total_abundance()
        composition.viral_fraction = new_viral_fraction
        composition.contamination_fraction = 1.0 - new_viral_fraction

        logger.info(f"  MDA amplification complete")
        logger.info(f"✓ MDA amplification applied successfully")

    def __repr__(self) -> str:
        return (f"MDAAmplification(time={self.amplification_time}h, "
                f"gc_bias={self.gc_bias_strength}, "
                f"stochasticity={self.stochasticity})")


class LinkerAmplification(AmplificationMethod):
    """
    Linker-based amplification (adapter ligation + PCR).

    Modern library prep protocols that ligate adapters to all fragments
    before amplification, reducing but not eliminating bias.

    Biases:
    - **Minimal length bias**: Adapters on all fragments = no length advantage
    - **Moderate GC bias**: Still present but weaker than RdAB
    - **Fewer cycles**: Typically 10-25 cycles vs 40 for RdAB

    Used by: Nextera, TruSeq, and other modern library prep kits

    Literature:
    - Marine et al. (2014) PeerJ 2:e868
    - Linker-based methods show 2-5x less bias than RdAB

    Attributes:
        cycles: PCR cycles (10-25 typical)
        gc_bias_strength: GC bias intensity (weaker than RdAB, ~0.3-0.7)
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.amplification import LinkerAmplification
        >>>
        >>> composition = create_mock_virome('gut', 'realistic')
        >>> linker = LinkerAmplification(cycles=20, gc_bias_strength=0.5)
        >>> linker.apply(composition)  # Minimal bias
    """

    def __init__(
        self,
        cycles: int = 20,
        gc_bias_strength: float = 0.5,
        random_seed: Optional[int] = None
    ):
        """
        Initialize linker-based amplification parameters.

        Args:
            cycles: PCR cycles (10-25 typical, fewer than RdAB)
            gc_bias_strength: GC bias intensity (0.3-0.7, weaker than RdAB)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if cycles < 0 or cycles > 40:
            raise ValueError(f"cycles must be 0-40, got {cycles}")

        if gc_bias_strength < 0:
            raise ValueError(f"gc_bias_strength must be >= 0, got {gc_bias_strength}")

        self.cycles = cycles
        self.gc_bias_strength = gc_bias_strength
        self.random_seed = random_seed

        # Create local random generator
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized Linker amplification: {self}")

    def calculate_gc_efficiency(self, gc_content: float) -> float:
        """
        Calculate moderate GC-dependent efficiency.

        Linker-based methods still have some GC bias, but weaker than
        random priming methods.

        Args:
            gc_content: GC content fraction (0-1)

        Returns:
            Amplification efficiency (0-1)
        """
        if self.gc_bias_strength == 0:
            return 1.0

        # Similar to RdAB but with weaker penalty and wider tolerance
        optimal_gc = 0.50
        gc_tolerance = 0.20  # Wider tolerance than RdAB (0.15)

        gc_deviation = abs(gc_content - optimal_gc)
        relative_deviation = gc_deviation / gc_tolerance

        # Use exponential model like RdAB (consistent implementation)
        penalty = (relative_deviation ** 2) * self.gc_bias_strength
        efficiency = np.exp(-penalty)

        return efficiency

    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply linker-based amplification bias.

        Modifies genome abundances based on:
        1. No length bias (all fragments have adapters)
        2. Moderate GC bias (weaker than RdAB)
        3. Fewer cycles (less bias accumulation)

        Args:
            composition: MockViromeComposition to modify in-place
        """
        logger.info(f"Applying Linker amplification to {composition.name}")
        logger.info(f"  Cycles: {self.cycles}")
        logger.info(f"  GC bias strength: {self.gc_bias_strength}")

        # Apply GC bias only (no length bias)
        for genome in composition.viral_community.genomes:
            gc_eff = self.calculate_gc_efficiency(genome.gc_content)
            amplification_factor = gc_eff ** self.cycles
            genome.abundance *= amplification_factor

        if composition.contamination_profile is not None:
            for contaminant in composition.contamination_profile.contaminants:
                gc_eff = self.calculate_gc_efficiency(contaminant.gc_content)
                amplification_factor = gc_eff ** self.cycles
                contaminant.abundance *= amplification_factor

        # Normalize
        total_abundance = composition.get_total_abundance()
        if total_abundance > 0:
            for genome in composition.viral_community.genomes:
                genome.abundance /= total_abundance

            if composition.contamination_profile is not None:
                for contaminant in composition.contamination_profile.contaminants:
                    contaminant.abundance /= total_abundance

        # Update viral_fraction
        new_viral_fraction = composition.viral_community.get_total_abundance()
        composition.viral_fraction = new_viral_fraction
        composition.contamination_fraction = 1.0 - new_viral_fraction

        logger.info(f"✓ Linker amplification applied successfully")

    def __repr__(self) -> str:
        return f"LinkerAmplification(cycles={self.cycles}, gc_bias={self.gc_bias_strength})"


class NoAmplification(AmplificationMethod):
    """
    No amplification (control for high-biomass samples).

    Some high-biomass virome samples have sufficient DNA for sequencing
    without amplification, eliminating amplification bias entirely.

    This is the ideal scenario for unbiased composition, but rare in
    practice due to low viral biomass in most samples.

    Biases:
    - None (preserves original composition)

    Example:
        >>> from viroforge.amplification import NoAmplification
        >>>
        >>> composition = create_mock_virome('gut', 'clean')
        >>> no_amp = NoAmplification()
        >>> no_amp.apply(composition)  # No change to abundances
    """

    def __init__(self):
        """Initialize no amplification (no parameters needed)."""
        logger.info("Initialized NoAmplification (control)")

    def apply(self, composition: MockViromeComposition) -> None:
        """
        Apply no amplification (no-op, preserves composition).

        This method does nothing - abundances are not modified.
        Included for completeness and pipeline consistency.

        Args:
            composition: MockViromeComposition (unchanged)
        """
        logger.info(f"Applying NoAmplification to {composition.name}")
        logger.info("  No bias applied (control)")
        logger.info("✓ NoAmplification complete (composition unchanged)")

    def __repr__(self) -> str:
        return "NoAmplification()"


# Pre-defined amplification protocols for convenience

def rdab_40_cycles() -> RdABAmplification:
    """
    Standard RdAB with 40 cycles (most common for viromes).

    Strong bias, typical for standard virome protocols.

    Returns:
        RdABAmplification with 40 cycles, standard bias

    Example:
        >>> from viroforge.amplification import rdab_40_cycles
        >>> amp = rdab_40_cycles()
        >>> amp.apply(composition)
    """
    return RdABAmplification(
        cycles=40,
        length_bias_strength=1.0,
        gc_bias_strength=1.0
    )


def rdab_30_cycles() -> RdABAmplification:
    """
    RdAB with 30 cycles (moderate bias).

    Weaker bias due to fewer cycles, used when input DNA is higher.

    Returns:
        RdABAmplification with 30 cycles, standard bias

    Example:
        >>> from viroforge.amplification import rdab_30_cycles
        >>> amp = rdab_30_cycles()
        >>> amp.apply(composition)
    """
    return RdABAmplification(
        cycles=30,
        length_bias_strength=1.0,
        gc_bias_strength=1.0
    )


def mda_standard() -> MDAAmplification:
    """
    Standard MDA with 4 hour amplification.

    Extreme bias, typical for low-biomass samples.

    Returns:
        MDAAmplification with 4h, strong bias

    Example:
        >>> from viroforge.amplification import mda_standard
        >>> amp = mda_standard()
        >>> amp.apply(composition)
    """
    return MDAAmplification(
        amplification_time_hours=4.0,
        gc_bias_strength=3.0,
        stochasticity=0.3
    )


def mda_overnight() -> MDAAmplification:
    """
    Overnight MDA with 16 hour amplification.

    Very extreme bias, maximum amplification for lowest biomass samples.

    Returns:
        MDAAmplification with 16h, very strong bias

    Example:
        >>> from viroforge.amplification import mda_overnight
        >>> amp = mda_overnight()
        >>> amp.apply(composition)
    """
    return MDAAmplification(
        amplification_time_hours=16.0,
        gc_bias_strength=3.5,
        stochasticity=0.4
    )


def linker_standard() -> LinkerAmplification:
    """
    Standard linker-based amplification with 20 cycles.

    Minimal bias, typical for modern library prep kits.

    Returns:
        LinkerAmplification with 20 cycles, weak bias

    Example:
        >>> from viroforge.amplification import linker_standard
        >>> amp = linker_standard()
        >>> amp.apply(composition)
    """
    return LinkerAmplification(
        cycles=20,
        gc_bias_strength=0.5
    )


def no_amplification() -> NoAmplification:
    """
    No amplification control.

    No bias, used for high-biomass samples or as control.

    Returns:
        NoAmplification instance

    Example:
        >>> from viroforge.amplification import no_amplification
        >>> amp = no_amplification()
        >>> amp.apply(composition)  # No change
    """
    return NoAmplification()
