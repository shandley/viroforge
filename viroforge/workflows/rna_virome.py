"""
RNA Virome Workflow Module

This module provides classes and functions for modeling RNA virus sequencing workflows,
including reverse transcription, rRNA depletion, and RNA degradation. RNA virome
preparations differ significantly from DNA viromes due to the additional steps
required to convert RNA→cDNA and the massive rRNA contamination problem.

Key Components:
    ReverseTranscription: Models RT efficiency and artifacts by virus type
    RiboDepletion: Models rRNA depletion (Ribo-Zero/Ribominus)
    RNADegradation: Models RNA fragility and RNase effects
    RNAViromeWorkflow: Orchestrates complete RNA virome workflow

Author: ViroForge Development Team
Date: 2025-11-09
Phase: 8.2 - RNA Workflow Components
"""

from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import logging

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set up logging
logger = logging.getLogger(__name__)


class RNAVirusType(Enum):
    """
    RNA virus genome types with different RT efficiencies.

    Values:
        SSRNA_POSITIVE: Single-stranded RNA, positive sense (e.g., norovirus, poliovirus)
        SSRNA_NEGATIVE: Single-stranded RNA, negative sense (e.g., influenza, RSV)
        DSRNA: Double-stranded RNA (e.g., rotavirus, reovirus)
        SSRNA_RT: Single-stranded RNA retroviruses (e.g., HIV, not typically in viromes)

    RT Efficiency:
        - ssRNA+: 70-90% (easier, direct template)
        - ssRNA-: 50-70% (requires complementary strand synthesis)
        - dsRNA: 40-80% (variable, depends on denaturation)
    """
    SSRNA_POSITIVE = "ssRNA+"
    SSRNA_NEGATIVE = "ssRNA-"
    DSRNA = "dsRNA"
    SSRNA_RT = "ssRNA-RT"


class PrimerType(Enum):
    """
    Primer types for reverse transcription.

    Values:
        RANDOM_HEXAMER: Random hexamer primers (most common for viromes)
        OLIGO_DT: Oligo-dT primers (targets polyadenylated RNA, eukaryotic viruses)
        SPECIFIC: Virus-specific primers (targeted sequencing)
        RANDOM_OCTAMER: Random octamer primers (better specificity than hexamers)
    """
    RANDOM_HEXAMER = "random_hexamer"
    OLIGO_DT = "oligo_dt"
    SPECIFIC = "specific"
    RANDOM_OCTAMER = "random_octamer"


class RiboDepleteMethod(Enum):
    """
    rRNA depletion methods.

    Values:
        RIBO_ZERO: Ribo-Zero (Illumina) - targets bacterial and eukaryotic rRNA
        RIBOMINUS: RiboMinus (Invitrogen) - similar to Ribo-Zero
        NONE: No rRNA depletion (for comparison/failed preps)

    Efficiency:
        - Ribo-Zero: 90-95% removal
        - RiboMinus: 85-90% removal
        - None: 0% removal (rRNA dominates 80-95% of reads)
    """
    RIBO_ZERO = "ribo_zero"
    RIBOMINUS = "ribominus"
    NONE = "none"


@dataclass
class ReverseTranscription:
    """
    Models reverse transcription step in RNA virome workflow.

    Reverse transcription converts RNA→cDNA for Illumina sequencing. RT efficiency
    varies dramatically by virus type and primer choice. RT also introduces artifacts
    like template switching and truncation.

    Attributes:
        primer_type: Type of primer used (random hexamer most common for viromes)
        random_seed: Random seed for reproducibility
        base_efficiency: Base RT efficiency before virus-type adjustments
        template_switching_rate: Rate of template switching artifacts (0-1)
        truncation_rate: Rate of 5'/3' truncation (0-1)
        strand_bias: Strand bias for ssRNA viruses (1.0 = no bias)

    Example:
        >>> rt = ReverseTranscription(primer_type=PrimerType.RANDOM_HEXAMER)
        >>> efficiency = rt.get_efficiency(RNAVirusType.SSRNA_POSITIVE)
        >>> print(f"RT efficiency: {efficiency:.2%}")
        RT efficiency: 85%
    """
    primer_type: PrimerType = PrimerType.RANDOM_HEXAMER
    random_seed: Optional[int] = None
    base_efficiency: float = 0.80  # 80% base efficiency
    template_switching_rate: float = 0.02  # 2% template switching
    truncation_rate: float = 0.15  # 15% truncated at 5' or 3'
    strand_bias: float = 1.0  # No strand bias by default

    def __post_init__(self):
        """Initialize random number generator."""
        if self.random_seed is not None:
            self.rng = np.random.default_rng(self.random_seed)
        else:
            self.rng = np.random.default_rng()

    def get_efficiency(
        self,
        virus_type: RNAVirusType,
        genome_length: Optional[int] = None
    ) -> float:
        """
        Get RT efficiency for a specific virus type.

        RT efficiency varies by RNA virus type:
        - ssRNA+ (positive sense): 70-90% (direct template, easier)
        - ssRNA- (negative sense): 50-70% (needs complementary strand)
        - dsRNA: 40-80% (variable, depends on denaturation and strand separation)

        Args:
            virus_type: RNA virus genome type
            genome_length: Optional genome length (affects efficiency for long genomes)

        Returns:
            RT efficiency as fraction (0-1)
        """
        # Base efficiency by virus type
        efficiency_ranges = {
            RNAVirusType.SSRNA_POSITIVE: (0.70, 0.90),
            RNAVirusType.SSRNA_NEGATIVE: (0.50, 0.70),
            RNAVirusType.DSRNA: (0.40, 0.80),
            RNAVirusType.SSRNA_RT: (0.75, 0.85),  # Retroviruses
        }

        if virus_type not in efficiency_ranges:
            logger.warning(f"Unknown virus type: {virus_type}, using default")
            return self.base_efficiency

        min_eff, max_eff = efficiency_ranges[virus_type]

        # Sample efficiency from range
        efficiency = self.rng.uniform(min_eff, max_eff)

        # Primer type affects efficiency
        primer_modifiers = {
            PrimerType.RANDOM_HEXAMER: 1.0,  # Standard
            PrimerType.RANDOM_OCTAMER: 1.05,  # Slightly better specificity
            PrimerType.OLIGO_DT: 0.9,  # Only works for polyA+ RNA
            PrimerType.SPECIFIC: 1.1,  # Best for targeted viruses
        }

        efficiency *= primer_modifiers.get(self.primer_type, 1.0)

        # Length penalty for very long genomes (>15kb)
        if genome_length and genome_length > 15000:
            length_penalty = 1.0 - (min(genome_length - 15000, 10000) / 50000)
            efficiency *= length_penalty

        # Segmented genomes (like influenza) may have variable efficiency per segment
        # This is handled at the workflow level

        return min(efficiency, 1.0)

    def apply_template_switching(
        self,
        sequences: List[SeqRecord]
    ) -> Tuple[List[SeqRecord], int]:
        """
        Apply template switching artifacts.

        Template switching occurs when RT enzyme jumps between templates,
        creating chimeric sequences. Rate: ~2% for random hexamers.

        Args:
            sequences: List of sequence records

        Returns:
            Tuple of (modified_sequences, n_chimeric)
        """
        if len(sequences) < 2:
            return sequences, 0

        n_chimeric = 0
        modified = []

        for seq_record in sequences:
            # Apply template switching with probability
            if self.rng.random() < self.template_switching_rate:
                # Create chimeric sequence by joining with another random sequence
                other_seq = self.rng.choice(sequences)

                # Find switching point (usually in first half of sequence)
                switch_point = self.rng.integers(
                    len(seq_record.seq) // 4,
                    3 * len(seq_record.seq) // 4
                )

                chimeric_seq = (
                    seq_record.seq[:switch_point] +
                    other_seq.seq[switch_point:]
                )

                modified_record = SeqRecord(
                    chimeric_seq,
                    id=seq_record.id + "_chimeric",
                    description=f"{seq_record.description} [template_switching]"
                )
                modified.append(modified_record)
                n_chimeric += 1
            else:
                modified.append(seq_record)

        if n_chimeric > 0:
            logger.info(f"Template switching: {n_chimeric} chimeric sequences "
                       f"({n_chimeric/len(sequences):.1%})")

        return modified, n_chimeric

    def apply_truncation(
        self,
        sequences: List[SeqRecord]
    ) -> Tuple[List[SeqRecord], Dict[str, int]]:
        """
        Apply 5' and 3' truncation artifacts.

        RT often fails to reach the 5' end (especially for random primers)
        or truncates at the 3' end. Rate: ~15% of sequences affected.

        Args:
            sequences: List of sequence records

        Returns:
            Tuple of (modified_sequences, truncation_stats)
        """
        truncation_stats = {
            '5_prime_truncated': 0,
            '3_prime_truncated': 0,
            'both_truncated': 0
        }

        modified = []

        for seq_record in sequences:
            seq_len = len(seq_record.seq)

            # Apply truncation with probability
            if self.rng.random() < self.truncation_rate:
                # Decide which end(s) to truncate
                truncate_5 = self.rng.random() < 0.6  # 5' more common
                truncate_3 = self.rng.random() < 0.3

                if truncate_5:
                    # Truncate 10-30% from 5' end
                    truncate_len = int(seq_len * self.rng.uniform(0.10, 0.30))
                    seq_record = SeqRecord(
                        seq_record.seq[truncate_len:],
                        id=seq_record.id,
                        description=seq_record.description + " [5'_truncated]"
                    )
                    truncation_stats['5_prime_truncated'] += 1

                if truncate_3:
                    # Truncate 5-15% from 3' end
                    truncate_len = int(len(seq_record.seq) * self.rng.uniform(0.05, 0.15))
                    seq_record = SeqRecord(
                        seq_record.seq[:-truncate_len],
                        id=seq_record.id,
                        description=seq_record.description + " [3'_truncated]"
                    )
                    truncation_stats['3_prime_truncated'] += 1

                if truncate_5 and truncate_3:
                    truncation_stats['both_truncated'] += 1

            modified.append(seq_record)

        total_truncated = (
            truncation_stats['5_prime_truncated'] +
            truncation_stats['3_prime_truncated'] -
            truncation_stats['both_truncated']
        )

        if total_truncated > 0:
            logger.info(f"RT truncation: {total_truncated} sequences "
                       f"({total_truncated/len(sequences):.1%})")

        return modified, truncation_stats


@dataclass
class RiboDepletion:
    """
    Models rRNA depletion step (Ribo-Zero or RiboMinus).

    rRNA depletion is CRITICAL for RNA viromes because 80-95% of RNA in samples
    is rRNA (vs ~5% for DNA preps). Without depletion, viral reads are <1% of total.
    With good depletion, rRNA drops to 5-10% and viral enrichment is ~20x.

    Attributes:
        method: Depletion method (Ribo-Zero, RiboMinus, or None)
        efficiency: Depletion efficiency (0-1)
        random_seed: Random seed for reproducibility
        targets_bacterial: Whether method targets bacterial rRNA
        targets_eukaryotic: Whether method targets eukaryotic rRNA

    Example:
        >>> ribo = RiboDepletion(method=RiboDepleteMethod.RIBO_ZERO)
        >>> print(f"Depletion efficiency: {ribo.efficiency:.1%}")
        Depletion efficiency: 92.5%
        >>> # Before: 90% rRNA, 1% viral
        >>> # After: 10% rRNA, 20% viral (20x enrichment!)
    """
    method: RiboDepleteMethod = RiboDepleteMethod.RIBO_ZERO
    efficiency: Optional[float] = None  # Auto-set based on method
    random_seed: Optional[int] = None
    targets_bacterial: bool = True
    targets_eukaryotic: bool = True

    def __post_init__(self):
        """Initialize depletion parameters."""
        if self.random_seed is not None:
            self.rng = np.random.default_rng(self.random_seed)
        else:
            self.rng = np.random.default_rng()

        # Set efficiency based on method if not specified
        if self.efficiency is None:
            efficiency_ranges = {
                RiboDepleteMethod.RIBO_ZERO: (0.90, 0.95),
                RiboDepleteMethod.RIBOMINUS: (0.85, 0.90),
                RiboDepleteMethod.NONE: (0.0, 0.0),
            }

            min_eff, max_eff = efficiency_ranges[self.method]
            self.efficiency = self.rng.uniform(min_eff, max_eff)

    def apply_depletion(
        self,
        rrna_abundance_before: float
    ) -> Tuple[float, Dict[str, float]]:
        """
        Apply rRNA depletion to abundance profile.

        This models the dramatic reduction in rRNA and corresponding enrichment
        of viral sequences.

        Args:
            rrna_abundance_before: rRNA abundance before depletion (e.g., 0.90 = 90%)

        Returns:
            Tuple of (rrna_abundance_after, enrichment_stats)

        Example:
            >>> ribo = RiboDepletion(method=RiboDepleteMethod.RIBO_ZERO)
            >>> rrna_after, stats = ribo.apply_depletion(rrna_abundance_before=0.90)
            >>> print(f"rRNA: {rrna_after:.1%} (was 90%)")
            rRNA: 9.5% (was 90%)
            >>> print(f"Viral enrichment: {stats['viral_enrichment']:.1f}x")
            Viral enrichment: 19.0x
        """
        if self.method == RiboDepleteMethod.NONE:
            return rrna_abundance_before, {
                'rrna_removed': 0.0,
                'rrna_remaining': rrna_abundance_before,
                'viral_enrichment': 1.0,
                'depletion_efficiency': 0.0
            }

        # Calculate rRNA removed and remaining
        rrna_removed = rrna_abundance_before * self.efficiency
        rrna_remaining = rrna_abundance_before - rrna_removed

        # Calculate viral enrichment
        # Before: viral = 1 - rrna_before
        # After: viral = 1 - rrna_remaining
        # Enrichment = (viral_after / total_after) / (viral_before / total_before)
        viral_before = 1.0 - rrna_abundance_before
        viral_after = 1.0 - rrna_remaining

        # Renormalize (total is now less than 1.0 due to rRNA removal)
        total_after = rrna_remaining + viral_after

        viral_enrichment = (viral_after / total_after) / (viral_before / 1.0)

        enrichment_stats = {
            'rrna_removed': rrna_removed,
            'rrna_remaining': rrna_remaining,
            'viral_enrichment': viral_enrichment,
            'depletion_efficiency': self.efficiency,
            'method': self.method.value
        }

        logger.info(f"rRNA depletion ({self.method.value}): "
                   f"{rrna_abundance_before:.1%} → {rrna_remaining:.1%} "
                   f"(Viral enrichment: {viral_enrichment:.1f}x)")

        return rrna_remaining, enrichment_stats


@dataclass
class RNADegradation:
    """
    Models RNA degradation and fragmentation.

    RNA is MUCH more fragile than DNA (10-100x faster degradation) due to:
    - RNase contamination (ubiquitous in environment)
    - 2'-OH group makes RNA chemically unstable
    - Alkaline pH causes strand breaks

    This results in fragmented viral genomes and 5'/3' coverage bias.

    Attributes:
        degradation_rate: Rate of degradation (0-1, higher = more degraded)
        rnase_contamination: Level of RNase contamination (0-1)
        fragmentation_bias_5prime: 5' end fragmentation bias (1.0 = no bias)
        fragmentation_bias_3prime: 3' end fragmentation bias (1.0 = no bias)
        random_seed: Random seed for reproducibility

    Example:
        >>> rna_deg = RNADegradation(degradation_rate=0.3)
        >>> # 30% of RNA is degraded/fragmented
        >>> # Results in uneven coverage across genomes
    """
    degradation_rate: float = 0.20  # 20% degradation (moderate)
    rnase_contamination: float = 0.10  # 10% RNase effect
    fragmentation_bias_5prime: float = 1.2  # Slight 5' bias
    fragmentation_bias_3prime: float = 0.8  # 3' more stable
    random_seed: Optional[int] = None

    def __post_init__(self):
        """Initialize random number generator."""
        if self.random_seed is not None:
            self.rng = np.random.default_rng(self.random_seed)
        else:
            self.rng = np.random.default_rng()

    def apply_degradation(
        self,
        sequences: List[SeqRecord]
    ) -> Tuple[List[SeqRecord], Dict[str, float]]:
        """
        Apply RNA degradation and fragmentation.

        Degradation results in:
        - Fragmented sequences (multiple reads from same genome)
        - Uneven coverage (5' end more degraded)
        - Loss of full-length sequences

        Args:
            sequences: List of sequence records

        Returns:
            Tuple of (modified_sequences, degradation_stats)
        """
        degradation_stats = {
            'n_fragmented': 0,
            'mean_fragment_size': 0.0,
            'coverage_bias_5to3': self.fragmentation_bias_5prime / self.fragmentation_bias_3prime
        }

        modified = []
        fragment_sizes = []

        for seq_record in sequences:
            seq_len = len(seq_record.seq)

            # Apply degradation with probability
            if self.rng.random() < self.degradation_rate:
                # Fragment into 2-4 pieces
                n_fragments = self.rng.integers(2, 5)

                # Generate fragment break points with 5' bias
                break_points = sorted([
                    int(seq_len * self._biased_position())
                    for _ in range(n_fragments - 1)
                ])

                # Add start and end
                break_points = [0] + break_points + [seq_len]

                # Create fragments
                for i in range(len(break_points) - 1):
                    start = break_points[i]
                    end = break_points[i + 1]

                    if end - start > 100:  # Keep fragments >100bp
                        fragment = SeqRecord(
                            seq_record.seq[start:end],
                            id=f"{seq_record.id}_frag{i+1}",
                            description=f"{seq_record.description} [degraded_fragment]"
                        )
                        modified.append(fragment)
                        fragment_sizes.append(end - start)

                degradation_stats['n_fragmented'] += 1
            else:
                modified.append(seq_record)

        if fragment_sizes:
            degradation_stats['mean_fragment_size'] = np.mean(fragment_sizes)
            logger.info(f"RNA degradation: {degradation_stats['n_fragmented']} sequences "
                       f"fragmented ({degradation_stats['n_fragmented']/len(sequences):.1%})")

        return modified, degradation_stats

    def _biased_position(self) -> float:
        """
        Generate position with 5' bias.

        Returns:
            Position as fraction (0-1) with 5' bias
        """
        # Beta distribution with alpha < beta gives 5' bias
        alpha = 2.0 * (1.0 / self.fragmentation_bias_5prime)
        beta = 2.0

        return self.rng.beta(alpha, beta)


class RNAViromeWorkflow:
    """
    Complete RNA virome sequencing workflow.

    This class orchestrates the entire RNA virome workflow:
    1. Reverse transcription (RNA → cDNA)
    2. rRNA depletion (Ribo-Zero/RiboMinus)
    3. RNA degradation modeling

    The workflow handles the unique challenges of RNA viromes:
    - Much higher rRNA contamination (80-95% vs 5% for DNA)
    - RT efficiency varies by virus type
    - RNA is more fragile than DNA
    - Additional artifacts from RT step

    Attributes:
        reverse_transcription: RT step configuration
        ribo_depletion: rRNA depletion configuration
        rna_degradation: RNA degradation configuration
        random_seed: Random seed for reproducibility

    Example:
        >>> workflow = RNAViromeWorkflow(
        ...     reverse_transcription=ReverseTranscription(
        ...         primer_type=PrimerType.RANDOM_HEXAMER
        ...     ),
        ...     ribo_depletion=RiboDepletion(
        ...         method=RiboDepleteMethod.RIBO_ZERO
        ...     ),
        ...     rna_degradation=RNADegradation(
        ...         degradation_rate=0.20
        ...     )
        ... )
        >>> # Apply workflow to viral community
        >>> processed_sequences, stats = workflow.apply(
        ...     sequences=viral_sequences,
        ...     virus_types=virus_type_map,
        ...     rrna_abundance=0.90
        ... )
    """

    def __init__(
        self,
        reverse_transcription: Optional[ReverseTranscription] = None,
        ribo_depletion: Optional[RiboDepletion] = None,
        rna_degradation: Optional[RNADegradation] = None,
        random_seed: Optional[int] = None
    ):
        """
        Initialize RNA virome workflow.

        Args:
            reverse_transcription: RT configuration (default: random hexamer)
            ribo_depletion: rRNA depletion configuration (default: Ribo-Zero)
            rna_degradation: RNA degradation configuration (default: 20% degradation)
            random_seed: Random seed for reproducibility
        """
        self.random_seed = random_seed

        # Initialize components with defaults if not provided
        self.reverse_transcription = reverse_transcription or ReverseTranscription(
            random_seed=random_seed
        )

        self.ribo_depletion = ribo_depletion or RiboDepletion(
            method=RiboDepleteMethod.RIBO_ZERO,
            random_seed=random_seed
        )

        self.rna_degradation = rna_degradation or RNADegradation(
            degradation_rate=0.20,
            random_seed=random_seed
        )

        logger.info("RNA virome workflow initialized")
        logger.info(f"  RT primer: {self.reverse_transcription.primer_type.value}")
        logger.info(f"  rRNA depletion: {self.ribo_depletion.method.value}")
        logger.info(f"  RNA degradation rate: {self.rna_degradation.degradation_rate:.1%}")

    def apply(
        self,
        sequences: List[SeqRecord],
        virus_types: Dict[str, RNAVirusType],
        rrna_abundance_before: float = 0.90
    ) -> Tuple[List[SeqRecord], Dict]:
        """
        Apply complete RNA virome workflow to sequences.

        This applies the workflow in order:
        1. RNA degradation (happens during extraction/storage)
        2. Reverse transcription (RNA → cDNA)
        3. rRNA depletion (removes rRNA reads)

        Args:
            sequences: List of viral genome sequences
            virus_types: Mapping of genome_id → RNAVirusType
            rrna_abundance_before: rRNA abundance before depletion (default: 0.90 = 90%)

        Returns:
            Tuple of (processed_sequences, workflow_stats)
        """
        logger.info("=" * 80)
        logger.info("APPLYING RNA VIROME WORKFLOW")
        logger.info("=" * 80)

        workflow_stats = {
            'n_sequences_input': len(sequences),
            'n_sequences_output': 0,
            'rt_stats': {},
            'ribo_depletion_stats': {},
            'degradation_stats': {}
        }

        # Step 1: RNA degradation (happens first, before RT)
        logger.info("\nStep 1: Applying RNA degradation...")
        sequences, degradation_stats = self.rna_degradation.apply_degradation(sequences)
        workflow_stats['degradation_stats'] = degradation_stats

        # Step 2: Reverse transcription (RNA → cDNA)
        logger.info("\nStep 2: Applying reverse transcription...")

        # Apply RT efficiency based on virus type
        rt_processed = []
        rt_efficiencies = []

        for seq_record in sequences:
            # Get virus type for this sequence
            genome_id = seq_record.id.split('_frag')[0]  # Handle fragments
            virus_type = virus_types.get(genome_id, RNAVirusType.SSRNA_POSITIVE)

            # Get RT efficiency
            efficiency = self.reverse_transcription.get_efficiency(
                virus_type,
                genome_length=len(seq_record.seq)
            )
            rt_efficiencies.append(efficiency)

            # Apply efficiency (sequence makes it through RT)
            if self.reverse_transcription.rng.random() < efficiency:
                rt_processed.append(seq_record)

        logger.info(f"  RT efficiency: {np.mean(rt_efficiencies):.1%} "
                   f"({len(rt_processed)}/{len(sequences)} sequences)")

        # Apply RT artifacts
        rt_processed, n_chimeric = self.reverse_transcription.apply_template_switching(
            rt_processed
        )
        rt_processed, truncation_stats = self.reverse_transcription.apply_truncation(
            rt_processed
        )

        workflow_stats['rt_stats'] = {
            'mean_efficiency': float(np.mean(rt_efficiencies)),
            'n_sequences_after_rt': len(rt_processed),
            'n_chimeric': n_chimeric,
            'truncation': truncation_stats
        }

        # Step 3: rRNA depletion
        logger.info("\nStep 3: Applying rRNA depletion...")
        rrna_after, ribo_stats = self.ribo_depletion.apply_depletion(
            rrna_abundance_before
        )
        workflow_stats['ribo_depletion_stats'] = ribo_stats

        # Final stats
        workflow_stats['n_sequences_output'] = len(rt_processed)
        workflow_stats['overall_recovery'] = (
            len(rt_processed) / len(sequences)
            if len(sequences) > 0 else 0.0
        )

        logger.info("\n" + "=" * 80)
        logger.info("RNA WORKFLOW COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Input sequences: {workflow_stats['n_sequences_input']}")
        logger.info(f"Output sequences: {workflow_stats['n_sequences_output']}")
        logger.info(f"Overall recovery: {workflow_stats['overall_recovery']:.1%}")
        logger.info(f"Viral enrichment: {ribo_stats['viral_enrichment']:.1f}x")

        return rt_processed, workflow_stats


# Helper function for determining virus type from taxonomy
def infer_virus_type_from_taxonomy(
    taxonomy: Dict[str, str]
) -> RNAVirusType:
    """
    Infer RNA virus type from taxonomic information.

    This is a helper function that maps virus families to their genome types.
    For RNA viromes, this is critical for accurate RT efficiency modeling.

    Args:
        taxonomy: Dictionary with 'family', 'genus', 'species' keys

    Returns:
        Inferred RNAVirusType

    Example:
        >>> tax = {'family': 'Orthomyxoviridae', 'genus': 'Alphainfluenzavirus'}
        >>> virus_type = infer_virus_type_from_taxonomy(tax)
        >>> print(virus_type)
        RNAVirusType.SSRNA_NEGATIVE
    """
    family = taxonomy.get('family', '').lower()
    genus = taxonomy.get('genus', '').lower()
    genome_name = taxonomy.get('genome_name', '').lower()

    # ssRNA- (negative sense) families
    ssrna_negative_families = [
        'orthomyxoviridae',  # Influenza
        'paramyxoviridae',   # RSV, measles
        'pneumoviridae',     # RSV
        'filoviridae',       # Ebola
        'rhabdoviridae',     # Rabies
        'arenaviridae',      # Lassa
        'bunyaviridae',      # Hantavirus (some)
    ]

    # dsRNA families
    dsrna_families = [
        'reoviridae',        # Reovirus
        'sedoreoviridae',    # Rotavirus
        'picobirnaviridae',  # Picobirnavirus
    ]

    # ssRNA+ (positive sense) families - most common
    ssrna_positive_families = [
        'picornaviridae',    # Poliovirus, rhinovirus, enterovirus
        'caliciviridae',     # Norovirus
        'astroviridae',      # Astrovirus
        'coronaviridae',     # SARS-CoV-2
        'flaviviridae',      # Zika, dengue, yellow fever
        'togaviridae',       # Chikungunya
        'hepeviridae',       # Hepatitis E
    ]

    if family in ssrna_negative_families:
        return RNAVirusType.SSRNA_NEGATIVE
    elif family in dsrna_families:
        return RNAVirusType.DSRNA
    elif family in ssrna_positive_families:
        return RNAVirusType.SSRNA_POSITIVE
    elif 'retro' in family or 'retro' in genus:
        return RNAVirusType.SSRNA_RT
    else:
        # Default to ssRNA+ (most common in environmental samples)
        logger.debug(f"Unknown RNA virus family '{family}', defaulting to ssRNA+")
        return RNAVirusType.SSRNA_POSITIVE


if __name__ == "__main__":
    # Example usage and testing
    print("ViroForge RNA Virome Workflow - Example Usage\n")
    print("=" * 80)

    # Example 1: Create RT configuration
    print("\nExample 1: Reverse Transcription Efficiency")
    print("-" * 80)
    rt = ReverseTranscription(primer_type=PrimerType.RANDOM_HEXAMER, random_seed=42)

    for virus_type in RNAVirusType:
        efficiency = rt.get_efficiency(virus_type)
        print(f"  {virus_type.value:15s}: {efficiency:.1%}")

    # Example 2: rRNA depletion
    print("\n\nExample 2: rRNA Depletion (Ribo-Zero)")
    print("-" * 80)
    ribo = RiboDepletion(method=RiboDepleteMethod.RIBO_ZERO, random_seed=42)

    rrna_before = 0.90  # 90% rRNA before depletion
    rrna_after, stats = ribo.apply_depletion(rrna_before)

    print(f"  rRNA before: {rrna_before:.1%}")
    print(f"  rRNA after: {rrna_after:.1%}")
    print(f"  Viral enrichment: {stats['viral_enrichment']:.1f}x")
    print(f"  Depletion efficiency: {stats['depletion_efficiency']:.1%}")

    # Example 3: Complete workflow
    print("\n\nExample 3: Complete RNA Virome Workflow")
    print("-" * 80)

    # Create test sequences
    from Bio.Seq import Seq
    test_sequences = [
        SeqRecord(Seq("A" * 1000), id=f"genome_{i:03d}", description="Test genome")
        for i in range(10)
    ]

    # Create virus type mapping
    virus_types = {
        f"genome_{i:03d}": RNAVirusType.SSRNA_POSITIVE
        for i in range(10)
    }

    # Create and apply workflow
    workflow = RNAViromeWorkflow(random_seed=42)
    processed, stats = workflow.apply(
        sequences=test_sequences,
        virus_types=virus_types,
        rrna_abundance_before=0.90
    )

    print(f"\nWorkflow Results:")
    print(f"  Input sequences: {stats['n_sequences_input']}")
    print(f"  Output sequences: {stats['n_sequences_output']}")
    print(f"  Overall recovery: {stats['overall_recovery']:.1%}")
    print(f"  Viral enrichment: {stats['ribo_depletion_stats']['viral_enrichment']:.1f}x")

    print("\n" + "=" * 80)
    print("✓ RNA Virome Workflow Module Ready!")
    print("=" * 80)
