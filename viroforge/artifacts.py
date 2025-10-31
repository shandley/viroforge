"""
Platform-Specific Sequencing Artifacts
=======================================

This module models platform-specific artifacts introduced during Illumina sequencing.
Different platforms have different artifact profiles that can significantly impact
virome analysis.

Key artifacts modeled:
- **PolyG tails**: Patterned flow cells (NovaSeq, NextSeq)
- **Optical duplicates**: Adjacent cluster duplicates
- **Index hopping**: Barcode misassignment

Literature basis:
- Costello et al. (2018) "Characterization and remediation of sample index swaps"
  BMC Genomics 19:332
- Chen et al. (2017) "NovaSeq 6000 System Sequencing Quality"
  Illumina Technical Note
- Sinha et al. (2017) "Index switching causes spurious variant calls"
  Genome Res 27:1962-1970

Author: ViroForge Development Team
Date: 2025-10-31
"""

import logging
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from collections import defaultdict

logger = logging.getLogger(__name__)


@dataclass
class ReadPair:
    """
    Represents a paired-end read with metadata.

    Attributes:
        read_id: Read identifier
        forward_seq: R1 sequence
        reverse_seq: R2 sequence (optional for single-end)
        forward_qual: R1 quality scores
        reverse_qual: R2 quality scores
        sample_index: Sample barcode/index
        genome_id: Source genome (ground truth)
        tile_x: Flow cell tile X coordinate
        tile_y: Flow cell tile Y coordinate
    """
    read_id: str
    forward_seq: str
    reverse_seq: Optional[str] = None
    forward_qual: Optional[str] = None
    reverse_qual: Optional[str] = None
    sample_index: str = "NNNNNNNN"
    genome_id: str = ""
    tile_x: int = 0
    tile_y: int = 0

    def __post_init__(self):
        """Set default quality scores if not provided."""
        if self.forward_qual is None:
            self.forward_qual = 'I' * len(self.forward_seq)
        if self.reverse_seq is not None and self.reverse_qual is None:
            self.reverse_qual = 'I' * len(self.reverse_seq)


class PlatformArtifact(ABC):
    """
    Abstract base class for platform-specific sequencing artifacts.

    All artifacts modify read pairs to simulate platform-specific biases
    and errors beyond the standard error models.
    """

    @abstractmethod
    def apply(self, reads: List[ReadPair], random_seed: Optional[int] = None) -> List[ReadPair]:
        """
        Apply artifact to a list of read pairs.

        Args:
            reads: List of ReadPair objects
            random_seed: Optional random seed for reproducibility

        Returns:
            Modified list of ReadPair objects

        Note:
            Some artifacts may add reads (e.g., optical duplicates),
            others modify existing reads (e.g., polyG tails),
            and some may remove reads (e.g., quality filtering).
        """
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """String representation of the artifact."""
        pass


class PolyGTailArtifact(PlatformArtifact):
    """
    Models polyG tail artifacts in patterned flow cells (NovaSeq, NextSeq).

    In patterned flow cells, incomplete quenching of fluorescent signal
    can cause spurious G-calls, resulting in long runs of G bases at the
    end of reads.

    This is NOT an issue in cluster-based flow cells (MiSeq, HiSeq).

    Mechanism:
    - Green laser excitation for G and T
    - Failed quenching → continued G signal
    - Manifests as polyG tails (rarely polyT)

    Literature:
    - Chen et al. (2017) NovaSeq technical note
    - More common in R2 than R1
    - Typical length: 10-50 bp
    - Frequency: 1-5% of reads depending on chemistry version

    Attributes:
        frequency: Fraction of reads affected (0.01-0.05)
        min_length: Minimum polyG tail length in bp
        max_length: Maximum polyG tail length in bp
        r1_rate: Relative rate for R1 (default 0.3, R1 less affected)
        r2_rate: Relative rate for R2 (default 1.0, R2 more affected)
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.artifacts import PolyGTailArtifact
        >>> artifact = PolyGTailArtifact(frequency=0.02, min_length=15, max_length=40)
        >>> reads_with_polyg = artifact.apply(reads)
    """

    def __init__(
        self,
        frequency: float = 0.02,
        min_length: int = 10,
        max_length: int = 50,
        r1_rate: float = 0.3,
        r2_rate: float = 1.0,
        random_seed: Optional[int] = None
    ):
        """
        Initialize polyG tail artifact parameters.

        Args:
            frequency: Fraction of reads affected (typical 0.01-0.05)
            min_length: Minimum polyG tail length
            max_length: Maximum polyG tail length
            r1_rate: Relative rate for R1 (R1 less affected)
            r2_rate: Relative rate for R2 (R2 more affected)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if not 0 <= frequency <= 1:
            raise ValueError(f"frequency must be 0-1, got {frequency}")

        if min_length < 1:
            raise ValueError(f"min_length must be >= 1, got {min_length}")

        if max_length < min_length:
            raise ValueError(f"max_length ({max_length}) must be >= min_length ({min_length})")

        if not 0 <= r1_rate <= 1:
            raise ValueError(f"r1_rate must be 0-1, got {r1_rate}")

        if not 0 <= r2_rate <= 1:
            raise ValueError(f"r2_rate must be 0-1, got {r2_rate}")

        self.frequency = frequency
        self.min_length = min_length
        self.max_length = max_length
        self.r1_rate = r1_rate
        self.r2_rate = r2_rate
        self.random_seed = random_seed

        # Create local random generator
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized PolyG tail artifact: {self}")

    def apply(self, reads: List[ReadPair], random_seed: Optional[int] = None) -> List[ReadPair]:
        """
        Apply polyG tail artifacts to reads.

        Args:
            reads: List of ReadPair objects
            random_seed: Optional override of instance random_seed

        Returns:
            List of ReadPair objects with polyG tails added
        """
        if random_seed is not None:
            rng = np.random.default_rng(random_seed)
        else:
            rng = self.rng

        modified_reads = []
        n_r1_affected = 0
        n_r2_affected = 0

        for read in reads:
            new_read = read  # Default: no modification

            # Check R1 for polyG tail
            if rng.random() < (self.frequency * self.r1_rate):
                tail_length = rng.integers(self.min_length, self.max_length + 1)
                polyg_tail = 'G' * tail_length
                polyg_qual = 'I' * tail_length  # High quality for G calls

                new_forward_seq = read.forward_seq + polyg_tail
                new_forward_qual = read.forward_qual + polyg_qual

                new_read = ReadPair(
                    read_id=read.read_id,
                    forward_seq=new_forward_seq,
                    reverse_seq=read.reverse_seq,
                    forward_qual=new_forward_qual,
                    reverse_qual=read.reverse_qual,
                    sample_index=read.sample_index,
                    genome_id=read.genome_id,
                    tile_x=read.tile_x,
                    tile_y=read.tile_y
                )
                n_r1_affected += 1

            # Check R2 for polyG tail (if paired-end)
            if read.reverse_seq is not None and rng.random() < (self.frequency * self.r2_rate):
                tail_length = rng.integers(self.min_length, self.max_length + 1)
                polyg_tail = 'G' * tail_length
                polyg_qual = 'I' * tail_length

                new_reverse_seq = new_read.reverse_seq + polyg_tail
                new_reverse_qual = new_read.reverse_qual + polyg_qual

                new_read = ReadPair(
                    read_id=new_read.read_id,
                    forward_seq=new_read.forward_seq,
                    reverse_seq=new_reverse_seq,
                    forward_qual=new_read.forward_qual,
                    reverse_qual=new_reverse_qual,
                    sample_index=new_read.sample_index,
                    genome_id=new_read.genome_id,
                    tile_x=new_read.tile_x,
                    tile_y=new_read.tile_y
                )
                n_r2_affected += 1

            modified_reads.append(new_read)

        logger.info(f"Applied polyG tails: {n_r1_affected} R1, {n_r2_affected} R2 "
                   f"({n_r1_affected/len(reads)*100:.2f}%, {n_r2_affected/len(reads)*100:.2f}%)")

        return modified_reads

    def __repr__(self) -> str:
        return (f"PolyGTailArtifact(freq={self.frequency:.3f}, "
                f"length={self.min_length}-{self.max_length}bp)")


class OpticalDuplicateArtifact(PlatformArtifact):
    """
    Models optical duplicates from adjacent clusters on flow cell.

    Optical duplicates occur when:
    1. Fluorescent signal from one cluster "bleeds" into adjacent area
    2. Image processing identifies it as separate cluster
    3. Results in duplicate read pairs with identical sequences

    Rate depends on:
    - Cluster density (higher density = more duplicates)
    - Flow cell type (patterned vs cluster-based)
    - Platform (NovaSeq > HiSeq > MiSeq)

    Literature:
    - Illumina recommends ~5-10% for optimal cluster density
    - >15% indicates over-clustering
    - Proximity threshold: typically 100 pixels apart

    Attributes:
        rate: Fraction of reads that are optical duplicates (0.01-0.15)
        proximity_threshold: Maximum distance for optical duplicates (pixels)
        tile_size: Flow cell tile size for coordinate generation
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.artifacts import OpticalDuplicateArtifact
        >>> artifact = OpticalDuplicateArtifact(rate=0.08, proximity_threshold=100)
        >>> reads_with_dups = artifact.apply(reads)
    """

    def __init__(
        self,
        rate: float = 0.05,
        proximity_threshold: int = 100,
        tile_size: int = 10000,
        random_seed: Optional[int] = None
    ):
        """
        Initialize optical duplicate artifact parameters.

        Args:
            rate: Fraction of reads that are optical duplicates
            proximity_threshold: Max pixel distance for optical dup
            tile_size: Flow cell tile dimensions (pixels)
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if not 0 <= rate <= 0.3:
            raise ValueError(f"rate must be 0-0.3, got {rate}")

        if proximity_threshold < 1:
            raise ValueError(f"proximity_threshold must be >= 1, got {proximity_threshold}")

        if tile_size < 100:
            raise ValueError(f"tile_size must be >= 100, got {tile_size}")

        self.rate = rate
        self.proximity_threshold = proximity_threshold
        self.tile_size = tile_size
        self.random_seed = random_seed

        # Create local random generator
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        logger.info(f"Initialized Optical duplicate artifact: {self}")

    def apply(self, reads: List[ReadPair], random_seed: Optional[int] = None) -> List[ReadPair]:
        """
        Apply optical duplicate artifacts to reads.

        Creates duplicate reads with nearby coordinates.

        Args:
            reads: List of ReadPair objects
            random_seed: Optional override of instance random_seed

        Returns:
            List of ReadPair objects including optical duplicates
        """
        if random_seed is not None:
            rng = np.random.default_rng(random_seed)
        else:
            rng = self.rng

        # Assign coordinates to reads if not already set
        for read in reads:
            if read.tile_x == 0 and read.tile_y == 0:
                read.tile_x = rng.integers(0, self.tile_size)
                read.tile_y = rng.integers(0, self.tile_size)

        # Create optical duplicates
        n_duplicates = int(len(reads) * self.rate)
        duplicate_sources = rng.choice(reads, size=n_duplicates, replace=False)

        duplicated_reads = list(reads)  # Start with original reads

        for i, source_read in enumerate(duplicate_sources):
            # Create duplicate with nearby coordinates
            offset_x = rng.integers(-self.proximity_threshold, self.proximity_threshold + 1)
            offset_y = rng.integers(-self.proximity_threshold, self.proximity_threshold + 1)

            dup_x = max(0, min(self.tile_size - 1, source_read.tile_x + offset_x))
            dup_y = max(0, min(self.tile_size - 1, source_read.tile_y + offset_y))

            duplicate_read = ReadPair(
                read_id=f"{source_read.read_id}_optdup{i}",
                forward_seq=source_read.forward_seq,
                reverse_seq=source_read.reverse_seq,
                forward_qual=source_read.forward_qual,
                reverse_qual=source_read.reverse_qual,
                sample_index=source_read.sample_index,
                genome_id=source_read.genome_id,
                tile_x=dup_x,
                tile_y=dup_y
            )

            duplicated_reads.append(duplicate_read)

        logger.info(f"Added {n_duplicates} optical duplicates "
                   f"({n_duplicates/len(reads)*100:.2f}% of original reads)")

        return duplicated_reads

    def __repr__(self) -> str:
        return f"OpticalDuplicateArtifact(rate={self.rate:.3f})"


class IndexHoppingArtifact(PlatformArtifact):
    """
    Models index hopping (barcode misassignment) in multiplexed libraries.

    Index hopping occurs when:
    1. Free adapters/indexes in library pool
    2. Template switching during bridge amplification
    3. Read assigned to wrong sample in demultiplexing

    More common in:
    - Patterned flow cells (NovaSeq: ~1-2%)
    - vs cluster-based flow cells (MiSeq: ~0.1%)

    Causes:
    - ExAmp chemistry in patterned flow cells
    - Over-clustering
    - Poor library QC (excess free adapters)

    Literature:
    - Costello et al. (2018) BMC Genomics
    - Sinha et al. (2017) Genome Res
    - Can cause false variant calls in low-frequency variants

    Attributes:
        rate: Fraction of reads that hop to wrong index (0.001-0.02)
        n_indexes: Number of unique indexes in multiplexing pool
        random_seed: Random seed for reproducibility

    Example:
        >>> from viroforge.artifacts import IndexHoppingArtifact
        >>> artifact = IndexHoppingArtifact(rate=0.01, n_indexes=96)
        >>> reads_with_hopping = artifact.apply(reads)
    """

    def __init__(
        self,
        rate: float = 0.01,
        n_indexes: int = 96,
        random_seed: Optional[int] = None
    ):
        """
        Initialize index hopping artifact parameters.

        Args:
            rate: Fraction of reads with index hopping
            n_indexes: Number of unique indexes in pool
            random_seed: Random seed for reproducibility

        Raises:
            ValueError: If parameters are out of valid ranges
        """
        if not 0 <= rate <= 0.1:
            raise ValueError(f"rate must be 0-0.1, got {rate}")

        if n_indexes < 2:
            raise ValueError(f"n_indexes must be >= 2, got {n_indexes}")

        self.rate = rate
        self.n_indexes = n_indexes
        self.random_seed = random_seed

        # Create local random generator
        if random_seed is not None:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

        # Generate pool of sample indexes (simple integer IDs)
        self.index_pool = [f"INDEX{i:04d}" for i in range(n_indexes)]

        logger.info(f"Initialized Index hopping artifact: {self}")

    def apply(self, reads: List[ReadPair], random_seed: Optional[int] = None) -> List[ReadPair]:
        """
        Apply index hopping artifacts to reads.

        Randomly reassigns sample indexes to simulate hopping.

        Args:
            reads: List of ReadPair objects
            random_seed: Optional override of instance random_seed

        Returns:
            List of ReadPair objects with some indexes reassigned
        """
        if random_seed is not None:
            rng = np.random.default_rng(random_seed)
        else:
            rng = self.rng

        modified_reads = []
        n_hopped = 0

        for read in reads:
            # Check if this read hops
            if rng.random() < self.rate:
                # Assign random index from pool
                new_index = rng.choice(self.index_pool)

                hopped_read = ReadPair(
                    read_id=read.read_id,
                    forward_seq=read.forward_seq,
                    reverse_seq=read.reverse_seq,
                    forward_qual=read.forward_qual,
                    reverse_qual=read.reverse_qual,
                    sample_index=new_index,  # Changed index
                    genome_id=read.genome_id,
                    tile_x=read.tile_x,
                    tile_y=read.tile_y
                )

                modified_reads.append(hopped_read)
                n_hopped += 1
            else:
                modified_reads.append(read)

        logger.info(f"Applied index hopping: {n_hopped} reads hopped "
                   f"({n_hopped/len(reads)*100:.2f}%)")

        return modified_reads

    def __repr__(self) -> str:
        return f"IndexHoppingArtifact(rate={self.rate:.3f}, n_indexes={self.n_indexes})"


class PlatformProfile:
    """
    Bundles platform-specific artifacts into a reusable profile.

    Different Illumina platforms have different artifact profiles based on:
    - Flow cell type (patterned vs cluster-based)
    - Chemistry version
    - Cluster density
    - Read length capabilities

    Attributes:
        name: Platform name (e.g., "NovaSeq 6000")
        flow_cell_type: "patterned" or "cluster"
        artifacts: List of PlatformArtifact objects to apply
        description: Optional description of platform characteristics

    Example:
        >>> from viroforge.artifacts import PlatformProfile, PolyGTailArtifact
        >>> profile = PlatformProfile(
        ...     name="NovaSeq 6000",
        ...     flow_cell_type="patterned",
        ...     artifacts=[PolyGTailArtifact(frequency=0.02)]
        ... )
        >>> reads = profile.apply(reads)
    """

    def __init__(
        self,
        name: str,
        flow_cell_type: str,
        artifacts: List[PlatformArtifact],
        description: str = ""
    ):
        """
        Initialize platform profile.

        Args:
            name: Platform name
            flow_cell_type: "patterned" or "cluster"
            artifacts: List of artifacts to apply
            description: Optional description

        Raises:
            ValueError: If flow_cell_type is invalid
        """
        if flow_cell_type not in ["patterned", "cluster"]:
            raise ValueError(f"flow_cell_type must be 'patterned' or 'cluster', got '{flow_cell_type}'")

        self.name = name
        self.flow_cell_type = flow_cell_type
        self.artifacts = artifacts
        self.description = description

        logger.info(f"Created platform profile: {self.name} ({flow_cell_type} flow cell)")

    def apply(self, reads: List[ReadPair], random_seed: Optional[int] = None) -> List[ReadPair]:
        """
        Apply all artifacts in the profile sequentially.

        Args:
            reads: List of ReadPair objects
            random_seed: Optional random seed for reproducibility

        Returns:
            List of ReadPair objects with all artifacts applied
        """
        logger.info(f"Applying {self.name} platform artifacts...")

        modified_reads = reads
        for artifact in self.artifacts:
            modified_reads = artifact.apply(modified_reads, random_seed=random_seed)

        logger.info(f"✓ {self.name} artifacts applied: {len(reads)} → {len(modified_reads)} reads")

        return modified_reads

    def get_summary(self) -> Dict[str, any]:
        """
        Get summary statistics for this platform profile.

        Returns:
            Dictionary with platform characteristics
        """
        return {
            'name': self.name,
            'flow_cell_type': self.flow_cell_type,
            'n_artifacts': len(self.artifacts),
            'artifacts': [repr(a) for a in self.artifacts],
            'description': self.description
        }

    def __repr__(self) -> str:
        return f"PlatformProfile('{self.name}', {len(self.artifacts)} artifacts)"


# ============================================================================
# Pre-defined Platform Profiles
# ============================================================================

def novaseq_6000() -> PlatformProfile:
    """
    NovaSeq 6000 platform profile (patterned flow cell).

    Characteristics:
    - Patterned flow cell (nanowell arrays)
    - ExAmp chemistry
    - High throughput (up to 6 Tb per run)
    - PolyG tails common (2-3% of reads)
    - Higher optical duplicate rate (8-10%)
    - Index hopping more common (1-2%)

    Use case: High-throughput virome studies, large cohorts

    Returns:
        PlatformProfile configured for NovaSeq 6000

    Example:
        >>> from viroforge.artifacts import novaseq_6000
        >>> platform = novaseq_6000()
        >>> reads = platform.apply(reads)
    """
    return PlatformProfile(
        name="NovaSeq 6000",
        flow_cell_type="patterned",
        artifacts=[
            PolyGTailArtifact(frequency=0.025, min_length=15, max_length=45),
            OpticalDuplicateArtifact(rate=0.09, proximity_threshold=100),
            IndexHoppingArtifact(rate=0.015, n_indexes=96)
        ],
        description="High-throughput patterned flow cell with ExAmp chemistry"
    )


def nextseq_2000() -> PlatformProfile:
    """
    NextSeq 2000 platform profile (patterned flow cell).

    Characteristics:
    - Patterned flow cell (similar to NovaSeq)
    - Two-color chemistry
    - Medium throughput (up to 360 Gb per run)
    - PolyG tails present (1.5-2%)
    - Moderate optical duplicates (6-8%)
    - Index hopping moderate (1%)

    Use case: Mid-scale virome studies, clinical samples

    Returns:
        PlatformProfile configured for NextSeq 2000

    Example:
        >>> from viroforge.artifacts import nextseq_2000
        >>> platform = nextseq_2000()
        >>> reads = platform.apply(reads)
    """
    return PlatformProfile(
        name="NextSeq 2000",
        flow_cell_type="patterned",
        artifacts=[
            PolyGTailArtifact(frequency=0.018, min_length=12, max_length=40),
            OpticalDuplicateArtifact(rate=0.07, proximity_threshold=100),
            IndexHoppingArtifact(rate=0.01, n_indexes=96)
        ],
        description="Mid-throughput patterned flow cell with two-color chemistry"
    )


def miseq() -> PlatformProfile:
    """
    MiSeq platform profile (cluster-based flow cell).

    Characteristics:
    - Cluster-based flow cell (bridge amplification)
    - Four-color chemistry
    - Low throughput (up to 15 Gb per run)
    - NO polyG tails (cluster-based flow cell)
    - Low optical duplicates (2-3%)
    - Minimal index hopping (0.1%)
    - Longer read lengths (2x300 bp)

    Use case: Small-scale virome studies, method development

    Returns:
        PlatformProfile configured for MiSeq

    Example:
        >>> from viroforge.artifacts import miseq
        >>> platform = miseq()
        >>> reads = platform.apply(reads)
    """
    return PlatformProfile(
        name="MiSeq",
        flow_cell_type="cluster",
        artifacts=[
            # No polyG tails in cluster-based flow cells
            OpticalDuplicateArtifact(rate=0.025, proximity_threshold=100),
            IndexHoppingArtifact(rate=0.001, n_indexes=96)
        ],
        description="Low-throughput cluster-based flow cell (no polyG tails)"
    )


def hiseq_2500() -> PlatformProfile:
    """
    HiSeq 2500 platform profile (cluster-based flow cell).

    Characteristics:
    - Cluster-based flow cell (bridge amplification)
    - Four-color chemistry
    - High throughput (up to 1 Tb per run)
    - NO polyG tails (cluster-based flow cell)
    - Moderate optical duplicates (4-5%)
    - Low index hopping (0.2%)

    Use case: Legacy platform, older virome datasets

    Returns:
        PlatformProfile configured for HiSeq 2500

    Example:
        >>> from viroforge.artifacts import hiseq_2500
        >>> platform = hiseq_2500()
        >>> reads = platform.apply(reads)
    """
    return PlatformProfile(
        name="HiSeq 2500",
        flow_cell_type="cluster",
        artifacts=[
            # No polyG tails in cluster-based flow cells
            OpticalDuplicateArtifact(rate=0.045, proximity_threshold=100),
            IndexHoppingArtifact(rate=0.002, n_indexes=96)
        ],
        description="High-throughput cluster-based flow cell (legacy platform)"
    )


def no_artifacts() -> PlatformProfile:
    """
    Control profile with no artifacts (ideal sequencing).

    Use case: Baseline comparison, method validation

    Returns:
        PlatformProfile with no artifacts

    Example:
        >>> from viroforge.artifacts import no_artifacts
        >>> platform = no_artifacts()
        >>> reads = platform.apply(reads)  # No changes
    """
    return PlatformProfile(
        name="Ideal (No Artifacts)",
        flow_cell_type="cluster",  # Default to cluster for compatibility
        artifacts=[],
        description="Control with no platform-specific artifacts"
    )
