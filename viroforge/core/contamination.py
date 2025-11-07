"""
Contamination profile generation module.

This module provides functionality for adding realistic contamination to
viral communities, including host DNA, rRNA, reagent bacteria, and PhiX
controls. It models the various sources of contamination commonly found
in virome sequencing experiments.

Classes:
    ContaminantGenome: Represents a contaminant sequence with metadata
    ContaminationProfile: Manages a collection of contaminants with abundances

Functions:
    create_contamination_profile: Create pre-defined contamination profiles
    add_host_contamination: Add host DNA contamination
    add_rrna_contamination: Add rRNA contamination
    add_reagent_contamination: Add reagent bacteria contamination
    add_phix_control: Add PhiX174 spike-in control
"""

from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from pathlib import Path
from enum import Enum
import random
import logging

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set up logging
logger = logging.getLogger(__name__)


class ContaminantType(Enum):
    """Enumeration of contamination types."""
    HOST_DNA = "host_dna"
    RRNA = "rrna"
    REAGENT_BACTERIA = "reagent_bacteria"
    PHIX = "phix"
    OTHER = "other"


@dataclass
class ContaminantGenome:
    """
    Represents a contaminant genome with metadata.

    Attributes:
        genome_id: Unique identifier for the contaminant
        sequence: Nucleotide sequence as BioPython Seq object
        contaminant_type: Type of contamination (host, rRNA, etc.)
        length: Genome/sequence length in base pairs
        description: Additional description
        organism: Organism name (e.g., "Homo sapiens", "Delftia acidovorans")
        source: Source of contamination (e.g., "GRCh38", "SILVA")
        abundance: Relative abundance (0 to 1)
        gc_content: GC content percentage
        metadata: Additional custom metadata
    """
    genome_id: str
    sequence: Union[Seq, str]
    contaminant_type: ContaminantType
    length: int = field(init=False)
    description: str = ""
    organism: str = "Unknown"
    source: str = "Unknown"
    abundance: float = 0.0
    gc_content: Optional[float] = None
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        """Calculate derived attributes after initialization."""
        # Convert string to Seq if needed
        if isinstance(self.sequence, str):
            self.sequence = Seq(self.sequence)

        # Calculate length
        self.length = len(self.sequence)

        # Calculate GC content if not provided
        if self.gc_content is None:
            self.gc_content = self._calculate_gc_content()

    def _calculate_gc_content(self) -> float:
        """
        Calculate GC content percentage, handling ambiguous bases.

        Ambiguous bases contribute fractionally:
        - S (G/C): 1.0 GC
        - R (A/G), K (G/T): 0.5 GC
        - Y (C/T), M (A/C): 0.5 GC
        - Other ambiguous: weighted by GC probability
        - N: excluded from both numerator and denominator

        Returns:
            GC content as percentage (0-100)
        """
        seq_str = str(self.sequence).upper()

        # IUPAC ambiguity code GC contributions
        gc_weights = {
            'G': 1.0, 'C': 1.0,
            'A': 0.0, 'T': 0.0, 'U': 0.0,
            'R': 0.5,  # A or G
            'Y': 0.5,  # C or T
            'M': 0.5,  # A or C
            'K': 0.5,  # G or T
            'S': 1.0,  # G or C
            'W': 0.0,  # A or T
            'B': 0.67, # C/G/T (not A)
            'D': 0.33, # A/G/T (not C)
            'H': 0.33, # A/C/T (not G)
            'V': 0.67, # A/C/G (not T)
            'N': None  # Exclude from calculation
        }

        gc_sum = 0.0
        valid_bases = 0

        for base in seq_str:
            weight = gc_weights.get(base, None)
            if weight is not None:
                gc_sum += weight
                valid_bases += 1

        if valid_bases == 0:
            return 0.0

        return (gc_sum / valid_bases) * 100

    def to_dict(self) -> Dict:
        """
        Convert contaminant to dictionary representation.

        Returns:
            Dictionary containing contaminant metadata
        """
        return {
            'genome_id': self.genome_id,
            'contaminant_type': self.contaminant_type.value,
            'organism': self.organism,
            'source': self.source,
            'length': self.length,
            'description': self.description,
            'abundance': self.abundance,
            'gc_content': self.gc_content,
            **self.metadata
        }

    def to_seqrecord(self) -> SeqRecord:
        """
        Convert to BioPython SeqRecord for writing to FASTA.

        Returns:
            SeqRecord object
        """
        return SeqRecord(
            self.sequence,
            id=self.genome_id,
            description=self.description
        )

    def validate(
        self,
        check_sequence: bool = True,
        check_length: bool = True,
        check_abundance: bool = True,
        check_gc: bool = False
    ) -> bool:
        """
        Validate contaminant genome data integrity.

        This method runs validation checks to ensure the contaminant data is
        consistent and valid. It's recommended to call this during development
        and testing, and can be made mandatory for production pipelines.

        Args:
            check_sequence: Validate sequence contains only valid DNA characters
            check_length: Validate length matches sequence length
            check_abundance: Validate abundance is in valid range [0, 1]
            check_gc: Validate GC content matches calculated value

        Returns:
            True if all checks pass

        Raises:
            ValueError: If any validation check fails

        Example:
            >>> from viroforge.core.contamination import ContaminantGenome, ContaminantType
            >>> contaminant = ContaminantGenome(
            ...     "test", "ATCG", ContaminantType.HOST_DNA
            ... )
            >>> contaminant.validate()  # Runs all default checks
            True
        """
        from ..utils.validation import (
            validate_sequence,
            validate_sequence_length,
            validate_abundance_range,
            validate_gc_content
        )

        if check_sequence:
            validate_sequence(self.sequence)

        if check_length:
            validate_sequence_length(self)

        if check_abundance:
            validate_abundance_range(self.abundance)

        if check_gc:
            validate_gc_content(self, warn_only=False)

        return True


class ContaminationProfile:
    """
    Manages a contamination profile with multiple contaminant sources.

    This class handles collections of contaminants with their abundances
    and provides methods for combining with viral communities.

    Attributes:
        name: Name/identifier for this contamination profile
        contaminants: List of ContaminantGenome objects
        total_contamination_pct: Total contamination as percentage (0-100)
    """

    def __init__(self, name: str = "contamination_profile"):
        """
        Initialize a contamination profile.

        Args:
            name: Name for this contamination profile
        """
        self.name = name
        self.contaminants: List[ContaminantGenome] = []
        self._abundance_normalized = False

    def add_contaminant(self, contaminant: ContaminantGenome) -> None:
        """
        Add a contaminant to the profile.

        Args:
            contaminant: ContaminantGenome object to add
        """
        self.contaminants.append(contaminant)
        self._abundance_normalized = False

    def add_contaminants(self, contaminants: List[ContaminantGenome]) -> None:
        """
        Add multiple contaminants to the profile.

        Args:
            contaminants: List of ContaminantGenome objects to add
        """
        self.contaminants.extend(contaminants)
        self._abundance_normalized = False

    def set_abundances(
        self,
        abundances: Union[List[float], np.ndarray],
        normalize: bool = True
    ) -> None:
        """
        Set abundances for all contaminants.

        Args:
            abundances: Array of abundance values (one per contaminant)
            normalize: Whether to normalize abundances to sum to 1.0

        Raises:
            ValueError: If number of abundances doesn't match contaminants
        """
        if len(abundances) != len(self.contaminants):
            raise ValueError(
                f"Number of abundances ({len(abundances)}) must match "
                f"number of contaminants ({len(self.contaminants)})"
            )

        abundances = np.array(abundances)

        if normalize:
            abundances = abundances / abundances.sum()
            self._abundance_normalized = True

        for contaminant, abundance in zip(self.contaminants, abundances):
            contaminant.abundance = float(abundance)

    def get_total_abundance(self) -> float:
        """
        Get the total abundance across all contaminants.

        Returns:
            Sum of all contaminant abundances
        """
        return sum(c.abundance for c in self.contaminants)

    def get_abundance_by_type(self) -> Dict[ContaminantType, float]:
        """
        Get total abundance for each contamination type.

        Returns:
            Dictionary mapping contamination type to total abundance
        """
        type_abundances = {}
        for contaminant in self.contaminants:
            ctype = contaminant.contaminant_type
            type_abundances[ctype] = type_abundances.get(ctype, 0.0) + contaminant.abundance

        return type_abundances

    def get_contamination_table(self) -> pd.DataFrame:
        """
        Get contamination table as pandas DataFrame.

        Returns:
            DataFrame with contaminant information
        """
        data = []
        for contaminant in self.contaminants:
            data.append({
                'genome_id': contaminant.genome_id,
                'type': contaminant.contaminant_type.value,
                'organism': contaminant.organism,
                'source': contaminant.source,
                'abundance': contaminant.abundance,
                'length': contaminant.length,
                'gc_content': contaminant.gc_content
            })

        return pd.DataFrame(data)

    def get_summary_stats(self) -> Dict:
        """
        Get summary statistics for the contamination profile.

        Returns:
            Dictionary containing summary statistics
        """
        if not self.contaminants:
            return {
                'n_contaminants': 0,
                'total_abundance': 0.0,
                'by_type': {}
            }

        type_abundances = self.get_abundance_by_type()

        return {
            'n_contaminants': len(self.contaminants),
            'total_abundance': self.get_total_abundance(),
            'by_type': {k.value: v for k, v in type_abundances.items()},
            'mean_length': np.mean([c.length for c in self.contaminants]),
            'total_length': sum(c.length for c in self.contaminants)
        }

    def export_contaminants_fasta(self, output_path: Path) -> None:
        """
        Export all contaminants to a FASTA file.

        Args:
            output_path: Path to output FASTA file
        """
        output_path = Path(output_path)

        with open(output_path, 'w') as f:
            SeqIO.write(
                [c.to_seqrecord() for c in self.contaminants],
                f,
                'fasta'
            )

        logger.info(f"Exported {len(self.contaminants)} contaminants to {output_path}")

    def export_contamination_table(self, output_path: Path) -> None:
        """
        Export contamination table to TSV file.

        Args:
            output_path: Path to output TSV file
        """
        output_path = Path(output_path)
        df = self.get_contamination_table()
        df.to_csv(output_path, sep='\t', index=False)

        logger.info(f"Exported contamination table to {output_path}")

    def __len__(self) -> int:
        """Return number of contaminants in the profile."""
        return len(self.contaminants)

    def __repr__(self) -> str:
        """String representation of the contamination profile."""
        stats = self.get_summary_stats()
        return (
            f"ContaminationProfile(name='{self.name}', "
            f"n_contaminants={stats['n_contaminants']}, "
            f"total_abundance={stats['total_abundance']:.3f})"
        )


def add_host_contamination(
    profile: ContaminationProfile,
    host_organism: str = "human",
    abundance_pct: float = 2.0,
    genome_path: Optional[Path] = None,
    n_fragments: int = 100,
    fragment_length: int = 10000,
    random_seed: Optional[int] = None
) -> ContaminationProfile:
    """
    Add host DNA contamination to a profile.

    This simulates host DNA contamination that commonly occurs in virome
    preparations when VLP enrichment is incomplete.

    Args:
        profile: ContaminationProfile to add host DNA to
        host_organism: Host organism ("human", "mouse", or other)
        abundance_pct: Percentage of total abundance (0-100)
        genome_path: Optional path to host genome FASTA
        n_fragments: Number of genome fragments to sample
        fragment_length: Length of each fragment (bp)
        random_seed: Random seed for reproducibility

    Returns:
        Updated ContaminationProfile

    Example:
        >>> profile = ContaminationProfile()
        >>> add_host_contamination(profile, "human", abundance_pct=5.0)
    """
    # Use local RNG instead of global state for thread safety
    if random_seed is not None:
        rng = np.random.default_rng(random_seed)
        random.seed(random_seed)  # Still needed for random.choice()
    else:
        rng = np.random.default_rng()

    # Map organism to genome info
    host_info = {
        'human': {
            'organism': 'Homo sapiens',
            'source': 'GRCh38',
            'gc_content': 41.0
        },
        'mouse': {
            'organism': 'Mus musculus',
            'source': 'GRCm39',
            'gc_content': 42.0
        },
        'rat': {
            'organism': 'Rattus norvegicus',
            'source': 'mRatBN7.2',
            'gc_content': 42.0
        }
    }

    organism_info = host_info.get(host_organism.lower(), {
        'organism': host_organism,
        'source': 'Unknown',
        'gc_content': 42.0
    })

    logger.info(
        f"Adding {abundance_pct}% {organism_info['organism']} host DNA contamination"
    )

    # Calculate abundance per fragment
    abundance_per_fragment = (abundance_pct / 100.0) / n_fragments

    if genome_path is not None and genome_path.exists():
        # Sample from actual genome
        logger.info(f"Sampling host DNA from {genome_path}")
        records = list(SeqIO.parse(genome_path, 'fasta'))

        for i in range(n_fragments):
            # Randomly select a chromosome/contig
            record = random.choice(records)

            # Randomly select a position
            if len(record.seq) > fragment_length:
                start = random.randint(0, len(record.seq) - fragment_length)
                fragment_seq = record.seq[start:start + fragment_length]
            else:
                fragment_seq = record.seq

            contaminant = ContaminantGenome(
                genome_id=f"host_{host_organism}_{i:04d}",
                sequence=fragment_seq,
                contaminant_type=ContaminantType.HOST_DNA,
                organism=organism_info['organism'],
                source=f"{organism_info['source']}_{record.id}",
                abundance=abundance_per_fragment,
                description=f"Host DNA fragment from {record.id}"
            )
            profile.add_contaminant(contaminant)

    else:
        # Create synthetic host DNA fragments
        logger.info("No genome path provided, creating synthetic host DNA fragments")

        for i in range(n_fragments):
            # Generate random sequence with appropriate GC content
            fragment_seq = _generate_sequence_with_gc(
                fragment_length,
                organism_info['gc_content']
            )

            contaminant = ContaminantGenome(
                genome_id=f"host_{host_organism}_{i:04d}",
                sequence=fragment_seq,
                contaminant_type=ContaminantType.HOST_DNA,
                organism=organism_info['organism'],
                source=organism_info['source'],
                abundance=abundance_per_fragment,
                description=f"Synthetic host DNA fragment"
            )
            profile.add_contaminant(contaminant)

    return profile


def add_rrna_contamination(
    profile: ContaminationProfile,
    abundance_pct: float = 5.0,
    rrna_database_path: Optional[Path] = None,
    domains: List[str] = ["bacteria", "archaea", "eukaryota"],
    n_sequences: int = 50,
    random_seed: Optional[int] = None
) -> ContaminationProfile:
    """
    Add rRNA contamination to a profile.

    This simulates rRNA contamination from incomplete nuclease treatment
    or DNase/RNase contamination.

    Args:
        profile: ContaminationProfile to add rRNA to
        abundance_pct: Percentage of total abundance (0-100)
        rrna_database_path: Optional path to rRNA database (e.g., SILVA)
        domains: List of domains to include (bacteria, archaea, eukaryota)
        n_sequences: Number of rRNA sequences to sample
        random_seed: Random seed for reproducibility

    Returns:
        Updated ContaminationProfile

    Example:
        >>> profile = ContaminationProfile()
        >>> add_rrna_contamination(profile, abundance_pct=3.0)
    """
    # Use local RNG instead of global state for thread safety
    if random_seed is not None:
        rng = np.random.default_rng(random_seed)
        random.seed(random_seed)  # Still needed for random.choice()
    else:
        rng = np.random.default_rng()

    logger.info(f"Adding {abundance_pct}% rRNA contamination")

    # Calculate abundance per sequence
    abundance_per_sequence = (abundance_pct / 100.0) / n_sequences

    if rrna_database_path is not None and rrna_database_path.exists():
        # Sample from actual database
        logger.info(f"Sampling rRNA from {rrna_database_path}")
        records = list(SeqIO.parse(rrna_database_path, 'fasta'))

        sampled_records = random.sample(records, min(n_sequences, len(records)))

        for i, record in enumerate(sampled_records):
            contaminant = ContaminantGenome(
                genome_id=f"rrna_{i:04d}",
                sequence=record.seq,
                contaminant_type=ContaminantType.RRNA,
                organism=record.id,
                source="SILVA",
                abundance=abundance_per_sequence,
                description=record.description
            )
            profile.add_contaminant(contaminant)

    else:
        # Create synthetic rRNA sequences
        logger.info("No database path provided, creating synthetic rRNA sequences")

        # Typical rRNA lengths
        rrna_types = {
            '16S': 1500,  # Bacterial/Archaeal small subunit
            '18S': 1800,  # Eukaryotic small subunit
            '23S': 2900,  # Bacterial/Archaeal large subunit
            '28S': 3400,  # Eukaryotic large subunit
        }

        sequences_per_domain = n_sequences // len(domains)

        for domain in domains:
            for i in range(sequences_per_domain):
                # Select random rRNA type
                rrna_type = random.choice(list(rrna_types.keys()))
                length = rrna_types[rrna_type]

                # rRNA typically has higher GC content
                gc_content = rng.normal(55.0, 5.0)
                gc_content = np.clip(gc_content, 40.0, 70.0)

                seq = _generate_sequence_with_gc(length, gc_content)

                contaminant = ContaminantGenome(
                    genome_id=f"rrna_{domain}_{rrna_type}_{i:04d}",
                    sequence=seq,
                    contaminant_type=ContaminantType.RRNA,
                    organism=f"{domain.capitalize()} {rrna_type} rRNA",
                    source="SILVA_synthetic",
                    abundance=abundance_per_sequence,
                    description=f"Synthetic {domain} {rrna_type} rRNA"
                )
                profile.add_contaminant(contaminant)

    return profile


def add_reagent_contamination(
    profile: ContaminationProfile,
    abundance_pct: float = 0.5,
    reagent_database_path: Optional[Path] = None,
    species: Optional[List[str]] = None,
    n_genomes: int = 5,
    random_seed: Optional[int] = None
) -> ContaminationProfile:
    """
    Add reagent bacteria contamination to a profile.

    This simulates contamination from reagents and kits, based on
    Salter et al. 2014 contamination survey.

    Args:
        profile: ContaminationProfile to add reagent bacteria to
        abundance_pct: Percentage of total abundance (0-100)
        reagent_database_path: Optional path to reagent bacteria genomes
        species: List of species to include (default: common contaminants)
        n_genomes: Number of genomes to include
        random_seed: Random seed for reproducibility

    Returns:
        Updated ContaminationProfile

    Example:
        >>> profile = ContaminationProfile()
        >>> add_reagent_contamination(profile, abundance_pct=0.5)
    """
    # Use local RNG instead of global state for thread safety
    if random_seed is not None:
        rng = np.random.default_rng(random_seed)
        random.seed(random_seed)  # Still needed for random.choice()
    else:
        rng = np.random.default_rng()

    # Common reagent contaminants from Salter et al. 2014
    if species is None:
        species = [
            'Delftia acidovorans',
            'Ralstonia pickettii',
            'Burkholderia cepacia',
            'Bradyrhizobium japonicum',
            'Pseudomonas fluorescens'
        ]

    logger.info(f"Adding {abundance_pct}% reagent bacteria contamination")

    # Calculate abundance per genome
    abundance_per_genome = (abundance_pct / 100.0) / n_genomes

    if reagent_database_path is not None and reagent_database_path.exists():
        # Sample from actual database
        logger.info(f"Sampling reagent bacteria from {reagent_database_path}")
        records = list(SeqIO.parse(reagent_database_path, 'fasta'))

        sampled_records = random.sample(records, min(n_genomes, len(records)))

        for i, record in enumerate(sampled_records):
            contaminant = ContaminantGenome(
                genome_id=f"reagent_{i:04d}",
                sequence=record.seq,
                contaminant_type=ContaminantType.REAGENT_BACTERIA,
                organism=record.id,
                source="Reagent_contamination",
                abundance=abundance_per_genome,
                description=record.description
            )
            profile.add_contaminant(contaminant)

    else:
        # Create synthetic reagent bacterial genomes
        logger.info("No database path provided, creating synthetic reagent bacteria")

        selected_species = random.sample(species, min(n_genomes, len(species)))

        for i, sp in enumerate(selected_species):
            # Typical bacterial genome: 3-6 Mbp
            genome_length = int(rng.uniform(3e6, 6e6))

            # Bacterial genomes typically 40-60% GC
            gc_content = rng.uniform(40.0, 60.0)

            seq = _generate_sequence_with_gc(genome_length, gc_content)

            contaminant = ContaminantGenome(
                genome_id=f"reagent_{sp.replace(' ', '_')}_{i:04d}",
                sequence=seq,
                contaminant_type=ContaminantType.REAGENT_BACTERIA,
                organism=sp,
                source="Reagent_kit_contamination",
                abundance=abundance_per_genome,
                description=f"Reagent contamination: {sp}"
            )
            profile.add_contaminant(contaminant)

    return profile


def add_phix_control(
    profile: ContaminationProfile,
    abundance_pct: float = 0.1,
    phix_sequence_path: Optional[Path] = None
) -> ContaminationProfile:
    """
    Add PhiX174 control spike-in to a profile.

    PhiX174 is commonly added as a sequencing control for Illumina runs.

    Args:
        profile: ContaminationProfile to add PhiX to
        abundance_pct: Percentage of total abundance (0-100)
        phix_sequence_path: Optional path to PhiX174 genome (NC_001422.1)

    Returns:
        Updated ContaminationProfile

    Example:
        >>> profile = ContaminationProfile()
        >>> add_phix_control(profile, abundance_pct=1.0)
    """
    logger.info(f"Adding {abundance_pct}% PhiX174 control")

    if phix_sequence_path is not None and phix_sequence_path.exists():
        # Read actual PhiX genome
        record = next(SeqIO.parse(phix_sequence_path, 'fasta'))
        phix_seq = record.seq
    else:
        # Use known PhiX174 genome (NC_001422.1) - 5,386 bp, 44.8% GC
        # For now, create synthetic sequence with correct properties
        logger.info("No PhiX path provided, creating synthetic PhiX sequence")
        phix_seq = _generate_sequence_with_gc(5386, 44.8)

    contaminant = ContaminantGenome(
        genome_id="NC_001422.1_PhiX174",
        sequence=phix_seq,
        contaminant_type=ContaminantType.PHIX,
        organism="Enterobacteria phage phiX174",
        source="Illumina_control",
        abundance=abundance_pct / 100.0,
        description="PhiX174 sequencing control"
    )

    profile.add_contaminant(contaminant)

    return profile


def create_contamination_profile(
    profile_type: str = "realistic",
    random_seed: Optional[int] = None,
    **kwargs
) -> ContaminationProfile:
    """
    Create a pre-defined contamination profile.

    Args:
        profile_type: Type of profile to create:
            - 'clean': Minimal contamination (VLP worked well)
            - 'realistic': Moderate contamination (typical virome)
            - 'heavy': High contamination (VLP enrichment failed)
            - 'failed': Very high contamination (complete failure)
        random_seed: Random seed for reproducibility
        **kwargs: Override default parameters for contamination levels

    Returns:
        ContaminationProfile with specified contamination levels

    Example:
        >>> profile = create_contamination_profile('realistic')
        >>> print(profile.get_summary_stats())
    """
    # Use local RNG instead of global state for thread safety
    if random_seed is not None:
        rng = np.random.default_rng(random_seed)
        random.seed(random_seed)  # Still needed for random.choice()
    else:
        rng = np.random.default_rng()

    # Pre-defined contamination profiles based on ViromeQC paper
    profile_configs = {
        'clean': {
            'name': 'clean_virome',
            'description': 'Clean VLP preparation with minimal contamination',
            'host_dna_pct': 0.1,
            'rrna_pct': 0.5,
            'reagent_pct': 0.01,
            'phix_pct': 0.1,
        },
        'realistic': {
            'name': 'realistic_virome',
            'description': 'Typical virome with moderate contamination',
            'host_dna_pct': 2.0,
            'rrna_pct': 5.0,
            'reagent_pct': 0.5,
            'phix_pct': 0.1,
        },
        'heavy': {
            'name': 'heavy_contamination',
            'description': 'Poor VLP enrichment with high contamination',
            'host_dna_pct': 10.0,
            'rrna_pct': 15.0,
            'reagent_pct': 2.0,
            'phix_pct': 0.1,
        },
        'failed': {
            'name': 'failed_vlp',
            'description': 'Failed VLP preparation, very high contamination',
            'host_dna_pct': 15.0,
            'rrna_pct': 20.0,
            'reagent_pct': 5.0,
            'phix_pct': 0.1,
        }
    }

    if profile_type not in profile_configs:
        raise ValueError(
            f"Unknown profile type: {profile_type}. "
            f"Choose from: {', '.join(profile_configs.keys())}"
        )

    config = profile_configs[profile_type]

    # Allow kwargs to override defaults
    host_dna_pct = kwargs.get('host_dna_pct', config['host_dna_pct'])
    rrna_pct = kwargs.get('rrna_pct', config['rrna_pct'])
    reagent_pct = kwargs.get('reagent_pct', config['reagent_pct'])
    phix_pct = kwargs.get('phix_pct', config['phix_pct'])

    # Create profile
    profile = ContaminationProfile(name=config['name'])

    # Add contaminants
    add_host_contamination(
        profile,
        host_organism=kwargs.get('host_organism', 'human'),
        abundance_pct=host_dna_pct,
        random_seed=random_seed
    )

    add_rrna_contamination(
        profile,
        abundance_pct=rrna_pct,
        random_seed=random_seed
    )

    add_reagent_contamination(
        profile,
        abundance_pct=reagent_pct,
        random_seed=random_seed
    )

    add_phix_control(
        profile,
        abundance_pct=phix_pct
    )

    logger.info(
        f"Created '{profile_type}' contamination profile: "
        f"Host={host_dna_pct}%, rRNA={rrna_pct}%, "
        f"Reagent={reagent_pct}%, PhiX={phix_pct}%"
    )

    return profile


def _generate_sequence_with_gc(length: int, gc_content: float) -> str:
    """
    Generate a random DNA sequence with specified GC content.

    Args:
        length: Length of sequence to generate
        gc_content: Desired GC content (0-100)

    Returns:
        Random DNA sequence string
    """
    gc_fraction = gc_content / 100.0

    # Calculate number of G/C and A/T bases
    n_gc = int(length * gc_fraction)
    n_at = length - n_gc

    # Split G/C and A/T evenly
    n_g = n_gc // 2
    n_c = n_gc - n_g
    n_a = n_at // 2
    n_t = n_at - n_a

    # Create sequence and shuffle
    seq_list = ['G'] * n_g + ['C'] * n_c + ['A'] * n_a + ['T'] * n_t
    random.shuffle(seq_list)

    return ''.join(seq_list)


if __name__ == "__main__":
    # Example usage and testing
    print("ViroForge Contamination Module - Example Usage\n")

    # Example 1: Create pre-defined contamination profiles
    print("Example 1: Creating pre-defined contamination profiles")
    print("=" * 60)

    for profile_type in ['clean', 'realistic', 'heavy']:
        profile = create_contamination_profile(profile_type, random_seed=42)
        print(f"\n{profile_type.upper()} Profile:")
        print(f"  {profile}")

        stats = profile.get_summary_stats()
        print(f"  Total contaminants: {stats['n_contaminants']}")
        print(f"  Total abundance: {stats['total_abundance']:.3f}")
        print(f"  By type:")
        for ctype, abundance in stats['by_type'].items():
            print(f"    {ctype}: {abundance*100:.2f}%")

    # Example 2: Custom contamination profile
    print("\n\nExample 2: Creating custom contamination profile")
    print("=" * 60)

    custom_profile = ContaminationProfile(name="custom_contamination")

    add_host_contamination(custom_profile, "human", abundance_pct=5.0, random_seed=42)
    add_rrna_contamination(custom_profile, abundance_pct=3.0, random_seed=42)
    add_phix_control(custom_profile, abundance_pct=1.0)

    print(f"\n{custom_profile}")
    print(f"Summary: {custom_profile.get_summary_stats()}")
