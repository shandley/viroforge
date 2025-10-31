"""
Viral community composition and sampling module.

This module provides functionality for creating realistic viral community
compositions with configurable abundance distributions and body-site specific
profiles. It handles genome sampling from viral databases and generates the
ground truth metadata required for validation.

Classes:
    ViralGenome: Represents a single viral genome with metadata
    ViralCommunity: Represents a complete viral community composition

Functions:
    sample_genomes_from_fasta: Sample genomes from FASTA file
    create_abundance_profile: Generate abundance distributions
    create_body_site_profile: Create body-site specific compositions
"""

from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from pathlib import Path
import random
import logging

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class ViralGenome:
    """
    Represents a viral genome with metadata and taxonomic information.

    Attributes:
        genome_id: Unique identifier for the genome (e.g., RefSeq accession)
        sequence: Nucleotide sequence as BioPython Seq object
        taxonomy: Taxonomic lineage (e.g., "Caudoviricetes;Siphoviridae;crAssphage")
        length: Genome length in base pairs
        description: Additional genome description/annotation
        family: Viral family name
        genus: Viral genus name (if available)
        species: Viral species name (if available)
        abundance: Relative abundance in the community (0 to 1)
        gc_content: GC content percentage (calculated if not provided)
        metadata: Additional custom metadata
    """
    genome_id: str
    sequence: Union[Seq, str]
    taxonomy: str
    length: int = field(init=False)
    description: str = ""
    family: str = "Unknown"
    genus: str = "Unknown"
    species: str = "Unknown"
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
        Calculate GC content percentage.

        Returns:
            GC content as percentage (0-100)
        """
        seq_str = str(self.sequence).upper()
        g_count = seq_str.count('G')
        c_count = seq_str.count('C')
        total = len(seq_str)

        if total == 0:
            return 0.0

        return ((g_count + c_count) / total) * 100

    def to_dict(self) -> Dict:
        """
        Convert genome to dictionary representation.

        Returns:
            Dictionary containing genome metadata
        """
        return {
            'genome_id': self.genome_id,
            'taxonomy': self.taxonomy,
            'length': self.length,
            'description': self.description,
            'family': self.family,
            'genus': self.genus,
            'species': self.species,
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
        Validate genome data integrity.

        This method runs validation checks to ensure the genome data is
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
            >>> genome = ViralGenome("test", "ATCG", "Test")
            >>> genome.validate()  # Runs all default checks
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


class ViralCommunity:
    """
    Represents a viral community composition with multiple genomes.

    This class manages a collection of viral genomes with their abundances
    and provides methods for sampling, abundance modeling, and exporting
    ground truth metadata.

    Attributes:
        genomes: List of ViralGenome objects
        name: Name/identifier for this community
        total_abundance: Total abundance (should sum to 1.0)
    """

    def __init__(self, name: str = "viral_community"):
        """
        Initialize a viral community.

        Args:
            name: Name for this community
        """
        self.name = name
        self.genomes: List[ViralGenome] = []
        self._abundance_normalized = False

    def add_genome(self, genome: ViralGenome) -> None:
        """
        Add a viral genome to the community.

        Args:
            genome: ViralGenome object to add
        """
        self.genomes.append(genome)
        self._abundance_normalized = False

    def add_genomes(self, genomes: List[ViralGenome]) -> None:
        """
        Add multiple viral genomes to the community.

        Args:
            genomes: List of ViralGenome objects to add
        """
        self.genomes.extend(genomes)
        self._abundance_normalized = False

    def set_abundances(
        self,
        abundances: Union[List[float], np.ndarray],
        normalize: bool = True
    ) -> None:
        """
        Set abundances for all genomes in the community.

        Args:
            abundances: Array of abundance values (one per genome)
            normalize: Whether to normalize abundances to sum to 1.0

        Raises:
            ValueError: If number of abundances doesn't match number of genomes
        """
        if len(abundances) != len(self.genomes):
            raise ValueError(
                f"Number of abundances ({len(abundances)}) must match "
                f"number of genomes ({len(self.genomes)})"
            )

        abundances = np.array(abundances)

        if normalize:
            abundances = abundances / abundances.sum()
            self._abundance_normalized = True

        for genome, abundance in zip(self.genomes, abundances):
            genome.abundance = float(abundance)

    def apply_abundance_distribution(
        self,
        distribution: str = "lognormal",
        **kwargs
    ) -> None:
        """
        Apply an abundance distribution to the community.

        Args:
            distribution: Distribution type ('lognormal', 'powerlaw', 'even')
            **kwargs: Additional parameters for the distribution
                     - For 'lognormal': mu (default=0), sigma (default=1.5)
                     - For 'powerlaw': alpha (default=1.5)
                     - For 'even': no parameters

        Raises:
            ValueError: If distribution type is not recognized
        """
        n_genomes = len(self.genomes)

        if n_genomes == 0:
            logger.warning("No genomes in community, cannot apply distribution")
            return

        abundances = create_abundance_profile(
            n_genomes=n_genomes,
            distribution=distribution,
            **kwargs
        )

        self.set_abundances(abundances, normalize=True)
        logger.info(
            f"Applied {distribution} distribution to {n_genomes} genomes"
        )

    def get_total_abundance(self) -> float:
        """
        Get the total abundance across all genomes.

        Returns:
            Sum of all genome abundances
        """
        return sum(genome.abundance for genome in self.genomes)

    def get_abundance_table(self) -> pd.DataFrame:
        """
        Get abundance table as pandas DataFrame.

        Returns:
            DataFrame with columns: genome_id, taxonomy, abundance, length
        """
        data = []
        for genome in self.genomes:
            data.append({
                'genome_id': genome.genome_id,
                'taxonomy': genome.taxonomy,
                'family': genome.family,
                'genus': genome.genus,
                'species': genome.species,
                'abundance': genome.abundance,
                'length': genome.length,
                'gc_content': genome.gc_content
            })

        return pd.DataFrame(data)

    def get_summary_stats(self) -> Dict:
        """
        Get summary statistics for the community.

        Returns:
            Dictionary containing summary statistics
        """
        if not self.genomes:
            return {
                'n_genomes': 0,
                'total_abundance': 0.0,
                'mean_abundance': 0.0,
                'median_abundance': 0.0,
                'min_abundance': 0.0,
                'max_abundance': 0.0,
                'mean_length': 0.0,
                'mean_gc': 0.0
            }

        abundances = [g.abundance for g in self.genomes]
        lengths = [g.length for g in self.genomes]
        gc_contents = [g.gc_content for g in self.genomes]

        return {
            'n_genomes': len(self.genomes),
            'total_abundance': sum(abundances),
            'mean_abundance': np.mean(abundances),
            'median_abundance': np.median(abundances),
            'min_abundance': min(abundances),
            'max_abundance': max(abundances),
            'mean_length': np.mean(lengths),
            'total_length': sum(lengths),
            'mean_gc': np.mean(gc_contents),
            'families': len(set(g.family for g in self.genomes)),
            'genera': len(set(g.genus for g in self.genomes))
        }

    def export_genomes_fasta(self, output_path: Path) -> None:
        """
        Export all genomes to a FASTA file.

        Args:
            output_path: Path to output FASTA file
        """
        output_path = Path(output_path)

        with open(output_path, 'w') as f:
            SeqIO.write(
                [genome.to_seqrecord() for genome in self.genomes],
                f,
                'fasta'
            )

        logger.info(f"Exported {len(self.genomes)} genomes to {output_path}")

    def export_abundance_table(self, output_path: Path) -> None:
        """
        Export abundance table to TSV file.

        Args:
            output_path: Path to output TSV file
        """
        output_path = Path(output_path)
        df = self.get_abundance_table()
        df.to_csv(output_path, sep='\t', index=False)

        logger.info(f"Exported abundance table to {output_path}")

    def __len__(self) -> int:
        """Return number of genomes in the community."""
        return len(self.genomes)

    def __repr__(self) -> str:
        """String representation of the community."""
        stats = self.get_summary_stats()
        return (
            f"ViralCommunity(name='{self.name}', "
            f"n_genomes={stats['n_genomes']}, "
            f"total_abundance={stats['total_abundance']:.3f})"
        )


def sample_genomes_from_fasta(
    fasta_path: Path,
    n_genomes: int,
    replace: bool = False,
    taxonomy_separator: str = "|",
    taxonomy_field: int = 1,
    random_seed: Optional[int] = None
) -> List[ViralGenome]:
    """
    Sample viral genomes from a FASTA file.

    Args:
        fasta_path: Path to FASTA file containing viral genomes
        n_genomes: Number of genomes to sample
        replace: Whether to sample with replacement
        taxonomy_separator: Separator for taxonomy in header
        taxonomy_field: Which field contains taxonomy (0-indexed)
        random_seed: Random seed for reproducibility

    Returns:
        List of ViralGenome objects

    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If n_genomes > available genomes and replace=False

    Example:
        >>> genomes = sample_genomes_from_fasta(
        ...     'viral_database.fasta',
        ...     n_genomes=50,
        ...     replace=False
        ... )
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    # Set random seed if provided
    if random_seed is not None:
        random.seed(random_seed)
        np.random.seed(random_seed)

    # Read all sequences
    logger.info(f"Reading genomes from {fasta_path}")
    all_records = list(SeqIO.parse(fasta_path, 'fasta'))

    if not replace and n_genomes > len(all_records):
        raise ValueError(
            f"Cannot sample {n_genomes} genomes without replacement "
            f"from {len(all_records)} available genomes"
        )

    # Sample records
    if replace:
        sampled_records = random.choices(all_records, k=n_genomes)
    else:
        sampled_records = random.sample(all_records, k=n_genomes)

    # Convert to ViralGenome objects
    genomes = []
    for record in sampled_records:
        # Parse taxonomy from header
        # Expected format: >genome_id|taxonomy|other_info
        header_parts = record.description.split(taxonomy_separator)

        if len(header_parts) > taxonomy_field:
            taxonomy = header_parts[taxonomy_field].strip()
        else:
            taxonomy = "Unknown"

        # Parse taxonomy levels
        tax_parts = taxonomy.split(';')
        family = tax_parts[0] if len(tax_parts) > 0 else "Unknown"
        genus = tax_parts[1] if len(tax_parts) > 1 else "Unknown"
        species = tax_parts[2] if len(tax_parts) > 2 else "Unknown"

        genome = ViralGenome(
            genome_id=record.id,
            sequence=record.seq,
            taxonomy=taxonomy,
            description=record.description,
            family=family,
            genus=genus,
            species=species
        )
        genomes.append(genome)

    logger.info(f"Sampled {len(genomes)} genomes")
    return genomes


def create_abundance_profile(
    n_genomes: int,
    distribution: str = "lognormal",
    **kwargs
) -> np.ndarray:
    """
    Create an abundance distribution for viral genomes.

    This function generates realistic abundance distributions following
    different statistical models commonly observed in virome studies.

    Args:
        n_genomes: Number of genomes to create abundances for
        distribution: Distribution type:
            - 'lognormal': Log-normal distribution (realistic)
            - 'powerlaw': Power-law distribution (highly uneven)
            - 'even': Even distribution (all equal abundance)
        **kwargs: Distribution-specific parameters:
            For 'lognormal':
                - mu: Mean of log-normal (default=0)
                - sigma: Standard deviation of log-normal (default=1.5)
            For 'powerlaw':
                - alpha: Power-law exponent (default=1.5)
            For 'even':
                - No parameters

    Returns:
        Array of abundances (normalized to sum to 1.0)

    Raises:
        ValueError: If distribution type is not recognized

    Examples:
        >>> # Create log-normal distribution (realistic virome)
        >>> abundances = create_abundance_profile(50, 'lognormal', sigma=2.0)

        >>> # Create power-law distribution (highly uneven)
        >>> abundances = create_abundance_profile(50, 'powerlaw', alpha=1.5)

        >>> # Create even distribution (equal abundances)
        >>> abundances = create_abundance_profile(50, 'even')
    """
    if distribution == 'lognormal':
        # Log-normal distribution - realistic for many viromes
        # Mean and sigma control the shape
        mu = kwargs.get('mu', 0)
        sigma = kwargs.get('sigma', 1.5)

        abundances = np.random.lognormal(mu, sigma, n_genomes)

    elif distribution == 'powerlaw':
        # Power-law distribution - highly uneven, few dominant species
        alpha = kwargs.get('alpha', 1.5)

        # Generate power-law distributed values
        # Using inverse transform sampling
        uniform = np.random.uniform(0, 1, n_genomes)
        abundances = (1 - uniform) ** (-1 / (alpha - 1))

    elif distribution == 'even':
        # Even distribution - all species equally abundant
        abundances = np.ones(n_genomes)

    else:
        raise ValueError(
            f"Unknown distribution: {distribution}. "
            f"Choose from: 'lognormal', 'powerlaw', 'even'"
        )

    # Normalize to sum to 1.0
    abundances = abundances / abundances.sum()

    return abundances


def create_body_site_profile(
    body_site: str,
    n_genomes: int = 50,
    database_path: Optional[Path] = None,
    random_seed: Optional[int] = None
) -> ViralCommunity:
    """
    Create a body-site specific viral community profile.

    This function creates realistic viral community compositions based on
    body site, using literature-derived taxonomic profiles and abundance
    distributions.

    Args:
        body_site: Body site name:
            - 'gut': Gut virome (crAssphage, Microviridae, etc.)
            - 'oral': Oral virome (Streptococcus/Actinomyces phages)
            - 'skin': Skin virome (Propionibacterium/Staphylococcus phages)
            - 'respiratory': Respiratory tract virome
            - 'environmental': Environmental phage communities
        n_genomes: Number of viral genomes to include
        database_path: Optional path to body-site specific database
        random_seed: Random seed for reproducibility

    Returns:
        ViralCommunity object with body-site specific composition

    Raises:
        ValueError: If body_site is not recognized

    Example:
        >>> community = create_body_site_profile(
        ...     'gut',
        ...     n_genomes=100,
        ...     random_seed=42
        ... )
    """
    if random_seed is not None:
        random.seed(random_seed)
        np.random.seed(random_seed)

    # Define body-site specific parameters
    body_site_configs = {
        'gut': {
            'name': 'gut_virome',
            'description': 'Gut virome with crAssphage and typical phages',
            'family_weights': {
                'Siphoviridae': 0.35,
                'Microviridae': 0.25,
                'Myoviridae': 0.20,
                'Podoviridae': 0.15,
                'Other': 0.05
            },
            'abundance_dist': 'lognormal',
            'abundance_params': {'mu': 0, 'sigma': 2.0}
        },
        'oral': {
            'name': 'oral_virome',
            'description': 'Oral virome with Streptococcus and Actinomyces phages',
            'family_weights': {
                'Siphoviridae': 0.40,
                'Myoviridae': 0.30,
                'Podoviridae': 0.20,
                'Other': 0.10
            },
            'abundance_dist': 'lognormal',
            'abundance_params': {'mu': 0, 'sigma': 1.5}
        },
        'skin': {
            'name': 'skin_virome',
            'description': 'Skin virome with Propionibacterium and Staphylococcus phages',
            'family_weights': {
                'Siphoviridae': 0.30,
                'Myoviridae': 0.25,
                'Podoviridae': 0.25,
                'Microviridae': 0.10,
                'Other': 0.10
            },
            'abundance_dist': 'powerlaw',
            'abundance_params': {'alpha': 1.6}
        },
        'respiratory': {
            'name': 'respiratory_virome',
            'description': 'Respiratory tract virome',
            'family_weights': {
                'Siphoviridae': 0.35,
                'Myoviridae': 0.30,
                'Podoviridae': 0.20,
                'Other': 0.15
            },
            'abundance_dist': 'lognormal',
            'abundance_params': {'mu': 0, 'sigma': 1.8}
        },
        'environmental': {
            'name': 'environmental_virome',
            'description': 'Environmental phage community',
            'family_weights': {
                'Siphoviridae': 0.30,
                'Myoviridae': 0.25,
                'Podoviridae': 0.20,
                'Microviridae': 0.15,
                'Other': 0.10
            },
            'abundance_dist': 'powerlaw',
            'abundance_params': {'alpha': 1.5}
        }
    }

    if body_site not in body_site_configs:
        raise ValueError(
            f"Unknown body site: {body_site}. "
            f"Choose from: {', '.join(body_site_configs.keys())}"
        )

    config = body_site_configs[body_site]

    # Create community
    community = ViralCommunity(name=config['name'])

    # If database path provided, sample from it
    if database_path is not None:
        genomes = sample_genomes_from_fasta(
            database_path,
            n_genomes=n_genomes,
            random_seed=random_seed
        )
        community.add_genomes(genomes)
    else:
        # Create synthetic genomes based on family weights
        logger.info(
            f"No database provided, creating synthetic {body_site} community"
        )

        families = list(config['family_weights'].keys())
        weights = list(config['family_weights'].values())

        for i in range(n_genomes):
            # Select family based on weights
            family = np.random.choice(families, p=weights)

            # Create synthetic genome
            # (In real usage, these would be sampled from a database)
            genome_length = int(np.random.lognormal(9.5, 0.8))  # ~30kb mean

            genome = ViralGenome(
                genome_id=f"{body_site}_genome_{i:04d}",
                sequence="N" * genome_length,  # Placeholder sequence
                taxonomy=f"{family};Unknown;Unknown",
                family=family,
                genus="Unknown",
                species=f"{body_site}_sp_{i:04d}",
                description=f"Synthetic {body_site} genome"
            )
            community.add_genome(genome)

    # Apply abundance distribution
    community.apply_abundance_distribution(
        distribution=config['abundance_dist'],
        **config['abundance_params']
    )

    logger.info(
        f"Created {body_site} community with {len(community)} genomes"
    )

    return community


if __name__ == "__main__":
    # Example usage and testing
    print("ViroForge Community Module - Example Usage\n")

    # Example 1: Create a synthetic gut virome community
    print("Example 1: Creating synthetic gut virome community")
    gut_community = create_body_site_profile(
        'gut',
        n_genomes=50,
        random_seed=42
    )
    print(gut_community)
    print(f"\nSummary stats: {gut_community.get_summary_stats()}\n")

    # Example 2: Create custom community with different distributions
    print("Example 2: Creating custom community with log-normal distribution")
    custom_community = ViralCommunity(name="custom_virome")

    # Create some synthetic genomes
    for i in range(20):
        genome = ViralGenome(
            genome_id=f"genome_{i:03d}",
            sequence="ATCG" * 1000,  # 4kb genome
            taxonomy=f"Family_{i%5};Genus_{i%10};Species_{i}",
            family=f"Family_{i%5}",
            genus=f"Genus_{i%10}",
            species=f"Species_{i}"
        )
        custom_community.add_genome(genome)

    # Apply log-normal distribution
    custom_community.apply_abundance_distribution('lognormal', sigma=2.0)
    print(custom_community)
    print(f"\nTop 5 most abundant genomes:")
    abundance_df = custom_community.get_abundance_table()
    print(abundance_df.nlargest(5, 'abundance')[['genome_id', 'abundance']])
