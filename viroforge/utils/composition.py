"""
Utilities for combining viral communities with contamination profiles.

This module provides functions to create complete mock virome compositions
by combining viral communities with contamination profiles and managing
the overall abundance distributions.
"""

from typing import Dict, List, Optional, Tuple
import logging

import numpy as np
import pandas as pd

from ..core.community import ViralCommunity, ViralGenome
from ..core.contamination import ContaminationProfile, ContaminantGenome

logger = logging.getLogger(__name__)


class MockViromeComposition:
    """
    Represents a complete mock virome composition.

    This class combines a viral community with a contamination profile,
    managing the relative abundances to create a realistic mock dataset.

    Attributes:
        name: Name for this composition
        viral_community: ViralCommunity object
        contamination_profile: ContaminationProfile object
        viral_fraction: Fraction of reads that are viral (0-1)
        contamination_fraction: Fraction of reads that are contamination (0-1)
    """

    def __init__(
        self,
        name: str,
        viral_community: ViralCommunity,
        contamination_profile: Optional[ContaminationProfile] = None,
        viral_fraction: float = 0.90
    ):
        """
        Initialize a mock virome composition.

        Args:
            name: Name for this composition
            viral_community: ViralCommunity to use
            contamination_profile: Optional ContaminationProfile
            viral_fraction: Fraction of total that should be viral (0-1)
        """
        self.name = name
        self.viral_community = viral_community
        self.contamination_profile = contamination_profile
        self.viral_fraction = viral_fraction
        self.contamination_fraction = 1.0 - viral_fraction

        # Normalize abundances
        self._normalize_abundances()

    def _normalize_abundances(self) -> None:
        """
        Normalize viral and contamination abundances to match target fractions.
        """
        # Scale viral abundances
        total_viral = self.viral_community.get_total_abundance()
        if total_viral > 0:
            for genome in self.viral_community.genomes:
                genome.abundance = (genome.abundance / total_viral) * self.viral_fraction

        # Scale contamination abundances if present
        if self.contamination_profile is not None:
            total_contam = self.contamination_profile.get_total_abundance()
            if total_contam > 0:
                for contaminant in self.contamination_profile.contaminants:
                    contaminant.abundance = (
                        contaminant.abundance / total_contam
                    ) * self.contamination_fraction

        logger.info(
            f"Normalized abundances: {self.viral_fraction*100:.1f}% viral, "
            f"{self.contamination_fraction*100:.1f}% contamination"
        )

    def get_total_abundance(self) -> float:
        """
        Get total abundance across all components.

        Returns:
            Total abundance (should be ~1.0 after normalization)
        """
        total = self.viral_community.get_total_abundance()

        if self.contamination_profile is not None:
            total += self.contamination_profile.get_total_abundance()

        return total

    def get_all_sequences(self) -> List[Tuple[str, object, float, str]]:
        """
        Get all sequences (viral + contaminants) with metadata.

        Returns:
            List of tuples (genome_id, sequence_obj, abundance, source_type)
        """
        sequences = []

        # Add viral genomes
        for genome in self.viral_community.genomes:
            sequences.append((
                genome.genome_id,
                genome,
                genome.abundance,
                'viral'
            ))

        # Add contaminants
        if self.contamination_profile is not None:
            for contaminant in self.contamination_profile.contaminants:
                sequences.append((
                    contaminant.genome_id,
                    contaminant,
                    contaminant.abundance,
                    contaminant.contaminant_type.value
                ))

        return sequences

    def get_composition_table(self) -> pd.DataFrame:
        """
        Get complete composition table with all sequences.

        Returns:
            DataFrame with all genomes and contaminants
        """
        data = []

        # Add viral genomes
        for genome in self.viral_community.genomes:
            data.append({
                'genome_id': genome.genome_id,
                'type': 'viral',
                'taxonomy': genome.taxonomy,
                'organism': f"{genome.family}; {genome.genus}",
                'abundance': genome.abundance,
                'length': genome.length,
                'gc_content': genome.gc_content
            })

        # Add contaminants
        if self.contamination_profile is not None:
            for contaminant in self.contamination_profile.contaminants:
                data.append({
                    'genome_id': contaminant.genome_id,
                    'type': contaminant.contaminant_type.value,
                    'taxonomy': 'N/A',
                    'organism': contaminant.organism,
                    'abundance': contaminant.abundance,
                    'length': contaminant.length,
                    'gc_content': contaminant.gc_content
                })

        return pd.DataFrame(data)

    def get_summary_stats(self) -> Dict:
        """
        Get summary statistics for the composition.

        Returns:
            Dictionary with composition statistics
        """
        viral_stats = self.viral_community.get_summary_stats()

        stats = {
            'name': self.name,
            'n_viral_genomes': viral_stats['n_genomes'],
            'viral_abundance': self.viral_community.get_total_abundance(),
            'viral_fraction': self.viral_fraction,
        }

        if self.contamination_profile is not None:
            contam_stats = self.contamination_profile.get_summary_stats()
            stats.update({
                'n_contaminants': contam_stats['n_contaminants'],
                'contamination_abundance': self.contamination_profile.get_total_abundance(),
                'contamination_fraction': self.contamination_fraction,
                'contamination_by_type': contam_stats['by_type']
            })
        else:
            stats.update({
                'n_contaminants': 0,
                'contamination_abundance': 0.0,
                'contamination_fraction': 0.0,
                'contamination_by_type': {}
            })

        stats['total_abundance'] = self.get_total_abundance()

        return stats

    def export_composition_table(self, output_path) -> None:
        """
        Export complete composition table to TSV.

        Args:
            output_path: Path to output file
        """
        df = self.get_composition_table()
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Exported composition table to {output_path}")

    def __repr__(self) -> str:
        """String representation."""
        stats = self.get_summary_stats()
        return (
            f"MockViromeComposition(name='{self.name}', "
            f"n_viral={stats['n_viral_genomes']}, "
            f"n_contaminants={stats['n_contaminants']}, "
            f"viral_frac={stats['viral_fraction']:.2f})"
        )


def create_mock_virome(
    name: str,
    body_site: str = "gut",
    contamination_level: str = "realistic",
    n_viral_genomes: int = 50,
    viral_fraction: float = 0.90,
    random_seed: Optional[int] = None,
    **kwargs
) -> MockViromeComposition:
    """
    Create a complete mock virome composition with viral genomes and contamination.

    This is a convenience function that creates both a viral community and
    contamination profile, then combines them into a MockViromeComposition.

    Args:
        name: Name for the mock virome
        body_site: Body site for viral community ('gut', 'oral', 'skin', etc.)
        contamination_level: Level of contamination ('clean', 'realistic', 'heavy', 'failed')
        n_viral_genomes: Number of viral genomes to include
        viral_fraction: Fraction of reads that should be viral (0-1)
        random_seed: Random seed for reproducibility
        **kwargs: Additional parameters for viral community or contamination

    Returns:
        MockViromeComposition object

    Example:
        >>> mock_virome = create_mock_virome(
        ...     'my_mock_dataset',
        ...     body_site='gut',
        ...     contamination_level='realistic',
        ...     n_viral_genomes=100,
        ...     viral_fraction=0.85
        ... )
    """
    from ..core.community import create_body_site_profile
    from ..core.contamination import create_contamination_profile

    logger.info(
        f"Creating mock virome '{name}': {body_site} site, "
        f"{contamination_level} contamination, {n_viral_genomes} viral genomes"
    )

    # Create viral community
    viral_community = create_body_site_profile(
        body_site=body_site,
        n_genomes=n_viral_genomes,
        random_seed=random_seed
    )

    # Create contamination profile
    contamination_profile = create_contamination_profile(
        profile_type=contamination_level,
        random_seed=random_seed,
        **kwargs
    )

    # Combine into composition
    composition = MockViromeComposition(
        name=name,
        viral_community=viral_community,
        contamination_profile=contamination_profile,
        viral_fraction=viral_fraction
    )

    logger.info(f"Created mock virome: {composition}")

    return composition


if __name__ == "__main__":
    # Example usage
    print("ViroForge Composition Utilities - Example Usage\n")

    # Example 1: Create a complete mock virome
    print("Example 1: Creating a complete mock virome")
    print("=" * 60)

    mock_virome = create_mock_virome(
        name="test_gut_virome",
        body_site="gut",
        contamination_level="realistic",
        n_viral_genomes=30,
        viral_fraction=0.85,
        random_seed=42
    )

    print(f"\n{mock_virome}")
    print("\nSummary Statistics:")
    stats = mock_virome.get_summary_stats()
    for key, value in stats.items():
        if key != 'contamination_by_type':
            print(f"  {key}: {value}")

    print("\nContamination breakdown:")
    for ctype, abundance in stats['contamination_by_type'].items():
        print(f"  {ctype}: {abundance*100:.2f}%")

    # Example 2: Get composition table
    print("\n\nExample 2: Composition table")
    print("=" * 60)
    comp_table = mock_virome.get_composition_table()
    print(f"\nTotal sequences: {len(comp_table)}")
    print(f"\nBreakdown by type:")
    print(comp_table.groupby('type')['abundance'].agg(['count', 'sum']))
