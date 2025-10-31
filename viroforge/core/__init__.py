"""
ViroForge core modules for viral community simulation.

This package contains the core functionality for creating realistic
viral community compositions and simulating virome-specific features.

Modules:
    community: Viral community composition and sampling
    contamination: Contamination profile generation
    enrichment: VLP enrichment modeling (to be implemented)
    artifacts: Sequencing and library prep artifacts (to be implemented)
"""

from .community import (
    ViralGenome,
    ViralCommunity,
    sample_genomes_from_fasta,
    create_abundance_profile,
    create_body_site_profile,
)

from .contamination import (
    ContaminantGenome,
    ContaminantType,
    ContaminationProfile,
    create_contamination_profile,
    add_host_contamination,
    add_rrna_contamination,
    add_reagent_contamination,
    add_phix_control,
)

__all__ = [
    # Community module
    'ViralGenome',
    'ViralCommunity',
    'sample_genomes_from_fasta',
    'create_abundance_profile',
    'create_body_site_profile',
    # Contamination module
    'ContaminantGenome',
    'ContaminantType',
    'ContaminationProfile',
    'create_contamination_profile',
    'add_host_contamination',
    'add_rrna_contamination',
    'add_reagent_contamination',
    'add_phix_control',
]
