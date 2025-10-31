"""
ViroForge utility modules.

This package contains utility functions for working with viral communities,
contamination profiles, genome databases, and creating complete mock virome compositions.

Modules:
    composition: Combine viral communities with contamination profiles
    validation: Quality control and data validation utilities
    genome_database: Download and manage viral genome databases from NCBI
    abundance: Abundance distribution modeling (to be implemented)
    metrics: Calculate ground truth metrics (to be implemented)
"""

from .composition import (
    MockViromeComposition,
    create_mock_virome,
)

from .validation import (
    # Sequence validation
    validate_sequence,
    validate_sequence_length,
    validate_sequence_not_empty,
    # Abundance validation
    validate_abundances,
    validate_abundance_range,
    # FASTQ validation
    validate_fastq_record,
    verify_fastq_file,
    # GC content validation
    validate_gc_content,
    # File validation
    validate_file_not_empty,
    validate_output_directory,
    # Batch validation
    validate_genome_collection,
)

from .genome_database import (
    get_genome_database,
    get_database_info,
    set_entrez_email,
    download_genome_set,
)

__all__ = [
    # Composition
    'MockViromeComposition',
    'create_mock_virome',
    # Validation
    'validate_sequence',
    'validate_sequence_length',
    'validate_sequence_not_empty',
    'validate_abundances',
    'validate_abundance_range',
    'validate_fastq_record',
    'verify_fastq_file',
    'validate_gc_content',
    'validate_file_not_empty',
    'validate_output_directory',
    'validate_genome_collection',
    # Genome Database
    'get_genome_database',
    'get_database_info',
    'set_entrez_email',
    'download_genome_set',
]
