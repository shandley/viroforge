"""
ViroForge sequencing simulators.

This package contains modules for simulating sequencing data from various
platforms and technologies.

Modules:
    illumina: Illumina sequencing simulation using InSilicoSeq
    longread: PacBio HiFi and Oxford Nanopore simulation using PBSIM3
"""

from .illumina import (
    generate_reads,
    quick_generate,
    estimate_file_size,
    check_insilicoseq_installed,
)

from .longread import (
    generate_long_reads,
    LongReadPlatform,
    PacBioHiFiConfig,
    NanoporeConfig,
    check_pbsim3_installed,
    check_pbccs_installed,
)

__all__ = [
    # Short-read (Illumina)
    'generate_reads',
    'quick_generate',
    'estimate_file_size',
    'check_insilicoseq_installed',

    # Long-read (PacBio HiFi, Nanopore)
    'generate_long_reads',
    'LongReadPlatform',
    'PacBioHiFiConfig',
    'NanoporeConfig',
    'check_pbsim3_installed',
    'check_pbccs_installed',
]
