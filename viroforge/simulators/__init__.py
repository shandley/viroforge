"""
ViroForge sequencing simulators.

This package contains modules for simulating sequencing data from various
platforms and technologies.

Modules:
    illumina: Illumina sequencing simulation using InSilicoSeq
    pacbio: PacBio long-read simulation (to be implemented)
    nanopore: Oxford Nanopore simulation (to be implemented)
"""

from .illumina import (
    generate_reads,
    quick_generate,
    estimate_file_size,
    check_insilicoseq_installed,
)

__all__ = [
    'generate_reads',
    'quick_generate',
    'estimate_file_size',
    'check_insilicoseq_installed',
]
