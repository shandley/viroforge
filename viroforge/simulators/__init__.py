"""
ViroForge sequencing simulators.

This package contains modules for simulating sequencing data from various
platforms and technologies.

Modules:
    illumina: Illumina sequencing simulation using InSilicoSeq
    longread: PacBio HiFi and Oxford Nanopore simulation using PBSIM3
"""

from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_paired_fastq(r1_path: Path, r2_path: Path) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Load paired FASTQ files and validate they have matching read counts.

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.

    Returns:
        Tuple of (r1_records, r2_records).

    Raises:
        ValueError: If R1 and R2 have different read counts.
    """
    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    return r1_records, r2_records


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
