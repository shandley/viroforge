"""
Adapter read-through contamination simulator.

Models the situation where DNA inserts are shorter than the read length,
causing the sequencer to read into the adapter sequence on the opposite
end. This is a common artifact in Illumina sequencing that QC tools
(fastp, Trimmomatic, BBDuk) are designed to detect and trim.

InSilicoSeq does not model adapter read-through, so this module
post-processes generated FASTQ files to inject realistic adapter
contamination at controlled rates.

Usage:
    from viroforge.simulators.adapters import add_adapter_readthrough

    stats = add_adapter_readthrough(
        r1_path=Path("reads_R1.fastq"),
        r2_path=Path("reads_R2.fastq"),
        output_r1=Path("reads_adapters_R1.fastq"),
        output_r2=Path("reads_adapters_R2.fastq"),
        adapter_rate=0.05,
        adapter_type="truseq",
    )
"""

import logging
import random
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# Standard Illumina adapter sequences that appear in read-through.
# When insert < read_length, the sequencer reads past the insert and into
# the adapter ligated to the opposite end of the fragment.
#
# For R1: reads into the reverse complement of the Read 2 adapter
# For R2: reads into the reverse complement of the Read 1 adapter
ADAPTER_SEQUENCES = {
    "truseq": {
        # TruSeq: R1 reads into this adapter at 3' end
        "r1_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        # TruSeq: R2 reads into this adapter at 3' end
        "r2_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    },
    "nextera": {
        "r1_adapter": "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        "r2_adapter": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    },
}


def _inject_adapter(
    seq: str,
    qual: str,
    adapter: str,
    insert_size: int,
    rng: random.Random,
) -> tuple[str, str]:
    """Replace the 3' end of a read with adapter sequence.

    Args:
        seq: Original read sequence.
        qual: Original quality string.
        adapter: Adapter sequence to inject.
        insert_size: Size of the actual DNA insert. Adapter starts at
            this position in the read.
        rng: Random number generator.

    Returns:
        Tuple of (modified_sequence, modified_quality).
    """
    read_len = len(seq)
    if insert_size >= read_len:
        return seq, qual

    # Biological sequence up to insert boundary
    bio_seq = seq[:insert_size]

    # Adapter portion that fills the rest of the read
    adapter_len_needed = read_len - insert_size
    # Extend adapter if needed (repeat it)
    extended_adapter = (adapter * ((adapter_len_needed // len(adapter)) + 1))
    adapter_portion = extended_adapter[:adapter_len_needed]

    new_seq = bio_seq + adapter_portion

    # Quality scores: adapter bases typically sequence at high quality
    # (Q30-37) since they're synthetic and amplify well
    bio_qual = qual[:insert_size]
    adapter_qual = "".join(
        chr(rng.randint(63, 70))  # Phred+33: Q30-Q37
        for _ in range(adapter_len_needed)
    )
    new_qual = bio_qual + adapter_qual

    return new_seq, new_qual


def add_adapter_readthrough(
    r1_path: Path,
    r2_path: Path,
    output_r1: Optional[Path] = None,
    output_r2: Optional[Path] = None,
    adapter_rate: float = 0.05,
    adapter_type: str = "truseq",
    min_insert: int = 50,
    max_insert: Optional[int] = None,
    random_seed: Optional[int] = None,
    in_place: bool = False,
) -> dict:
    """Post-process FASTQ files to add adapter read-through contamination.

    Simulates adapter contamination by replacing the 3' end of a fraction
    of reads with adapter sequences. The amount of adapter in each
    affected read varies based on a randomly sampled insert size.

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_r1: Output path for modified R1 (default: overwrite input).
        output_r2: Output path for modified R2 (default: overwrite input).
        adapter_rate: Fraction of read pairs to contaminate (0.0 to 1.0).
        adapter_type: Adapter type ("truseq" or "nextera").
        min_insert: Minimum insert size for affected reads (bp).
        max_insert: Maximum insert size for affected reads. Defaults to
            read_length - 1 (at least 1 bp of adapter).
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place (overwrite originals).

    Returns:
        Dict with statistics:
            - reads_total: Total read pairs processed
            - reads_modified: Number of read pairs with adapter injected
            - adapter_type: Adapter type used
            - adapter_lengths: List of adapter contamination lengths
            - mean_adapter_length: Mean adapter bases per modified read
    """
    if adapter_type not in ADAPTER_SEQUENCES:
        raise ValueError(
            f"Unknown adapter type: {adapter_type}. "
            f"Choose from: {', '.join(ADAPTER_SEQUENCES.keys())}"
        )

    if not (0.0 <= adapter_rate <= 1.0):
        raise ValueError(f"adapter_rate must be between 0 and 1, got {adapter_rate}")

    if adapter_rate == 0.0:
        logger.info("Adapter rate is 0, skipping adapter injection")
        return {
            "reads_total": 0,
            "reads_modified": 0,
            "adapter_type": adapter_type,
            "adapter_lengths": [],
            "mean_adapter_length": 0.0,
        }

    adapters = ADAPTER_SEQUENCES[adapter_type]
    rng = random.Random(random_seed)

    # Determine output paths
    if in_place:
        output_r1 = r1_path
        output_r2 = r2_path
    elif output_r1 is None or output_r2 is None:
        output_r1 = output_r1 or r1_path.parent / f"{r1_path.stem}_adapters{r1_path.suffix}"
        output_r2 = output_r2 or r2_path.parent / f"{r2_path.stem}_adapters{r2_path.suffix}"

    logger.info(
        f"Adding {adapter_rate:.1%} adapter read-through "
        f"({adapter_type}) to {r1_path.name}"
    )

    # Read all records
    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    total_reads = len(r1_records)
    n_to_modify = int(total_reads * adapter_rate)

    # Determine read length from first record
    read_length = len(r1_records[0].seq) if r1_records else 150
    if max_insert is None:
        max_insert = read_length - 1

    # Select which read pairs to modify
    modify_indices = set(rng.sample(range(total_reads), min(n_to_modify, total_reads)))

    adapter_lengths = []
    modified_r1 = []
    modified_r2 = []

    for i, (r1, r2) in enumerate(zip(r1_records, r2_records)):
        if i in modify_indices:
            # Random insert size determines how much adapter appears
            insert_size = rng.randint(min_insert, max_insert)
            adapter_len = read_length - insert_size

            r1_seq = str(r1.seq)
            r1_qual = r1.letter_annotations["phred_quality"]
            r1_qual_str = "".join(chr(q + 33) for q in r1_qual)

            r2_seq = str(r2.seq)
            r2_qual = r2.letter_annotations["phred_quality"]
            r2_qual_str = "".join(chr(q + 33) for q in r2_qual)

            new_r1_seq, new_r1_qual = _inject_adapter(
                r1_seq, r1_qual_str, adapters["r1_adapter"], insert_size, rng
            )
            new_r2_seq, new_r2_qual = _inject_adapter(
                r2_seq, r2_qual_str, adapters["r2_adapter"], insert_size, rng
            )

            # Create new records
            new_r1 = SeqRecord(
                Seq(new_r1_seq),
                id=r1.id,
                name=r1.name,
                description=r1.description,
            )
            new_r1.letter_annotations["phred_quality"] = [
                ord(c) - 33 for c in new_r1_qual
            ]

            new_r2 = SeqRecord(
                Seq(new_r2_seq),
                id=r2.id,
                name=r2.name,
                description=r2.description,
            )
            new_r2.letter_annotations["phred_quality"] = [
                ord(c) - 33 for c in new_r2_qual
            ]

            modified_r1.append(new_r1)
            modified_r2.append(new_r2)
            adapter_lengths.append(adapter_len)
        else:
            modified_r1.append(r1)
            modified_r2.append(r2)

    # Write output
    SeqIO.write(modified_r1, output_r1, "fastq")
    SeqIO.write(modified_r2, output_r2, "fastq")

    mean_adapter_len = (
        sum(adapter_lengths) / len(adapter_lengths) if adapter_lengths else 0.0
    )

    logger.info(
        f"Adapter injection complete: {len(adapter_lengths)}/{total_reads} reads modified "
        f"(mean adapter length: {mean_adapter_len:.1f} bp)"
    )

    return {
        "reads_total": total_reads,
        "reads_modified": len(adapter_lengths),
        "adapter_type": adapter_type,
        "adapter_lengths": adapter_lengths,
        "mean_adapter_length": mean_adapter_len,
    }
