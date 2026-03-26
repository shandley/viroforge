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
    mean_insert: Optional[int] = None,
    std_insert: int = 30,
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Post-process FASTQ files to add adapter read-through contamination.

    Simulates adapter contamination by replacing the 3' end of a fraction
    of reads with adapter sequences. Insert sizes are drawn from a normal
    distribution centered on mean_insert, producing realistic adapter
    length patterns where most reads have small adapter tails and few
    have large adapter content.

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
        mean_insert: Mean insert size for normal distribution. Defaults to
            (min_insert + max_insert) / 2. Most reads will have inserts near
            this value, producing small adapter tails.
        std_insert: Standard deviation of insert size distribution (default: 30).
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place (overwrite originals).
        manifest_path: If set, write a TSV manifest of modified reads with
            columns: read_id, adapter_type, adapter_length, insert_size.

    Returns:
        Dict with statistics:
            - reads_total: Total read pairs processed
            - reads_modified: Number of read pairs with adapter injected
            - adapter_type: Adapter type used
            - adapter_lengths: List of adapter contamination lengths
            - mean_adapter_length: Mean adapter bases per modified read
            - modified_read_ids: List of modified read IDs
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
    if mean_insert is None:
        mean_insert = (min_insert + max_insert) // 2

    # Select which read pairs to modify
    modify_indices = set(rng.sample(range(total_reads), min(n_to_modify, total_reads)))

    adapter_lengths = []
    modified_read_ids = []
    manifest_rows = []
    modified_r1 = []
    modified_r2 = []

    for i, (r1, r2) in enumerate(zip(r1_records, r2_records)):
        if i in modify_indices:
            # Normal distribution for insert size (most reads have small adapter tails)
            insert_size = int(rng.gauss(mean_insert, std_insert))
            insert_size = max(min_insert, min(max_insert, insert_size))
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

            # Create new records with adapter tag in description
            adapter_tag = f"adapter_injected={adapter_type}:{adapter_len}bp"
            new_r1 = SeqRecord(
                Seq(new_r1_seq),
                id=r1.id,
                name=r1.name,
                description=f"{r1.description} {adapter_tag}".strip(),
            )
            new_r1.letter_annotations["phred_quality"] = [
                ord(c) - 33 for c in new_r1_qual
            ]

            new_r2 = SeqRecord(
                Seq(new_r2_seq),
                id=r2.id,
                name=r2.name,
                description=f"{r2.description} {adapter_tag}".strip(),
            )
            new_r2.letter_annotations["phred_quality"] = [
                ord(c) - 33 for c in new_r2_qual
            ]

            modified_r1.append(new_r1)
            modified_r2.append(new_r2)
            adapter_lengths.append(adapter_len)
            modified_read_ids.append(r1.id)
            manifest_rows.append({
                "read_id": r1.id,
                "adapter_type": adapter_type,
                "adapter_length": adapter_len,
                "insert_size": insert_size,
            })
        else:
            modified_r1.append(r1)
            modified_r2.append(r2)

    # Write output
    SeqIO.write(modified_r1, output_r1, "fastq")
    SeqIO.write(modified_r2, output_r2, "fastq")

    # Write manifest if requested
    if manifest_path and manifest_rows:
        with open(manifest_path, "w") as f:
            f.write("read_id\tadapter_type\tadapter_length\tinsert_size\n")
            for row in manifest_rows:
                f.write(
                    f"{row['read_id']}\t{row['adapter_type']}\t"
                    f"{row['adapter_length']}\t{row['insert_size']}\n"
                )
        logger.info(f"Wrote adapter manifest: {manifest_path}")

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
        "modified_read_ids": modified_read_ids,
    }


def add_insert_size_adapters(
    r1_path: Path,
    r2_path: Path,
    output_r1: Optional[Path] = None,
    output_r2: Optional[Path] = None,
    adapter_type: str = "truseq",
    mean_insert: int = 350,
    std_insert: int = 50,
    chimera_rate: float = 0.0,
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Insert-size-driven adapter contamination.

    Models adapter read-through as a natural consequence of the insert
    size distribution. Every read pair is assigned an insert size from
    a normal distribution. When insert < read_length, the read naturally
    contains adapter sequence at the 3' end. This produces realistic
    adapter contamination rates that depend on the relationship between
    insert size and read length, matching real library behavior.

    Optionally injects chimeric adapter reads (internal adapter from
    chimeric ligation events, not just 3' read-through).

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_r1: Output path for modified R1.
        output_r2: Output path for modified R2.
        adapter_type: Adapter type ("truseq" or "nextera").
        mean_insert: Mean insert size in bp (default: 350).
        std_insert: Standard deviation of insert size (default: 50).
        chimera_rate: Fraction of reads with internal adapter chimeras
            (default: 0.0). Models chimeric ligation events.
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place.
        manifest_path: If set, write manifest of adapter-affected reads.

    Returns:
        Dict with statistics including emergent adapter rate.
    """
    if adapter_type not in ADAPTER_SEQUENCES:
        raise ValueError(
            f"Unknown adapter type: {adapter_type}. "
            f"Choose from: {', '.join(ADAPTER_SEQUENCES.keys())}"
        )

    adapters = ADAPTER_SEQUENCES[adapter_type]
    rng = random.Random(random_seed)

    if in_place:
        output_r1 = r1_path
        output_r2 = r2_path
    elif output_r1 is None or output_r2 is None:
        output_r1 = output_r1 or r1_path.parent / f"{r1_path.stem}_adapters{r1_path.suffix}"
        output_r2 = output_r2 or r2_path.parent / f"{r2_path.stem}_adapters{r2_path.suffix}"

    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    total_reads = len(r1_records)
    read_length = len(r1_records[0].seq) if r1_records else 150

    logger.info(
        f"Insert-size-driven adapter contamination: mean={mean_insert}bp, "
        f"sd={std_insert}bp, read_length={read_length}bp ({adapter_type})"
    )

    readthrough_lengths = []
    chimera_count = 0
    modified_read_ids = []
    manifest_rows = []
    modified_r1 = []
    modified_r2 = []
    insert_sizes = []

    for i, (r1, r2) in enumerate(zip(r1_records, r2_records)):
        # Sample insert size for this fragment
        insert_size = max(30, int(rng.gauss(mean_insert, std_insert)))
        insert_sizes.append(insert_size)

        modified = False

        # Read-through: when insert < read_length, adapter appears at 3' end
        if insert_size < read_length:
            adapter_len = read_length - insert_size

            r1_seq = str(r1.seq)
            r1_qual_str = "".join(chr(q + 33) for q in r1.letter_annotations["phred_quality"])
            r2_seq = str(r2.seq)
            r2_qual_str = "".join(chr(q + 33) for q in r2.letter_annotations["phred_quality"])

            new_r1_seq, new_r1_qual = _inject_adapter(
                r1_seq, r1_qual_str, adapters["r1_adapter"], insert_size, rng
            )
            new_r2_seq, new_r2_qual = _inject_adapter(
                r2_seq, r2_qual_str, adapters["r2_adapter"], insert_size, rng
            )

            tag = f"adapter_readthrough={adapter_type}:{adapter_len}bp insert_size={insert_size}"
            r1 = SeqRecord(
                Seq(new_r1_seq), id=r1.id, name=r1.name,
                description=f"{r1.description} {tag}".strip(),
            )
            r1.letter_annotations["phred_quality"] = [ord(c) - 33 for c in new_r1_qual]
            r2 = SeqRecord(
                Seq(new_r2_seq), id=r2.id, name=r2.name,
                description=f"{r2.description} {tag}".strip(),
            )
            r2.letter_annotations["phred_quality"] = [ord(c) - 33 for c in new_r2_qual]

            readthrough_lengths.append(adapter_len)
            modified = True
            manifest_rows.append({
                "read_id": r1.id, "type": "readthrough",
                "adapter_type": adapter_type, "adapter_length": adapter_len,
                "insert_size": insert_size,
            })

        # Chimeric adapter: internal adapter sequence within the read
        elif chimera_rate > 0 and rng.random() < chimera_rate:
            # Insert adapter at a random internal position
            pos = rng.randint(20, read_length - 20)
            adapter_seq = adapters["r1_adapter"][:rng.randint(10, 30)]

            r1_seq = str(r1.seq)
            chimeric_seq = r1_seq[:pos] + adapter_seq + r1_seq[pos + len(adapter_seq):]
            chimeric_seq = chimeric_seq[:read_length]

            r1_qual = r1.letter_annotations["phred_quality"][:]
            tag = f"adapter_chimera={adapter_type}:{len(adapter_seq)}bp pos={pos}"
            r1 = SeqRecord(
                Seq(chimeric_seq), id=r1.id, name=r1.name,
                description=f"{r1.description} {tag}".strip(),
            )
            r1.letter_annotations["phred_quality"] = r1_qual

            chimera_count += 1
            modified = True
            manifest_rows.append({
                "read_id": r1.id, "type": "chimera",
                "adapter_type": adapter_type, "adapter_length": len(adapter_seq),
                "insert_size": insert_size,
            })

        if modified:
            modified_read_ids.append(r1.id)

        modified_r1.append(r1)
        modified_r2.append(r2)

    SeqIO.write(modified_r1, output_r1, "fastq")
    SeqIO.write(modified_r2, output_r2, "fastq")

    if manifest_path and manifest_rows:
        with open(manifest_path, "w") as f:
            f.write("read_id\ttype\tadapter_type\tadapter_length\tinsert_size\n")
            for row in manifest_rows:
                f.write(
                    f"{row['read_id']}\t{row['type']}\t{row['adapter_type']}\t"
                    f"{row['adapter_length']}\t{row['insert_size']}\n"
                )
        logger.info(f"Wrote adapter manifest: {manifest_path}")

    n_readthrough = len(readthrough_lengths)
    emergent_rate = n_readthrough / total_reads if total_reads > 0 else 0.0
    mean_adapter_len = (
        sum(readthrough_lengths) / len(readthrough_lengths)
        if readthrough_lengths else 0.0
    )

    logger.info(
        f"Insert-size adapter injection complete: "
        f"{n_readthrough} read-through ({emergent_rate:.1%}), "
        f"{chimera_count} chimeras, "
        f"mean adapter length: {mean_adapter_len:.1f}bp"
    )

    return {
        "reads_total": total_reads,
        "reads_readthrough": n_readthrough,
        "reads_chimera": chimera_count,
        "reads_modified": n_readthrough + chimera_count,
        "emergent_adapter_rate": emergent_rate,
        "adapter_type": adapter_type,
        "adapter_lengths": readthrough_lengths,
        "mean_adapter_length": mean_adapter_len,
        "mean_insert_size": sum(insert_sizes) / len(insert_sizes) if insert_sizes else 0,
        "insert_size_sd": std_insert,
        "modified_read_ids": modified_read_ids,
    }
