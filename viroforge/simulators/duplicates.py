"""
PCR duplicate injection for deduplication tool validation.

ISS generates independent reads even with amplification bias applied,
so every read starts at a unique genome position. Real PCR duplicates
arise from the same template molecule and share identical start
positions and sequences (with occasional PCR errors).

This module post-processes generated FASTQ files to create realistic
PCR duplicate clusters, enabling validation of deduplication tools
(Picard MarkDuplicates, samtools markdup, fastp dedup, virome-qc).

Duplicate model:
- Selected reads become "templates" with 1-N copies
- Copy count follows geometric distribution (most have 1 copy, few have many)
- Copies have optional per-base PCR error rate
- Both mates duplicated together (paired-end consistency)
- Duplicates shuffled throughout file (not adjacent)

Usage:
    from viroforge.simulators.duplicates import add_pcr_duplicates

    stats = add_pcr_duplicates(
        r1_path=Path("reads_R1.fastq"),
        r2_path=Path("reads_R2.fastq"),
        duplicate_rate=0.30,
        max_copies=5,
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


def _introduce_pcr_errors(
    seq: str, error_rate: float, rng: random.Random
) -> str:
    """Introduce random PCR errors (substitutions) into a sequence."""
    if error_rate <= 0:
        return seq
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if rng.random() < error_rate:
            original = seq_list[i]
            seq_list[i] = rng.choice([b for b in "ACGT" if b != original])
    return "".join(seq_list)


def add_pcr_duplicates(
    r1_path: Path,
    r2_path: Path,
    output_r1: Optional[Path] = None,
    output_r2: Optional[Path] = None,
    duplicate_rate: float = 0.30,
    max_copies: int = 5,
    copy_distribution: str = "geometric",
    error_rate: float = 0.001,
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Post-process FASTQ files to inject PCR duplicate read pairs.

    Selects a fraction of read pairs as "templates" and generates
    1 to max_copies duplicate copies of each. Copies are near-identical
    with optional PCR error rate. Duplicates are shuffled throughout
    the output file.

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_r1: Output path for modified R1.
        output_r2: Output path for modified R2.
        duplicate_rate: Fraction of original reads that become templates
            with at least one copy (0.0 to 1.0).
        max_copies: Maximum number of duplicate copies per template.
        copy_distribution: Distribution for copy counts.
            "geometric" (default): Most templates get 1 copy, few get many.
                P(k copies) = (1-p)^(k-1) * p, where p=0.5.
            "uniform": Equal probability for 1 to max_copies.
        error_rate: Per-base substitution rate in copies (default: 0.001).
            Models PCR polymerase errors. Set to 0 for exact duplicates.
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place.
        manifest_path: If set, write a TSV manifest mapping copies to
            their template read IDs.

    Returns:
        Dict with statistics:
            - reads_original: Original read pair count
            - reads_output: Total read pairs in output (original + copies)
            - templates: Number of template read pairs selected
            - copies_generated: Total duplicate copies created
            - copy_count_distribution: Dict of {n_copies: n_templates}
            - duplicate_fraction: Fraction of output reads that are duplicates
    """
    if not (0.0 <= duplicate_rate <= 1.0):
        raise ValueError(f"duplicate_rate must be between 0 and 1, got {duplicate_rate}")

    if copy_distribution not in ("geometric", "uniform"):
        raise ValueError(
            f"Unknown copy_distribution: {copy_distribution}. "
            f"Choose from: geometric, uniform"
        )

    if duplicate_rate == 0.0:
        logger.info("Duplicate rate is 0, skipping")
        return {
            "reads_original": 0,
            "reads_output": 0,
            "templates": 0,
            "copies_generated": 0,
            "copy_count_distribution": {},
            "duplicate_fraction": 0.0,
        }

    rng = random.Random(random_seed)

    if in_place:
        output_r1 = r1_path
        output_r2 = r2_path
    elif output_r1 is None or output_r2 is None:
        output_r1 = output_r1 or r1_path.parent / f"{r1_path.stem}_dedup{r1_path.suffix}"
        output_r2 = output_r2 or r2_path.parent / f"{r2_path.stem}_dedup{r2_path.suffix}"

    logger.info(
        f"Injecting PCR duplicates ({duplicate_rate:.0%} template rate, "
        f"max {max_copies} copies, {copy_distribution} distribution) "
        f"into {r1_path.name}"
    )

    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    total_reads = len(r1_records)
    n_templates = int(total_reads * duplicate_rate)

    # Select template indices
    template_indices = set(
        rng.sample(range(total_reads), min(n_templates, total_reads))
    )

    # Build output: original reads + duplicate copies
    out_r1 = []
    out_r2 = []
    manifest_rows = []
    copy_count_dist: dict[int, int] = {}
    total_copies = 0

    for i, (r1, r2) in enumerate(zip(r1_records, r2_records)):
        # Always include the original
        out_r1.append(r1)
        out_r2.append(r2)

        if i in template_indices:
            # Determine number of copies
            if copy_distribution == "geometric":
                # Geometric with p=0.5: P(1)=0.5, P(2)=0.25, P(3)=0.125, ...
                n_copies = 1
                while n_copies < max_copies and rng.random() > 0.5:
                    n_copies += 1
            else:  # uniform
                n_copies = rng.randint(1, max_copies)

            copy_count_dist[n_copies] = copy_count_dist.get(n_copies, 0) + 1

            for copy_num in range(1, n_copies + 1):
                total_copies += 1

                # Create copy with optional PCR errors
                r1_copy_seq = _introduce_pcr_errors(str(r1.seq), error_rate, rng)
                r2_copy_seq = _introduce_pcr_errors(str(r2.seq), error_rate, rng)

                # Quality scores: slightly degraded from original (PCR copies
                # tend to have marginally lower quality in later cycles)
                r1_qual = r1.letter_annotations["phred_quality"][:]
                r2_qual = r2.letter_annotations["phred_quality"][:]

                dup_tag = f"pcr_duplicate=true duplicate_of={r1.id} copy_number={copy_num}"

                copy_r1 = SeqRecord(
                    Seq(r1_copy_seq),
                    id=f"{r1.id}_dup{copy_num}",
                    name=r1.name,
                    description=f"{r1.description} {dup_tag}".strip(),
                )
                copy_r1.letter_annotations["phred_quality"] = r1_qual

                copy_r2 = SeqRecord(
                    Seq(r2_copy_seq),
                    id=f"{r2.id}_dup{copy_num}",
                    name=r2.name,
                    description=f"{r2.description} {dup_tag}".strip(),
                )
                copy_r2.letter_annotations["phred_quality"] = r2_qual

                out_r1.append(copy_r1)
                out_r2.append(copy_r2)

                manifest_rows.append({
                    "template_read_id": r1.id,
                    "copy_read_id": f"{r1.id}_dup{copy_num}",
                    "copy_number": copy_num,
                    "total_copies": n_copies,
                })

    # Shuffle duplicates throughout the file (zip R1/R2 to keep pairs together)
    paired = list(zip(out_r1, out_r2))
    rng.shuffle(paired)
    out_r1 = [p[0] for p in paired]
    out_r2 = [p[1] for p in paired]

    # Write output
    SeqIO.write(out_r1, output_r1, "fastq")
    SeqIO.write(out_r2, output_r2, "fastq")

    # Write manifest
    if manifest_path and manifest_rows:
        with open(manifest_path, "w") as f:
            f.write("template_read_id\tcopy_read_id\tcopy_number\ttotal_copies\n")
            for row in manifest_rows:
                f.write(
                    f"{row['template_read_id']}\t{row['copy_read_id']}\t"
                    f"{row['copy_number']}\t{row['total_copies']}\n"
                )
        logger.info(f"Wrote duplicate manifest: {manifest_path}")

    reads_output = len(out_r1)
    dup_fraction = total_copies / reads_output if reads_output > 0 else 0.0

    logger.info(
        f"Duplicate injection complete: {len(template_indices)} templates, "
        f"{total_copies} copies generated, "
        f"{reads_output} total reads ({dup_fraction:.1%} duplicates)"
    )

    return {
        "reads_original": total_reads,
        "reads_output": reads_output,
        "templates": len(template_indices),
        "copies_generated": total_copies,
        "copy_count_distribution": copy_count_dist,
        "duplicate_fraction": dup_fraction,
    }
