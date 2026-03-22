"""
Low-complexity read injection for sequencing artifact simulation.

ISS generates reads exclusively from reference genomes, which are never
low-complexity. Real sequencing data contains low-complexity artifacts
from adapter dimers, PCR failures, and cluster initialization problems.

This module post-processes generated FASTQ files to replace a configurable
fraction of reads with synthetic low-complexity sequences, enabling
validation of complexity filters (fastp, prinseq, BBDuk).

Artifact types modeled:
- Homopolymer runs (AAAA..., GGGG...) - flow cell signal saturation
- Dinucleotide repeats (ATATAT..., GCGCGC...) - PCR stutter artifacts
- Simple sequence repeats (AATAATAAT...) - polymerase slippage
- Low-entropy random sequence (biased base composition) - failed clusters

Usage:
    from viroforge.simulators.low_complexity import add_low_complexity_reads

    stats = add_low_complexity_reads(
        r1_path=Path("reads_R1.fastq"),
        r2_path=Path("reads_R2.fastq"),
        rate=0.01,
    )
"""

import logging
import math
import random
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# Artifact type definitions with relative frequencies based on
# empirical observations from Illumina sequencing QC reports.
ARTIFACT_TYPES = {
    "homopolymer": {
        "weight": 0.30,
        "motifs": ["A", "T", "G", "C"],
    },
    "dinucleotide_repeat": {
        "weight": 0.25,
        "motifs": ["AT", "TA", "GC", "CG", "AC", "GT", "AG", "CT"],
    },
    "simple_repeat": {
        "weight": 0.20,
        "motifs": ["AAT", "AAC", "AAG", "TTA", "GGA", "CCT", "AAAT", "TTAC"],
    },
    "low_entropy": {
        "weight": 0.25,
        "motifs": None,  # generated procedurally
    },
}


def _generate_homopolymer(length: int, rng: random.Random) -> str:
    """Generate a homopolymer run with minor noise."""
    base = rng.choice(["A", "T", "G", "C"])
    # ~95% dominant base, ~5% noise (real homopolymer artifacts aren't perfect)
    seq = []
    for _ in range(length):
        if rng.random() < 0.95:
            seq.append(base)
        else:
            seq.append(rng.choice([b for b in "ACGT" if b != base]))
    return "".join(seq)


def _generate_repeat(length: int, motif: str, rng: random.Random) -> str:
    """Generate a tandem repeat with occasional errors."""
    repeats = (motif * ((length // len(motif)) + 1))[:length]
    # ~3% error rate (sequencer errors in repetitive regions)
    seq = list(repeats)
    for i in range(length):
        if rng.random() < 0.03:
            seq[i] = rng.choice("ACGT")
    return "".join(seq)


def _generate_low_entropy(length: int, rng: random.Random) -> str:
    """Generate low-entropy sequence with heavily biased base composition."""
    # Pick one or two dominant bases (70-90% of sequence)
    dominant = rng.choice(["A", "T"])
    secondary = rng.choice([b for b in "ACGT" if b != dominant])
    dom_freq = rng.uniform(0.70, 0.90)
    sec_freq = rng.uniform(0.05, 1.0 - dom_freq)
    other_freq = (1.0 - dom_freq - sec_freq) / 2

    others = [b for b in "ACGT" if b not in (dominant, secondary)]
    weights = {
        dominant: dom_freq,
        secondary: sec_freq,
        others[0]: other_freq,
        others[1]: other_freq,
    }

    bases = list(weights.keys())
    probs = list(weights.values())
    return "".join(rng.choices(bases, weights=probs, k=length))


def _shannon_entropy(seq: str) -> float:
    """Compute Shannon entropy of a DNA sequence (bits, max 2.0 for 4 bases)."""
    n = len(seq)
    if n == 0:
        return 0.0
    freqs = [seq.count(b) / n for b in "ACGT"]
    return -sum(f * math.log2(f) for f in freqs if f > 0)


def _generate_target_entropy(
    length: int, target_entropy: float, rng: random.Random
) -> str:
    """Generate a sequence with a specific target Shannon entropy.

    Uses base frequency manipulation to achieve target entropy.
    Shannon entropy H = -sum(p_i * log2(p_i)) for base frequencies p_i.
    By controlling how skewed the base distribution is, we control entropy.

    Args:
        length: Sequence length.
        target_entropy: Target Shannon entropy in bits (0.0 to 2.0).
        rng: Random number generator.

    Returns:
        DNA sequence with entropy close to target.
    """
    target_entropy = max(0.05, min(1.95, target_entropy))

    # Binary search for dominant base frequency that produces target entropy.
    # With one dominant base at frequency p and the rest sharing (1-p)/3:
    #   H = -p*log2(p) - 3*((1-p)/3)*log2((1-p)/3)
    # Higher p = lower entropy.
    lo, hi = 0.25, 0.99  # p range (0.25 = uniform = max entropy)
    best_p = 0.25
    best_diff = float("inf")

    for _ in range(30):
        mid = (lo + hi) / 2
        other = (1.0 - mid) / 3.0
        h = 0.0
        for f in [mid, other, other, other]:
            if f > 0:
                h -= f * math.log2(f)
        diff = abs(h - target_entropy)
        if diff < best_diff:
            best_diff = diff
            best_p = mid
        if h > target_entropy:
            lo = mid  # need more skew (higher p) to lower entropy
        else:
            hi = mid
        if best_diff < 0.005:
            break

    # Generate sequence with the computed base frequencies
    dominant = rng.choice("ACGT")
    others = [b for b in "ACGT" if b != dominant]
    other_freq = (1.0 - best_p) / 3.0
    bases = [dominant] + others
    weights = [best_p, other_freq, other_freq, other_freq]

    return "".join(rng.choices(bases, weights=weights, k=length))


def _generate_artifact_sequence(
    length: int,
    rng: random.Random,
    target_entropy: Optional[float] = None,
) -> tuple[str, str]:
    """Generate a low-complexity artifact sequence.

    Args:
        length: Read length.
        rng: Random number generator.
        target_entropy: If set, generate a sequence with this specific
            Shannon entropy (0.0-2.0) instead of using artifact type weights.

    Returns:
        Tuple of (sequence, artifact_subtype).
    """
    if target_entropy is not None:
        return _generate_target_entropy(length, target_entropy, rng), "controlled_entropy"

    # Select artifact type by weight
    types = list(ARTIFACT_TYPES.keys())
    weights = [ARTIFACT_TYPES[t]["weight"] for t in types]
    artifact_type = rng.choices(types, weights=weights, k=1)[0]

    if artifact_type == "homopolymer":
        return _generate_homopolymer(length, rng), "homopolymer"
    elif artifact_type == "low_entropy":
        return _generate_low_entropy(length, rng), "low_entropy"
    else:
        motifs = ARTIFACT_TYPES[artifact_type]["motifs"]
        motif = rng.choice(motifs)
        return _generate_repeat(length, motif, rng), artifact_type


def add_low_complexity_reads(
    r1_path: Path,
    r2_path: Path,
    output_r1: Optional[Path] = None,
    output_r2: Optional[Path] = None,
    rate: float = 0.01,
    entropy_range: Optional[tuple[float, float]] = None,
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Post-process FASTQ files to inject low-complexity artifact reads.

    Replaces a fraction of read pairs with synthetic low-complexity
    sequences. Both R1 and R2 of selected pairs are replaced (artifact
    reads are never properly paired).

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_r1: Output path for modified R1.
        output_r2: Output path for modified R2.
        rate: Fraction of read pairs to replace (0.0 to 1.0).
        entropy_range: If set, generate reads with Shannon entropy uniformly
            distributed across this range instead of using predefined artifact
            types. Tuple of (min_entropy, max_entropy) in bits (0.0-2.0).
            Useful for testing complexity filter threshold sensitivity.
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place.
        manifest_path: If set, write a TSV manifest of replaced reads.

    Returns:
        Dict with statistics.
    """
    if not (0.0 <= rate <= 1.0):
        raise ValueError(f"rate must be between 0 and 1, got {rate}")

    if rate == 0.0:
        logger.info("Low-complexity rate is 0, skipping")
        return {
            "reads_total": 0,
            "reads_modified": 0,
            "artifact_counts": {},
            "modified_read_ids": [],
        }

    rng = random.Random(random_seed)

    if in_place:
        output_r1 = r1_path
        output_r2 = r2_path
    elif output_r1 is None or output_r2 is None:
        output_r1 = output_r1 or r1_path.parent / f"{r1_path.stem}_lc{r1_path.suffix}"
        output_r2 = output_r2 or r2_path.parent / f"{r2_path.stem}_lc{r2_path.suffix}"

    logger.info(
        f"Injecting {rate:.1%} low-complexity artifact reads into {r1_path.name}"
    )

    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    total_reads = len(r1_records)
    n_to_modify = int(total_reads * rate)
    read_length = len(r1_records[0].seq) if r1_records else 150

    modify_indices = set(rng.sample(range(total_reads), min(n_to_modify, total_reads)))

    modified_r1 = []
    modified_r2 = []
    modified_read_ids = []
    manifest_rows = []
    artifact_counts: dict[str, int] = {}

    for i, (r1, r2) in enumerate(zip(r1_records, r2_records)):
        if i in modify_indices:
            # Generate artifact sequences for both R1 and R2
            if entropy_range is not None:
                target = rng.uniform(entropy_range[0], entropy_range[1])
                r1_artifact_seq, artifact_type = _generate_artifact_sequence(
                    read_length, rng, target_entropy=target
                )
                r2_artifact_seq, _ = _generate_artifact_sequence(
                    read_length, rng, target_entropy=target
                )
            else:
                r1_artifact_seq, artifact_type = _generate_artifact_sequence(
                    read_length, rng
                )
                r2_artifact_seq, _ = _generate_artifact_sequence(read_length, rng)

            # Low quality scores for artifact reads (Q10-Q20 typical)
            artifact_qual = [rng.randint(10, 20) for _ in range(read_length)]

            actual_entropy = _shannon_entropy(r1_artifact_seq)

            source_tag = "source=artifact_low_complexity"
            type_tag = f"artifact_type={artifact_type}"
            entropy_tag = f"entropy={actual_entropy:.3f}"

            new_r1 = SeqRecord(
                Seq(r1_artifact_seq),
                id=r1.id,
                name=r1.name,
                description=f"{r1.description} {source_tag} {type_tag} {entropy_tag}".strip(),
            )
            new_r1.letter_annotations["phred_quality"] = artifact_qual

            new_r2 = SeqRecord(
                Seq(r2_artifact_seq),
                id=r2.id,
                name=r2.name,
                description=f"{r2.description} {source_tag} {type_tag} {entropy_tag}".strip(),
            )
            new_r2.letter_annotations["phred_quality"] = artifact_qual[:]

            modified_r1.append(new_r1)
            modified_r2.append(new_r2)
            modified_read_ids.append(r1.id)
            artifact_counts[artifact_type] = artifact_counts.get(artifact_type, 0) + 1
            manifest_rows.append({
                "read_id": r1.id,
                "artifact_type": artifact_type,
                "entropy": round(actual_entropy, 3),
            })
        else:
            modified_r1.append(r1)
            modified_r2.append(r2)

    SeqIO.write(modified_r1, output_r1, "fastq")
    SeqIO.write(modified_r2, output_r2, "fastq")

    if manifest_path and manifest_rows:
        with open(manifest_path, "w") as f:
            f.write("read_id\tartifact_type\tentropy\n")
            for row in manifest_rows:
                f.write(f"{row['read_id']}\t{row['artifact_type']}\t{row['entropy']}\n")
        logger.info(f"Wrote low-complexity manifest: {manifest_path}")

    logger.info(
        f"Low-complexity injection complete: {len(modified_read_ids)}/{total_reads} "
        f"reads replaced ({', '.join(f'{k}={v}' for k, v in sorted(artifact_counts.items()))})"
    )

    return {
        "reads_total": total_reads,
        "reads_modified": len(modified_read_ids),
        "artifact_counts": artifact_counts,
        "modified_read_ids": modified_read_ids,
    }
