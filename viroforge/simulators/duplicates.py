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


def _compute_gc(seq: str) -> float:
    """Compute GC fraction of a sequence."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if len(seq) > 0 else 0.5


def _extract_genome_and_index(record_id: str):
    """Extract genome ID and read index from ISS read header.

    ISS format: {genome_id}_{read_index}_{pair}/{read}
    e.g., GCF_000882635.1_9203_0/1 -> ("GCF_000882635.1", 9203)
    """
    import re
    m = re.match(r'((?:GCF_\d+\.\d+|[a-zA-Z_]+\d+))_(\d+)_\d+/\d+', record_id)
    if m:
        return m.group(1), int(m.group(2))
    # Fallback: try splitting
    parts = record_id.rsplit('_', 2)
    if len(parts) >= 3:
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return record_id, 0


def _compute_mda_hotspot_weights(
    records: list[SeqRecord], rng: random.Random,
    n_hotspots_per_genome: int = 8,
    hotspot_width: float = 0.05,
) -> list[float]:
    """Compute MDA priming hotspot weights for within-genome clustering.

    In real MDA, phi29 uses random hexamer primers. Some genomic regions
    have more favorable priming sites, creating amplification hotspots
    where coverage is 10-100x higher than cold spots.

    This function:
    1. Groups reads by genome
    2. Places random priming hotspots along each genome
    3. Weights reads by proximity to nearest hotspot
    4. Combines with GC bias and stochastic noise

    Args:
        records: List of SeqRecord objects.
        rng: Random number generator.
        n_hotspots_per_genome: Number of priming hotspots per genome
            (default 8). More hotspots = more even coverage.
        hotspot_width: Width of hotspot influence as fraction of genome
            index range (default 0.05 = 5%). Narrower = sharper peaks.

    Returns:
        List of weights (higher = more likely to be a template).
    """
    # Group reads by genome and extract positional indices
    genome_reads = {}  # genome_id -> [(list_index, read_index)]
    for i, record in enumerate(records):
        genome_id, read_idx = _extract_genome_and_index(record.id)
        if genome_id not in genome_reads:
            genome_reads[genome_id] = []
        genome_reads[genome_id].append((i, read_idx))

    weights = [0.0] * len(records)

    for genome_id, read_list in genome_reads.items():
        if not read_list:
            continue

        # Get the range of read indices for this genome
        indices = [idx for _, idx in read_list]
        min_idx = min(indices)
        max_idx = max(indices)
        idx_range = max(max_idx - min_idx, 1)

        # Place random priming hotspots along the genome
        # Scale hotspots to number based on genome size (more reads = larger genome)
        n_hotspots = max(3, min(n_hotspots_per_genome, len(read_list) // 50))
        hotspot_positions = [rng.uniform(0, 1) for _ in range(n_hotspots)]

        # Hotspot intensities: not all hotspots are equal
        # Some are strong priming sites, others weak
        hotspot_intensities = [rng.lognormvariate(0, 1.0) for _ in range(n_hotspots)]

        # Width in index space
        width = hotspot_width * idx_range

        for list_idx, read_idx in read_list:
            # Normalize read position to [0, 1]
            norm_pos = (read_idx - min_idx) / idx_range if idx_range > 0 else 0.5

            # Compute hotspot proximity weight: sum of Gaussian contributions
            # from all hotspots
            hotspot_weight = 0.0
            for hp, intensity in zip(hotspot_positions, hotspot_intensities):
                dist = abs(norm_pos - hp)
                # Wrap-around for circular genomes
                dist = min(dist, 1.0 - dist)
                hotspot_weight += intensity * 2.718 ** (-(dist ** 2) / (2 * (hotspot_width ** 2)))

            # Add baseline weight (even cold regions get some amplification)
            hotspot_weight = max(0.1, hotspot_weight)

            # GC bias (phi29 optimal at 40%)
            gc = _compute_gc(str(records[list_idx].seq))
            gc_weight = 2.718 ** (-((gc - 0.40) ** 2) / (2 * 0.10 ** 2))

            # Stochastic noise
            stochastic = rng.lognormvariate(0, 0.5)

            weight = hotspot_weight * gc_weight * stochastic
            weights[list_idx] = max(0.001, weight)

    return weights


def _compute_pcr_bias_weights(
    records: list[SeqRecord], rng: random.Random,
    amplification_method: str = "pcr",
) -> list[float]:
    """Compute per-read amplification bias weights.

    Different amplification methods have different bias profiles:
    - PCR/Linker/RdAB: Moderate GC bias (optimal ~50%), mild length bias
    - MDA (phi29): Extreme GC bias (optimal ~40%), strong stochasticity,
      within-genome hotspot clustering (reads near priming hotspots are
      much more likely to be duplicated)

    Args:
        records: List of SeqRecord objects.
        rng: Random number generator.
        amplification_method: "pcr" (default) or "mda".

    Returns:
        List of weights (higher = more likely to be a template).
    """
    if amplification_method == "mda":
        return _compute_mda_hotspot_weights(records, rng)

    # PCR/Linker/RdAB: standard PCR bias
    weights = []
    for record in records:
        seq = str(record.seq)
        read_len = len(seq)
        gc = _compute_gc(seq)

        # GC bias: peak efficiency at ~50% GC, drops at extremes
        gc_weight = 2.718 ** (-((gc - 0.50) ** 2) / (2 * 0.15 ** 2))

        # Length bias: shorter fragments amplify better
        len_weight = max(0.3, 1.0 - (read_len - 100) / 500)

        # Combined weight with some noise
        weight = gc_weight * len_weight * rng.uniform(0.7, 1.3)
        weights.append(max(0.01, weight))

    return weights


def add_pcr_duplicates(
    r1_path: Path,
    r2_path: Path,
    output_r1: Optional[Path] = None,
    output_r2: Optional[Path] = None,
    duplicate_rate: float = 0.30,
    max_copies: int = 5,
    copy_distribution: str = "geometric",
    error_rate: float = 0.001,
    gc_length_bias: bool = True,
    amplification_method: str = "pcr",
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Post-process FASTQ files to inject duplicate read pairs.

    Selects a fraction of read pairs as "templates" and generates
    1 to max_copies duplicate copies of each. Copies are near-identical
    with optional polymerase error rate. Duplicates are shuffled
    throughout the output file.

    Behavior changes based on amplification_method:
    - "pcr" (default): Geometric copy distribution, moderate GC/length
      bias (Taq optimal ~50% GC), 0.001 error rate.
    - "mda": Power-law copy distribution (heavy tail — few templates
      massively over-amplified), extreme GC bias (phi29 optimal ~40%),
      no length bias, 0.0001 error rate (phi29 high fidelity), higher
      max_copies (up to 20).

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
            "power_law": Heavy-tailed (MDA-like), few templates get
                very many copies. P(k) proportional to k^(-alpha).
        error_rate: Per-base substitution rate in copies (default: 0.001).
            Models polymerase errors. phi29 (MDA) has ~10x lower error
            rate than Taq (PCR).
        gc_length_bias: If True (default), weight template selection by
            GC content and fragment length, adjusted for amplification
            method.
        amplification_method: "pcr" (default) or "mda". Adjusts bias
            weights and copy distribution to match the amplification
            chemistry.
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

    if copy_distribution not in ("geometric", "uniform", "power_law"):
        raise ValueError(
            f"Unknown copy_distribution: {copy_distribution}. "
            f"Choose from: geometric, uniform, power_law"
        )

    # Apply MDA-specific defaults if not explicitly overridden
    if amplification_method == "mda":
        if max_copies == 5:  # default, not user-set
            max_copies = 20
        if error_rate == 0.001:  # default, not user-set
            error_rate = 0.0001  # phi29 has ~10x lower error rate
        if copy_distribution == "geometric":  # default, not user-set
            copy_distribution = "power_law"

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

    amp_label = "MDA (phi29)" if amplification_method == "mda" else "PCR"
    logger.info(
        f"Injecting {amp_label} duplicates ({duplicate_rate:.0%} template rate, "
        f"max {max_copies} copies, {copy_distribution} distribution, "
        f"error rate {error_rate}) into {r1_path.name}"
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

    # Select template indices (weighted by GC/length bias or uniform)
    if gc_length_bias and total_reads > 0:
        weights = _compute_pcr_bias_weights(r1_records, rng, amplification_method)
        # Weighted sampling without replacement
        indices = list(range(total_reads))
        selected = []
        remaining_weights = weights[:]
        remaining_indices = indices[:]
        for _ in range(min(n_templates, total_reads)):
            chosen = rng.choices(remaining_indices, weights=remaining_weights, k=1)[0]
            selected.append(chosen)
            idx = remaining_indices.index(chosen)
            remaining_indices.pop(idx)
            remaining_weights.pop(idx)
            if not remaining_indices:
                break
        template_indices = set(selected)
    else:
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
            elif copy_distribution == "power_law":
                # Power-law (MDA-like): heavy tail, few templates get
                # very many copies. Uses Zipf-like distribution.
                # alpha=1.5 gives a heavy tail: most get 1-2, some get 10+
                raw = int(rng.paretovariate(1.5))
                n_copies = min(max(1, raw), max_copies)
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
        "max_copies_used": max_copies,
        "error_rate_used": error_rate,
        "amplification_method": amplification_method,
        "copy_distribution": copy_distribution,
    }


def add_mda_chimeras(
    r1_path: Path,
    r2_path: Path,
    chimera_rate: float = 0.15,
    random_seed: Optional[int] = None,
    in_place: bool = False,
    manifest_path: Optional[Path] = None,
) -> dict:
    """Inject MDA chimeric reads into FASTQ files.

    MDA (phi29) creates chimeric artifacts when the polymerase
    displaces a growing strand and the displaced strand re-primes
    on a different template or region. The result is a read where
    the first portion comes from one genomic location and the
    second portion comes from another.

    Types of chimeras modeled:
    - Inter-genomic: Two halves from different genomes (most common)
    - Intra-genomic: Two halves from different regions of the same
      genome (inverted or direct)

    Chimeric reads are created by joining the first half of one read
    with the second half of another read. The junction point varies
    randomly (not always at the midpoint).

    Args:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        chimera_rate: Fraction of reads to replace with chimeras
            (0.0-0.3, default 0.15). Literature: 5-30% for MDA.
        random_seed: Random seed for reproducibility.
        in_place: If True, modify files in place.
        manifest_path: If set, write a TSV manifest of chimeric reads.

    Returns:
        Dict with statistics:
            - reads_total: Total reads in file
            - chimeras_created: Number of chimeric reads injected
            - chimera_fraction: Fraction of output that are chimeras
            - inter_genomic: Count of inter-genomic chimeras
            - intra_genomic: Count of intra-genomic chimeras
    """
    if not (0.0 <= chimera_rate <= 1.0):
        raise ValueError(f"chimera_rate must be 0-1, got {chimera_rate}")

    if chimera_rate == 0.0:
        return {
            "reads_total": 0, "chimeras_created": 0,
            "chimera_fraction": 0.0,
            "inter_genomic": 0, "intra_genomic": 0,
        }

    rng = random.Random(random_seed)

    logger.info(
        f"Injecting MDA chimeras ({chimera_rate:.0%} rate) into {r1_path.name}"
    )

    r1_records = list(SeqIO.parse(r1_path, "fastq"))
    r2_records = list(SeqIO.parse(r2_path, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 and R2 have different read counts: "
            f"{len(r1_records)} vs {len(r2_records)}"
        )

    total_reads = len(r1_records)
    n_chimeras = int(total_reads * chimera_rate)

    if n_chimeras == 0 or total_reads < 2:
        return {
            "reads_total": total_reads, "chimeras_created": 0,
            "chimera_fraction": 0.0,
            "inter_genomic": 0, "intra_genomic": 0,
        }

    # Parse genome IDs from read headers to distinguish inter vs intra
    def _get_genome_id(record):
        """Extract genome ID from ISS read header."""
        rid = record.id
        parts = rid.rsplit('_', 2)
        return parts[0] if len(parts) >= 3 else rid

    # Select reads to replace with chimeras
    chimera_indices = set(rng.sample(range(total_reads), min(n_chimeras, total_reads)))

    manifest_rows = []
    inter_count = 0
    intra_count = 0

    for idx in chimera_indices:
        # Pick a random donor read (different from this one)
        donor_idx = idx
        while donor_idx == idx:
            donor_idx = rng.randint(0, total_reads - 1)

        r1_target = r1_records[idx]
        r1_donor = r1_records[donor_idx]
        r2_target = r2_records[idx]
        r2_donor = r2_records[donor_idx]

        read_len = len(r1_target.seq)

        # Junction point: varies randomly, weighted toward center
        # (branch migration more likely in middle of displaced strand)
        junction = int(rng.gauss(read_len * 0.5, read_len * 0.15))
        junction = max(10, min(junction, read_len - 10))

        # Create chimeric R1: first part from target, second from donor
        chimera_seq = str(r1_target.seq[:junction]) + str(r1_donor.seq[junction:])
        chimera_qual = (
            r1_target.letter_annotations["phred_quality"][:junction]
            + r1_donor.letter_annotations["phred_quality"][junction:]
        )

        target_genome = _get_genome_id(r1_target)
        donor_genome = _get_genome_id(r1_donor)
        is_inter = target_genome != donor_genome

        if is_inter:
            inter_count += 1
        else:
            intra_count += 1

        chimera_tag = (
            f"mda_chimera=true junction={junction} "
            f"donor={r1_donor.id} type={'inter' if is_inter else 'intra'}"
        )

        # Replace R1 with chimera
        r1_records[idx] = SeqRecord(
            Seq(chimera_seq),
            id=r1_target.id,
            name=r1_target.name,
            description=f"{r1_target.description} {chimera_tag}".strip(),
        )
        r1_records[idx].letter_annotations["phred_quality"] = chimera_qual

        # R2 stays as-is (chimera affects one end of the insert)
        # But tag it for ground truth
        r2_desc = f"{r2_target.description} mda_chimera_mate=true".strip()
        r2_records[idx] = SeqRecord(
            r2_target.seq,
            id=r2_target.id,
            name=r2_target.name,
            description=r2_desc,
        )
        r2_records[idx].letter_annotations["phred_quality"] = (
            r2_target.letter_annotations["phred_quality"]
        )

        manifest_rows.append({
            "chimera_read_id": r1_target.id,
            "target_genome": target_genome,
            "donor_genome": donor_genome,
            "junction_pos": junction,
            "chimera_type": "inter_genomic" if is_inter else "intra_genomic",
        })

    # Write output
    output_r1 = r1_path if in_place else r1_path.parent / f"{r1_path.stem}_chimera{r1_path.suffix}"
    output_r2 = r2_path if in_place else r2_path.parent / f"{r2_path.stem}_chimera{r2_path.suffix}"

    SeqIO.write(r1_records, output_r1, "fastq")
    SeqIO.write(r2_records, output_r2, "fastq")

    # Write manifest
    if manifest_path and manifest_rows:
        with open(manifest_path, "w") as f:
            f.write("chimera_read_id\ttarget_genome\tdonor_genome\t"
                    "junction_pos\tchimera_type\n")
            for row in manifest_rows:
                f.write(
                    f"{row['chimera_read_id']}\t{row['target_genome']}\t"
                    f"{row['donor_genome']}\t{row['junction_pos']}\t"
                    f"{row['chimera_type']}\n"
                )
        logger.info(f"Wrote chimera manifest: {manifest_path}")

    chimera_fraction = len(chimera_indices) / total_reads

    logger.info(
        f"MDA chimera injection complete: {len(chimera_indices)} chimeras "
        f"({chimera_fraction:.1%}), {inter_count} inter-genomic, "
        f"{intra_count} intra-genomic"
    )

    return {
        "reads_total": total_reads,
        "chimeras_created": len(chimera_indices),
        "chimera_fraction": chimera_fraction,
        "inter_genomic": inter_count,
        "intra_genomic": intra_count,
    }
