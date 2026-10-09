"""
Long-read sequencing simulation using PBSIM3.

This module provides functions to generate realistic Oxford Nanopore
long-read sequencing data from a mock virome composition using PBSIM3
as the underlying simulator.

Key Features:
- Oxford Nanopore: Ultra-long reads with characteristic homopolymer errors
- Compatible with ViroForge VLP enrichment and contamination workflows
- Complete ground truth tracking (read-to-genome mappings)
- Reproducible with random seeds

Dependencies:
- PBSIM3: Install with `conda install -c bioconda pbsim3`

Example:
    from viroforge.utils import create_mock_virome
    from viroforge.simulators.longread import generate_long_reads, LongReadPlatform

    # Create composition
    composition = create_mock_virome(
        name='gut_virome',
        body_site='gut',
        contamination_level='realistic'
    )

    # Generate Nanopore reads
    output = generate_long_reads(
        composition=composition,
        output_prefix='my_dataset',
        platform=LongReadPlatform.NANOPORE,
        depth=10.0,
        random_seed=42
    )

Author: ViroForge Development Team
Date: 2025-11-10
Phase: 10 - Long-Read Sequencing Support
"""

import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Optional, Dict
from dataclasses import dataclass
from enum import Enum
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils.validation import validate_output_directory

logger = logging.getLogger(__name__)


class LongReadPlatform(Enum):
    """
    Long-read sequencing platforms.

    Values:
        NANOPORE: Oxford Nanopore Technologies
                 - Accuracy: ~95% (R10.4 chemistry)
                 - Read length: 10kb-2Mb (ultra-long possible)
                 - Applications: Structural variants, phasing, assembly
    """
    NANOPORE = "nanopore"


@dataclass
class NanoporeConfig:
    """
    Configuration for Oxford Nanopore simulation.

    Oxford Nanopore generates ultra-long reads with characteristic
    homopolymer errors and quality-length relationships.

    Attributes:
        chemistry: Nanopore pore chemistry version (default: "R10.4")
                  Options: "R9.4" (older), "R10.4" (current)
        read_length_mean: Mean read length in bp (default: 20000)
                         Nanopore can generate much longer (100kb-2Mb)
        read_length_sd: Read length standard deviation (default: 10000)
        error_rate: Base error rate (default: 0.05, i.e., 5%)
                   R10.4 chemistry: ~5%, R9.4: ~10%
        hp_del_bias: Homopolymer deletion bias parameter (default: 6)
                    Higher values → stronger homopolymer deletion bias
                    This is a characteristic Nanopore error mode
        quality_mean: Mean quality score (default: 10)

    Example:
        >>> config = NanoporeConfig(chemistry="R10.4", read_length_mean=50000)
        >>> # Ultra-long reads configuration
    """
    chemistry: str = "R10.4"
    read_length_mean: int = 20000
    read_length_sd: int = 10000
    error_rate: float = 0.05
    hp_del_bias: int = 6
    quality_mean: int = 10


def check_pbsim3_installed() -> bool:
    """
    Check if PBSIM3 is installed and available.

    Returns:
        bool: True if PBSIM3 is installed, False otherwise

    Example:
        >>> if not check_pbsim3_installed():
        ...     print("Please install PBSIM3: conda install -c bioconda pbsim3")
    """
    try:
        result = subprocess.run(
            ['pbsim', '--version'],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False




def _write_genome_fasta(composition, output_path: Path) -> Dict[str, str]:
    """
    Write all genomes to a multi-FASTA file for PBSIM3.

    Args:
        composition: MockViromeComposition object
        output_path: Path to output FASTA file

    Returns:
        Dict mapping FASTA IDs to genome types ('viral' or contaminant type)

    Note:
        PBSIM3 requires unique sequence IDs. We use the genome_id
        from the composition, ensuring uniqueness.

        This is the same format as used by InSilicoSeq (illumina.py),
        allowing code reuse.
    """
    genome_types = {}
    records = []

    # Add viral genomes
    for genome in composition.viral_community.genomes:
        record = SeqRecord(
            Seq(str(genome.sequence)),
            id=genome.genome_id,
            description=f"viral|{genome.taxonomy}"
        )
        records.append(record)
        genome_types[genome.genome_id] = 'viral'

    # Add contaminant genomes
    if composition.contamination_profile:
        for contaminant in composition.contamination_profile.contaminants:
            record = SeqRecord(
                Seq(str(contaminant.sequence)),
                id=contaminant.genome_id,
                description=f"{contaminant.contaminant_type.value}|{contaminant.organism}"
            )
            records.append(record)
            genome_types[contaminant.genome_id] = contaminant.contaminant_type.value

    # Write to file
    SeqIO.write(records, output_path, "fasta")
    logger.info(f"Wrote {len(records)} genomes to {output_path}")

    return genome_types


def _calculate_depth_per_genome(
    composition,
    total_depth: float
) -> Dict[str, float]:
    """
    Calculate sequencing depth for each genome based on abundances.

    Args:
        composition: MockViromeComposition object
        total_depth: Total sequencing depth (coverage)

    Returns:
        Dict mapping genome_id to depth (coverage)

    Note:
        PBSIM3 requires per-genome depth specification.
        depth_i = total_depth * abundance_i

    Example:
        If total_depth=10 and genome has abundance=0.2,
        then genome depth = 10 * 0.2 = 2x coverage
    """
    depths = {}

    # Viral genomes
    for genome in composition.viral_community.genomes:
        depths[genome.genome_id] = total_depth * genome.abundance

    # Contaminants
    if composition.contamination_profile:
        for contaminant in composition.contamination_profile.contaminants:
            depths[contaminant.genome_id] = total_depth * contaminant.abundance

    return depths


def _run_pbsim3_nanopore(
    genomes_fasta: Path,
    output_prefix: str,
    depths: Dict[str, float],
    config: NanoporeConfig,
    seed: Optional[int]
) -> Path:
    """
    Run PBSIM3 to generate Nanopore reads.

    Single-step process for Nanopore simulation using PBSIM3's
    ONT error models with homopolymer deletion bias.

    Runs PBSIM3 separately per genome with abundance-weighted depth,
    then merges all FASTQ outputs. This ensures the community composition
    is preserved (PBSIM3 only accepts a single --depth value).

    Args:
        genomes_fasta: Path to multi-FASTA with all genomes
        output_prefix: Output file prefix
        depths: Dict mapping genome_id to sequencing depth
        config: NanoporeConfig object
        seed: Random seed for reproducibility

    Returns:
        Path to FASTQ file containing Nanopore reads

    Raises:
        RuntimeError: If PBSIM3 fails
    """
    import glob

    logger.info("Running PBSIM3 to generate Nanopore reads...")
    logger.info(f"Chemistry: {config.chemistry}")
    logger.info(f"Read length: {config.read_length_mean} ± {config.read_length_sd} bp")
    logger.info(f"Homopolymer deletion bias: {config.hp_del_bias}")

    # Parse multi-FASTA into individual genome records
    genomes = list(SeqIO.parse(genomes_fasta, "fasta"))
    logger.info(f"Simulating {len(genomes)} genomes with per-genome depth")

    # Skip genomes with negligible depth (< 0.01x) to avoid empty PBSIM3 runs
    min_depth = 0.01
    skipped = 0

    merged_fastq = Path(f"{output_prefix}_merged.fastq")
    total_reads = 0

    with open(merged_fastq, 'w') as outf:
        for i, genome in enumerate(genomes):
            genome_depth = depths.get(genome.id, 0.0)
            if genome_depth < min_depth:
                skipped += 1
                continue

            # Write single-genome FASTA
            genome_fasta = Path(f"{output_prefix}_genome_{i:04d}.fasta")
            SeqIO.write([genome], str(genome_fasta), "fasta")

            genome_prefix = f"{output_prefix}_genome_{i:04d}"

            cmd = [
                'pbsim',
                '--strategy', 'wgs',
                '--method', 'errhmm',
                '--errhmm', 'ERRHMM-ONT',
                '--depth', str(genome_depth),
                '--genome', str(genome_fasta),
                '--length-mean', str(config.read_length_mean),
                '--length-sd', str(config.read_length_sd),
                '--accuracy-mean', str(1.0 - config.error_rate),
                '--hp-del-bias', str(config.hp_del_bias),
                '--prefix', genome_prefix,
            ]

            if seed is not None:
                cmd.extend(['--seed', str(seed + i)])

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                logger.warning(f"PBSIM3 failed for genome {genome.id}: {e.stderr}")
                continue

            # Collect FASTQ output for this genome
            # PBSIM3 outputs .fq.gz (gzipped) files
            import gzip
            fastq_files = glob.glob(f"{genome_prefix}*.fq.gz")
            if not fastq_files:
                # Fall back to uncompressed .fastq
                fastq_files = glob.glob(f"{genome_prefix}*.fastq")
            for fq in sorted(fastq_files):
                opener = gzip.open if fq.endswith('.gz') else open
                with opener(fq, 'rt') as inf:
                    for line in inf:
                        outf.write(line)
                        if line.startswith('@'):
                            total_reads += 1
                Path(fq).unlink()

            # Clean up per-genome temp files
            genome_fasta.unlink(missing_ok=True)
            for f in glob.glob(f"{genome_prefix}*.ref"):
                Path(f).unlink(missing_ok=True)
            for f in glob.glob(f"{genome_prefix}*.maf*"):
                Path(f).unlink(missing_ok=True)

            if (i + 1) % 50 == 0:
                logger.info(f"  Processed {i + 1}/{len(genomes)} genomes ({total_reads} reads so far)")

    if skipped > 0:
        logger.info(f"Skipped {skipped} genomes with depth < {min_depth}x")
    logger.info(f"Nanopore reads generated: {merged_fastq} ({total_reads} reads)")
    return merged_fastq


def _create_ground_truth_mapping(
    composition,
    output_dir: Path,
    genome_types: Dict[str, str],
    platform: LongReadPlatform
) -> Path:
    """
    Create ground truth file mapping genomes to their properties.

    Args:
        composition: MockViromeComposition object
        output_dir: Output directory
        genome_types: Dict mapping genome_id to type
        platform: Long-read platform used

    Returns:
        Path to ground truth TSV file

    Note:
        Extended from short-read ground truth to include platform
        and read_type fields for downstream analysis.
    """
    ground_truth = []

    # Add viral genomes
    for genome in composition.viral_community.genomes:
        ground_truth.append({
            'genome_id': genome.genome_id,
            'genome_type': 'viral',
            'taxonomy': genome.taxonomy,
            'length': genome.length,
            'gc_content': genome.gc_content,
            'abundance': genome.abundance,
            'source': 'viral_community',
            'platform': platform.value,
            'read_type': 'long'
        })

    # Add contaminants
    if composition.contamination_profile:
        for contaminant in composition.contamination_profile.contaminants:
            ground_truth.append({
                'genome_id': contaminant.genome_id,
                'genome_type': contaminant.contaminant_type.value,
                'taxonomy': contaminant.organism,
                'length': contaminant.length,
                'gc_content': contaminant.gc_content,
                'abundance': contaminant.abundance,
                'source': 'contamination',
                'platform': platform.value,
                'read_type': 'long'
            })

    # Create DataFrame
    df = pd.DataFrame(ground_truth)

    # Save to file
    output_path = output_dir / 'ground_truth_genomes.tsv'
    df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote ground truth for {len(df)} genomes to {output_path}")

    return output_path


def generate_long_reads(
    composition,
    output_prefix: str,
    platform: LongReadPlatform,
    depth: float = 10.0,
    platform_config: Optional[NanoporeConfig] = None,
    validate_output: bool = True,
    random_seed: Optional[int] = None,
    keep_temp_files: bool = False
) -> Dict[str, Path]:
    """
    Generate long reads from a mock virome composition.

    This function:
    1. Writes genomes to temporary FASTA file
    2. Calls PBSIM3 to generate Nanopore reads
    3. Creates ground truth metadata
    4. Cleans up temporary files

    Args:
        composition: MockViromeComposition object containing viral genomes
                    and contamination profile
        output_prefix: Prefix for output files (e.g., 'my_dataset')
                      Files created: my_dataset.fastq,
                      my_dataset_ground_truth_genomes.tsv
        platform: Long-read platform (NANOPORE)
        depth: Sequencing depth (coverage) (default: 10.0)
        platform_config: NanoporeConfig object (default: None uses defaults)
        validate_output: Validate output files after generation (default: True)
        random_seed: Random seed for reproducibility (default: None)
        keep_temp_files: Keep temporary files for debugging (default: False)

    Returns:
        Dict containing paths to generated files:
        {
            'reads': Path to FASTQ file,
            'ground_truth': Path to ground truth TSV file,
            'temp_fasta': Path to temp FASTA (if keep_temp_files=True)
        }

    Raises:
        RuntimeError: If PBSIM3 is not installed or fails
        ValueError: If depth <= 0

    See Also:
        - PBSIM3 documentation: https://github.com/yukiteruono/pbsim3
    """
    # Validate inputs
    if depth <= 0:
        raise ValueError(f"Depth must be positive, got {depth}")

    # Check PBSIM3 is installed
    if not check_pbsim3_installed():
        raise RuntimeError(
            "PBSIM3 is not installed. Install with: "
            "conda install -c bioconda pbsim3"
        )

    # Use default config if not provided
    if platform_config is None:
        platform_config = NanoporeConfig()

    logger.info(f"Generating {platform.value} reads at {depth}x depth")
    logger.info(f"Output prefix: {output_prefix}")

    # Validate output directory
    output_path = Path(output_prefix)
    output_dir = output_path.parent if output_path.parent != Path('.') else Path.cwd()
    validate_output_directory(output_dir, create=True)

    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix='viroforge_longread_')
    temp_dir_path = Path(temp_dir)
    logger.debug(f"Created temporary directory: {temp_dir}")

    try:
        # Write genomes to FASTA
        genomes_fasta = temp_dir_path / 'genomes.fasta'
        genome_types = _write_genome_fasta(composition, genomes_fasta)

        # Calculate per-genome depths
        depths = _calculate_depth_per_genome(composition, depth)

        # Generate Nanopore reads
        reads_path = _run_pbsim3_nanopore(
            genomes_fasta=genomes_fasta,
            output_prefix=str(temp_dir_path / 'nanopore'),
            depths=depths,
            config=platform_config,
            seed=random_seed
        )

        # Move to output location
        final_reads = Path(f"{output_prefix}.fastq")
        shutil.move(str(reads_path), str(final_reads))
        reads_path = final_reads

        # Create ground truth mapping
        ground_truth_path = _create_ground_truth_mapping(
            composition,
            output_dir,
            genome_types,
            platform
        )

        # Prepare return dictionary
        result = {
            'reads': reads_path,
            'ground_truth': ground_truth_path,
        }

        # Keep temp files if requested
        if keep_temp_files:
            kept_fasta = output_dir / 'input_genomes.fasta'
            shutil.copy(genomes_fasta, kept_fasta)
            result['temp_fasta'] = kept_fasta
            logger.info(f"Kept temporary files in {output_dir}")

        logger.info("✓ Long-read generation complete!")
        logger.info(f"  Reads: {reads_path}")
        logger.info(f"  Ground truth: {ground_truth_path}")

        return result

    finally:
        # Clean up temporary directory
        if not keep_temp_files:
            shutil.rmtree(temp_dir)
            logger.debug(f"Removed temporary directory: {temp_dir}")
