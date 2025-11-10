"""
Long-read sequencing simulation using PBSIM3.

This module provides functions to generate realistic long-read sequencing data
(PacBio HiFi and Oxford Nanopore) from a mock virome composition using PBSIM3
as the underlying simulator.

Key Features:
- PacBio HiFi: High-accuracy reads (>99.9%, QV20+) via multi-pass CCS
- Oxford Nanopore: Ultra-long reads with characteristic homopolymer errors
- Compatible with ViroForge VLP enrichment and contamination workflows
- Complete ground truth tracking (read-to-genome mappings)
- Reproducible with random seeds

Dependencies:
- PBSIM3: Install with `conda install -c bioconda pbsim3`
- PacBio ccs (for HiFi only): `conda install -c bioconda pbccs`
- SAMtools: `conda install -c bioconda samtools`

Example:
    from viroforge.utils import create_mock_virome
    from viroforge.simulators.longread import generate_long_reads, LongReadPlatform

    # Create composition
    composition = create_mock_virome(
        name='gut_virome',
        body_site='gut',
        contamination_level='realistic'
    )

    # Generate PacBio HiFi reads
    output = generate_long_reads(
        composition=composition,
        output_prefix='my_dataset',
        platform=LongReadPlatform.PACBIO_HIFI,
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
from typing import Optional, Dict, Union
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
        PACBIO_HIFI: PacBio HiFi (Circular Consensus Sequencing)
                    - Accuracy: >99.9% (QV20+)
                    - Read length: 10-30kb typical
                    - Applications: Complete genome assembly, variant calling

        NANOPORE: Oxford Nanopore Technologies
                 - Accuracy: ~95% (R10.4 chemistry)
                 - Read length: 10kb-2Mb (ultra-long possible)
                 - Applications: Structural variants, phasing, assembly
    """
    PACBIO_HIFI = "pacbio-hifi"
    NANOPORE = "nanopore"


@dataclass
class PacBioHiFiConfig:
    """
    Configuration for PacBio HiFi simulation.

    PacBio HiFi uses circular consensus sequencing (CCS) to generate
    high-accuracy reads (QV20+, >99.9% accuracy) from multiple passes
    of the same ZMW (zero-mode waveguide).

    Attributes:
        passes: Number of passes for multi-pass sequencing (default: 10)
               More passes → higher accuracy but shorter read length
        min_passes: Minimum passes required for CCS calling (default: 3)
        accuracy_model: PBSIM3 error model to use (default: "QSHMM-RSII")
                       Options: "QSHMM-RSII", "QSHMM-SEQUEL"
        read_length_mean: Mean read length in bp (default: 15000)
                         Typical range: 10-25kb
        read_length_sd: Read length standard deviation (default: 5000)
        clr_error_rate: CLR (continuous long read) error rate before CCS (default: 0.15)
                       15% is typical for PacBio CLR

    Example:
        >>> config = PacBioHiFiConfig(passes=15, read_length_mean=20000)
        >>> # Higher passes and longer reads (lower throughput)
    """
    passes: int = 10
    min_passes: int = 3
    accuracy_model: str = "QSHMM-RSII"
    read_length_mean: int = 15000
    read_length_sd: int = 5000
    clr_error_rate: float = 0.15


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


def check_pbccs_installed() -> bool:
    """
    Check if PacBio ccs software is installed (required for HiFi simulation).

    Returns:
        bool: True if ccs is installed, False otherwise

    Note:
        Only required for PacBio HiFi simulation, not for Nanopore.

    Example:
        >>> if not check_pbccs_installed():
        ...     print("Please install ccs: conda install -c bioconda pbccs")
    """
    try:
        result = subprocess.run(
            ['ccs', '--version'],
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


def _run_pbsim3_clr(
    genomes_fasta: Path,
    output_prefix: str,
    depths: Dict[str, float],
    config: PacBioHiFiConfig,
    seed: Optional[int]
) -> Path:
    """
    Run PBSIM3 to generate CLR (Continuous Long Reads) with multi-pass sequencing.

    This is Step 1 of PacBio HiFi simulation. The output CLR BAM file
    will be input to PacBio ccs for consensus calling.

    Args:
        genomes_fasta: Path to multi-FASTA with all genomes
        output_prefix: Output file prefix
        depths: Dict mapping genome_id to sequencing depth
        config: PacBioHiFiConfig object
        seed: Random seed for reproducibility

    Returns:
        Path to output BAM file containing CLR with multiple passes

    Raises:
        RuntimeError: If PBSIM3 fails

    Note:
        PBSIM3 with --pass-num ≥2 enables multi-pass sequencing simulation.
        Output is in BAM format containing subreads from multiple passes.
    """
    # Build PBSIM3 command
    # PBSIM3 simulates per-genome, so we need to run it for each genome
    # and merge results

    logger.info(f"Running PBSIM3 to generate CLR with {config.passes} passes...")
    logger.info(f"Read length: {config.read_length_mean} ± {config.read_length_sd} bp")

    # For simplicity in initial implementation, simulate all genomes together
    # with average depth
    avg_depth = sum(depths.values()) / len(depths) if depths else 10.0

    cmd = [
        'pbsim',
        '--strategy', 'wgs',
        '--method', 'qshmm',
        '--qshmm', config.accuracy_model,
        '--depth', str(avg_depth),
        '--genome', str(genomes_fasta),
        '--pass-num', str(config.passes),
        '--length-mean', str(config.read_length_mean),
        '--length-sd', str(config.read_length_sd),
        '--accuracy-mean', str(1.0 - config.clr_error_rate),
        '--prefix', output_prefix,
    ]

    if seed is not None:
        cmd.extend(['--seed', str(seed)])

    logger.info(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info("PBSIM3 CLR generation completed successfully")
        logger.debug(f"PBSIM3 stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"PBSIM3 failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise RuntimeError(f"PBSIM3 failed: {e.stderr}")

    # PBSIM3 output file (BAM format for multi-pass)
    bam_path = Path(f"{output_prefix}.sam")  # PBSIM3 outputs SAM by default

    if not bam_path.exists():
        raise RuntimeError(f"PBSIM3 did not create expected output: {bam_path}")

    logger.info(f"CLR generated: {bam_path}")
    return bam_path


def _run_pacbio_ccs(
    clr_sam: Path,
    output_prefix: str,
    config: PacBioHiFiConfig
) -> Path:
    """
    Run PacBio ccs to generate HiFi consensus reads from CLR.

    This is Step 2 of PacBio HiFi simulation. Takes multi-pass CLR
    and calls consensus to generate high-accuracy HiFi reads.

    Args:
        clr_sam: Path to CLR SAM file from PBSIM3
        output_prefix: Output file prefix
        config: PacBioHiFiConfig object

    Returns:
        Path to HiFi FASTQ file

    Raises:
        RuntimeError: If ccs fails

    Note:
        PacBio ccs requires BAM input. We'll convert SAM→BAM first.
        Output is FASTQ.GZ format.
    """
    logger.info("Converting SAM to BAM for ccs...")

    # Convert SAM to BAM using samtools
    bam_path = Path(f"{output_prefix}_clr.bam")

    sam_to_bam_cmd = [
        'samtools', 'view',
        '-b',  # Output BAM
        '-o', str(bam_path),
        str(clr_sam)
    ]

    try:
        result = subprocess.run(
            sam_to_bam_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info(f"Converted to BAM: {bam_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"SAM to BAM conversion failed: {e.stderr}")
        raise RuntimeError(f"samtools view failed: {e.stderr}")

    # Run PacBio ccs
    logger.info(f"Running PacBio ccs (min {config.min_passes} passes)...")

    hifi_fastq = Path(f"{output_prefix}_hifi.fastq.gz")

    ccs_cmd = [
        'ccs',
        str(bam_path),
        str(hifi_fastq),
        '--min-passes', str(config.min_passes),
        '--min-rq', '0.99',  # Minimum read quality (Q20)
        '--log-level', 'INFO'
    ]

    logger.info(f"Running: {' '.join(ccs_cmd)}")

    try:
        result = subprocess.run(
            ccs_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info("PacBio ccs completed successfully")
        logger.debug(f"ccs stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"PacBio ccs failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise RuntimeError(f"PacBio ccs failed: {e.stderr}")

    if not hifi_fastq.exists():
        raise RuntimeError(f"ccs did not create expected output: {hifi_fastq}")

    logger.info(f"HiFi reads generated: {hifi_fastq}")
    return hifi_fastq


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
    logger.info("Running PBSIM3 to generate Nanopore reads...")
    logger.info(f"Chemistry: {config.chemistry}")
    logger.info(f"Read length: {config.read_length_mean} ± {config.read_length_sd} bp")
    logger.info(f"Homopolymer deletion bias: {config.hp_del_bias}")

    # Average depth across all genomes
    avg_depth = sum(depths.values()) / len(depths) if depths else 10.0

    cmd = [
        'pbsim',
        '--strategy', 'wgs',
        '--method', 'errhmm',
        '--errhmm', f'ERRHMM-ONT',  # ONT error model
        '--depth', str(avg_depth),
        '--genome', str(genomes_fasta),
        '--length-mean', str(config.read_length_mean),
        '--length-sd', str(config.read_length_sd),
        '--accuracy-mean', str(1.0 - config.error_rate),
        '--hp-del-bias', str(config.hp_del_bias),
        '--prefix', output_prefix,
    ]

    if seed is not None:
        cmd.extend(['--seed', str(seed)])

    logger.info(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info("PBSIM3 Nanopore generation completed successfully")
        logger.debug(f"PBSIM3 stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"PBSIM3 failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise RuntimeError(f"PBSIM3 failed: {e.stderr}")

    # PBSIM3 creates .fastq files (possibly multiple for different chromosomes)
    # Look for output FASTQ
    fastq_pattern = f"{output_prefix}*.fastq"
    import glob
    fastq_files = glob.glob(fastq_pattern)

    if not fastq_files:
        raise RuntimeError(f"PBSIM3 did not create expected FASTQ output: {fastq_pattern}")

    # If multiple files, merge them (typical for multi-chromosome genomes)
    if len(fastq_files) == 1:
        fastq_path = Path(fastq_files[0])
    else:
        # Merge multiple FASTQ files
        merged_fastq = Path(f"{output_prefix}_merged.fastq")
        with open(merged_fastq, 'w') as outf:
            for fq in sorted(fastq_files):
                with open(fq) as inf:
                    outf.write(inf.read())
        fastq_path = merged_fastq
        logger.info(f"Merged {len(fastq_files)} FASTQ files into {fastq_path}")

    logger.info(f"Nanopore reads generated: {fastq_path}")
    return fastq_path


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
    platform_config: Optional[Union[PacBioHiFiConfig, NanoporeConfig]] = None,
    validate_output: bool = True,
    random_seed: Optional[int] = None,
    keep_temp_files: bool = False
) -> Dict[str, Path]:
    """
    Generate long reads from a mock virome composition.

    This function:
    1. Writes genomes to temporary FASTA file
    2. Calls PBSIM3 to generate reads (±PacBio ccs for HiFi)
    3. Creates ground truth metadata
    4. Cleans up temporary files

    Args:
        composition: MockViromeComposition object containing viral genomes
                    and contamination profile
        output_prefix: Prefix for output files (e.g., 'my_dataset')
                      Files created:
                      - PacBio HiFi: my_dataset_hifi.fastq.gz
                      - Nanopore: my_dataset.fastq
                      - Both: my_dataset_ground_truth_genomes.tsv
        platform: Long-read platform (PACBIO_HIFI or NANOPORE)
        depth: Sequencing depth (coverage) (default: 10.0)
              Higher depth → more reads per genome
        platform_config: Platform-specific configuration
                        PacBioHiFiConfig for HiFi, NanoporeConfig for Nanopore
                        If None, uses default configuration
        validate_output: Validate output files after generation (default: True)
        random_seed: Random seed for reproducibility (default: None)
        keep_temp_files: Keep temporary files for debugging (default: False)

    Returns:
        Dict containing paths to generated files:
        {
            'reads': Path to FASTQ file (HiFi: .fastq.gz, Nanopore: .fastq),
            'ground_truth': Path to ground truth TSV file,
            'temp_fasta': Path to temp FASTA (if keep_temp_files=True)
        }

    Raises:
        RuntimeError: If PBSIM3/ccs is not installed or fails
        ValueError: If depth ≤ 0 or invalid platform

    Example:
        >>> from viroforge.utils import create_mock_virome
        >>> from viroforge.simulators.longread import (
        ...     generate_long_reads,
        ...     LongReadPlatform
        ... )
        >>>
        >>> # Create composition
        >>> composition = create_mock_virome(
        ...     name='gut_virome',
        ...     body_site='gut',
        ...     contamination_level='realistic'
        ... )
        >>>
        >>> # Generate PacBio HiFi reads at 10x depth
        >>> output = generate_long_reads(
        ...     composition=composition,
        ...     output_prefix='output/gut_hifi',
        ...     platform=LongReadPlatform.PACBIO_HIFI,
        ...     depth=10.0,
        ...     random_seed=42
        ... )
        >>>
        >>> print(f"HiFi reads: {output['reads']}")

    Note:
        PacBio HiFi Workflow (two-step):
        1. PBSIM3 generates CLR with multi-pass sequencing
        2. PacBio ccs calls consensus to generate HiFi reads

        Nanopore Workflow (single-step):
        1. PBSIM3 generates ONT reads with homopolymer errors

    See Also:
        - PBSIM3 documentation: https://github.com/yukiteruono/pbsim3
        - PacBio ccs: https://ccs.how/
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

    # Check ccs for PacBio HiFi
    if platform == LongReadPlatform.PACBIO_HIFI and not check_pbccs_installed():
        raise RuntimeError(
            "PacBio ccs is not installed (required for HiFi). Install with: "
            "conda install -c bioconda pbccs"
        )

    # Use default config if not provided
    if platform_config is None:
        if platform == LongReadPlatform.PACBIO_HIFI:
            platform_config = PacBioHiFiConfig()
        else:
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

        # Generate reads based on platform
        if platform == LongReadPlatform.PACBIO_HIFI:
            # Two-step: PBSIM3 CLR → ccs HiFi
            clr_sam = _run_pbsim3_clr(
                genomes_fasta=genomes_fasta,
                output_prefix=str(temp_dir_path / 'clr'),
                depths=depths,
                config=platform_config,
                seed=random_seed
            )

            reads_path = _run_pacbio_ccs(
                clr_sam=clr_sam,
                output_prefix=output_prefix,
                config=platform_config
            )

        else:  # NANOPORE
            # Single-step: PBSIM3 ONT
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
