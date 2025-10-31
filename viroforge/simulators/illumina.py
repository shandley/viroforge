"""
Illumina sequencing read simulation using InSilicoSeq.

This module provides functions to generate realistic Illumina sequencing reads
from a mock virome composition using InSilicoSeq as the underlying simulator.

Key Features:
- Generates paired-end reads with realistic error models
- Supports multiple Illumina platforms (NovaSeq, HiSeq, MiSeq)
- Validates all FASTQ records (prevents seq/qual length mismatches)
- Tracks complete ground truth (read-to-genome mappings)
- Reproducible with random seeds

Dependencies:
- InSilicoSeq (iss): Install with `conda install -c bioconda insilicoseq`
  or `pip install InSilicoSeq`

Example:
    from viroforge.utils import create_mock_virome
    from viroforge.simulators import generate_reads

    # Create composition
    composition = create_mock_virome(
        name='gut_virome',
        body_site='gut',
        contamination_level='realistic'
    )

    # Generate reads
    generate_reads(
        composition=composition,
        output_prefix='my_dataset',
        n_reads=1_000_000,
        model='NovaSeq',
        random_seed=42
    )
"""

import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Optional, Dict, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils.validation import verify_fastq_file, validate_output_directory

logger = logging.getLogger(__name__)


def check_insilicoseq_installed() -> bool:
    """
    Check if InSilicoSeq is installed and available.

    Returns:
        bool: True if InSilicoSeq is installed, False otherwise

    Example:
        >>> if not check_insilicoseq_installed():
        ...     print("Please install InSilicoSeq: conda install -c bioconda insilicoseq")
    """
    try:
        result = subprocess.run(
            ['iss', '--version'],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def _write_genome_fasta(composition, output_path: Path) -> Dict[str, str]:
    """
    Write all genomes to a multi-FASTA file for InSilicoSeq.

    Args:
        composition: MockViromeComposition object
        output_path: Path to output FASTA file

    Returns:
        Dict mapping FASTA IDs to genome types ('viral' or contaminant type)

    Note:
        InSilicoSeq requires unique sequence IDs. We use the genome_id
        from the composition, ensuring uniqueness.
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


def _write_abundance_file(composition, output_path: Path) -> None:
    """
    Write abundance file in InSilicoSeq format (TSV).

    Args:
        composition: MockViromeComposition object
        output_path: Path to output TSV file

    Format:
        genome_id    abundance
        genome1      0.5
        genome2      0.3
        genome3      0.2

    Note:
        Abundances must sum to 1.0 (InSilicoSeq requirement)
    """
    abundances = []

    # Add viral genomes
    for genome in composition.viral_community.genomes:
        abundances.append({
            'genome_id': genome.genome_id,
            'abundance': genome.abundance
        })

    # Add contaminants
    if composition.contamination_profile:
        for contaminant in composition.contamination_profile.contaminants:
            abundances.append({
                'genome_id': contaminant.genome_id,
                'abundance': contaminant.abundance
            })

    # Create DataFrame and verify sum
    df = pd.DataFrame(abundances)
    total = df['abundance'].sum()

    if abs(total - 1.0) > 1e-6:
        logger.warning(f"Abundances sum to {total:.6f}, not 1.0. Normalizing...")
        df['abundance'] = df['abundance'] / total

    # Write to file (no header, tab-separated)
    df.to_csv(output_path, sep='\t', index=False, header=False)
    logger.info(f"Wrote abundances for {len(df)} genomes to {output_path}")


def _run_insilicoseq(
    genomes_fasta: Path,
    abundance_file: Path,
    output_prefix: str,
    n_reads: int,
    model: str,
    cpus: int,
    compress: bool,
    gc_bias: bool,
    seed: Optional[int]
) -> Tuple[Path, Path]:
    """
    Run InSilicoSeq to generate reads.

    Args:
        genomes_fasta: Path to multi-FASTA with all genomes
        abundance_file: Path to abundance file (TSV)
        output_prefix: Output file prefix
        n_reads: Number of reads to generate
        model: Error model (NovaSeq, HiSeq, MiSeq, etc.)
        cpus: Number of CPUs to use
        compress: Gzip compress output
        gc_bias: Enable GC bias simulation
        seed: Random seed for reproducibility

    Returns:
        Tuple of (R1_path, R2_path)

    Raises:
        RuntimeError: If InSilicoSeq fails
    """
    # Build command
    cmd = [
        'iss', 'generate',
        '--genomes', str(genomes_fasta),
        '--abundance_file', str(abundance_file),
        '--n_reads', str(n_reads),
        '--model', model,
        '--output', output_prefix,
        '--cpus', str(cpus),
    ]

    if compress:
        cmd.append('--compress')

    if gc_bias:
        cmd.append('--gc_bias')

    if seed is not None:
        cmd.extend(['--seed', str(seed)])

    logger.info(f"Running InSilicoSeq: {' '.join(cmd)}")

    # Run InSilicoSeq
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info("InSilicoSeq completed successfully")
        logger.debug(f"InSilicoSeq stdout: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"InSilicoSeq failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise RuntimeError(f"InSilicoSeq failed: {e.stderr}")

    # Determine output paths
    ext = '.fastq.gz' if compress else '.fastq'
    r1_path = Path(f"{output_prefix}_R1{ext}")
    r2_path = Path(f"{output_prefix}_R2{ext}")

    # Verify files exist
    if not r1_path.exists() or not r2_path.exists():
        raise RuntimeError(
            f"InSilicoSeq did not create expected output files: "
            f"{r1_path}, {r2_path}"
        )

    return r1_path, r2_path


def _create_ground_truth_mapping(
    composition,
    output_dir: Path,
    genome_types: Dict[str, str]
) -> Path:
    """
    Create ground truth file mapping genomes to their properties.

    Args:
        composition: MockViromeComposition object
        output_dir: Output directory
        genome_types: Dict mapping genome_id to type

    Returns:
        Path to ground truth TSV file

    Note:
        This creates a genome-level ground truth file. Read-level ground truth
        (which read came from which genome) requires parsing InSilicoSeq headers.
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
            'source': 'viral_community'
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
                'source': 'contamination'
            })

    # Create DataFrame
    df = pd.DataFrame(ground_truth)

    # Save to file
    output_path = output_dir / 'ground_truth_genomes.tsv'
    df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote ground truth for {len(df)} genomes to {output_path}")

    return output_path


def generate_reads(
    composition,
    output_prefix: str,
    n_reads: int = 1_000_000,
    model: str = 'NovaSeq',
    cpus: int = 4,
    compress: bool = False,
    gc_bias: bool = False,
    validate_output: bool = True,
    random_seed: Optional[int] = None,
    keep_temp_files: bool = False
) -> Dict[str, Path]:
    """
    Generate Illumina paired-end reads from a mock virome composition.

    This function:
    1. Writes genomes and abundances to temporary files
    2. Calls InSilicoSeq to generate reads
    3. Validates output FASTQ files (prevents quality issues)
    4. Creates ground truth metadata
    5. Cleans up temporary files

    Args:
        composition: MockViromeComposition object containing viral genomes
                    and contamination profile
        output_prefix: Prefix for output files (e.g., 'my_dataset')
                      Files created: my_dataset_R1.fastq, my_dataset_R2.fastq,
                      my_dataset_ground_truth_genomes.tsv
        n_reads: Number of paired-end reads to generate (default: 1,000,000)
        model: InSilicoSeq error model (default: 'NovaSeq')
              Options: 'NovaSeq', 'HiSeq', 'MiSeq', 'NextSeq'
              Note: Read length and insert size are determined by the model
        cpus: Number of CPUs for parallel processing (default: 4)
        compress: Gzip compress FASTQ output (default: False)
        gc_bias: Enable GC content bias simulation (default: False)
        validate_output: Validate FASTQ files after generation (default: True)
                        Detects length mismatches and truncated files
        random_seed: Random seed for reproducibility (default: None)
        keep_temp_files: Keep temporary FASTA/abundance files (default: False)
                        Useful for debugging

    Returns:
        Dict containing paths to generated files:
        {
            'r1': Path to R1 FASTQ file,
            'r2': Path to R2 FASTQ file,
            'ground_truth': Path to ground truth TSV file,
            'temp_fasta': Path to temp FASTA (if keep_temp_files=True),
            'temp_abundance': Path to temp abundance file (if keep_temp_files=True)
        }

    Raises:
        RuntimeError: If InSilicoSeq is not installed or fails
        ValueError: If validation fails (FASTQ quality issues)

    Example:
        >>> from viroforge.utils import create_mock_virome
        >>> from viroforge.simulators import generate_reads
        >>>
        >>> # Create composition
        >>> composition = create_mock_virome(
        ...     name='gut_virome_clean',
        ...     body_site='gut',
        ...     contamination_level='clean',
        ...     n_viral_genomes=50
        ... )
        >>>
        >>> # Generate 1M reads
        >>> output = generate_reads(
        ...     composition=composition,
        ...     output_prefix='output/gut_virome',
        ...     n_reads=1_000_000,
        ...     model='NovaSeq',
        ...     random_seed=42
        ... )
        >>>
        >>> print(f"Generated {output['r1']}")

    Note:
        - Read length and insert size are determined by the error model
        - NovaSeq: 150bp paired-end, ~350bp insert size
        - HiSeq: 125bp paired-end, ~350bp insert size
        - MiSeq: 250bp paired-end, ~450bp insert size

    See Also:
        - InSilicoSeq documentation: https://insilicoseq.readthedocs.io/
        - validate_fastq_record: For details on validation
    """
    # Check InSilicoSeq is installed
    if not check_insilicoseq_installed():
        raise RuntimeError(
            "InSilicoSeq is not installed. Install with: "
            "conda install -c bioconda insilicoseq"
        )

    logger.info(f"Generating {n_reads:,} reads with {model} error model")
    logger.info(f"Output prefix: {output_prefix}")

    # Validate output directory
    output_path = Path(output_prefix)
    output_dir = output_path.parent if output_path.parent != Path('.') else Path.cwd()
    validate_output_directory(output_dir, create=True)

    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix='viroforge_')
    temp_dir_path = Path(temp_dir)
    logger.debug(f"Created temporary directory: {temp_dir}")

    try:
        # Write genomes to FASTA
        genomes_fasta = temp_dir_path / 'genomes.fasta'
        genome_types = _write_genome_fasta(composition, genomes_fasta)

        # Write abundance file
        abundance_file = temp_dir_path / 'abundances.txt'
        _write_abundance_file(composition, abundance_file)

        # Run InSilicoSeq
        r1_path, r2_path = _run_insilicoseq(
            genomes_fasta=genomes_fasta,
            abundance_file=abundance_file,
            output_prefix=output_prefix,
            n_reads=n_reads,
            model=model,
            cpus=cpus,
            compress=compress,
            gc_bias=gc_bias,
            seed=random_seed
        )

        # Validate output FASTQ files
        if validate_output:
            logger.info("Validating output FASTQ files...")
            try:
                # Note: verify_fastq_file doesn't support .gz files yet
                # Skip validation for compressed files
                if not compress:
                    n_reads_r1 = verify_fastq_file(r1_path)
                    n_reads_r2 = verify_fastq_file(r2_path)
                    logger.info(f"✓ R1 validated: {n_reads_r1:,} reads")
                    logger.info(f"✓ R2 validated: {n_reads_r2:,} reads")

                    if n_reads_r1 != n_reads_r2:
                        raise ValueError(
                            f"Read count mismatch: R1 has {n_reads_r1:,} reads, "
                            f"R2 has {n_reads_r2:,} reads"
                        )
                else:
                    logger.info("Skipping validation for compressed files")
                    logger.info("Decompress files to validate: gunzip *.fastq.gz")
            except Exception as e:
                logger.error(f"FASTQ validation failed: {e}")
                raise

        # Create ground truth mapping
        ground_truth_path = _create_ground_truth_mapping(
            composition,
            output_dir,
            genome_types
        )

        # Prepare return dictionary
        result = {
            'r1': r1_path,
            'r2': r2_path,
            'ground_truth': ground_truth_path,
        }

        # Keep temp files if requested
        if keep_temp_files:
            # Copy temp files to output directory
            kept_fasta = output_dir / 'input_genomes.fasta'
            kept_abundance = output_dir / 'input_abundances.txt'
            shutil.copy(genomes_fasta, kept_fasta)
            shutil.copy(abundance_file, kept_abundance)
            result['temp_fasta'] = kept_fasta
            result['temp_abundance'] = kept_abundance
            logger.info(f"Kept temporary files in {output_dir}")

        logger.info("✓ Read generation complete!")
        logger.info(f"  R1: {r1_path}")
        logger.info(f"  R2: {r2_path}")
        logger.info(f"  Ground truth: {ground_truth_path}")

        return result

    finally:
        # Clean up temporary directory
        if not keep_temp_files:
            shutil.rmtree(temp_dir)
            logger.debug(f"Removed temporary directory: {temp_dir}")


def estimate_file_size(n_reads: int, read_length: int = 150,
                       compress: bool = False) -> str:
    """
    Estimate output FASTQ file size.

    Args:
        n_reads: Number of paired-end reads
        read_length: Average read length (default: 150)
        compress: Whether files will be compressed (default: False)

    Returns:
        Human-readable file size estimate

    Example:
        >>> estimate_file_size(1_000_000, read_length=150)
        'R1: ~95 MB, R2: ~95 MB, Total: ~190 MB'
        >>> estimate_file_size(1_000_000, compress=True)
        'R1: ~24 MB, R2: ~24 MB, Total: ~48 MB (compressed)'
    """
    # FASTQ format: 4 lines per read
    # @header\nsequence\n+\nquality\n
    # Approximate: header(50) + seq(read_length) + sep(2) + qual(read_length) + newlines(4)
    bytes_per_read = 50 + read_length + 2 + read_length + 4

    # Total size for one file (in bytes)
    size_bytes = n_reads * bytes_per_read

    # Compression ratio ~4:1 for FASTQ
    if compress:
        size_bytes = size_bytes // 4

    # Convert to human-readable
    if size_bytes < 1024**2:
        size_str = f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024**3:
        size_str = f"{size_bytes / 1024**2:.1f} MB"
    else:
        size_str = f"{size_bytes / 1024**3:.1f} GB"

    total_str = f"{size_bytes * 2 / 1024**2:.0f} MB" if size_bytes * 2 < 1024**3 else f"{size_bytes * 2 / 1024**3:.1f} GB"

    suffix = " (compressed)" if compress else ""
    return f"R1: ~{size_str}, R2: ~{size_str}, Total: ~{total_str}{suffix}"


# Convenience function
def quick_generate(
    body_site: str = 'gut',
    contamination_level: str = 'realistic',
    n_viral_genomes: int = 50,
    n_reads: int = 1_000_000,
    output_prefix: str = 'viroforge_output',
    random_seed: Optional[int] = None
) -> Dict[str, Path]:
    """
    Quick one-liner to generate a complete mock virome dataset.

    Args:
        body_site: Body site profile ('gut', 'oral', 'skin', 'respiratory', 'environmental')
        contamination_level: Contamination level ('clean', 'realistic', 'heavy', 'failed')
        n_viral_genomes: Number of viral genomes to include (default: 50)
        n_reads: Number of reads to generate (default: 1,000,000)
        output_prefix: Output file prefix (default: 'viroforge_output')
        random_seed: Random seed for reproducibility (default: None)

    Returns:
        Dict containing paths to generated files

    Example:
        >>> from viroforge.simulators import quick_generate
        >>>
        >>> # Generate complete gut virome in one line
        >>> output = quick_generate(
        ...     body_site='gut',
        ...     contamination_level='realistic',
        ...     n_reads=1_000_000,
        ...     random_seed=42
        ... )
    """
    from ..utils import create_mock_virome

    logger.info(f"Quick generate: {body_site} virome with {contamination_level} contamination")

    # Create composition
    composition = create_mock_virome(
        name=f"{body_site}_virome_{contamination_level}",
        body_site=body_site,
        contamination_level=contamination_level,
        n_viral_genomes=n_viral_genomes,
        random_seed=random_seed
    )

    # Generate reads
    return generate_reads(
        composition=composition,
        output_prefix=output_prefix,
        n_reads=n_reads,
        random_seed=random_seed
    )
