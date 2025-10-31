"""
Quality control and validation utilities for ViroForge.

This module provides validation functions to ensure data integrity throughout
the simulation pipeline. These validators help prevent common errors like:
- FASTQ length mismatches (sequence length ≠ quality score length)
- Invalid DNA characters
- File truncation/corruption
- Abundance calculation errors
- Inconsistent metadata

All validators can be used in two modes:
- Strict mode (raises ValueError on failure)
- Warning mode (logs warning and continues)

Example:
    >>> from viroforge.utils.validation import validate_sequence
    >>> validate_sequence("ATCGN")  # Valid
    True
    >>> validate_sequence("ATCGX")  # Invalid
    ValueError: Invalid DNA characters found: {'X'}
"""

from typing import Union, List, Optional, Any
from pathlib import Path
import logging

from Bio.Seq import Seq

logger = logging.getLogger(__name__)


# ============================================================================
# Sequence Validation
# ============================================================================

def validate_sequence(
    seq: Union[str, Seq],
    allow_n: bool = True,
    allow_lowercase: bool = True
) -> bool:
    """
    Validate DNA sequence contains only valid characters.

    Args:
        seq: DNA sequence to validate
        allow_n: Whether to allow 'N' characters (default True)
        allow_lowercase: Whether to allow lowercase letters (default True)

    Returns:
        True if valid

    Raises:
        ValueError: If invalid characters found

    Examples:
        >>> validate_sequence("ATCG")
        True
        >>> validate_sequence("ATCGN")
        True
        >>> validate_sequence("atcg")
        True
        >>> validate_sequence("ATCGX")
        ValueError: Invalid DNA characters found: {'X'}
    """
    seq_str = str(seq)

    # Convert to uppercase for checking if lowercase allowed
    if allow_lowercase:
        seq_str = seq_str.upper()

    # Define valid character set
    valid_chars = set('ATCGN') if allow_n else set('ATCG')

    # Find invalid characters
    invalid = set(seq_str) - valid_chars

    if invalid:
        raise ValueError(
            f"Invalid DNA characters found: {invalid}\n"
            f"Valid characters: {valid_chars}"
        )

    return True


def validate_sequence_length(
    genome: Any,
    attribute_name: str = 'length',
    genome_id_attr: str = 'genome_id'
) -> bool:
    """
    Validate that genome.length matches len(genome.sequence).

    This prevents inconsistencies between stored length metadata and
    actual sequence length.

    Args:
        genome: Genome object with .sequence and .length attributes
        attribute_name: Name of length attribute (default 'length')
        genome_id_attr: Name of ID attribute for error messages

    Returns:
        True if lengths match

    Raises:
        ValueError: If lengths don't match

    Examples:
        >>> from viroforge.core.community import ViralGenome
        >>> genome = ViralGenome(genome_id="test", sequence="ATCG", taxonomy="Test")
        >>> validate_sequence_length(genome)
        True
    """
    expected_length = getattr(genome, attribute_name)
    actual_length = len(genome.sequence)
    genome_id = getattr(genome, genome_id_attr, "unknown")

    if expected_length != actual_length:
        raise ValueError(
            f"Length mismatch for {genome_id}: "
            f"{attribute_name}={expected_length} but "
            f"sequence is {actual_length} bp"
        )

    return True


def validate_sequence_not_empty(
    seq: Union[str, Seq],
    name: str = "sequence"
) -> bool:
    """
    Validate that sequence is not empty.

    Args:
        seq: Sequence to validate
        name: Name for error message

    Returns:
        True if not empty

    Raises:
        ValueError: If sequence is empty
    """
    if len(seq) == 0:
        raise ValueError(f"{name} is empty (length 0)")

    return True


# ============================================================================
# Abundance Validation
# ============================================================================

def validate_abundances(
    objects: List[Any],
    tolerance: float = 1e-6,
    warn_only: bool = True,
    abundance_attr: str = 'abundance'
) -> bool:
    """
    Validate that abundances sum to approximately 1.0.

    Due to floating point precision, abundances may not sum to exactly 1.0.
    This validator checks that the sum is within a tolerance.

    Args:
        objects: List of objects with .abundance attribute
        tolerance: Maximum allowed deviation from 1.0 (default 1e-6)
        warn_only: If True, warn instead of raising error (default True)
        abundance_attr: Name of abundance attribute (default 'abundance')

    Returns:
        True if within tolerance, False otherwise

    Raises:
        ValueError: If outside tolerance and warn_only=False

    Examples:
        >>> class MockGenome:
        ...     def __init__(self, abundance):
        ...         self.abundance = abundance
        >>> genomes = [MockGenome(0.5), MockGenome(0.5)]
        >>> validate_abundances(genomes)
        True
        >>> genomes = [MockGenome(0.5), MockGenome(0.4)]
        >>> validate_abundances(genomes, warn_only=True)
        False  # Logs warning
    """
    total = sum(getattr(obj, abundance_attr) for obj in objects)
    deviation = abs(total - 1.0)

    if deviation > tolerance:
        msg = (
            f"Abundances sum to {total:.10f}, not 1.0 "
            f"(deviation: {deviation:.2e}, tolerance: {tolerance:.2e})"
        )
        if warn_only:
            logger.warning(msg)
        else:
            raise ValueError(msg)
        return False

    return True


def validate_abundance_range(
    abundance: float,
    min_val: float = 0.0,
    max_val: float = 1.0,
    name: str = "abundance"
) -> bool:
    """
    Validate that abundance is within valid range.

    Args:
        abundance: Abundance value to check
        min_val: Minimum valid value (default 0.0)
        max_val: Maximum valid value (default 1.0)
        name: Name for error message

    Returns:
        True if in range

    Raises:
        ValueError: If outside range
    """
    if not (min_val <= abundance <= max_val):
        raise ValueError(
            f"{name} ({abundance}) must be between {min_val} and {max_val}"
        )

    return True


# ============================================================================
# FASTQ Validation
# ============================================================================

def validate_fastq_record(
    read_id: str,
    sequence: str,
    quality: str,
    allow_n: bool = True
) -> bool:
    """
    Validate a FASTQ record before writing.

    **CRITICAL:** This MUST be called before writing every FASTQ record
    to prevent the most common FASTQ errors:
    - Length mismatch between sequence and quality scores
    - Invalid DNA characters
    - Invalid quality score characters

    Args:
        read_id: Read identifier (for error messages)
        sequence: DNA sequence string
        quality: Quality score string (Phred33 encoding)
        allow_n: Whether to allow 'N' in sequence (default True)

    Returns:
        True if valid

    Raises:
        ValueError: If validation fails with detailed error message

    Examples:
        >>> validate_fastq_record("read1", "ATCG", "IIII")
        True
        >>> validate_fastq_record("read1", "ATCG", "III")
        ValueError: FASTQ LENGTH MISMATCH for read1...
        >>> validate_fastq_record("read1", "ATCGX", "IIIII")
        ValueError: Invalid DNA characters found...
    """
    # Length check (MOST CRITICAL)
    if len(sequence) != len(quality):
        raise ValueError(
            f"FASTQ LENGTH MISMATCH for {read_id}:\n"
            f"  Sequence: {len(sequence)} bp\n"
            f"  Quality:  {len(quality)} chars\n"
            f"These MUST be equal!\n"
            f"This is the most common FASTQ error and will cause downstream failures."
        )

    # Empty check
    if len(sequence) == 0:
        raise ValueError(f"Empty FASTQ record for {read_id}")

    # Sequence character validation
    try:
        validate_sequence(sequence, allow_n=allow_n, allow_lowercase=False)
    except ValueError as e:
        raise ValueError(f"Invalid sequence in {read_id}: {e}")

    # Quality score validation (Phred33: ASCII 33-126)
    for i, char in enumerate(quality):
        ascii_val = ord(char)
        if not (33 <= ascii_val <= 126):
            raise ValueError(
                f"Invalid quality score at position {i} in {read_id}:\n"
                f"  Character: '{char}' (ASCII {ascii_val})\n"
                f"  Valid range: ASCII 33-126 (Phred33 encoding)\n"
                f"  Common quality chars: !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
            )

    return True


def verify_fastq_file(
    file_path: Union[str, Path],
    max_records_to_check: Optional[int] = None,
    verbose: bool = False
) -> int:
    """
    Verify FASTQ file integrity after writing.

    This function checks that a FASTQ file is properly formatted and
    not truncated. It's recommended to call this after writing large files.

    Checks:
    - Proper 4-line format (@header, sequence, +, quality)
    - Sequence and quality lengths match
    - No truncation (file ends properly)
    - No empty records

    Args:
        file_path: Path to FASTQ file to verify
        max_records_to_check: Optional limit on records to check (for large files)
        verbose: If True, log progress periodically

    Returns:
        Number of valid records in file

    Raises:
        ValueError: If file is malformed with specific error details
        FileNotFoundError: If file doesn't exist

    Examples:
        >>> # After writing FASTQ file
        >>> n_reads = verify_fastq_file("output.fastq")
        >>> print(f"Verified {n_reads} reads")
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {file_path}")

    n_records = 0

    with open(file_path) as f:
        while True:
            # Check if we've hit the limit
            if max_records_to_check and n_records >= max_records_to_check:
                logger.info(f"Verified {n_records} records (limit reached)")
                break

            # Read 4 lines for one FASTQ record
            header = f.readline()

            # Check for EOF
            if not header:
                break  # Normal end of file

            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            # Check for truncation (incomplete record at EOF)
            if not qual:
                raise ValueError(
                    f"File appears truncated at record {n_records + 1}:\n"
                    f"  Found header but missing quality line\n"
                    f"  This usually means the file was not written completely"
                )

            # Validate format
            if not header.startswith('@'):
                raise ValueError(
                    f"Invalid header at record {n_records + 1}:\n"
                    f"  Expected line starting with '@', got: {header[:50]}\n"
                    f"  This may indicate file corruption"
                )

            if not plus.startswith('+'):
                raise ValueError(
                    f"Invalid separator at record {n_records + 1}:\n"
                    f"  Expected line starting with '+', got: {plus[:50]}\n"
                    f"  This may indicate file corruption"
                )

            # Validate lengths match
            seq = seq.strip()
            qual = qual.strip()

            if len(seq) != len(qual):
                raise ValueError(
                    f"Length mismatch at record {n_records + 1}:\n"
                    f"  Sequence: {len(seq)} bp\n"
                    f"  Quality:  {len(qual)} chars\n"
                    f"  Header: {header.strip()}"
                )

            n_records += 1

            # Verbose progress
            if verbose and n_records % 100000 == 0:
                logger.info(f"Verified {n_records:,} records...")

    logger.info(f"Verified {n_records:,} FASTQ records in {file_path}")

    return n_records


# ============================================================================
# GC Content Validation
# ============================================================================

def validate_gc_content(
    genome: Any,
    tolerance: float = 5.0,
    warn_only: bool = True,
    sequence_attr: str = 'sequence',
    gc_attr: str = 'gc_content',
    id_attr: str = 'genome_id'
) -> bool:
    """
    Verify that actual GC content matches expected value.

    This helps catch errors in sequence generation or metadata.

    Args:
        genome: Genome object with sequence and gc_content attributes
        tolerance: Maximum allowed deviation in percentage points (default 5.0)
        warn_only: If True, warn instead of raising error (default True)
        sequence_attr: Name of sequence attribute
        gc_attr: Name of GC content attribute
        id_attr: Name of ID attribute

    Returns:
        True if within tolerance, False otherwise

    Raises:
        ValueError: If outside tolerance and warn_only=False
    """
    seq = getattr(genome, sequence_attr)
    expected_gc = getattr(genome, gc_attr)
    genome_id = getattr(genome, id_attr, "unknown")

    # Calculate actual GC content
    seq_str = str(seq).upper()
    if len(seq_str) == 0:
        actual_gc = 0.0
    else:
        g_count = seq_str.count('G')
        c_count = seq_str.count('C')
        actual_gc = ((g_count + c_count) / len(seq_str)) * 100

    deviation = abs(actual_gc - expected_gc)

    if deviation > tolerance:
        msg = (
            f"GC content mismatch for {genome_id}:\n"
            f"  Expected: {expected_gc:.2f}%\n"
            f"  Actual:   {actual_gc:.2f}%\n"
            f"  Deviation: {deviation:.2f}% (tolerance: {tolerance:.2f}%)"
        )
        if warn_only:
            logger.warning(msg)
        else:
            raise ValueError(msg)
        return False

    return True


# ============================================================================
# File Validation
# ============================================================================

def validate_file_not_empty(file_path: Union[str, Path]) -> bool:
    """
    Validate that file exists and is not empty.

    Args:
        file_path: Path to file

    Returns:
        True if file exists and has content

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is empty
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    if file_path.stat().st_size == 0:
        raise ValueError(f"File is empty: {file_path}")

    return True


def validate_output_directory(dir_path: Union[str, Path], create: bool = False) -> bool:
    """
    Validate that output directory exists and is writable.

    Args:
        dir_path: Path to directory
        create: If True, create directory if it doesn't exist

    Returns:
        True if directory is valid

    Raises:
        FileNotFoundError: If directory doesn't exist and create=False
        PermissionError: If directory is not writable
    """
    dir_path = Path(dir_path)

    if not dir_path.exists():
        if create:
            dir_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created output directory: {dir_path}")
        else:
            raise FileNotFoundError(f"Directory not found: {dir_path}")

    if not dir_path.is_dir():
        raise ValueError(f"Not a directory: {dir_path}")

    # Check if writable by trying to create a temp file
    import tempfile
    try:
        with tempfile.NamedTemporaryFile(dir=dir_path, delete=True):
            pass
    except PermissionError:
        raise PermissionError(f"Directory is not writable: {dir_path}")

    return True


# ============================================================================
# Batch Validation
# ============================================================================

def validate_genome_collection(
    genomes: List[Any],
    check_sequences: bool = True,
    check_lengths: bool = True,
    check_abundances: bool = True,
    check_gc: bool = False,
    abundance_tolerance: float = 1e-6,
    gc_tolerance: float = 5.0
) -> bool:
    """
    Validate a collection of genomes (viral or contaminant).

    This is a convenience function that runs multiple validators on
    a collection of genomes.

    Args:
        genomes: List of genome objects
        check_sequences: Validate sequence characters
        check_lengths: Validate length consistency
        check_abundances: Validate abundance sum
        check_gc: Validate GC content
        abundance_tolerance: Tolerance for abundance sum
        gc_tolerance: Tolerance for GC content

    Returns:
        True if all checks pass

    Raises:
        ValueError: If any validation fails
    """
    if not genomes:
        logger.warning("Empty genome collection")
        return True

    # Check sequences
    if check_sequences:
        for genome in genomes:
            try:
                validate_sequence(genome.sequence)
            except ValueError as e:
                raise ValueError(f"Genome {genome.genome_id}: {e}")

    # Check lengths
    if check_lengths:
        for genome in genomes:
            validate_sequence_length(genome)

    # Check abundances
    if check_abundances:
        validate_abundances(genomes, tolerance=abundance_tolerance, warn_only=False)

    # Check GC content
    if check_gc:
        for genome in genomes:
            validate_gc_content(genome, tolerance=gc_tolerance, warn_only=False)

    logger.info(f"Validated {len(genomes)} genomes successfully")

    return True


if __name__ == "__main__":
    # Example usage and testing
    print("ViroForge Validation Module - Example Usage\n")
    print("=" * 60)

    # Example 1: Sequence validation
    print("Example 1: Sequence Validation")
    try:
        validate_sequence("ATCGN")
        print("  ✅ Valid DNA sequence")
    except ValueError as e:
        print(f"  ❌ {e}")

    try:
        validate_sequence("ATCGX")
        print("  ✅ Valid DNA sequence")
    except ValueError as e:
        print(f"  ❌ Invalid: {e}")

    # Example 2: FASTQ validation
    print("\nExample 2: FASTQ Record Validation")
    try:
        validate_fastq_record("read1", "ATCG", "IIII")
        print("  ✅ Valid FASTQ record")
    except ValueError as e:
        print(f"  ❌ {e}")

    try:
        validate_fastq_record("read2", "ATCG", "III")  # Length mismatch
        print("  ✅ Valid FASTQ record")
    except ValueError as e:
        print(f"  ❌ Length mismatch detected: {str(e)[:80]}...")

    # Example 3: Abundance validation
    print("\nExample 3: Abundance Validation")

    class MockGenome:
        def __init__(self, abundance):
            self.abundance = abundance

    genomes = [MockGenome(0.5), MockGenome(0.5)]
    if validate_abundances(genomes, warn_only=True):
        print("  ✅ Abundances sum correctly")

    genomes_bad = [MockGenome(0.5), MockGenome(0.4)]
    if not validate_abundances(genomes_bad, warn_only=True):
        print("  ⚠️  Abundances don't sum to 1.0 (warning logged)")

    print("\n" + "=" * 60)
    print("Validation module ready for use!")
