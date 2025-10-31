"""
Unit tests for validation module.

Tests all validation functions to ensure they correctly identify valid and
invalid inputs, and raise appropriate errors.
"""

import pytest
import tempfile
from pathlib import Path

from viroforge.utils.validation import (
    validate_sequence,
    validate_sequence_length,
    validate_sequence_not_empty,
    validate_abundances,
    validate_abundance_range,
    validate_fastq_record,
    verify_fastq_file,
    validate_gc_content,
    validate_file_not_empty,
    validate_output_directory,
    validate_genome_collection,
)
from viroforge.core.community import ViralGenome
from viroforge.core.contamination import ContaminantGenome, ContaminantType


class TestSequenceValidation:
    """Test sequence validation functions."""

    def test_validate_sequence_valid(self):
        """Test validation of valid DNA sequences."""
        assert validate_sequence("ATCG")
        assert validate_sequence("ATCGN")
        assert validate_sequence("atcg")  # Lowercase
        assert validate_sequence("AAAA")  # Homopolymer
        assert validate_sequence("N" * 100)  # All N

    def test_validate_sequence_invalid_characters(self):
        """Test rejection of invalid characters."""
        with pytest.raises(ValueError, match="Invalid DNA characters"):
            validate_sequence("ATCGX")
        with pytest.raises(ValueError, match="Invalid DNA characters"):
            validate_sequence("ATCG123")
        with pytest.raises(ValueError, match="Invalid DNA characters"):
            validate_sequence("ATCG-N")

    def test_validate_sequence_no_n(self):
        """Test validation with N not allowed."""
        assert validate_sequence("ATCG", allow_n=False)
        with pytest.raises(ValueError, match="Invalid DNA characters"):
            validate_sequence("ATCGN", allow_n=False)

    def test_validate_sequence_not_empty(self):
        """Test empty sequence detection."""
        with pytest.raises(ValueError, match="empty"):
            validate_sequence_not_empty("")
        assert validate_sequence_not_empty("A")
        assert validate_sequence_not_empty("ATCG")


class TestLengthValidation:
    """Test length validation functions."""

    def test_validate_sequence_length_matching(self):
        """Test validation when lengths match."""
        genome = ViralGenome(
            genome_id="test",
            sequence="ATCG",
            taxonomy="Test"
        )
        assert validate_sequence_length(genome)

    def test_validate_sequence_length_mismatch(self):
        """Test detection of length mismatches."""
        genome = ViralGenome(
            genome_id="test",
            sequence="ATCG",
            taxonomy="Test"
        )
        # Manually break the length
        genome.length = 10
        with pytest.raises(ValueError, match="Length mismatch"):
            validate_sequence_length(genome)


class TestAbundanceValidation:
    """Test abundance validation functions."""

    def test_validate_abundances_correct_sum(self):
        """Test validation when abundances sum to 1.0."""
        class MockGenome:
            def __init__(self, abundance):
                self.abundance = abundance

        genomes = [MockGenome(0.5), MockGenome(0.5)]
        assert validate_abundances(genomes, warn_only=False)

    def test_validate_abundances_floating_point_precision(self):
        """Test tolerance for floating point precision errors."""
        class MockGenome:
            def __init__(self, abundance):
                self.abundance = abundance

        # Sum = 0.99999999 (close to 1.0)
        genomes = [MockGenome(0.3333333), MockGenome(0.3333333), MockGenome(0.3333333)]
        total = sum(g.abundance for g in genomes)
        assert abs(total - 1.0) < 0.01  # Not exactly 1.0
        assert validate_abundances(genomes, tolerance=1e-2, warn_only=False)

    def test_validate_abundances_incorrect_sum(self):
        """Test detection when abundances don't sum to 1.0."""
        class MockGenome:
            def __init__(self, abundance):
                self.abundance = abundance

        genomes = [MockGenome(0.5), MockGenome(0.4)]
        assert not validate_abundances(genomes, warn_only=True)
        with pytest.raises(ValueError, match="sum to"):
            validate_abundances(genomes, warn_only=False)

    def test_validate_abundance_range_valid(self):
        """Test validation of abundance values in range."""
        assert validate_abundance_range(0.0)
        assert validate_abundance_range(0.5)
        assert validate_abundance_range(1.0)

    def test_validate_abundance_range_invalid(self):
        """Test detection of out-of-range abundances."""
        with pytest.raises(ValueError, match="must be between"):
            validate_abundance_range(-0.1)
        with pytest.raises(ValueError, match="must be between"):
            validate_abundance_range(1.1)


class TestFASTQValidation:
    """Test FASTQ validation functions."""

    def test_validate_fastq_record_valid(self):
        """Test validation of valid FASTQ records."""
        assert validate_fastq_record("read1", "ATCG", "IIII")
        assert validate_fastq_record("read2", "ATCGN", "IIIII")
        assert validate_fastq_record("read3", "A" * 150, "I" * 150)

    def test_validate_fastq_record_length_mismatch(self):
        """Test detection of sequence/quality length mismatch."""
        with pytest.raises(ValueError, match="LENGTH MISMATCH"):
            validate_fastq_record("read1", "ATCG", "III")
        with pytest.raises(ValueError, match="LENGTH MISMATCH"):
            validate_fastq_record("read1", "ATCG", "IIIII")

    def test_validate_fastq_record_invalid_sequence(self):
        """Test detection of invalid sequence characters."""
        with pytest.raises(ValueError, match="Invalid"):
            validate_fastq_record("read1", "ATCGX", "IIIII")

    def test_validate_fastq_record_invalid_quality(self):
        """Test detection of invalid quality scores."""
        # Null character (ASCII 0)
        with pytest.raises(ValueError, match="Invalid quality score"):
            validate_fastq_record("read1", "ATCG", "III\x00")
        # Non-printable character
        with pytest.raises(ValueError, match="Invalid quality score"):
            validate_fastq_record("read1", "ATCG", "III\x1f")

    def test_validate_fastq_record_empty(self):
        """Test detection of empty records."""
        with pytest.raises(ValueError, match="Empty"):
            validate_fastq_record("read1", "", "")

    def test_verify_fastq_file_valid(self):
        """Test verification of valid FASTQ file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            temp_path = f.name
            # Write valid FASTQ records
            f.write("@read1\n")
            f.write("ATCG\n")
            f.write("+\n")
            f.write("IIII\n")
            f.write("@read2\n")
            f.write("GCTA\n")
            f.write("+\n")
            f.write("JJJJ\n")

        try:
            n_records = verify_fastq_file(temp_path)
            assert n_records == 2
        finally:
            Path(temp_path).unlink()

    def test_verify_fastq_file_length_mismatch(self):
        """Test detection of length mismatch in file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            temp_path = f.name
            # Write invalid FASTQ record (length mismatch)
            f.write("@read1\n")
            f.write("ATCG\n")
            f.write("+\n")
            f.write("III\n")  # Only 3 quality scores for 4 bases

        try:
            with pytest.raises(ValueError, match="Length mismatch"):
                verify_fastq_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_verify_fastq_file_truncated(self):
        """Test detection of truncated file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            temp_path = f.name
            # Write incomplete FASTQ record (missing quality line)
            f.write("@read1\n")
            f.write("ATCG\n")
            f.write("+\n")
            # Missing quality line!

        try:
            with pytest.raises(ValueError, match="truncated"):
                verify_fastq_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_verify_fastq_file_malformed_header(self):
        """Test detection of malformed headers."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            temp_path = f.name
            # Write record with bad header
            f.write("read1\n")  # Missing @ symbol
            f.write("ATCG\n")
            f.write("+\n")
            f.write("IIII\n")

        try:
            with pytest.raises(ValueError, match="Invalid header"):
                verify_fastq_file(temp_path)
        finally:
            Path(temp_path).unlink()


class TestGCContentValidation:
    """Test GC content validation."""

    def test_validate_gc_content_matching(self):
        """Test validation when GC content matches."""
        genome = ViralGenome(
            genome_id="test",
            sequence="ATCG",  # 50% GC
            taxonomy="Test"
        )
        # GC content is calculated automatically
        assert validate_gc_content(genome, tolerance=5.0, warn_only=False)

    def test_validate_gc_content_mismatch(self):
        """Test detection when GC content doesn't match."""
        genome = ViralGenome(
            genome_id="test",
            sequence="AAAA",  # 0% GC
            taxonomy="Test"
        )
        # Manually set wrong GC content
        genome.gc_content = 50.0

        # Should fail with strict checking
        with pytest.raises(ValueError, match="GC content mismatch"):
            validate_gc_content(genome, tolerance=5.0, warn_only=False)

        # Should warn with lenient checking
        assert not validate_gc_content(genome, tolerance=5.0, warn_only=True)


class TestFileValidation:
    """Test file validation functions."""

    def test_validate_file_not_empty_valid(self):
        """Test validation of non-empty file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            temp_path = f.name
            f.write("content")

        try:
            assert validate_file_not_empty(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_validate_file_not_empty_missing(self):
        """Test detection of missing file."""
        with pytest.raises(FileNotFoundError):
            validate_file_not_empty("/tmp/nonexistent_file_12345.txt")

    def test_validate_file_not_empty_empty(self):
        """Test detection of empty file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            temp_path = f.name
            # Don't write anything

        try:
            with pytest.raises(ValueError, match="empty"):
                validate_file_not_empty(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_validate_output_directory_exists(self):
        """Test validation of existing directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            assert validate_output_directory(temp_dir)

    def test_validate_output_directory_create(self):
        """Test creation of directory when it doesn't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            new_dir = Path(temp_dir) / "subdir"
            assert not new_dir.exists()
            assert validate_output_directory(new_dir, create=True)
            assert new_dir.exists()

    def test_validate_output_directory_not_exists(self):
        """Test detection of non-existent directory without create."""
        with pytest.raises(FileNotFoundError):
            validate_output_directory("/tmp/nonexistent_dir_12345", create=False)


class TestBatchValidation:
    """Test batch validation of genome collections."""

    def test_validate_genome_collection_valid(self):
        """Test validation of valid genome collection."""
        genomes = [
            ViralGenome(
                genome_id=f"genome{i}",
                sequence="ATCG" * 25,  # 100 bp
                taxonomy=f"Family{i};Genus{i};Species{i}",
                abundance=0.5
            )
            for i in range(2)
        ]

        assert validate_genome_collection(
            genomes,
            check_sequences=True,
            check_lengths=True,
            check_abundances=True,
            check_gc=True
        )

    def test_validate_genome_collection_invalid_sequence(self):
        """Test detection of invalid sequences in collection."""
        genome1 = ViralGenome(
            genome_id="genome1",
            sequence="ATCG",
            taxonomy="Test",
            abundance=0.5
        )
        genome2 = ViralGenome(
            genome_id="genome2",
            sequence="ATCGX",  # Invalid
            taxonomy="Test",
            abundance=0.5
        )

        with pytest.raises(ValueError, match="Invalid DNA characters"):
            validate_genome_collection([genome1, genome2], check_sequences=True)

    def test_validate_genome_collection_empty(self):
        """Test validation of empty collection."""
        assert validate_genome_collection([])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
