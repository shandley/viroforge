"""
Tests for genome database management utilities.

Tests downloading, validation, and caching of viral genomes from NCBI.
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viroforge.utils.genome_database import (
    get_database_info,
    validate_genome_record,
    extract_genome_metadata,
    set_entrez_email,
    download_genome_from_ncbi,
    get_genome_database,
    GENOMES_DIR,
    METADATA_DIR
)


class TestDatabaseInfo:
    """Test database information retrieval."""

    def test_get_minimal_database_info(self):
        """Test getting information about minimal database."""
        info = get_database_info('minimal')

        assert info['dataset'] == 'minimal'
        assert info['total_genomes'] == 20
        assert info['families'] == 6
        assert 'Microviridae' in info['family_counts']
        assert 'Siphoviridae' in info['family_counts']

    def test_get_full_database_info(self):
        """Test getting information about full database."""
        info = get_database_info('full')

        assert info['dataset'] == 'full'
        assert info['total_genomes'] > 100  # Should be ~105
        assert info['families'] > 5
        assert 'accessions' in info

    def test_invalid_dataset_raises_error(self):
        """Test that invalid dataset name raises error."""
        with pytest.raises(ValueError, match="Invalid dataset"):
            get_database_info('invalid_name')


class TestGenomeValidation:
    """Test genome validation functions."""

    def test_validate_valid_genome(self):
        """Test validation of a valid genome."""
        # Create a valid genome record
        record = SeqRecord(
            Seq("ATCGATCGATCG" * 100),  # 1200 bp
            id="TEST_001",
            description="Test genome"
        )
        record.annotations['organism'] = "Test virus"

        is_valid, errors = validate_genome_record(record)

        assert is_valid is True
        assert len(errors) == 0

    def test_validate_short_genome(self):
        """Test that short genomes are flagged."""
        record = SeqRecord(
            Seq("ATCG" * 100),  # 400 bp (too short)
            id="TEST_002"
        )
        record.annotations['organism'] = "Test virus"

        is_valid, errors = validate_genome_record(record, min_length=1000)

        assert is_valid is False
        assert any("too short" in err.lower() for err in errors)

    def test_validate_excessive_n_bases(self):
        """Test that excessive N bases are flagged."""
        record = SeqRecord(
            Seq("N" * 600 + "ATCG" * 100),  # 60% N bases
            id="TEST_003"
        )
        record.annotations['organism'] = "Test virus"

        is_valid, errors = validate_genome_record(record, max_n_percent=5.0)

        assert is_valid is False
        assert any("excessive n" in err.lower() for err in errors)

    def test_validate_extreme_gc_content(self):
        """Test that extreme GC content is flagged."""
        # Create very low GC sequence
        record = SeqRecord(
            Seq("A" * 1500),  # 0% GC
            id="TEST_004"
        )
        record.annotations['organism'] = "Test virus"

        is_valid, errors = validate_genome_record(record, min_gc=20.0)

        assert is_valid is False
        assert any("gc content" in err.lower() for err in errors)

    def test_validate_invalid_characters(self):
        """Test that invalid DNA characters are flagged."""
        record = SeqRecord(
            Seq("ATCGATCGXYZ" * 100),  # Invalid chars XYZ
            id="TEST_005"
        )
        record.annotations['organism'] = "Test virus"

        is_valid, errors = validate_genome_record(record)

        assert is_valid is False
        assert any("invalid" in err.lower() for err in errors)

    def test_validate_missing_organism(self):
        """Test that missing organism annotation is flagged."""
        record = SeqRecord(
            Seq("ATCGATCG" * 150),
            id="TEST_006"
        )
        # No organism annotation

        is_valid, errors = validate_genome_record(record)

        assert is_valid is False
        assert any("organism" in err.lower() for err in errors)


class TestMetadataExtraction:
    """Test metadata extraction from GenBank records."""

    def test_extract_basic_metadata(self):
        """Test extracting basic metadata from a record."""
        record = SeqRecord(
            Seq("ATCGATCG" * 625),  # 5000 bp, 50% GC
            id="NC_001422.1",
            description="Enterobacteria phage phiX174"
        )
        record.annotations['organism'] = "Escherichia virus phiX174"
        record.annotations['taxonomy'] = ['Viruses', 'Microviridae']
        record.annotations['source'] = "Enterobacteria phage"

        metadata = extract_genome_metadata(record)

        assert metadata['accession'] == "NC_001422.1"
        assert metadata['organism'] == "Escherichia virus phiX174"
        assert metadata['length'] == 5000
        assert 49.0 <= metadata['gc_content'] <= 51.0  # Should be ~50%
        assert 'Microviridae' in metadata['family']

    def test_extract_metadata_with_high_gc(self):
        """Test metadata extraction for high-GC genome."""
        record = SeqRecord(
            Seq("GCGCGC" * 250),  # 1500 bp, 100% GC
            id="TEST_GC",
            description="High GC test"
        )
        record.annotations['organism'] = "High GC virus"
        record.annotations['taxonomy'] = ['Viruses']

        metadata = extract_genome_metadata(record)

        assert metadata['gc_content'] == 100.0
        assert metadata['length'] == 1500

    def test_extract_metadata_unknown_taxonomy(self):
        """Test metadata extraction when taxonomy is missing."""
        record = SeqRecord(
            Seq("ATCG" * 300),
            id="TEST_UNKNOWN",
            description="Unknown taxonomy"
        )
        record.annotations['organism'] = "Unknown virus"
        # No taxonomy annotation

        metadata = extract_genome_metadata(record)

        assert metadata['family'] == 'Unknown'
        assert metadata['taxonomy'] == 'Unknown'


class TestEntrezEmail:
    """Test Entrez email setting."""

    def test_set_entrez_email(self):
        """Test setting Entrez email."""
        from Bio import Entrez

        original_email = Entrez.email
        test_email = "test@example.com"

        set_entrez_email(test_email)
        assert Entrez.email == test_email

        # Restore original
        Entrez.email = original_email


class TestGenomeDatabaseCache:
    """Test genome database caching behavior."""

    def test_database_directories_exist(self):
        """Test that database directories are created."""
        assert GENOMES_DIR.exists()
        assert METADATA_DIR.exists()
        assert GENOMES_DIR.is_dir()
        assert METADATA_DIR.is_dir()

    @pytest.mark.skipif(
        not (GENOMES_DIR / "minimal_test_set.fasta").exists(),
        reason="Minimal database not downloaded yet"
    )
    def test_load_cached_minimal_database(self):
        """Test loading minimal database from cache (if downloaded)."""
        genomes = get_genome_database('minimal')

        assert len(genomes) == 20
        assert all(hasattr(g, 'genome_id') for g in genomes)
        assert all(hasattr(g, 'sequence') for g in genomes)
        assert all(hasattr(g, 'family') for g in genomes)


class TestDownloadMocking:
    """Test download functionality with mocking (no actual NCBI calls)."""

    @patch('viroforge.utils.genome_database.Entrez.efetch')
    @patch('viroforge.utils.genome_database.SeqIO.read')
    def test_download_genome_success(self, mock_seqio, mock_efetch):
        """Test successful genome download (mocked)."""
        # Create mock record
        mock_record = SeqRecord(
            Seq("ATCGATCG" * 625),
            id="NC_001422.1",
            description="phiX174"
        )
        mock_record.annotations['organism'] = "phiX174"
        mock_record.annotations['taxonomy'] = ['Viruses', 'Microviridae']

        # Setup mocks
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_seqio.return_value = mock_record

        # Test download
        record = download_genome_from_ncbi("NC_001422.1")

        assert record.id == "NC_001422.1"
        assert len(record.seq) == 5000
        mock_efetch.assert_called_once()
        mock_handle.close.assert_called_once()

    @patch('viroforge.utils.genome_database.Entrez.efetch')
    def test_download_genome_retry_on_failure(self, mock_efetch):
        """Test that download retries on failure."""
        mock_efetch.side_effect = Exception("Network error")

        with pytest.raises(Exception):  # Should raise after retries
            download_genome_from_ncbi("NC_001422.1", retries=2, delay=0.1)

        # Should have tried multiple times
        assert mock_efetch.call_count == 2


def test_curated_genome_lists_loaded():
    """Test that curated genome lists are properly loaded."""
    from viroforge.data.curated_genomes import (
        MINIMAL_TEST_SET,
        FULL_PRODUCTION_SET,
        get_database_summary
    )

    # Check minimal set
    assert len(MINIMAL_TEST_SET) > 0
    assert 'Microviridae' in MINIMAL_TEST_SET
    assert len(MINIMAL_TEST_SET['Microviridae']) == 3

    # Check full set
    assert len(FULL_PRODUCTION_SET) > 0
    assert all(isinstance(accs, list) for accs in FULL_PRODUCTION_SET.values())

    # Check summary function
    summary = get_database_summary()
    assert summary['minimal_set']['total_genomes'] == 20
    assert summary['full_set']['total_genomes'] > 100


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
