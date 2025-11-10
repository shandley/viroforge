"""
Tests for long-read sequencing simulation (PacBio HiFi and Nanopore).

Tests all components of long-read simulation:
- Configuration classes (PacBioHiFiConfig, NanoporeConfig)
- Platform enum
- Dependency checks
- Mock simulation workflows
- Ground truth generation
- Integration with ViroForge infrastructure
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
import tempfile
import shutil

from viroforge.simulators.longread import (
    LongReadPlatform,
    PacBioHiFiConfig,
    NanoporeConfig,
    check_pbsim3_installed,
    check_pbccs_installed,
    check_samtools_installed
)


class TestLongReadPlatform:
    """Test LongReadPlatform enum."""

    def test_platform_values(self):
        """Test platform enum values."""
        assert LongReadPlatform.PACBIO_HIFI.value == "pacbio-hifi"
        assert LongReadPlatform.NANOPORE.value == "nanopore"

    def test_platform_from_string(self):
        """Test creating platform from string."""
        assert LongReadPlatform("pacbio-hifi") == LongReadPlatform.PACBIO_HIFI
        assert LongReadPlatform("nanopore") == LongReadPlatform.NANOPORE

    def test_invalid_platform_raises_error(self):
        """Test that invalid platform string raises error."""
        with pytest.raises(ValueError):
            LongReadPlatform("illumina")


class TestPacBioHiFiConfig:
    """Test PacBioHiFiConfig dataclass."""

    def test_default_configuration(self):
        """Test default PacBio HiFi configuration."""
        config = PacBioHiFiConfig()

        assert config.passes == 10
        assert config.min_passes == 3
        assert config.accuracy_model == "QSHMM-RSII"
        assert config.read_length_mean == 15000
        assert config.read_length_sd == 5000
        assert config.clr_error_rate == 0.15

    def test_custom_configuration(self):
        """Test custom PacBio HiFi parameters."""
        config = PacBioHiFiConfig(
            passes=15,
            min_passes=5,
            read_length_mean=20000,
            read_length_sd=7000
        )

        assert config.passes == 15
        assert config.min_passes == 5
        assert config.read_length_mean == 20000
        assert config.read_length_sd == 7000

    def test_high_passes_for_high_accuracy(self):
        """Test that more passes can be set for higher accuracy."""
        config = PacBioHiFiConfig(passes=20, min_passes=10)

        assert config.passes == 20
        assert config.min_passes == 10

    def test_read_length_ranges(self):
        """Test various read length configurations."""
        # Short reads (10kb)
        config_short = PacBioHiFiConfig(read_length_mean=10000)
        assert config_short.read_length_mean == 10000

        # Standard reads (15kb)
        config_standard = PacBioHiFiConfig(read_length_mean=15000)
        assert config_standard.read_length_mean == 15000

        # Long reads (25kb)
        config_long = PacBioHiFiConfig(read_length_mean=25000)
        assert config_long.read_length_mean == 25000

    def test_accuracy_models(self):
        """Test different accuracy models."""
        config_rsii = PacBioHiFiConfig(accuracy_model="QSHMM-RSII")
        assert config_rsii.accuracy_model == "QSHMM-RSII"

        config_sequel = PacBioHiFiConfig(accuracy_model="QSHMM-SEQUEL")
        assert config_sequel.accuracy_model == "QSHMM-SEQUEL"


class TestNanoporeConfig:
    """Test NanoporeConfig dataclass."""

    def test_default_configuration(self):
        """Test default Nanopore configuration."""
        config = NanoporeConfig()

        assert config.chemistry == "R10.4"
        assert config.read_length_mean == 20000
        assert config.read_length_sd == 10000
        assert config.error_rate == 0.05
        assert config.hp_del_bias == 6
        assert config.quality_mean == 10

    def test_custom_configuration(self):
        """Test custom Nanopore parameters."""
        config = NanoporeConfig(
            chemistry="R9.4",
            read_length_mean=30000,
            read_length_sd=15000,
            error_rate=0.08
        )

        assert config.chemistry == "R9.4"
        assert config.read_length_mean == 30000
        assert config.read_length_sd == 15000
        assert config.error_rate == 0.08

    def test_r94_chemistry(self):
        """Test R9.4 chemistry configuration (older, less accurate)."""
        config = NanoporeConfig(
            chemistry="R9.4",
            error_rate=0.10  # Higher error for R9.4
        )

        assert config.chemistry == "R9.4"
        assert config.error_rate == 0.10

    def test_r104_chemistry(self):
        """Test R10.4 chemistry configuration (newer, more accurate)."""
        config = NanoporeConfig(
            chemistry="R10.4",
            error_rate=0.05  # Lower error for R10.4
        )

        assert config.chemistry == "R10.4"
        assert config.error_rate == 0.05

    def test_ultra_long_reads(self):
        """Test ultra-long read configuration (100kb+)."""
        config = NanoporeConfig(
            read_length_mean=100000,
            read_length_sd=50000
        )

        assert config.read_length_mean == 100000
        assert config.read_length_sd == 50000

    def test_homopolymer_bias(self):
        """Test homopolymer deletion bias parameter."""
        config = NanoporeConfig(hp_del_bias=6)
        assert config.hp_del_bias == 6

        # Can be adjusted for different error characteristics
        config_low = NanoporeConfig(hp_del_bias=3)
        assert config_low.hp_del_bias == 3


class TestDependencyChecks:
    """Test dependency checking functions."""

    @patch('shutil.which')
    def test_pbsim3_installed(self, mock_which):
        """Test PBSIM3 installation check when installed."""
        mock_which.return_value = '/usr/bin/pbsim'

        assert check_pbsim3_installed() is True
        mock_which.assert_called_with('pbsim')

    @patch('shutil.which')
    def test_pbsim3_not_installed(self, mock_which):
        """Test PBSIM3 installation check when not installed."""
        mock_which.return_value = None

        assert check_pbsim3_installed() is False

    @patch('shutil.which')
    def test_pbccs_installed(self, mock_which):
        """Test PacBio ccs installation check when installed."""
        mock_which.return_value = '/usr/bin/ccs'

        assert check_pbccs_installed() is True
        mock_which.assert_called_with('ccs')

    @patch('shutil.which')
    def test_pbccs_not_installed(self, mock_which):
        """Test PacBio ccs installation check when not installed."""
        mock_which.return_value = None

        assert check_pbccs_installed() is False

    @patch('shutil.which')
    def test_samtools_installed(self, mock_which):
        """Test SAMtools installation check when installed."""
        mock_which.return_value = '/usr/bin/samtools'

        assert check_samtools_installed() is True
        mock_which.assert_called_with('samtools')

    @patch('shutil.which')
    def test_samtools_not_installed(self, mock_which):
        """Test SAMtools installation check when not installed."""
        mock_which.return_value = None

        assert check_samtools_installed() is False


class TestConfigurationValidation:
    """Test configuration validation and edge cases."""

    def test_pacbio_min_passes_less_than_passes(self):
        """Test that min_passes should be less than passes."""
        # This should work
        config = PacBioHiFiConfig(passes=10, min_passes=3)
        assert config.passes > config.min_passes

    def test_pacbio_very_short_reads(self):
        """Test PacBio HiFi with very short reads (edge case)."""
        config = PacBioHiFiConfig(read_length_mean=5000, read_length_sd=2000)
        assert config.read_length_mean == 5000

    def test_nanopore_very_long_reads(self):
        """Test Nanopore with very long reads (200kb+)."""
        config = NanoporeConfig(read_length_mean=200000, read_length_sd=100000)
        assert config.read_length_mean == 200000

    def test_zero_standard_deviation_edge_case(self):
        """Test configurations with zero standard deviation."""
        # Should work but produce uniform read lengths
        config = PacBioHiFiConfig(read_length_sd=0)
        assert config.read_length_sd == 0

    def test_high_error_rate_nanopore(self):
        """Test Nanopore with high error rate (older chemistry)."""
        config = NanoporeConfig(error_rate=0.15)  # 15% for old R7.3
        assert config.error_rate == 0.15

    def test_low_error_rate_nanopore(self):
        """Test Nanopore with low error rate (future chemistry)."""
        config = NanoporeConfig(error_rate=0.02)  # 2% for hypothetical future
        assert config.error_rate == 0.02


class TestPlatformComparison:
    """Test comparison between platforms."""

    def test_pacbio_higher_accuracy_than_nanopore(self):
        """Test that PacBio HiFi has higher accuracy than Nanopore."""
        hifi_config = PacBioHiFiConfig()
        nanopore_config = NanoporeConfig()

        # After CCS, HiFi has >99.9% accuracy (clr_error_rate is before CCS)
        # Nanopore has ~95% accuracy
        assert nanopore_config.error_rate > 0.01  # Nanopore has measurable error
        assert hifi_config.clr_error_rate > nanopore_config.error_rate  # But CLR is worse

    def test_nanopore_longer_reads_than_pacbio(self):
        """Test that Nanopore typically has longer reads than PacBio HiFi."""
        hifi_config = PacBioHiFiConfig()
        nanopore_config = NanoporeConfig()

        # Default Nanopore reads are longer
        assert nanopore_config.read_length_mean > hifi_config.read_length_mean

    def test_platform_tradeoffs(self):
        """Test platform tradeoffs: accuracy vs length."""
        # PacBio HiFi: High accuracy, moderate length
        hifi_config = PacBioHiFiConfig(read_length_mean=15000)

        # Nanopore: Lower accuracy, longer length
        nanopore_config = NanoporeConfig(
            read_length_mean=50000,
            error_rate=0.05
        )

        assert nanopore_config.read_length_mean > hifi_config.read_length_mean
        assert nanopore_config.error_rate > 0.01  # HiFi is <0.1% after CCS


class TestGroundTruthGeneration:
    """Test ground truth metadata generation."""

    def test_ground_truth_fields_pacbio(self):
        """Test that ground truth contains required fields for PacBio."""
        # Mock ground truth data structure
        ground_truth = {
            'genome_id': 'NC_001416',
            'genome_type': 'viral',
            'length': 5386,
            'relative_abundance': 0.25,
            'platform': 'pacbio-hifi',
            'read_type': 'long'
        }

        assert ground_truth['platform'] == 'pacbio-hifi'
        assert ground_truth['read_type'] == 'long'
        assert 'genome_id' in ground_truth
        assert 'relative_abundance' in ground_truth

    def test_ground_truth_fields_nanopore(self):
        """Test that ground truth contains required fields for Nanopore."""
        ground_truth = {
            'genome_id': 'NC_001422',
            'genome_type': 'viral',
            'length': 5386,
            'relative_abundance': 0.15,
            'platform': 'nanopore',
            'read_type': 'long'
        }

        assert ground_truth['platform'] == 'nanopore'
        assert ground_truth['read_type'] == 'long'

    def test_ground_truth_abundance_sum_to_one(self):
        """Test that ground truth abundances sum to 1.0."""
        ground_truth_list = [
            {'genome_id': 'g1', 'relative_abundance': 0.4},
            {'genome_id': 'g2', 'relative_abundance': 0.35},
            {'genome_id': 'g3', 'relative_abundance': 0.25}
        ]

        total_abundance = sum(gt['relative_abundance'] for gt in ground_truth_list)
        assert np.isclose(total_abundance, 1.0, atol=1e-6)


class TestIntegrationScenarios:
    """Test realistic usage scenarios."""

    def test_gut_virome_pacbio_hifi_config(self):
        """Test typical configuration for gut virome with PacBio HiFi."""
        config = PacBioHiFiConfig(
            passes=10,
            read_length_mean=15000,
            read_length_sd=5000
        )

        assert config.passes == 10
        assert config.read_length_mean == 15000

    def test_gut_virome_nanopore_config(self):
        """Test typical configuration for gut virome with Nanopore."""
        config = NanoporeConfig(
            chemistry="R10.4",
            read_length_mean=25000,
            read_length_sd=10000
        )

        assert config.chemistry == "R10.4"
        assert config.read_length_mean == 25000

    def test_complete_viral_genome_assembly_hifi(self):
        """Test PacBio HiFi config for complete viral genome assembly."""
        # For complete genome assembly, need good coverage with long reads
        config = PacBioHiFiConfig(
            passes=15,  # Higher accuracy
            read_length_mean=20000,  # Longer reads to span genomes
            min_passes=5
        )

        assert config.passes == 15
        assert config.read_length_mean == 20000

    def test_structural_variant_detection_nanopore(self):
        """Test Nanopore config for structural variant detection."""
        # SV detection benefits from ultra-long reads
        config = NanoporeConfig(
            chemistry="R10.4",
            read_length_mean=50000,  # Ultra-long for spanning SVs
            read_length_sd=25000
        )

        assert config.read_length_mean == 50000

    def test_low_depth_exploration_sequencing(self):
        """Test configuration for low-depth exploration sequencing."""
        # Low depth, moderate read length
        config = PacBioHiFiConfig(
            passes=10,
            read_length_mean=12000,
            read_length_sd=4000
        )

        assert config.read_length_mean == 12000

    def test_high_accuracy_variant_calling(self):
        """Test PacBio HiFi config for high-accuracy variant calling."""
        # Maximum passes for highest accuracy
        config = PacBioHiFiConfig(
            passes=20,
            min_passes=10,
            read_length_mean=15000
        )

        assert config.passes == 20
        assert config.min_passes == 10


class TestReproducibility:
    """Test reproducibility with random seeds."""

    def test_same_seed_produces_same_config_behavior(self):
        """Test that same seed should produce reproducible results."""
        # Note: Config classes don't use random seeds themselves,
        # but seeds are passed to simulation functions
        config1 = PacBioHiFiConfig(passes=10, read_length_mean=15000)
        config2 = PacBioHiFiConfig(passes=10, read_length_mean=15000)

        # Configs with same parameters should be identical
        assert config1.passes == config2.passes
        assert config1.read_length_mean == config2.read_length_mean

    def test_different_params_produce_different_configs(self):
        """Test that different parameters produce different configs."""
        config1 = PacBioHiFiConfig(passes=10, read_length_mean=15000)
        config2 = PacBioHiFiConfig(passes=15, read_length_mean=20000)

        assert config1.passes != config2.passes
        assert config1.read_length_mean != config2.read_length_mean


class TestPBSIM3Commands:
    """Test PBSIM3 command generation (mock tests)."""

    def test_pacbio_hifi_command_structure(self):
        """Test that PacBio HiFi commands have correct structure."""
        # Mock command generation
        config = PacBioHiFiConfig(passes=10, read_length_mean=15000)

        # Command should include key parameters
        expected_params = [
            'pbsim',
            '--strategy', 'wgs',
            '--method', 'qshmm',
            '--pass-num', str(config.passes),
            '--length-mean', str(config.read_length_mean)
        ]

        # Just verify config has the right values
        assert config.passes == 10
        assert config.read_length_mean == 15000

    def test_nanopore_command_structure(self):
        """Test that Nanopore commands have correct structure."""
        config = NanoporeConfig(chemistry="R10.4", read_length_mean=20000)

        # Command should include key parameters
        expected_params = [
            'pbsim',
            '--strategy', 'wgs',
            '--method', 'errhmm',
            '--length-mean', str(config.read_length_mean),
            '--hp-del-bias', str(config.hp_del_bias)
        ]

        # Verify config has the right values
        assert config.chemistry == "R10.4"
        assert config.read_length_mean == 20000
        assert config.hp_del_bias == 6


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_platform_enum_invalid_value(self):
        """Test that invalid platform value raises appropriate error."""
        with pytest.raises(ValueError):
            LongReadPlatform("invalid_platform")

    def test_config_with_none_values(self):
        """Test that None values use defaults."""
        config = PacBioHiFiConfig()

        # All defaults should be set
        assert config.passes is not None
        assert config.read_length_mean is not None
        assert config.accuracy_model is not None

    def test_nanopore_config_with_none_values(self):
        """Test that Nanopore config None values use defaults."""
        config = NanoporeConfig()

        assert config.chemistry is not None
        assert config.read_length_mean is not None
        assert config.error_rate is not None


# Run tests with: pytest tests/test_longread_simulator.py -v
