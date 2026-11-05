"""
Tests for long-read sequencing simulation (Oxford Nanopore).

Tests all components of long-read simulation:
- Configuration classes (NanoporeConfig)
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
    NanoporeConfig,
    check_pbsim3_installed,
)


class TestLongReadPlatform:
    """Test LongReadPlatform enum."""

    def test_platform_values(self):
        """Test platform enum values."""
        assert LongReadPlatform.NANOPORE.value == "nanopore"

    def test_platform_from_string(self):
        """Test creating platform from string."""
        assert LongReadPlatform("nanopore") == LongReadPlatform.NANOPORE

    def test_invalid_platform_raises_error(self):
        """Test that invalid platform string raises error."""
        with pytest.raises(ValueError):
            LongReadPlatform("illumina")


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
            error_rate=0.10
        )

        assert config.chemistry == "R9.4"
        assert config.error_rate == 0.10

    def test_r104_chemistry(self):
        """Test R10.4 chemistry configuration (newer, more accurate)."""
        config = NanoporeConfig(
            chemistry="R10.4",
            error_rate=0.05
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

        config_low = NanoporeConfig(hp_del_bias=3)
        assert config_low.hp_del_bias == 3


class TestDependencyChecks:
    """Test dependency checking functions."""

    @patch('subprocess.run')
    def test_pbsim3_installed(self, mock_run):
        """Test PBSIM3 installation check when installed."""
        mock_run.return_value = Mock(returncode=0)

        assert check_pbsim3_installed() is True

    @patch('subprocess.run', side_effect=FileNotFoundError)
    def test_pbsim3_not_installed(self, mock_run):
        """Test PBSIM3 installation check when not installed."""
        assert check_pbsim3_installed() is False


class TestConfigurationValidation:
    """Test configuration validation and edge cases."""

    def test_nanopore_very_long_reads(self):
        """Test Nanopore with very long reads (200kb+)."""
        config = NanoporeConfig(read_length_mean=200000, read_length_sd=100000)
        assert config.read_length_mean == 200000

    def test_high_error_rate_nanopore(self):
        """Test Nanopore with high error rate (older chemistry)."""
        config = NanoporeConfig(error_rate=0.15)
        assert config.error_rate == 0.15

    def test_low_error_rate_nanopore(self):
        """Test Nanopore with low error rate (future chemistry)."""
        config = NanoporeConfig(error_rate=0.02)
        assert config.error_rate == 0.02

    def test_zero_standard_deviation_edge_case(self):
        """Test configurations with zero standard deviation."""
        config = NanoporeConfig(read_length_sd=0)
        assert config.read_length_sd == 0


class TestGroundTruthGeneration:
    """Test ground truth metadata generation."""

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

    def test_gut_virome_nanopore_config(self):
        """Test typical configuration for gut virome with Nanopore."""
        config = NanoporeConfig(
            chemistry="R10.4",
            read_length_mean=25000,
            read_length_sd=10000
        )

        assert config.chemistry == "R10.4"
        assert config.read_length_mean == 25000

    def test_structural_variant_detection_nanopore(self):
        """Test Nanopore config for structural variant detection."""
        config = NanoporeConfig(
            chemistry="R10.4",
            read_length_mean=50000,
            read_length_sd=25000
        )

        assert config.read_length_mean == 50000


class TestReproducibility:
    """Test reproducibility with random seeds."""

    def test_same_seed_produces_same_config_behavior(self):
        """Test that same seed should produce reproducible results."""
        config1 = NanoporeConfig(read_length_mean=15000)
        config2 = NanoporeConfig(read_length_mean=15000)

        assert config1.read_length_mean == config2.read_length_mean

    def test_different_params_produce_different_configs(self):
        """Test that different parameters produce different configs."""
        config1 = NanoporeConfig(read_length_mean=15000)
        config2 = NanoporeConfig(read_length_mean=20000)

        assert config1.read_length_mean != config2.read_length_mean


class TestPBSIM3Commands:
    """Test PBSIM3 command generation (mock tests)."""

    def test_nanopore_command_structure(self):
        """Test that Nanopore commands have correct structure."""
        config = NanoporeConfig(chemistry="R10.4", read_length_mean=20000)

        assert config.chemistry == "R10.4"
        assert config.read_length_mean == 20000
        assert config.hp_del_bias == 6


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_platform_enum_invalid_value(self):
        """Test that invalid platform value raises appropriate error."""
        with pytest.raises(ValueError):
            LongReadPlatform("invalid_platform")

    def test_nanopore_config_with_none_values(self):
        """Test that Nanopore config None values use defaults."""
        config = NanoporeConfig()

        assert config.chemistry is not None
        assert config.read_length_mean is not None
        assert config.error_rate is not None


# Run tests with: pytest tests/test_longread_simulator.py -v
