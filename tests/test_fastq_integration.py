"""
Integration tests for FASTQ generation with enhanced VLP enrichment

Tests the complete workflow from collection loading through FASTQ generation
with VLP enrichment and contamination reduction.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import json
import subprocess
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestFASTQGeneration:
    """Test FASTQ generation workflow"""

    @pytest.fixture
    def temp_output(self):
        """Create temporary output directory"""
        temp_dir = Path(tempfile.mkdtemp())
        yield temp_dir
        # Cleanup
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    @pytest.fixture
    def generation_script(self):
        """Path to generation script"""
        script = Path(__file__).parent.parent / 'scripts' / 'generate_fastq_dataset.py'
        assert script.exists(), f"Generation script not found: {script}"
        return script

    def test_vlp_tangential_flow_generation(self, generation_script, temp_output):
        """Test FASTQ generation with tangential flow VLP protocol"""
        cmd = [
            'python', str(generation_script),
            '--collection-id', '16',  # Mouse gut (smallest)
            '--output', str(temp_output),
            '--n-reads', '100',  # Very small for speed
            '--vlp-protocol', 'tangential_flow',
            '--contamination-level', 'clean',
            '--seed', '42'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check success
        assert result.returncode == 0, f"Generation failed: {result.stderr}"

        # Check output files exist
        fasta_dir = temp_output / 'fasta'
        fastq_dir = temp_output / 'fastq'
        metadata_dir = temp_output / 'metadata'

        assert fasta_dir.exists(), "FASTA directory not created"
        assert fastq_dir.exists(), "FASTQ directory not created"
        assert metadata_dir.exists(), "Metadata directory not created"

        # Check FASTQ files
        fastq_files = list(fastq_dir.glob('*_R1.fastq'))
        assert len(fastq_files) == 1, "R1 FASTQ not found"

        fastq_files = list(fastq_dir.glob('*_R2.fastq'))
        assert len(fastq_files) == 1, "R2 FASTQ not found"

        # Check metadata
        metadata_files = list(metadata_dir.glob('*_metadata.json'))
        assert len(metadata_files) == 1, "Metadata JSON not found"

        # Validate metadata content
        with open(metadata_files[0]) as f:
            metadata = json.load(f)

        assert 'enrichment_stats' in metadata, "Enrichment stats missing"
        assert metadata['enrichment_stats']['vlp_protocol'] == 'tangential_flow'
        assert 'viral_enrichment' in metadata['enrichment_stats']
        assert 'contamination_reduction' in metadata['enrichment_stats']

    def test_bulk_metagenome_generation(self, generation_script, temp_output):
        """Test FASTQ generation without VLP (bulk metagenome)"""
        cmd = [
            'python', str(generation_script),
            '--collection-id', '16',
            '--output', str(temp_output),
            '--n-reads', '100',
            '--no-vlp',
            '--contamination-level', 'realistic',
            '--seed', '42'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Generation failed: {result.stderr}"

        # Load metadata
        metadata_file = list((temp_output / 'metadata').glob('*_metadata.json'))[0]
        with open(metadata_file) as f:
            metadata = json.load(f)

        # Check bulk mode
        assert metadata['enrichment_stats']['vlp_protocol'] == 'none'

        # Contamination should NOT be reduced
        contam_fraction = metadata['enrichment_stats']['contamination_fraction']
        assert contam_fraction > 0.05, f"Contamination too low for bulk: {contam_fraction}"

    def test_all_vlp_protocols(self, generation_script, temp_output):
        """Test that all VLP protocols generate successfully"""
        protocols = ['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen']

        for protocol in protocols:
            protocol_output = temp_output / protocol
            cmd = [
                'python', str(generation_script),
                '--collection-id', '16',
                '--output', str(protocol_output),
                '--n-reads', '100',
                '--vlp-protocol', protocol,
                '--contamination-level', 'clean',
                '--seed', '42'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Protocol {protocol} failed: {result.stderr}"

            # Check metadata
            metadata_file = list((protocol_output / 'metadata').glob('*_metadata.json'))[0]
            with open(metadata_file) as f:
                metadata = json.load(f)

            assert metadata['enrichment_stats']['vlp_protocol'] == protocol

    def test_contamination_levels(self, generation_script, temp_output):
        """Test different contamination levels"""
        levels = ['clean', 'realistic', 'heavy']
        # After VLP enrichment with tangential flow (~85-95% reduction):
        # - Clean (0.7% initial) → ~0.05-0.15% final
        # - Realistic (7.4% initial) → ~0.4-1.5% final
        # - Heavy (27% initial) → ~1.5-5% final
        expected_ranges = {
            'clean': (0.0005, 0.0020),      # 0.05-0.2% after VLP
            'realistic': (0.004, 0.020),    # 0.4-2% after VLP
            'heavy': (0.015, 0.050)         # 1.5-5% after VLP
        }

        for level in levels:
            level_output = temp_output / level
            cmd = [
                'python', str(generation_script),
                '--collection-id', '16',
                '--output', str(level_output),
                '--n-reads', '100',
                '--vlp-protocol', 'tangential_flow',
                '--contamination-level', level,
                '--seed', '42'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Contamination level {level} failed"

            # Check contamination fraction
            metadata_file = list((level_output / 'metadata').glob('*_metadata.json'))[0]
            with open(metadata_file) as f:
                metadata = json.load(f)

            contam_fraction = metadata['enrichment_stats']['contamination_fraction']
            min_expected, max_expected = expected_ranges[level]

            assert min_expected <= contam_fraction <= max_expected, \
                f"Contamination {contam_fraction} outside expected range {expected_ranges[level]} for level {level}"

    def test_dry_run_mode(self, generation_script, temp_output):
        """Test dry-run mode (no FASTQ generation)"""
        cmd = [
            'python', str(generation_script),
            '--collection-id', '16',
            '--output', str(temp_output),
            '--coverage', '10',
            '--vlp-protocol', 'tangential_flow',
            '--dry-run'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Dry run failed: {result.stderr}"

        # FASTA and metadata should exist
        assert (temp_output / 'fasta').exists()
        assert (temp_output / 'metadata').exists()

        # FASTQ should NOT be generated
        fastq_files = list((temp_output / 'fastq').glob('*.fastq'))
        assert len(fastq_files) == 0, "FASTQ files generated in dry-run mode"

    def test_amplification_rdab(self, generation_script, temp_output):
        """Test FASTQ generation with RdAB amplification"""
        cmd = [
            'python', str(generation_script),
            '--collection-id', '16',
            '--output', str(temp_output),
            '--n-reads', '100',
            '--vlp-protocol', 'tangential_flow',
            '--contamination-level', 'clean',
            '--amplification', 'rdab',
            '--seed', '42'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"RdAB generation failed: {result.stderr}"

        # Check metadata includes amplification stats
        metadata_file = list((temp_output / 'metadata').glob('*_metadata.json'))[0]
        with open(metadata_file) as f:
            metadata = json.load(f)

        assert 'amplification_stats' in metadata, "Amplification stats missing"
        assert metadata['amplification_stats']['method'] == 'rdab'
        assert metadata['amplification_stats']['bias_applied'] == True
        assert 'mean_fold_change' in metadata['amplification_stats']

    def test_all_amplification_methods(self, generation_script, temp_output):
        """Test that all amplification methods generate successfully"""
        methods = ['none', 'rdab', 'rdab-30', 'mda', 'mda-long', 'linker']

        for method in methods:
            method_output = temp_output / method
            cmd = [
                'python', str(generation_script),
                '--collection-id', '16',
                '--output', str(method_output),
                '--n-reads', '100',
                '--vlp-protocol', 'tangential_flow',
                '--contamination-level', 'clean',
                '--amplification', method,
                '--seed', '42'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Amplification method {method} failed: {result.stderr}"

            # Check metadata
            metadata_file = list((method_output / 'metadata').glob('*_metadata.json'))[0]
            with open(metadata_file) as f:
                metadata = json.load(f)

            assert metadata['amplification_stats']['method'] == method

            if method == 'none':
                assert metadata['amplification_stats']['bias_applied'] == False
            else:
                assert metadata['amplification_stats']['bias_applied'] == True
                # Should have non-trivial fold changes for biased methods
                assert metadata['amplification_stats']['mean_fold_change'] > 0


class TestVLPEnrichmentStats:
    """Test VLP enrichment statistics"""

    @pytest.fixture
    def temp_output(self):
        """Create temporary output directory"""
        temp_dir = Path(tempfile.mkdtemp())
        yield temp_dir
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    @pytest.fixture
    def generation_script(self):
        """Path to generation script"""
        script = Path(__file__).parent.parent / 'scripts' / 'generate_fastq_dataset.py'
        return script

    def _generate_and_load_metadata(self, generation_script, temp_output, vlp_protocol, contam_level):
        """Helper to generate dataset and load metadata"""
        cmd = [
            'python', str(generation_script),
            '--collection-id', '16',
            '--output', str(temp_output),
            '--n-reads', '100',
            '--vlp-protocol', vlp_protocol,
            '--contamination-level', contam_level,
            '--seed', '42'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Generation failed: {result.stderr}"

        metadata_file = list((temp_output / 'metadata').glob('*_metadata.json'))[0]
        with open(metadata_file) as f:
            return json.load(f)

    def test_viral_fraction_high_with_vlp(self, generation_script, temp_output):
        """Test that VLP enrichment results in high viral fraction"""
        metadata = self._generate_and_load_metadata(
            generation_script, temp_output, 'tangential_flow', 'realistic'
        )

        viral_fraction = metadata['enrichment_stats']['viral_fraction']

        # Should be >95% viral after VLP
        assert viral_fraction > 0.95, f"Viral fraction too low: {viral_fraction}"

    def test_contamination_reduction_documented(self, generation_script, temp_output):
        """Test that contamination reduction statistics are complete"""
        metadata = self._generate_and_load_metadata(
            generation_script, temp_output, 'tangential_flow', 'realistic'
        )

        contam_stats = metadata['enrichment_stats']['contamination_reduction']

        # Check required fields
        assert 'original_total_contamination' in contam_stats
        assert 'reduced_total_contamination' in contam_stats
        assert 'overall_reduction_factor' in contam_stats
        assert 'reduction_by_type' in contam_stats

        # Check type-specific reduction
        reduction_by_type = contam_stats['reduction_by_type']
        expected_types = ['host_dna', 'rrna', 'reagent_bacteria', 'phix']

        for ctype in expected_types:
            assert ctype in reduction_by_type, f"Missing contamination type: {ctype}"
            assert 'reduction_factor' in reduction_by_type[ctype]
            assert 'reduction_pct' in reduction_by_type[ctype]

    def test_size_bias_correlation(self, generation_script, temp_output):
        """Test that size bias correlation is strong (>0.9)"""
        metadata = self._generate_and_load_metadata(
            generation_script, temp_output, 'tangential_flow', 'realistic'
        )

        viral_stats = metadata['enrichment_stats']['viral_enrichment']
        size_bias = viral_stats['size_bias']

        correlation = size_bias['size_enrichment_correlation']

        # Should have strong positive correlation (larger viruses enriched)
        assert correlation > 0.9, f"Size bias correlation too weak: {correlation}"

    def test_host_dna_highly_reduced(self, generation_script, temp_output):
        """Test that host DNA is highly reduced by nuclease (>85%)"""
        metadata = self._generate_and_load_metadata(
            generation_script, temp_output, 'tangential_flow', 'realistic'
        )

        reduction_by_type = metadata['enrichment_stats']['contamination_reduction']['reduction_by_type']
        host_dna_reduction = reduction_by_type['host_dna']['reduction_factor']

        # Nuclease should remove >85% host DNA
        assert host_dna_reduction > 0.85, f"Host DNA reduction too low: {host_dna_reduction}"

    def test_bacteria_highly_reduced(self, generation_script, temp_output):
        """Test that reagent bacteria are highly reduced by filtration (>95%)"""
        metadata = self._generate_and_load_metadata(
            generation_script, temp_output, 'tangential_flow', 'realistic'
        )

        reduction_by_type = metadata['enrichment_stats']['contamination_reduction']['reduction_by_type']
        bacteria_reduction = reduction_by_type['reagent_bacteria']['reduction_factor']

        # Filtration should remove >95% bacteria (1-5 μm >> 0.2 μm pore)
        assert bacteria_reduction > 0.95, f"Bacteria reduction too low: {bacteria_reduction}"


class TestBatchGeneration:
    """Test batch FASTQ generation"""

    @pytest.fixture
    def temp_output(self):
        """Create temporary output directory"""
        temp_dir = Path(tempfile.mkdtemp())
        yield temp_dir
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    @pytest.fixture
    def batch_script(self):
        """Path to batch generation script"""
        script = Path(__file__).parent.parent / 'scripts' / 'batch_generate_fastq.py'
        assert script.exists(), f"Batch script not found: {script}"
        return script

    def test_quick_test_preset(self, batch_script, temp_output):
        """Test quick-test preset"""
        cmd = [
            'python', str(batch_script),
            '--preset', 'quick-test',
            '--output', str(temp_output),
            '--dry-run'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Batch generation failed: {result.stderr}"

        # Should generate 3 datasets (collections 16, 11, 10)
        datasets = list(temp_output.glob('collection_*'))
        assert len(datasets) == 3, f"Expected 3 datasets, got {len(datasets)}"

    def test_vlp_protocol_comparison_preset(self, batch_script, temp_output):
        """Test VLP protocol comparison preset"""
        cmd = [
            'python', str(batch_script),
            '--preset', 'vlp-protocol-comparison',
            '--output', str(temp_output),
            '--dry-run'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Batch generation failed: {result.stderr}"

        # Should generate 5 datasets (4 VLP protocols + bulk)
        datasets = list(temp_output.glob('collection_*'))
        assert len(datasets) == 5, f"Expected 5 datasets, got {len(datasets)}"

        # Check protocol names in dataset folders
        protocol_names = [d.name for d in datasets]
        expected_protocols = ['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen', 'bulk']

        for protocol in expected_protocols:
            matches = [name for name in protocol_names if protocol in name]
            assert len(matches) >= 1, f"Protocol {protocol} not found in generated datasets"

    def test_amplification_comparison_preset(self, batch_script, temp_output):
        """Test amplification comparison preset"""
        cmd = [
            'python', str(batch_script),
            '--preset', 'amplification-comparison',
            '--output', str(temp_output),
            '--dry-run'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Batch generation failed: {result.stderr}"

        # Should generate 6 datasets (6 amplification methods)
        datasets = list(temp_output.glob('collection_*'))
        assert len(datasets) == 6, f"Expected 6 datasets, got {len(datasets)}"

        # Check amplification method names in dataset folders
        dataset_names = [d.name for d in datasets]
        expected_methods = ['rdab', 'rdab-30', 'mda', 'mda-long', 'linker']

        # Check that we have one baseline (no suffix for 'none')
        baseline = [name for name in dataset_names if not any(method in name for method in expected_methods)]
        assert len(baseline) == 1, f"Expected 1 baseline dataset (none), got {len(baseline)}"

        # Check that all amplification methods are present
        for method in expected_methods:
            matches = [name for name in dataset_names if method in name]
            assert len(matches) >= 1, f"Amplification method {method} not found in generated datasets"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
