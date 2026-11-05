"""Tests for the auto-generated HTML diagnostic report."""

import json
from pathlib import Path

import pytest

from viroforge.cli.auto_report import generate_auto_report


@pytest.fixture
def sample_dataset(tmp_path):
    """Create a minimal dataset directory with metadata and composition."""
    ds = tmp_path / "gut_10x"
    ds.mkdir()

    metadata_dir = ds / "metadata"
    metadata_dir.mkdir()
    fastq_dir = ds / "fastq"
    fastq_dir.mkdir()
    fasta_dir = ds / "fasta"
    fasta_dir.mkdir()

    metadata = {
        "metadata_version": "1.1",
        "generation_info": {
            "timestamp": "2026-09-21T12:00:00",
            "viroforge_version": "0.22.0",
            "random_seed": 42,
        },
        "collection": {
            "id": 1,
            "name": "Gut Virome - Adult Healthy",
            "description": "Western diet adult gut",
            "n_viral_genomes": 133,
            "n_contaminants": 66,
            "total_sequences": 199,
        },
        "configuration": {
            "coverage": 10,
            "read_length": 150,
            "insert_size": 350,
            "platform": "novaseq",
            "molecule_type": "dna",
            "vlp_protocol": "tangential_flow",
            "contamination_level": "realistic",
            "amplification": "linker",
            "dark_matter": {
                "dark_matter_fraction": 0.30,
                "dark_matter_genomes": 66,
                "known_viral_genomes": 133,
                "total_genomes": 199,
            },
        },
        "enrichment_stats": {
            "viral_fraction": 0.82,
            "contamination_fraction": 0.18,
            "n_viral_genomes": 133,
            "n_contaminants": 66,
            "viral_enrichment": {"mean_enrichment_factor": 4.2},
            "contamination_reduction": {
                "original_total_contamination": 0.80,
                "reduced_total_contamination": 0.18,
                "reduction_by_type": {
                    "host_dna": {
                        "original_abundance": 0.05,
                        "reduced_abundance": 0.002,
                        "reduction_factor": 0.95,
                    },
                    "bacterial_background": {
                        "original_abundance": 0.70,
                        "reduced_abundance": 0.015,
                        "reduction_factor": 0.98,
                    },
                    "reagent_bacteria": {
                        "original_abundance": 0.005,
                        "reduced_abundance": 0.0003,
                        "reduction_factor": 0.90,
                    },
                    "rrna": {
                        "original_abundance": 0.03,
                        "reduced_abundance": 0.003,
                        "reduction_factor": 0.90,
                    },
                    "phix": {
                        "original_abundance": 0.001,
                        "reduced_abundance": 0.0006,
                        "reduction_factor": 0.36,
                    },
                }
            },
        },
        "amplification_stats": {
            "method": "linker",
            "gc_bias_range": "0.7-1.3",
        },
        "adapter_stats": {
            "emergent_adapter_rate": 0.03,
            "reads_readthrough": 450,
            "adapter_type": "truseq",
        },
        "low_complexity_stats": {
            "rate": 0.01,
            "reads_modified": 150,
            "artifact_counts": {"homopolymer": 80, "dinucleotide": 70},
        },
        "duplicate_stats": {
            "rate": 0.30,
            "duplicates_created": 4500,
        },
        "benchmarking": {
            "expected_coverage": {
                "total_genome_length": 500000,
                "per_genome": [
                    {"genome_id": "GCF_001", "sequence_type": "viral",
                     "length": 100000, "relative_abundance": 0.15,
                     "expected_coverage": 1500.0, "coverage_category": "complete"},
                    {"genome_id": "GCF_002", "sequence_type": "viral",
                     "length": 5000, "relative_abundance": 0.08,
                     "expected_coverage": 800.0, "coverage_category": "complete"},
                    {"genome_id": "GCF_003", "sequence_type": "viral",
                     "length": 45000, "relative_abundance": 0.05,
                     "expected_coverage": 50.0, "coverage_category": "high"},
                    {"genome_id": "GCF_004", "sequence_type": "viral",
                     "length": 3000, "relative_abundance": 0.001,
                     "expected_coverage": 0.5, "coverage_category": "fragmented"},
                ],
            },
        },
    }
    (metadata_dir / "gut_10x_metadata.json").write_text(
        json.dumps(metadata, indent=2)
    )

    comp_lines = [
        "genome_id\tgenome_name\tsequence_type\trelative_abundance\tlength",
        "GCF_001\tCrAssphage\tviral\t0.15\t100000",
        "GCF_002\tBacteroides phage\tviral\t0.08\t5000",
        "GCF_003\tLactobacillus phage\tviral\t0.05\t45000",
        "GCF_004\tSmall phage\tviral\t0.001\t3000",
    ]
    (metadata_dir / "gut_10x_composition.tsv").write_text("\n".join(comp_lines))

    (fastq_dir / "gut_10x_R1.fastq").write_text("@read1\nACGT\n+\nIIII\n")
    (fastq_dir / "gut_10x_R2.fastq").write_text("@read1\nACGT\n+\nIIII\n")
    (fasta_dir / "gut_10x.fasta").write_text(">GCF_001\nACGT\n")

    return ds


def test_report_generated(sample_dataset):
    report = generate_auto_report(sample_dataset)
    assert report is not None
    assert report.exists()
    assert report.name == "gut_10x_report.html"


def test_report_contains_key_sections(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "ViroForge Dataset Report" in html
    assert "Generation Summary" in html
    assert "Genome Composition" in html
    assert "VLP Enrichment" in html
    assert "Artifact Injection" in html
    assert "Output Files" in html


def test_report_contains_metadata_values(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Gut Virome - Adult Healthy" in html
    assert "NOVASEQ" in html
    assert "10x" in html
    assert "0.22.0" in html
    assert "42" in html


def test_report_contains_composition(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "CrAssphage" in html
    assert "15.0000%" in html


def test_report_contains_artifacts(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Adapter Read-Through" in html
    assert "Low-Complexity" in html
    assert "PCR Duplicates" in html


def test_report_contains_dark_matter(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Dark Matter" in html
    assert "66" in html


def test_report_contains_file_inventory(sample_dataset):
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "gut_10x_R1.fastq" in html
    assert "gut_10x.fasta" in html


def test_report_self_contained(sample_dataset):
    """Report should not reference external CDNs."""
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "cdn.jsdelivr.net" not in html
    assert "cdnjs.cloudflare.com" not in html


def test_report_returns_none_for_missing_metadata(tmp_path):
    ds = tmp_path / "empty_dataset"
    ds.mkdir()
    assert generate_auto_report(ds) is None


def test_report_html_escapes_special_chars(tmp_path):
    ds = tmp_path / "test_escape"
    ds.mkdir()
    md = ds / "metadata"
    md.mkdir()
    metadata = {
        "generation_info": {"timestamp": "2026-01-01", "viroforge_version": "0.22.0"},
        "collection": {
            "name": 'Test <script>alert("xss")</script>',
            "id": 99,
            "n_viral_genomes": 1,
            "n_contaminants": 0,
            "total_sequences": 1,
        },
        "configuration": {"platform": "novaseq"},
        "enrichment_stats": {},
    }
    (md / "test_escape_metadata.json").write_text(json.dumps(metadata))
    report = generate_auto_report(ds)
    html = report.read_text()
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_report_contains_vlp_composition_chart(sample_dataset):
    """Pre/post VLP composition chart should be present with SVG bars."""
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Pre-VLP" in html
    assert "Post-VLP" in html
    assert "<svg" in html
    assert "Bacterial Background" in html
    assert "Host DNA" in html


def test_report_contains_rank_abundance_chart(sample_dataset):
    """Rank-abundance curve should appear when composition has multiple genomes."""
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Rank-Abundance Curve" in html
    assert "polyline" in html
    assert "Rank" in html


def test_report_contains_size_vs_abundance_chart(sample_dataset):
    """Genome size vs abundance scatter should appear with length data."""
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Genome Size vs Abundance" in html
    assert "Genome Length" in html


def test_report_contains_coverage_distribution(sample_dataset):
    """Expected coverage bar chart should appear with benchmarking data."""
    report = generate_auto_report(sample_dataset)
    html = report.read_text()
    assert "Expected Genome Recovery" in html
    assert "Complete" in html
    assert "Fragmented" in html
