"""Tests for Module 8: Novel Virus Discovery Benchmarking."""

import json
import textwrap
from pathlib import Path

import pytest

from viroforge.benchmarking.discovery import benchmark_discovery


TAXONOMY_GT = {
    "NC_KNOWN1": {"ncbi_taxid": 100, "is_known": True, "family": "Siphoviridae"},
    "NC_KNOWN2": {"ncbi_taxid": 200, "is_known": True, "family": "Myoviridae"},
    "NC_DARK1": {"ncbi_taxid": 300, "is_known": False, "family": "Unknown"},
    "NC_DARK2": {"ncbi_taxid": 400, "is_known": False, "family": "Unknown"},
    "NC_DARK3": {"ncbi_taxid": 500, "is_known": False, "family": "Unknown"},
}


class TestBenchmarkDiscovery:

    def test_perfect_detection(self):
        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_KNOWN2_1_1/1": 200,
            "NC_DARK1_1_1/1": 999,
            "NC_DARK2_1_1/1": 998,
            "NC_DARK3_1_1/1": 997,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        assert m["reliable"]
        assert m["n_known_reads"] == 2
        assert m["n_dark_matter_reads"] == 3
        assert m["known_detection_rate"] == 1.0
        assert m["dark_matter_detection_rate"] == 1.0
        assert m["n_dark_genomes"] == 3
        assert m["n_dark_genomes_detected"] == 3
        assert m["genome_detection_rate"] == 1.0

    def test_no_dark_matter_detected(self):
        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_DARK1_1_1/1": None,
            "NC_DARK2_1_1/1": None,
            "NC_DARK3_1_1/1": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        assert m["dark_matter_detection_rate"] == 0.0
        assert m["n_dark_genomes_detected"] == 0

    def test_partial_detection(self):
        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_DARK1_1_1/1": 999,
            "NC_DARK1_2_1/1": 999,
            "NC_DARK2_1_1/1": None,
            "NC_DARK2_2_1/1": None,
            "NC_DARK3_1_1/1": 997,
            "NC_DARK3_2_1/1": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        assert m["n_dark_matter_reads"] == 6
        assert m["n_dark_detected"] == 3
        assert m["dark_matter_detection_rate"] == 0.5
        assert m["n_dark_genomes_detected"] == 2  # DARK1 and DARK3

    def test_non_viral_reads_excluded(self):
        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_DARK1_1_1/1": None,
            "CONTAMINANT_1_1/1": 777,
            "HOST_DNA_2_1/1": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        assert m["n_non_viral_reads"] == 2
        assert m["n_viral_reads"] == 2

    def test_empty_assignments(self):
        m = benchmark_discovery({}, TAXONOMY_GT)
        assert not m["reliable"]
        assert m["n_viral_reads"] == 0
        assert m["known_detection_rate"] is None
        assert m["dark_matter_detection_rate"] is None

    def test_no_dark_matter_genomes(self):
        gt_known_only = {
            "NC_A": {"ncbi_taxid": 1, "is_known": True, "family": "Siphoviridae"},
        }
        assignments = {"NC_A_1_1/1": 1}
        m = benchmark_discovery(assignments, gt_known_only)
        assert m["n_dark_matter_reads"] == 0
        assert m["dark_matter_detection_rate"] is None
        assert m["genome_detection_rate"] is None

    def test_per_genome_details(self):
        assignments = {
            "NC_DARK1_1_1/1": 999,
            "NC_DARK1_2_1/1": 999,
            "NC_DARK2_1_1/1": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        pg = {g["genome_id"]: g for g in m["per_genome"]}
        assert pg["NC_DARK1"]["total_reads"] == 2
        assert pg["NC_DARK1"]["detected_reads"] == 2
        assert pg["NC_DARK1"]["detection_rate"] == 1.0
        assert pg["NC_DARK2"]["total_reads"] == 1
        assert pg["NC_DARK2"]["detected_reads"] == 0
        assert pg["NC_DARK2"]["detection_rate"] == 0.0

    def test_detection_by_depth_bins(self):
        assignments = {}
        for i in range(5):
            assignments[f"NC_DARK1_{i}_1/1"] = None
        for i in range(15):
            assignments[f"NC_DARK2_{i}_1/1"] = 998
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        bins = {b["bin"]: b for b in m["detection_by_depth"]}
        assert "1-10" in bins
        assert bins["1-10"]["n_genomes"] == 1  # DARK1 has 5 reads
        assert bins["1-10"]["n_detected"] == 0
        assert "11-50" in bins
        assert bins["11-50"]["n_genomes"] == 1  # DARK2 has 15 reads
        assert bins["11-50"]["n_detected"] == 1

    def test_duplicate_reads_parsed(self):
        assignments = {
            "NC_DARK1_1_1/1": 999,
            "NC_DARK1_1_1_dup1": 999,
            "NC_DARK1_1_1_dup2": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        pg = {g["genome_id"]: g for g in m["per_genome"]}
        assert pg["NC_DARK1"]["total_reads"] == 3
        assert pg["NC_DARK1"]["detected_reads"] == 2


class TestDiscoveryReport:

    def test_markdown_renders(self):
        from viroforge.benchmarking.report import discovery_to_markdown

        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_DARK1_1_1/1": 999,
            "NC_DARK2_1_1/1": None,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        md = discovery_to_markdown(m)
        assert "# Novel Discovery Benchmark" in md
        assert "Dark matter detection" in md
        assert "Known virus detection" in md

    def test_unreliable_warning(self):
        from viroforge.benchmarking.report import discovery_to_markdown

        m = benchmark_discovery({}, TAXONOMY_GT)
        md = discovery_to_markdown(m)
        assert "WARNING" in md


class TestDiscoveryWriteReports:

    def test_json_and_markdown_written(self, tmp_path):
        from viroforge.benchmarking.report import write_reports

        assignments = {
            "NC_KNOWN1_1_1/1": 100,
            "NC_DARK1_1_1/1": 999,
        }
        m = benchmark_discovery(assignments, TAXONOMY_GT)
        json_path = tmp_path / "discovery.json"
        md_path = tmp_path / "discovery.md"
        write_reports(m, json_path=str(json_path), md_path=str(md_path), kind="discovery")

        assert json_path.exists()
        assert md_path.exists()
        data = json.loads(json_path.read_text())
        assert data["n_known_reads"] == 1
        assert "Novel Discovery" in md_path.read_text()


class TestDiscoveryCLIIntegration:

    def test_discovery_subcommand_registered(self):
        from viroforge.cli import main
        import sys
        old_argv = sys.argv
        sys.argv = ["viroforge", "benchmark", "discovery", "--help"]
        try:
            main()
        except SystemExit as e:
            assert e.code == 0
        finally:
            sys.argv = old_argv

    def test_discovery_runs_with_kraken2_file(self, tmp_path):
        gt = {
            "benchmarking": {
                "taxonomy": TAXONOMY_GT,
            }
        }
        gt_path = tmp_path / "metadata.json"
        gt_path.write_text(json.dumps(gt))

        kraken = tmp_path / "kraken2.out"
        kraken.write_text(textwrap.dedent("""\
            C\tNC_KNOWN1_1_1/1\t100\t150\t100:150
            U\tNC_DARK1_1_1/1\t0\t150\t0:150
            C\tNC_DARK2_1_1/1\t998\t150\t998:150
        """))

        output_json = tmp_path / "out.json"
        output_md = tmp_path / "out.md"

        from viroforge.cli.benchmark import _run_discovery
        import argparse
        args = argparse.Namespace(
            pipeline_output=str(kraken),
            ground_truth=str(gt_path),
            format="kraken2",
            read_id_column=1,
            taxid_column=2,
            output=str(output_json),
            markdown=str(output_md),
        )
        rc = _run_discovery(args)
        assert rc == 0
        data = json.loads(output_json.read_text())
        assert data["n_dark_detected"] == 1  # DARK2 classified
        assert data["n_dark_matter_reads"] == 2
