"""Taxonomy benchmark tests with an independently-counted oracle.

Assignments with a fixed, recorded mix of correct / unclassified / misclassified
across the known-virus and dark-matter strata; the engine must recover the
hand-computed per-stratum metrics.
"""

import pytest

from viroforge.benchmarking.taxonomy import (
    benchmark_taxonomy,
    parse_generic,
    parse_genome_id,
    parse_kraken2,
)

TAX_GT = {
    "K0": {"ncbi_taxid": 100, "is_known": True},
    "K1": {"ncbi_taxid": 101, "is_known": True},
    "D0": {"ncbi_taxid": 200, "is_known": False},
    "D1": {"ncbi_taxid": 201, "is_known": False},
}


def _assignments():
    a = {}
    # K0 (taxid 100): 7 correct, 2 unclassified, 1 misclassified
    for i in range(7):
        a[f"K0_{i}_0/1"] = 100
    for i in range(7, 9):
        a[f"K0_{i}_0/1"] = None
    a["K0_9_0/1"] = 999
    # K1 (taxid 101): 10 correct
    for i in range(10):
        a[f"K1_{i}_0/1"] = 101
    # D0 (dark, taxid 200): 1 correct, 5 unclassified
    a["D0_0_0/1"] = 200
    for i in range(1, 6):
        a[f"D0_{i}_0/1"] = None
    # D1 (dark, taxid 201): 3 unclassified, 1 misclassified
    for i in range(3):
        a[f"D1_{i}_0/1"] = None
    a["D1_3_0/1"] = 888
    # 3 non-viral (contaminant) reads
    for i in range(3):
        a[f"host_human_000{i}_0_0/1"] = 9606
    return a


def test_parse_genome_id():
    assert parse_genome_id("GCF_000819615.1_991_0/1") == "GCF_000819615.1"
    assert parse_genome_id("GCF_000819615.1_991_0") == "GCF_000819615.1"
    assert parse_genome_id("host_human_0000_0_0/1") == "host_human_0000"
    assert parse_genome_id("GCF_1.1_5_0/1_dup3") == "GCF_1.1"


def test_parse_kraken2(tmp_path):
    p = tmp_path / "k2.out"
    p.write_text(
        "C\tK0_0_0/1\t100\t150\t100:1\n"
        "U\tK0_7_0/1\t0\t150\t0:1\n"
        "C\tK1_0_0/1\tSinsheimervirus (taxid 101)\t150\t101:1\n"
    )
    a = parse_kraken2(p)
    assert a["K0_0_0/1"] == 100
    assert a["K0_7_0/1"] is None
    assert a["K1_0_0/1"] == 101  # --use-names form


def test_parse_generic(tmp_path):
    p = tmp_path / "gen.tsv"
    p.write_text("# read_id\ttaxid\nK0_0_0/1\t100\nK0_7_0/1\t0\n")
    a = parse_generic(p)
    assert a["K0_0_0/1"] == 100
    assert a["K0_7_0/1"] is None


def test_oracle_metrics():
    r = benchmark_taxonomy(_assignments(), TAX_GT)

    assert r["reliable"] is True
    assert r["n_assignments"] == 33
    assert r["n_viral_reads"] == 30
    assert r["n_non_viral_reads"] == 3

    k = r["known_viruses"]
    assert (k["n"], k["correct"], k["unclassified"], k["misclassified"]) == (20, 17, 2, 1)
    assert k["sensitivity"] == pytest.approx(0.85)
    assert k["precision"] == pytest.approx(17 / 18)
    assert k["unclassified_rate"] == pytest.approx(0.10)

    d = r["dark_matter"]
    assert (d["n"], d["correct"], d["unclassified"], d["misclassified"]) == (10, 1, 8, 1)
    assert d["sensitivity"] == pytest.approx(0.10)
    assert d["precision"] == pytest.approx(0.5)
    assert d["unclassified_rate"] == pytest.approx(0.80)

    ab = r["abundance_profile"]
    # union of true {100,101,200,201} and observed {100,101,200,888,999}
    assert ab["n_taxa"] == 6
    assert 0.0 <= ab["bray_curtis"] <= 1.0
    assert ab["pearson"] is not None


def test_unreliable_when_no_viral():
    # read ids that map to no known viral genome
    r = benchmark_taxonomy({"foo_0_0/1": 100, "bar_0_0/1": None}, TAX_GT)
    assert r["reliable"] is False
    assert r["n_viral_reads"] == 0
