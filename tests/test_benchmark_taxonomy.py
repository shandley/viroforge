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


# ---- higher-rank (genus/family) metrics ------------------------------------

def _write_taxdump(tmp_path):
    # family 10,11 ; genus 20,21 (in fam10), 22 (in fam11) ; species 100,101 (gen20),
    # 102 (gen21), 103 (gen22)
    nodes = "\n".join([
        "1\t|\t1\t|\tno rank\t|",
        "10\t|\t1\t|\tfamily\t|",
        "11\t|\t1\t|\tfamily\t|",
        "20\t|\t10\t|\tgenus\t|",
        "21\t|\t10\t|\tgenus\t|",
        "22\t|\t11\t|\tgenus\t|",
        "100\t|\t20\t|\tspecies\t|",
        "101\t|\t20\t|\tspecies\t|",
        "102\t|\t21\t|\tspecies\t|",
        "103\t|\t22\t|\tspecies\t|",
    ]) + "\n"
    p = tmp_path / "nodes.dmp"
    p.write_text(nodes)
    return p


def test_ncbi_tree_rank_resolution(tmp_path):
    from viroforge.benchmarking.ncbi_tree import NcbiTree
    tree = NcbiTree(_write_taxdump(tmp_path))
    assert tree.rank_taxid(100, "species") == 100
    assert tree.rank_taxid(100, "genus") == 20
    assert tree.rank_taxid(100, "family") == 10
    assert tree.rank_taxid(20, "family") == 10
    assert tree.rank_taxid(20, "species") is None  # genus has no species ancestor


def test_per_rank_metrics(tmp_path):
    from viroforge.benchmarking.ncbi_tree import NcbiTree
    tree = NcbiTree(_write_taxdump(tmp_path))
    gt = {
        "KA": {"ncbi_taxid": 100, "is_known": True},   # species 100, genus 20, family 10
        "DA": {"ncbi_taxid": 103, "is_known": False},  # dark, excluded
    }
    a = {
        "KA_0_0/1": 100,   # correct species/genus/family
        "KA_1_0/1": 101,   # wrong species, same genus/family
        "KA_2_0/1": 102,   # wrong species/genus, same family
        "KA_3_0/1": 103,   # wrong at all ranks
        "KA_4_0/1": 20,    # genus-level call: species too shallow, genus/family correct
        "KA_5_0/1": None,  # unclassified everywhere
        "DA_0_0/1": None,  # dark, must be excluded
    }
    r = benchmark_taxonomy(a, gt, ncbi_tree=tree)
    pr = r["per_rank"]

    assert (pr["species"]["correct"], pr["species"]["misclassified"],
            pr["species"]["unclassified"], pr["species"]["n"]) == (1, 3, 2, 6)
    assert pr["species"]["precision"] == pytest.approx(0.25)
    assert pr["species"]["recall"] == pytest.approx(1 / 6)

    assert (pr["genus"]["correct"], pr["genus"]["misclassified"],
            pr["genus"]["unclassified"]) == (3, 2, 1)
    assert pr["genus"]["precision"] == pytest.approx(0.6)
    assert pr["genus"]["recall"] == pytest.approx(0.5)

    assert (pr["family"]["correct"], pr["family"]["misclassified"],
            pr["family"]["unclassified"]) == (4, 1, 1)
    assert pr["family"]["precision"] == pytest.approx(0.8)
    assert pr["family"]["recall"] == pytest.approx(4 / 6)


def test_per_rank_absent_without_tree():
    r = benchmark_taxonomy(_assignments(), TAX_GT)
    assert r["per_rank"] is None
