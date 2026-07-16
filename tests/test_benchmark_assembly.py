"""Assembly benchmark tests with mock contigs of known recovery.

Synthetic genomes are built, then contigs are carved from them with fixed,
known completeness (complete / partial / fragmented / missing) plus a deliberate
two-genome chimera. The oracle is the construction here; the benchmark must
recover those figures.
"""

import random

import pytest

from viroforge.benchmarking.align import parse_coverage
from viroforge.benchmarking.assembly import benchmark_assembly

GLEN = 10000


def _seq(seed, n=GLEN):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def genomes():
    # viral A-F plus a host contaminant
    return {
        "A": _seq(1), "B": _seq(2), "C": _seq(3),
        "D": _seq(4), "E": _seq(5), "F": _seq(6),
        "HOST": _seq(7, 5000),
    }


def _write_fasta(path, records):
    with open(path, "w") as f:
        for name, header, seq in records:
            f.write(f">{header}\n{seq}\n")


def _metadata(genomes):
    types = {"HOST": "host_dna"}
    ab = {"A": 0.35, "B": 0.20, "C": 0.12, "D": 0.08, "E": 0.10, "F": 0.15, "HOST": 0.001}
    per = []
    for g, s in genomes.items():
        per.append({
            "genome_id": g, "sequence_type": types.get(g, "viral"),
            "length": len(s), "relative_abundance": ab[g],
            "expected_coverage": 30.0,
            "expected_completeness": 1.0 if g != "D" else 0.4,
        })
    return {"benchmarking": {"expected_coverage": {"per_genome": per}}}


@pytest.fixture
def dataset(tmp_path, genomes):
    g = genomes
    gfa = tmp_path / "genomes.fasta"
    _write_fasta(gfa, [(n, f"{n} | {'host_dna' if n == 'HOST' else 'Unknown'} | abundance=0.1", s)
                       for n, s in g.items()])
    # contigs with known recovery
    contigs = [
        ("cA", "NODE_1_length_10000_cov_100.0", g["A"]),          # A complete (1.0)
        ("cB", "NODE_2_length_6000_cov_40.0", g["B"][:6000]),     # B partial (0.60)
        ("cC", "NODE_3_length_3000_cov_20.0", g["C"][:3000]),     # C fragmented (0.30)
        ("cChim", "NODE_4_length_6000_cov_30.0", g["E"][:3000] + g["F"][:3000]),  # chimera
        ("cHost", "NODE_5_length_5000_cov_10.0", g["HOST"]),      # contaminant
        # D: no contig -> missing
    ]
    cfa = tmp_path / "contigs.fasta"
    _write_fasta(cfa, contigs)
    return cfa, gfa, _metadata(g)


def test_parse_coverage():
    assert parse_coverage("NODE_1_length_10000_cov_45.2") == pytest.approx(45.2)
    assert parse_coverage("k141_1 flag=1 multi=12.0 len=5000") == pytest.approx(12.0)
    assert parse_coverage("plain_contig_name") is None


def test_recovery_categories(dataset):
    cfa, gfa, meta = dataset
    r = benchmark_assembly(cfa, gfa, meta)

    pg = r["per_genome"]
    assert pg["A"]["completeness"] == pytest.approx(1.0, abs=0.02)
    assert pg["A"]["category"] == "complete"
    assert pg["B"]["completeness"] == pytest.approx(0.60, abs=0.02)
    assert pg["B"]["category"] == "partial"
    assert pg["C"]["completeness"] == pytest.approx(0.30, abs=0.02)
    assert pg["C"]["category"] == "fragmented"
    assert pg["D"]["completeness"] == 0.0
    assert pg["D"]["category"] == "missing"
    # E and F each recovered ~0.30 via the chimera
    assert pg["E"]["completeness"] == pytest.approx(0.30, abs=0.02)
    assert pg["F"]["completeness"] == pytest.approx(0.30, abs=0.02)

    gr = r["genome_recovery"]
    assert gr["n_viral_genomes"] == 6
    assert gr["recovered"] == 5
    assert gr["recovery_rate"] == pytest.approx(5 / 6)
    assert gr["categories"] == {
        "complete": 1, "high_quality": 0, "partial": 1,
        "fragmented": 3, "missing": 1,
    }
    # identity should be ~1.0 for exact carve-outs
    assert pg["A"]["identity"] == pytest.approx(1.0, abs=0.02)


def test_chimera_detected(dataset):
    cfa, gfa, meta = dataset
    r = benchmark_assembly(cfa, gfa, meta)
    assert r["chimeras"]["n"] == 1
    ex = r["chimeras"]["examples"][0]
    assert ex["contig"].startswith("NODE_4")  # the chimeric contig's header name
    assert ex["length"] == 6000
    assert {s["genome"] for s in ex["segments"]} == {"E", "F"}


def test_contig_classification(dataset):
    cfa, gfa, meta = dataset
    r = benchmark_assembly(cfa, gfa, meta)
    c = r["contigs"]
    assert c["total"] == 5
    assert c["matched"] == 5
    assert c["unmatched"] == 0
    assert c["viral"] == 4  # cA, cB, cC, cChim; cHost is contaminant


def test_assembly_stats(dataset):
    cfa, gfa, meta = dataset
    r = benchmark_assembly(cfa, gfa, meta)
    s = r["assembly_stats"]
    assert s["n_contigs"] == 5
    assert s["longest"] == 10000
    assert s["total_bp"] == 10000 + 6000 + 3000 + 6000 + 5000


def test_abundance_accuracy_available(dataset):
    cfa, gfa, meta = dataset
    r = benchmark_assembly(cfa, gfa, meta)
    a = r["abundance_accuracy"]
    assert a["available"] is True
    assert a["n_genomes"] == 6
    assert a["spearman"] is not None


def test_unmatched_contig(dataset, tmp_path, genomes):
    cfa, gfa, meta = dataset
    # a contig of pure novel sequence should be unmatched
    novel = _seq(999, 4000)
    with open(cfa, "a") as f:
        f.write(f">NODE_9_length_4000_cov_5.0\n{novel}\n")
    r = benchmark_assembly(cfa, gfa, meta)
    assert r["contigs"]["unmatched"] == 1
