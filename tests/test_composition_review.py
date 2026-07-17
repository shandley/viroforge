"""Independent-oracle tests for the biological-accuracy composition machinery.

The composition tools live in scripts/ (not the package), so they are imported via
sys.path. Every expected value here is hand-computed from the constructed inputs,
not read back from the tool, so the tests check the logic rather than restate it.
"""
import sqlite3
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import evaluate_composition as ec
import fix_collection_composition as fc
import enrich_taxonomy_from_ncbi as enr
import renormalize_abundances as rn
import build_composition_reference as bcr


# --------------------------------------------------------------------------- #
# observe_collection: host balance, name-based phage, unclassified, NA mix
# --------------------------------------------------------------------------- #
def _row(family, name, ab, cls, prev=1.0):
    return {"family": family, "genome_name": name, "relative_abundance": ab,
            "prevalence": prev, "class": cls}


PROPS = {
    "Suoliviridae": {"host_type": "phage", "na_class": "DNA"},
    "Microviridae": {"host_type": "phage", "na_class": "DNA"},
    "Papillomaviridae": {"host_type": "eukaryotic", "na_class": "DNA"},
}


def test_observe_host_balance_and_inference():
    rows = [
        _row("Suoliviridae", "Suoliviridae sp.", 0.4, "Caudoviricetes"),
        _row("Microviridae", "Microviridae sp.", 0.2, "Malgrandaviricetes"),
        _row("Papillomaviridae", "Human papillomavirus", 0.1, "Papovaviricetes"),
        _row("Unknown", "Escherichia phage T4", 0.2, "Caudoviricetes"),  # name-inferred phage, has class
        _row("Unknown", "Novel marine virus", 0.1, None),                # no class -> unclassified
    ]
    obs = ec.observe_collection(rows, PROPS)

    hb = obs["host_balance"]["abundance"]
    assert hb["phage"] == pytest.approx(0.8)      # 0.4 + 0.2 + 0.2 (inferred)
    assert hb["eukaryotic"] == pytest.approx(0.1)
    assert hb["unknown"] == pytest.approx(0.1)

    na = obs["na_class_mix"]["abundance"]
    assert na["DNA"] == pytest.approx(0.7)        # only mapped families get NA
    assert na["unknown"] == pytest.approx(0.3)

    assert obs["inferred_phage"]["abundance"] == pytest.approx(0.2)
    assert obs["unknown_family"]["abundance"] == pytest.approx(0.3)   # both family=Unknown rows
    assert obs["unclassified"]["abundance"] == pytest.approx(0.1)     # only the no-class row
    assert obs["unclassified"]["count"] == 1
    # top family by abundance is the dominant one
    assert obs["top_families"][0]["family"] == "Suoliviridae"


# --------------------------------------------------------------------------- #
# band status + severity + strictness dial
# --------------------------------------------------------------------------- #
def test_band_status_within_below_above():
    band = {"low": 0.4, "expected": 0.6, "high": 0.8}
    assert ec._band_status(0.5, band, 1.0)["status"] == "within"
    assert ec._band_status(0.3, band, 1.0)["status"] == "below"
    assert ec._band_status(0.9, band, 1.0)["status"] == "above"


def test_strictness_narrows_band():
    band = {"low": 0.4, "expected": 0.6, "high": 0.8}
    # strictness 2.0 halves each margin: effective band [0.5, 0.7]
    st = ec._band_status(0.45, band, 2.0)
    assert st["effective_band"] == [0.5, 0.7]
    assert st["status"] == "below"          # 0.45 < 0.5 now, was within at strictness 1.0
    assert ec._band_status(0.45, band, 1.0)["status"] == "within"


def test_severity_thresholds():
    assert ec._severity(0.0) == "ok"
    assert ec._severity(0.3) == "minor"
    assert ec._severity(1.0) == "moderate"
    assert ec._severity(2.0) == "major"


# --------------------------------------------------------------------------- #
# evaluate_site: metrics, signature roles, biology vs data-quality split
# --------------------------------------------------------------------------- #
def _obs_for_eval():
    rows = [
        _row("Suoliviridae", "Suoliviridae sp.", 0.4, "Caudoviricetes"),
        _row("Microviridae", "Microviridae sp.", 0.2, "Malgrandaviricetes"),
        _row("Papillomaviridae", "Human papillomavirus", 0.1, "Papovaviricetes"),
        _row("Unknown", "Escherichia phage T4", 0.2, "Caudoviricetes"),
        _row("Unknown", "Novel marine virus", 0.1, None),
    ]
    return ec.observe_collection(rows, PROPS)


def test_evaluate_biology_ok_dataquality_flag_split():
    obs = _obs_for_eval()
    expect = {
        # phage/(phage+euk) = 0.8/0.9 = 0.889 -> within [0.7, 1.0]
        "phage_fraction": {"low": 0.7, "expected": 0.9, "high": 1.0},
        # unclassified 0.1 > ceiling 0.05 -> data-quality flag, not biology
        "max_unclassified": {"high": 0.05},
        "signature_taxa": [{"family": "Suoliviridae", "role": "dominant"}],  # rank 1 -> ok
    }
    res = ec.evaluate_site(obs, expect, 1.0)
    assert res["site_verdict"] == "ok"                 # biology clean
    assert res["data_quality_verdict"] != "ok"         # unclassified flagged separately
    dq = [f for f in res["findings"] if f.get("kind") == "data_quality"]
    assert dq and dq[0]["metric"] == "unclassified_fraction"


def test_signature_dominant_requires_rank_one():
    obs = _obs_for_eval()  # top family is Suoliviridae; Papillomaviridae is not rank 1
    # dominant but not rank 1 -> present_not_dominant (moderate)
    res = ec.evaluate_site(obs, {"signature_taxa": [{"family": "Papillomaviridae", "role": "dominant"}]}, 1.0)
    sig = [f for f in res["findings"] if f["metric"].startswith("signature:")][0]
    assert sig["status"] == "present_not_dominant"
    assert res["site_verdict"] == "moderate"
    # missing signature -> major for dominant role
    res2 = ec.evaluate_site(obs, {"signature_taxa": [{"family": "Anelloviridae", "role": "dominant"}]}, 1.0)
    assert res2["site_verdict"] == "major"


def test_phage_fraction_below_band_flags_biology():
    # eukaryote-heavy community; expect high phage -> below band
    rows = [_row("Papillomaviridae", "HPV", 0.9, "Papovaviricetes"),
            _row("Suoliviridae", "phage", 0.1, "Caudoviricetes")]
    obs = ec.observe_collection(rows, PROPS)
    res = ec.evaluate_site(obs, {"phage_fraction": {"low": 0.8, "expected": 0.9, "high": 1.0}}, 1.0)
    pf = [f for f in res["findings"] if f["metric"] == "phage_fraction_of_classified"][0]
    assert pf["observed"] == pytest.approx(0.1)
    assert pf["status"] == "below"
    assert res["site_verdict"] != "ok"


# --------------------------------------------------------------------------- #
# property-map helpers (build_composition_reference)
# --------------------------------------------------------------------------- #
def test_host_type_classification():
    assert bcr._host_type("bacteria") == "phage"
    assert bcr._host_type("archaea") == "phage"
    assert bcr._host_type("vertebrates") == "eukaryotic"
    assert bcr._host_type("invertebrates, vertebrates") == "eukaryotic"
    assert bcr._host_type("soil (S)") == "unknown"
    assert bcr._host_type(None) == "unknown"


def test_na_class_and_coarse():
    assert bcr._na_class("dsDNA") == "DNA"
    assert bcr._na_class("ssRNA(+)") == "RNA"
    assert bcr._na_class("dsDNA-RT") == "RT"
    assert bcr._na_coarse("ssRNA(-); ssRNA(+/-)") == "ssRNA"   # mixed strand, coarse is unambiguous
    assert bcr._na_coarse("dsDNA-RT") == "dsDNA"               # RT folded to packaged form
    assert bcr._na_coarse("ssDNA(+/-)") == "ssDNA"


# --------------------------------------------------------------------------- #
# collection reweighting math (fix_collection_composition)
# --------------------------------------------------------------------------- #
def test_apply_family_shares_hits_targets_and_normalizes():
    rows = [("g1", "A", "n1", 0.5), ("g2", "A", "n2", 0.3), ("g3", "B", "n3", 0.2)]
    new = fc.apply_family_shares(rows, {"A": 0.6})   # B fills remainder 0.4
    assert sum(new.values()) == pytest.approx(1.0)
    assert new["g1"] + new["g2"] == pytest.approx(0.6)
    assert new["g3"] == pytest.approx(0.4)
    # within family A the split stays proportional to the original weights
    assert new["g1"] == pytest.approx(0.6 * 0.5 / 0.8)


def test_host_of_and_host_target():
    hmap = {"Suoliviridae": "phage", "Papillomaviridae": "eukaryotic"}
    assert fc.host_of("Suoliviridae", "x", hmap) == "phage"
    assert fc.host_of("Unknown", "Escherichia phage T4", hmap) == "phage"   # name fallback
    assert fc.host_of("Unknown", "novel virus", hmap) == "unknown"

    rows = [("g1", "Suoliviridae", "n", 0.3), ("g2", "Papillomaviridae", "n", 0.7)]
    new = fc.apply_host_target(rows, {"phage": 0.85, "eukaryotic": 0.15}, hmap)
    assert new["g1"] == pytest.approx(0.85)
    assert new["g2"] == pytest.approx(0.15)


# --------------------------------------------------------------------------- #
# NCBI efetch XML parsing (enrich_taxonomy_from_ncbi)
# --------------------------------------------------------------------------- #
EFETCH_XML = """<?xml version="1.0"?>
<TaxaSet>
  <Taxon>
    <TaxId>1229790</TaxId>
    <ScientificName>Pahexavirus P100A</ScientificName>
    <AkaTaxIds><TaxId>999999</TaxId></AkaTaxIds>
    <Rank>species</Rank>
    <LineageEx>
      <Taxon><TaxId>2731341</TaxId><ScientificName>Duplodnaviria</ScientificName><Rank>realm</Rank></Taxon>
      <Taxon><TaxId>2731619</TaxId><ScientificName>Caudoviricetes</ScientificName><Rank>class</Rank></Taxon>
      <Taxon><TaxId>3123456</TaxId><ScientificName>Pahexavirus</ScientificName><Rank>genus</Rank></Taxon>
    </LineageEx>
  </Taxon>
  <Taxon>
    <TaxId>1511883</TaxId>
    <ScientificName>Human bocavirus 4</ScientificName>
    <Rank>species</Rank>
    <LineageEx>
      <Taxon><TaxId>10239</TaxId><ScientificName>Monodnaviria</ScientificName><Rank>realm</Rank></Taxon>
      <Taxon><TaxId>2732090</TaxId><ScientificName>Parvoviridae</ScientificName><Rank>family</Rank></Taxon>
      <Taxon><TaxId>40119</TaxId><ScientificName>Bocaparvovirus</ScientificName><Rank>genus</Rank></Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>"""


def test_parse_efetch_phage_no_family_and_aka():
    out = enr._parse_efetch(EFETCH_XML)
    phage = out["1229790"]
    assert phage["class"] == "Caudoviricetes"
    assert phage["genus"] == "Pahexavirus"
    assert "family" not in phage                # ICTV assigns no family -> stays absent
    assert out["999999"] == phage               # AkaTaxIds map to the same lineage


def test_parse_efetch_virus_with_family():
    out = enr._parse_efetch(EFETCH_XML)
    boca = out["1511883"]
    assert boca["family"] == "Parvoviridae"
    assert boca["genus"] == "Bocaparvovirus"


# --------------------------------------------------------------------------- #
# renormalization audit against a synthetic DB
# --------------------------------------------------------------------------- #
def test_renormalize_audit_detects_off_sums(tmp_path):
    db = tmp_path / "t.db"
    conn = sqlite3.connect(str(db))
    conn.executescript(
        """CREATE TABLE body_site_collections(collection_id INTEGER PRIMARY KEY, collection_name TEXT);
           CREATE TABLE collection_genomes(collection_id INTEGER, genome_id TEXT, relative_abundance REAL);"""
    )
    conn.execute("INSERT INTO body_site_collections VALUES (1,'norm'), (2,'off')")
    conn.executemany("INSERT INTO collection_genomes VALUES (?,?,?)",
                     [(1, "a", 0.6), (1, "b", 0.4), (2, "c", 0.5), (2, "d", 0.3)])
    conn.commit()
    rows = {cid: (s, n) for cid, _name, s, n in rn.audit(conn)}
    conn.close()
    assert rows[1][0] == pytest.approx(1.0)   # collection 1 already normalized
    assert rows[2][0] == pytest.approx(0.8)   # collection 2 sums to 0.8
    assert rows[2][1] == 2


def test_renormalize_apply_end_to_end(tmp_path):
    """Exercise the real --apply path (as folded into setup-db) via subprocess."""
    import subprocess

    db = tmp_path / "t.db"
    conn = sqlite3.connect(str(db))
    conn.executescript(
        """CREATE TABLE body_site_collections(collection_id INTEGER PRIMARY KEY, collection_name TEXT);
           CREATE TABLE collection_genomes(collection_id INTEGER, genome_id TEXT, relative_abundance REAL);"""
    )
    conn.execute("INSERT INTO body_site_collections VALUES (1,'off')")
    conn.executemany("INSERT INTO collection_genomes VALUES (?,?,?)", [(1, "c", 0.5), (1, "d", 0.3)])
    conn.commit()
    conn.close()

    proc = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "renormalize_abundances.py"), "--db", str(db), "--apply"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    conn = sqlite3.connect(str(db))
    total = conn.execute("SELECT SUM(relative_abundance) FROM collection_genomes WHERE collection_id=1").fetchone()[0]
    a, b = conn.execute("SELECT relative_abundance FROM collection_genomes WHERE collection_id=1 ORDER BY genome_id").fetchall()
    conn.close()
    assert total == pytest.approx(1.0)          # rescaled to sum 1.0
    assert a[0] == pytest.approx(0.5 / 0.8)     # proportions preserved (0.625)
    assert b[0] == pytest.approx(0.3 / 0.8)     # 0.375
