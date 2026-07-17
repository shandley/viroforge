"""Independent-oracle tests for collection-specific contamination profiles.

The level acts as a multiplier on a per-collection baseline; expected fractions are
hand-computed from the baseline x multiplier (with clamping), not read back from
the tool. No network, no real reference files (use_real_references=False).
"""
import logging
import sqlite3
import sys
from pathlib import Path

import pytest

logging.disable(logging.CRITICAL)  # silence the module's INFO logging

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from viroforge.core.contamination import create_contamination_profile, LEVEL_MULTIPLIERS
import populate_contamination_defaults as pcd

BLOOD = {"collection_name": "Blood", "default_host_pct": 40.0, "default_rrna_pct": 0.5,
         "default_reagent_pct": 1.0, "default_phix_pct": 0.1, "host_organism": "human"}


def _host_frac(profile):
    return sum(c.abundance or 0 for c in profile.contaminants
               if c.contaminant_type.value == "host_dna")


def _profile(level, **kw):
    return create_contamination_profile(level, use_real_references=False, random_seed=1, **kw)


def test_level_scales_baseline():
    # baseline host 40% -> clean 0.25x=10%, realistic 1.0x=40%
    assert _host_frac(_profile("clean", collection_defaults=BLOOD)) == pytest.approx(0.10, abs=0.01)
    assert _host_frac(_profile("realistic", collection_defaults=BLOOD)) == pytest.approx(0.40, abs=0.01)


def test_high_multiplier_clamped():
    # heavy 2.5x=100% and failed 4.0x=160% both clamp to the 95% ceiling
    assert _host_frac(_profile("heavy", collection_defaults=BLOOD)) == pytest.approx(0.95, abs=0.01)
    assert _host_frac(_profile("failed", collection_defaults=BLOOD)) == pytest.approx(0.95, abs=0.01)


def test_multipliers_are_ordered():
    assert LEVEL_MULTIPLIERS["clean"] < LEVEL_MULTIPLIERS["realistic"] < LEVEL_MULTIPLIERS["heavy"] < LEVEL_MULTIPLIERS["failed"]


def test_null_column_falls_back_to_preset():
    # only host set; realistic multiplier is 1.0 so host stays at its baseline
    partial = {"collection_name": "X", "default_host_pct": 5.0}  # other columns NULL/absent
    assert _host_frac(_profile("realistic", collection_defaults=partial)) == pytest.approx(0.05, abs=0.01)


def test_backwards_compatible_without_defaults():
    # no collection_defaults -> original fixed preset (realistic host = 2%)
    assert _host_frac(_profile("realistic")) == pytest.approx(0.02, abs=0.01)
    # and NOT scaled: heavy preset host = 10% (a fixed value, not 2% x multiplier)
    assert _host_frac(_profile("heavy")) == pytest.approx(0.10, abs=0.02)


def test_host_free_collection_stays_near_zero():
    marine = {"collection_name": "Marine", "default_host_pct": 0.05, "default_rrna_pct": 5.0,
              "host_organism": "none"}
    assert _host_frac(_profile("realistic", collection_defaults=marine)) < 0.005


# ---- populate script -------------------------------------------------------- #
def test_load_defaults_parses_tsv():
    rows = pcd.load_defaults(REPO / "data/reference_profiles/contamination_defaults.tsv")
    assert len(rows) == 20
    blood = next(r for r in rows if r["collection_id"] == 17)
    assert blood["default_host_pct"] == 40.0 and blood["host_organism"] == "human"
    marine = next(r for r in rows if r["collection_id"] == 5)
    assert marine["host_organism"] == "none" and marine["default_host_pct"] == 0.05


def test_ensure_columns_adds_missing(tmp_path):
    db = tmp_path / "t.db"
    conn = sqlite3.connect(str(db))
    conn.execute("CREATE TABLE body_site_collections(collection_id INTEGER PRIMARY KEY, collection_name TEXT)")
    added = pcd.ensure_columns(conn)
    assert set(added) == set(pcd.COLUMNS)
    # idempotent: second call adds nothing
    assert pcd.ensure_columns(conn) == []
    cols = {r[1] for r in conn.execute("PRAGMA table_info(body_site_collections)")}
    assert "default_host_pct" in cols and "host_organism" in cols
    conn.close()
