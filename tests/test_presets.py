"""Tests for the bundled generation presets.

Every preset shipped a stale collection_id until 0.19.0: the 9-28 -> 1-20
renumbering in PR #6 never updated them, so `--preset gut-standard` generated a
WASTEWATER virome and `--preset respiratory-rna` pointed at collection 21, which
does not exist. Nothing caught it because no test loaded a preset and checked
where it pointed.
"""

import sqlite3
from pathlib import Path

import pytest
import yaml

PRESET_DIR = Path(__file__).parent.parent / "viroforge" / "presets"
DB = Path(__file__).parent.parent / "viroforge" / "data" / "viral_genomes.db"

PRESETS = sorted(PRESET_DIR.glob("*.yaml"))


def _collections():
    con = sqlite3.connect(DB)
    try:
        return dict(con.execute(
            "SELECT collection_id, collection_name FROM body_site_collections"))
    finally:
        con.close()


def _load(path):
    return yaml.safe_load(path.read_text())


def _requires_script(preset):
    """Presets that configure a standalone script, not `viroforge generate`."""
    return (preset.get("metadata") or {}).get("requires_script")


def test_presets_exist():
    assert PRESETS, "no presets found"


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_preset_is_wellformed(path):
    preset = _load(path)
    assert preset["name"] == path.stem, "name must match the filename"
    assert preset.get("description")
    assert isinstance(preset.get("parameters"), dict)


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_collection_id_exists(path):
    """The renumbering bug: an id that no longer resolves."""
    cid = _load(path)["parameters"].get("collection_id")
    assert cid in _collections(), (
        f"{path.stem} points at collection {cid}, which is not in the database"
    )


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_preset_name_matches_its_collection(path):
    """The subtler half: an id that resolves, to the wrong thing.

    gut-standard pointed at collection 9, which exists but is wastewater, so it
    generated successfully and silently produced the wrong dataset.
    """
    preset = _load(path)
    cid = preset["parameters"]["collection_id"]
    actual = _collections()[cid].lower()

    # Only checked for presets whose name names a body site.
    site_words = {
        "gut": "gut", "marine": "marine", "respiratory": "respiratory",
        "stool": "gut", "oral": "oral", "skin": "skin", "soil": "soil",
    }
    for word, expected in site_words.items():
        if word in path.stem:
            assert expected in actual, (
                f"{path.stem} points at '{actual}', which does not look like a "
                f"'{expected}' collection"
            )
            return


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_parameters_are_real_generator_arguments(path):
    """A typo'd key would otherwise be dropped or crash at generation time."""
    if _requires_script(_load(path)):
        pytest.skip("configures a standalone script, not `viroforge generate`")
    from viroforge.generator import build_parser

    known = vars(build_parser().parse_args([]))
    unknown = sorted(k for k in _load(path)["parameters"] if k not in known)
    assert not unknown, (
        f"{path.stem} sets {unknown}, which build_parser() does not define"
    )


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_preset_converts_to_a_namespace(path):
    """End to end: the preset must survive _params_to_namespace()."""
    if _requires_script(_load(path)):
        pytest.skip("configures a standalone script, not `viroforge generate`")
    from viroforge.cli.generate import _params_to_namespace

    params = dict(_load(path)["parameters"])
    params["output"] = "/tmp/unused"
    args = _params_to_namespace(params)
    assert args.collection_id == _load(path)["parameters"]["collection_id"]


class TestSampleTypePresets:
    """The trio from issue #37: bulk, VLP-enriched, and host-dominated."""

    def test_all_three_exist(self):
        names = {p.stem for p in PRESETS}
        assert {"gut-bulk", "gut-vlp", "stool-clinical"} <= names

    def test_bulk_and_vlp_differ_in_enrichment(self):
        bulk = _load(PRESET_DIR / "gut-bulk.yaml")["parameters"]
        vlp = _load(PRESET_DIR / "gut-vlp.yaml")["parameters"]
        assert bulk.get("no_vlp") is True
        assert vlp.get("vlp_protocol") and not vlp.get("no_vlp")
        assert bulk["collection_id"] == vlp["collection_id"], (
            "the pair is meant to isolate enrichment, so the collection must match"
        )

    def test_stool_clinical_pins_host_absolutely(self):
        """40% host is unreachable from the gut collection's 5% baseline via
        --contamination-level, so it needs the absolute override."""
        params = _load(PRESET_DIR / "stool-clinical.yaml")["parameters"]
        assert params["host_fraction"] == pytest.approx(0.40)
        assert params.get("no_vlp") is True


@pytest.mark.parametrize("path", PRESETS, ids=lambda p: p.stem)
def test_script_only_presets_fail_clearly(path):
    """hybrid-standard cannot run through `viroforge generate`. It must say so,
    not die on parameters build_parser() has never heard of."""
    preset = _load(path)
    script = _requires_script(preset)
    if not script:
        pytest.skip("runnable through `viroforge generate`")

    from viroforge.cli.generate import generate_with_preset

    class Args:
        preset = path.stem
        output = None
        seed = None
        coverage = None
        depth = None
        verbose = False

    assert generate_with_preset(Args()) == 1
    assert Path(script).exists(), f"{path.stem} points at a missing script"
