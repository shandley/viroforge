"""Tests for the in-process generation path used by `viroforge generate`.

`viroforge/cli/generate.py` calls `run_generation()` directly rather than
shelling out to `scripts/generate_fastq_dataset.py`, so it must hand over a
Namespace carrying every argument the generator reads. An earlier version
restated the parser's defaults in a second dict, which drifted to 6 missing
attributes and 11 wrong values and made `viroforge generate` raise
AttributeError. These tests pin the contract instead of the values.
"""

import ast
from pathlib import Path

import pytest

from viroforge.cli.generate import _params_to_namespace
from viroforge.generator import build_parser

GENERATOR_SRC = Path(__file__).parent.parent / "viroforge" / "generator.py"


def _args_read_by_generator() -> set[str]:
    """Every `args.<name>` the generator module reads."""
    tree = ast.parse(GENERATOR_SRC.read_text())
    return {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "args"
        and isinstance(node.ctx, ast.Load)
    }


def test_namespace_supplies_every_arg_the_generator_reads():
    ns = _params_to_namespace({})
    missing = sorted(a for a in _args_read_by_generator() if not hasattr(ns, a))
    assert not missing, (
        f"run_generation() reads {missing} but the CLI namespace does not "
        "provide them; `viroforge generate` would raise AttributeError"
    )


def test_defaults_come_from_the_real_parser():
    """No second copy of the defaults: the namespace must match parse_args([])."""
    parser_defaults = vars(build_parser().parse_args([]))
    ns = vars(_params_to_namespace({}))
    mismatched = {
        k: (ns[k], v) for k, v in parser_defaults.items()
        if k in ns and ns[k] != v
    }
    assert not mismatched, f"CLI defaults drifted from the parser: {mismatched}"


def test_params_override_defaults():
    ns = _params_to_namespace({"collection_id": 8, "n_reads": 500, "platform": "miseq"})
    assert (ns.collection_id, ns.n_reads, ns.platform) == (8, 500, "miseq")


def test_no_vlp_maps_to_protocol_none():
    assert _params_to_namespace({"no_vlp": True}).vlp_protocol == "none"


def test_unknown_parameter_is_rejected():
    """A typo used to be dropped silently, generating with the wrong settings."""
    with pytest.raises(ValueError, match="unknown generation parameter"):
        _params_to_namespace({"collection_id": 8, "adaptor_rate": 0.03})
