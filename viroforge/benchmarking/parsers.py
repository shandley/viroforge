"""FASTQ parsing for QC benchmarking.

Reads ViroForge raw reads to recover per-read ground truth (the `source=` label
and the PCR-duplicate tag), and reads a QC tool's cleaned reads to recover the
set of surviving read names. Plain and gzipped FASTQ are both accepted.

Ground truth lives in the raw FASTQ headers (validated accurate in
docs/DATA_QUALITY_EVALUATION.md); the metadata JSON is an optional cross-check,
not a required input.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

_MATE = re.compile(r"/[12]$")


def _open(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def _as_paths(paths) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(p) for p in paths]


def strip_mate(name: str) -> str:
    """Drop a trailing /1 or /2 mate suffix. Leaves .version and _dupN intact."""
    return _MATE.sub("", name)


def read_labels(paths) -> dict[str, tuple[str, bool]]:
    """Raw R1 FASTQ(s) -> {read_name: (source_label, is_pcr_duplicate)}.

    read_name is the first whitespace token of the header (without the leading @).
    """
    out: dict[str, tuple[str, bool]] = {}
    for p in _as_paths(paths):
        with _open(p) as fh:
            for i, line in enumerate(fh):
                if i % 4 != 0:
                    continue
                toks = line[1:].rstrip("\n").split()
                if not toks:
                    continue
                name = toks[0]
                source = "unknown"
                dup = False
                for t in toks[1:]:
                    if t.startswith("source="):
                        # Low-complexity artifact reads carry both the original
                        # label and the artifact label (source=viral
                        # source=artifact_low_complexity); the last one is the
                        # read's current nature, which is what QC should act on.
                        source = t.split("=", 1)[1]
                    elif t.startswith("pcr_duplicate="):
                        dup = t.split("=", 1)[1] == "true"
                out[name] = (source, dup)
    return out


def read_names(paths) -> set[str]:
    """Cleaned R1 FASTQ(s) -> set of surviving read names (first token)."""
    out: set[str] = set()
    for p in _as_paths(paths):
        with _open(p) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    toks = line[1:].rstrip("\n").split()
                    if toks:
                        out.add(toks[0])
    return out
