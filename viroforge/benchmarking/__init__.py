"""ViroForge benchmarking framework.

Modules that validate virome analysis pipelines against ViroForge ground truth.
Module 1 (QC) validates contamination removal and quality filtering.
Module 8 (Discovery) validates novel virus detection.
"""

from .qc import DEFAULT_KEEP_REMOVE, benchmark_qc
from .discovery import benchmark_discovery
from .parsers import read_labels, read_names

__all__ = [
    "benchmark_qc",
    "benchmark_discovery",
    "DEFAULT_KEEP_REMOVE",
    "read_labels",
    "read_names",
]
