"""ViroForge benchmarking framework.

Modules that validate virome analysis pipelines against ViroForge ground truth.
Module 1 (QC) validates contamination removal and quality filtering.
"""

from .qc import DEFAULT_KEEP_REMOVE, benchmark_qc
from .parsers import read_labels, read_names

__all__ = ["benchmark_qc", "DEFAULT_KEEP_REMOVE", "read_labels", "read_names"]
