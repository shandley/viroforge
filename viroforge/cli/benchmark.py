"""`viroforge benchmark` command: validate analysis pipelines against ground truth.

v1 implements Module 1 (QC): contamination removal and viral retention.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ..benchmarking import benchmark_qc, read_labels, read_names
from ..benchmarking.report import to_markdown, write_reports


def run_benchmark(args) -> int:
    if args.benchmark_command == "qc":
        return _run_qc(args)
    print("Usage: viroforge benchmark qc --raw-reads R1 --cleaned-reads R1", file=sys.stderr)
    return 2


def _run_qc(args) -> int:
    raw_paths = [Path(p) for p in args.raw_reads]
    cleaned_paths = [Path(p) for p in args.cleaned_reads]
    for p in raw_paths + cleaned_paths:
        if not p.exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 2

    keep_remove = None
    if args.keep_remove:
        keep_remove = {}
        for item in args.keep_remove:
            src, _, decision = item.partition(":")
            if decision not in ("keep", "remove"):
                print(f"ERROR: --keep-remove expects source:keep|remove, got {item}", file=sys.stderr)
                return 2
            keep_remove[src] = decision

    raw = read_labels(raw_paths)
    kept = read_names(cleaned_paths)
    metrics = benchmark_qc(raw, kept, keep_remove=keep_remove)

    write_reports(metrics, json_path=args.output, md_path=args.markdown)
    print(to_markdown(metrics))

    if not metrics["reliable"]:
        return 1  # unreliable match rate: signal via exit code
    return 0
