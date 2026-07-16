"""Render QC benchmark metrics as JSON and a human-readable summary.

The summary leads with viral retention and over-filtering, the failure mode that
matters most for low-biomass viromes, then per-type removal rates, then dedup and
the aggregate confusion matrix.
"""

from __future__ import annotations

import json
from pathlib import Path


def _pct(x) -> str:
    return "n/a" if x is None else f"{x * 100:.1f}%"


def to_markdown(m: dict) -> str:
    lines = ["# QC Benchmark", ""]

    if not m.get("reliable", True):
        lines += [
            f"> WARNING: only {_pct(m['match_rate'])} of cleaned reads mapped back "
            "to the raw reads (need 99%). The cleaned reads may be renamed, "
            "reordered, or paired with the wrong raw file. Metrics below are "
            "unreliable.",
            "",
        ]

    lines += [
        f"Raw reads: {m['n_raw_reads']:,}   Cleaned reads: {m['n_cleaned_reads']:,}   "
        f"Name match: {_pct(m['match_rate'])}"
        + ("  (matched after dropping mate suffix)" if m.get("match_rate_after_mate_strip") else ""),
        "",
        "## Viral retention (critical)",
        "",
        f"- Viral retention: **{_pct(m['viral_retention'])}**",
        f"- Over-filtering (viral reads wrongly removed): **{_pct(m['over_filtering_rate'])}**",
        "",
        "## Contamination removal by type",
        "",
        "| type | reads | removal rate |",
        "|---|---:|---:|",
    ]
    for src, d in sorted(m["removal_rate_by_type"].items()):
        lines.append(f"| {src} | {d['n']:,} | {_pct(d['removal_rate'])} |")

    c = m["contamination"]
    dd = m["dedup"]
    un = m["unknown_source"]
    lines += [
        "",
        "## Aggregate (non-duplicate reads, positive class = should-be-removed)",
        "",
        f"- Precision: {_pct(c['precision'])}   Recall: {_pct(c['recall'])}   F1: {_pct(c['f1'])}",
        f"- TP {c['tp']:,}  FP {c['fp']:,}  FN {c['fn']:,}  TN {c['tn']:,}",
        "",
        "## Deduplication (orthogonal to contamination)",
        "",
        f"- PCR-duplicate reads: {dd['n_duplicates']:,}   removed: {_pct(dd['dedup_rate'])}",
    ]
    if un["n"]:
        lines += [
            "",
            f"## Unclassified: {un['n']:,} reads with an unknown source label "
            f"(removal {_pct(un['removal_rate'])})",
        ]
    return "\n".join(lines) + "\n"


def write_reports(m: dict, json_path=None, md_path=None) -> None:
    if json_path:
        Path(json_path).write_text(json.dumps(m, indent=2))
    if md_path:
        Path(md_path).write_text(to_markdown(m))
