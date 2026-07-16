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


def assembly_to_markdown(m: dict) -> str:
    gr = m["genome_recovery"]
    st = m["assembly_stats"]
    c = m["contigs"]
    ch = m["chimeras"]
    lines = [
        "# Assembly Benchmark",
        "",
        "## Genome recovery",
        "",
        f"- Viral genomes recovered: **{gr['recovered']}/{gr['n_viral_genomes']}** "
        f"({_pct(gr['recovery_rate'])})",
        f"- Mean completeness: {_pct(gr['mean_completeness'])}",
        "",
        "| category | genomes |",
        "|---|---:|",
    ]
    for cat in ("complete", "high_quality", "partial", "fragmented", "missing"):
        lines.append(f"| {cat} | {gr['categories'][cat]:,} |")

    cve = m.get("completeness_vs_expected")
    if cve:
        lines += [
            "",
            f"Observed vs expected completeness: {cve['n_met_expectation']}/"
            f"{cve['n_scored']} genomes met expectation "
            f"(mean delta {cve['mean_delta']:+.3f}).",
        ]

    lines += [
        "",
        "## Contigs",
        "",
        f"- Total {c['total']:,}   matched {c['matched']:,}   viral {c['viral']:,}   "
        f"unmatched {c['unmatched']:,}",
        f"- Chimeric contigs: **{ch['n']:,}** ({_pct(ch['rate'])})",
        "",
        "## Contiguity",
        "",
        f"- Contigs {st['n_contigs']:,}   total {st['total_bp']:,} bp   "
        f"longest {st['longest']:,}   N50 {st['n50']:,}   L50 {st['l50']:,}",
    ]

    ab = m.get("abundance_accuracy", {})
    if ab.get("available"):
        lines += [
            "",
            "## Abundance accuracy (estimated vs true relative abundance)",
            "",
            f"- Genomes scored: {ab['n_genomes']:,}   Spearman: "
            f"{ab['spearman']:.3f}   mean abs error: {ab['mean_abs_error']:.4f}"
            if ab["spearman"] is not None else "- (insufficient genomes to correlate)",
        ]
    else:
        lines += [
            "",
            f"## Abundance accuracy: not computed ({ab.get('reason', 'no coverage')})",
        ]
    return "\n".join(lines) + "\n"


def write_reports(m: dict, json_path=None, md_path=None, kind="qc") -> None:
    render = assembly_to_markdown if kind == "assembly" else to_markdown
    if json_path:
        Path(json_path).write_text(json.dumps(m, indent=2))
    if md_path:
        Path(md_path).write_text(render(m))
