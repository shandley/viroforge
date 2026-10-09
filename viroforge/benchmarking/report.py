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


def _num(x, nd=3):
    return "n/a" if x is None else f"{x:.{nd}f}"


def taxonomy_to_markdown(m: dict) -> str:
    lines = ["# Taxonomy Benchmark (read-based, taxid-exact)", ""]
    if not m.get("reliable", True):
        lines += [
            "> WARNING: no read IDs mapped to known viral genomes. Check that the "
            "classifier output uses the ViroForge read names and that the "
            "ground-truth metadata matches this dataset.",
            "",
        ]
    lines += [
        f"Reads: {m['n_assignments']:,}   viral {m['n_viral_reads']:,}   "
        f"non-viral {m['n_non_viral_reads']:,}",
        "",
        "## Known viruses (ICTV family assigned)",
        "",
    ]
    k = m["known_viruses"]
    lines += [
        f"- Reads: {k['n']:,}",
        f"- Sensitivity (correct / all): **{_pct(k['sensitivity'])}**",
        f"- Precision (correct / classified): {_pct(k['precision'])}",
        f"- Left unclassified: {_pct(k['unclassified_rate'])}   "
        f"misclassified: {k['misclassified']:,}",
        "",
        "## Dark matter (family Unknown / novel)",
        "",
    ]
    d = m["dark_matter"]
    lines += [
        f"- Reads: {d['n']:,}",
        f"- Correctly left unclassified: **{_pct(d['unclassified_rate'])}** "
        "(expected high; classifiers should not force a call on novel content)",
        f"- Classified anyway: {d['correct'] + d['misclassified']:,} "
        f"(of which exactly correct: {d['correct']:,})",
        "",
    ]
    pr = m.get("per_rank")
    if pr:
        lines += [
            "## Per-rank accuracy (known viruses, via NCBI lineage)",
            "",
            "| rank | reads | precision | recall | F1 |",
            "|---|---:|---:|---:|---:|",
        ]
        for rank in ("species", "genus", "family"):
            if rank in pr:
                s = pr[rank]
                lines.append(
                    f"| {rank} | {s['n']:,} | {_pct(s['precision'])} | "
                    f"{_pct(s['recall'])} | {_pct(s['f1'])} |")
        lines.append("")
    lines += [
        "## Abundance profile (classifier vs true, over NCBI taxids)",
        "",
    ]
    ab = m["abundance_profile"]
    lines += [
        f"- Taxa: {ab['n_taxa']:,}   Bray-Curtis dissimilarity: {_num(ab['bray_curtis'])}",
        f"- Pearson: {_num(ab['pearson'])}   Spearman: {_num(ab['spearman'])}   "
        f"mean abs error: {_num(ab['mean_abs_error'], 4)}",
    ]
    prov = m.get("provenance")
    if prov:
        lines += [
            "",
            "## By genome provenance",
            "",
            "| provenance | reads | sensitivity | precision |",
            "|---|---:|---:|---:|",
        ]
        for cat, s in prov.items():
            k = s.get("known_viruses", {})
            lines.append(
                f"| {cat} | {s['n_reads']:,} | "
                f"{_pct(k.get('sensitivity'))} | {_pct(k.get('precision'))} |"
            )
        lines += [
            "",
            "> Prophage and endogenous reads originate from integrated sequences. "
            "A classifier that calls them 'viral' is sequence-correct but "
            "context-wrong in a bulk metagenome. Use VIBRANT, CheckV, or "
            "RepeatMasker to distinguish integrated from free viral particles.",
        ]
    return "\n".join(lines) + "\n"


def _strata_and_rank_md(m: dict) -> list:
    """Shared known/dark/per-rank/abundance rendering for taxonomy reports."""
    lines = ["## Known viruses (ICTV family assigned)", ""]
    k = m["known_viruses"]
    lines += [
        f"- Units: {k['n']:,}",
        f"- Sensitivity (correct / all): **{_pct(k['sensitivity'])}**",
        f"- Precision (correct / classified): {_pct(k['precision'])}",
        f"- Left unclassified: {_pct(k['unclassified_rate'])}   "
        f"misclassified: {k['misclassified']:,}",
        "",
        "## Dark matter (family Unknown / novel)",
        "",
    ]
    d = m["dark_matter"]
    lines += [
        f"- Units: {d['n']:,}",
        f"- Correctly left unclassified: **{_pct(d['unclassified_rate'])}** "
        "(expected high; classifiers should not force a call on novel content)",
        "",
    ]
    pr = m.get("per_rank")
    if pr:
        lines += [
            "## Per-rank accuracy (known viruses, via NCBI lineage)",
            "",
            "| rank | units | precision | recall | F1 |",
            "|---|---:|---:|---:|---:|",
        ]
        for rank in ("species", "genus", "family"):
            if rank in pr:
                s = pr[rank]
                lines.append(f"| {rank} | {s['n']:,} | {_pct(s['precision'])} | "
                             f"{_pct(s['recall'])} | {_pct(s['f1'])} |")
        lines.append("")
    ab = m["abundance_profile"]
    lines += [
        "## Abundance profile (classifier vs true, over NCBI taxids)",
        "",
        f"- Taxa: {ab['n_taxa']:,}   Bray-Curtis dissimilarity: {_num(ab['bray_curtis'])}",
        f"- Pearson: {_num(ab['pearson'])}   Spearman: {_num(ab['spearman'])}   "
        f"mean abs error: {_num(ab['mean_abs_error'], 4)}",
    ]
    return lines


def taxonomy_contig_to_markdown(m: dict) -> str:
    lines = ["# Taxonomy Benchmark (contig-based, taxid-exact)", ""]
    if not m.get("reliable", True):
        lines += ["> WARNING: no contigs mapped to known viral genomes.", ""]
    lines += [
        f"Contigs: {m['n_contigs']:,}   viral {m['n_viral_contigs']:,}   "
        f"novel {m['n_novel_contigs']:,}   contaminant {m['n_contaminant_contigs']:,}   "
        f"chimeric {m['n_chimeric_contigs']:,} ({m['chimera_handling']})",
        "",
    ]
    lines += _strata_and_rank_md(m)
    return "\n".join(lines) + "\n"


def discovery_to_markdown(m: dict) -> str:
    lines = ["# Novel Discovery Benchmark", ""]
    if not m.get("reliable", True):
        lines += [
            "> WARNING: no read IDs mapped to known viral genomes. Check that the "
            "classifier output uses the ViroForge read names and that the "
            "ground-truth metadata matches this dataset.",
            "",
        ]
    lines += [
        f"Reads: {m['n_assignments']:,}   viral {m['n_viral_reads']:,}   "
        f"non-viral {m['n_non_viral_reads']:,}",
        "",
        "## Known virus detection (baseline)",
        "",
        f"- Known reads: {m['n_known_reads']:,}   detected: {m['n_known_detected']:,}   "
        f"rate: **{_pct(m['known_detection_rate'])}**",
        "",
        "## Dark matter detection (the novel-discovery question)",
        "",
        f"- Dark matter reads: {m['n_dark_matter_reads']:,}   detected: {m['n_dark_detected']:,}   "
        f"rate: **{_pct(m['dark_matter_detection_rate'])}**",
        f"- Dark matter genomes: {m['n_dark_genomes']:,}   detected: {m['n_dark_genomes_detected']:,}   "
        f"rate: {_pct(m['genome_detection_rate'])}",
        "",
    ]
    bins = m.get("detection_by_depth", [])
    if bins:
        lines += [
            "## Detection by read depth (dark matter genomes)",
            "",
            "| reads/genome | genomes | detected | rate |",
            "|---|---:|---:|---:|",
        ]
        for b in bins:
            lines.append(
                f"| {b['bin']} | {b['n_genomes']:,} | {b['n_detected']:,} | "
                f"{_pct(b['detection_rate'])} |"
            )
        lines.append("")

    per_genome = m.get("per_genome", [])
    if per_genome:
        lines += ["## Per-genome detection (dark matter)", ""]
        detected = [g for g in per_genome if g["detection_rate"] > 0]
        missed = [g for g in per_genome if g["detection_rate"] == 0]
        if detected:
            lines += [
                "### Detected",
                "",
                "| genome | reads | detected | rate |",
                "|---|---:|---:|---:|",
            ]
            for g in sorted(detected, key=lambda x: -x["detection_rate"]):
                lines.append(
                    f"| {g['genome_id']} | {g['total_reads']:,} | "
                    f"{g['detected_reads']:,} | {_pct(g['detection_rate'])} |"
                )
            lines.append("")
        if missed:
            lines += [
                "### Missed (0% detection)",
                "",
                "| genome | reads |",
                "|---|---:|",
            ]
            for g in sorted(missed, key=lambda x: -x["total_reads"]):
                lines.append(f"| {g['genome_id']} | {g['total_reads']:,} |")
            lines.append("")
    return "\n".join(lines) + "\n"


def write_reports(m: dict, json_path=None, md_path=None, kind="qc") -> None:
    render = {"assembly": assembly_to_markdown,
              "taxonomy": taxonomy_to_markdown,
              "taxonomy_contig": taxonomy_contig_to_markdown,
              "discovery": discovery_to_markdown}.get(kind, to_markdown)
    if json_path:
        Path(json_path).write_text(json.dumps(m, indent=2))
    if md_path:
        Path(md_path).write_text(render(m))
