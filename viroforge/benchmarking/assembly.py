"""Assembly benchmarking metrics engine (Module 2).

Aligns a user's assembled contigs to the true ViroForge genomes (via minimap2
through mappy) and reports per-genome recovery/completeness/identity, chimeric
contigs, assembly contiguity (N50/L50), and, when per-contig coverage is
available, abundance-estimation accuracy against the true relative abundances.

Ground truth: the source genome FASTA plus the per-genome expected values in the
metadata (`benchmarking.expected_coverage.per_genome`), which already account for
what is recoverable at the dataset's coverage (Lander-Waterman).

Scope (v1): genome recovery, chimera detection, assembly stats, and abundance
accuracy. Coverage-uniformity/GC-bias plots and HTML visualizations are deferred.
"""

from __future__ import annotations

import re
from collections import defaultdict

from .align import align_contigs, merge_len, read_fasta

# genome recovery category cutoffs (fraction of the genome covered by contigs)
COMPLETE = 0.95
HIGH = 0.75
PARTIAL = 0.50


def _assembly_stats(lengths: list[int]) -> dict:
    if not lengths:
        return {"n_contigs": 0, "total_bp": 0, "longest": 0, "n50": 0, "l50": 0}
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2
    acc = 0
    n50 = l50 = 0
    for i, x in enumerate(lengths, 1):
        acc += x
        if acc >= half:
            n50, l50 = x, i
            break
    return {"n_contigs": len(lengths), "total_bp": total,
            "longest": lengths[0], "n50": n50, "l50": l50}


def _spearman(xs, ys):
    if len(xs) < 3:
        return None
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


def load_ground_truth(genome_fasta, metadata: dict | None):
    """Return {genome_id: {length, type, abundance}} from metadata (preferred) or FASTA."""
    gt: dict[str, dict] = {}
    if metadata:
        for g in metadata.get("benchmarking", {}).get("expected_coverage", {}).get("per_genome", []):
            gt[g["genome_id"]] = {
                "length": g["length"],
                "type": g.get("sequence_type", "viral"),
                "abundance": g.get("relative_abundance"),
                "expected_completeness": g.get("expected_completeness"),
            }
    if not gt:  # fall back to the FASTA headers
        for name, header, seq in read_fasta(genome_fasta):
            mid = header.split("|")[1].strip() if "|" in header else ""
            typ = mid if mid in ("host_dna", "rrna", "phix", "reagent_bacteria") else "viral"
            m = re.search(r"abundance=([0-9.eE+-]+)", header)
            gt[name] = {"length": len(seq), "type": typ,
                        "abundance": float(m.group(1)) if m else None,
                        "expected_completeness": None}
    return gt


def benchmark_assembly(contigs_fasta, genome_fasta, metadata=None) -> dict:
    gt = load_ground_truth(genome_fasta, metadata)
    viral_ids = {g for g, v in gt.items() if v["type"] == "viral"}

    al = align_contigs(contigs_fasta, genome_fasta)
    contigs = al["contigs"]
    genome_intervals = al["genome_intervals"]
    genome_match_bp = al["genome_match_bp"]
    genome_block_bp = al["genome_block_bp"]

    contig_lengths = [c["length"] for c in contigs]
    contig_primary_genome = {c["name"]: c["primary_genome"]
                             for c in contigs if c["primary_genome"] is not None}
    n_unmatched = sum(1 for c in contigs if c["primary_genome"] is None)
    chimeras = [{"contig": c["name"], "length": c["length"], "segments": c["chimera_segments"]}
                for c in contigs if c["is_chimera"]]

    # per-genome completeness (viral only for recovery)
    per_genome = {}
    cats = {"complete": 0, "high_quality": 0, "partial": 0, "fragmented": 0, "missing": 0}
    completeness_vals = []
    for gid in viral_ids:
        length = gt[gid]["length"]
        covered = merge_len(genome_intervals.get(gid, []))
        comp = covered / length if length else 0.0
        ident = (genome_match_bp[gid] / genome_block_bp[gid]) if genome_block_bp[gid] else None
        if comp >= COMPLETE:
            cat = "complete"
        elif comp >= HIGH:
            cat = "high_quality"
        elif comp >= PARTIAL:
            cat = "partial"
        elif comp > 0:
            cat = "fragmented"
        else:
            cat = "missing"
        cats[cat] += 1
        completeness_vals.append(comp)
        per_genome[gid] = {
            "completeness": comp, "identity": ident, "category": cat,
            "expected_completeness": gt[gid].get("expected_completeness"),
        }

    n_viral = len(viral_ids)
    recovered = n_viral - cats["missing"]

    # abundance accuracy (viral genomes with per-contig coverage available)
    abundance = _abundance_accuracy(contigs, gt, viral_ids)

    # observed vs expected completeness
    exp_pairs = [(per_genome[g]["completeness"], gt[g]["expected_completeness"])
                 for g in viral_ids if gt[g].get("expected_completeness") is not None]
    completeness_vs_expected = None
    if exp_pairs:
        deltas = [o - e for o, e in exp_pairs]
        met = sum(1 for o, e in exp_pairs if o >= e - 0.05)
        completeness_vs_expected = {
            "n_scored": len(exp_pairs),
            "mean_delta": sum(deltas) / len(deltas),
            "n_met_expectation": met,
        }

    return {
        "assembly_stats": _assembly_stats(contig_lengths),
        "genome_recovery": {
            "n_viral_genomes": n_viral,
            "recovered": recovered,
            "recovery_rate": recovered / n_viral if n_viral else None,
            "categories": cats,
            "mean_completeness": (sum(completeness_vals) / len(completeness_vals)
                                  if completeness_vals else None),
        },
        "contigs": {
            "total": len(contig_lengths),
            "matched": len(contig_primary_genome),
            "unmatched": n_unmatched,
            "viral": sum(1 for g in contig_primary_genome.values() if g in viral_ids),
        },
        "chimeras": {"n": len(chimeras), "rate": (len(chimeras) / len(contig_lengths)
                     if contig_lengths else None), "examples": chimeras[:10]},
        "abundance_accuracy": abundance,
        "completeness_vs_expected": completeness_vs_expected,
        "per_genome": per_genome,
    }


def _abundance_accuracy(contigs, gt, viral_ids):
    viral_with_cov = [c for c in contigs
                      if c["coverage"] is not None and c["primary_genome"] in viral_ids]
    if len(viral_with_cov) < 3:
        return {"available": False, "reason": "per-contig coverage not found in headers"}
    est: dict[str, float] = defaultdict(float)
    for c in viral_with_cov:
        est[c["primary_genome"]] += c["coverage"] * c["length"]
    genomes = [g for g in viral_ids if gt[g]["abundance"] is not None]
    total_est = sum(est.values()) or 1.0
    total_true = sum(gt[g]["abundance"] for g in genomes) or 1.0
    xs = [est.get(g, 0.0) / total_est for g in genomes]
    ys = [gt[g]["abundance"] / total_true for g in genomes]
    mae = sum(abs(a - b) for a, b in zip(xs, ys)) / len(genomes) if genomes else None
    return {
        "available": True,
        "n_genomes": len(genomes),
        "spearman": _spearman(xs, ys),
        "mean_abs_error": mae,
    }
