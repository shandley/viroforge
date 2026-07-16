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

import mappy

# ---- fixed thresholds --------------------------------------------------------
ASM_PRESET = "asm5"          # contigs vs their source genome: near-identical
MIN_IDENTITY = 0.90          # drop alignments below this (spurious)
COMPLETE = 0.95              # genome recovery category cutoffs (fraction covered)
HIGH = 0.75
PARTIAL = 0.50
CHIMERA_MIN_SEGMENT = 0.20   # each genome must cover >=20% of the contig
CHIMERA_MAX_SEGMENT_OVERLAP = 0.50  # segments must be largely non-overlapping

_SPADES_COV = re.compile(r"_cov_([0-9]+\.?[0-9]*)")
_MEGAHIT_COV = re.compile(r"multi=([0-9]+\.?[0-9]*)")


def parse_coverage(header: str) -> float | None:
    """Per-contig coverage from a SPAdes or MEGAHIT contig header, else None."""
    m = _SPADES_COV.search(header) or _MEGAHIT_COV.search(header)
    return float(m.group(1)) if m else None


def _read_fasta(path):
    name = None
    header = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(buf)
                header = line[1:].rstrip("\n")
                name = header.split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, header, "".join(buf)


def _merge(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by a set of [start,end) intervals after merging."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    total += cur_e - cur_s
    return total


def _overlap_frac(a: tuple[int, int], b: tuple[int, int]) -> float:
    """Fraction of the shorter interval that overlaps the other."""
    ov = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    shorter = min(a[1] - a[0], b[1] - b[0]) or 1
    return ov / shorter


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
        for name, header, seq in _read_fasta(genome_fasta):
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

    aligner = mappy.Aligner(str(genome_fasta), preset=ASM_PRESET)
    if not aligner:
        raise RuntimeError(f"failed to index genome FASTA: {genome_fasta}")

    genome_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    genome_match_bp: dict[str, int] = defaultdict(int)
    genome_block_bp: dict[str, int] = defaultdict(int)
    contig_lengths: list[int] = []
    contig_primary_genome: dict[str, str] = {}
    contig_coverage: dict[str, float | None] = {}
    chimeras: list[dict] = []
    n_unmatched = 0

    for name, header, seq in _read_fasta(contigs_fasta):
        contig_lengths.append(len(seq))
        contig_coverage[name] = parse_coverage(header)
        # gather per-genome contig-coordinate spans for this contig
        per_genome_qspans: dict[str, list[tuple[int, int]]] = defaultdict(list)
        best_blen = 0
        for hit in aligner.map(seq):
            if hit.blen == 0 or hit.mlen / hit.blen < MIN_IDENTITY:
                continue
            genome_intervals[hit.ctg].append((hit.r_st, hit.r_en))
            genome_match_bp[hit.ctg] += hit.mlen
            genome_block_bp[hit.ctg] += hit.blen
            per_genome_qspans[hit.ctg].append((hit.q_st, hit.q_en))
            if hit.blen > best_blen:
                best_blen = hit.blen
                contig_primary_genome[name] = hit.ctg
        if name not in contig_primary_genome:
            n_unmatched += 1
            continue
        # chimera: >=2 genomes each covering a distinct large span of the contig
        clen = len(seq)
        spans = {}
        for g, qs in per_genome_qspans.items():
            covered = _merge(qs)
            if covered / clen >= CHIMERA_MIN_SEGMENT:
                spans[g] = (min(s for s, _ in qs), max(e for _, e in qs), covered)
        if len(spans) >= 2:
            items = sorted(spans.items(), key=lambda kv: -kv[1][2])[:2]
            (g1, s1), (g2, s2) = items
            if _overlap_frac(s1[:2], s2[:2]) < CHIMERA_MAX_SEGMENT_OVERLAP:
                chimeras.append({
                    "contig": name, "length": clen,
                    "segments": [
                        {"genome": g1, "contig_span": [s1[0], s1[1]]},
                        {"genome": g2, "contig_span": [s2[0], s2[1]]},
                    ],
                })

    # per-genome completeness (viral only for recovery)
    per_genome = {}
    cats = {"complete": 0, "high_quality": 0, "partial": 0, "fragmented": 0, "missing": 0}
    completeness_vals = []
    for gid in viral_ids:
        length = gt[gid]["length"]
        covered = _merge(genome_intervals.get(gid, []))
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
    abundance = _abundance_accuracy(
        contig_primary_genome, contig_coverage, contig_lengths,
        list(_read_fasta(contigs_fasta)), gt, viral_ids)

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


def _abundance_accuracy(primary, coverage, _lengths, contig_records, gt, viral_ids):
    length_by_contig = {name: len(seq) for name, _h, seq in contig_records}
    have_cov = [c for c in primary if coverage.get(c) is not None]
    if len(have_cov) < 3:
        return {"available": False, "reason": "per-contig coverage not found in headers"}
    est: dict[str, float] = defaultdict(float)
    for contig, genome in primary.items():
        cov = coverage.get(contig)
        if cov is None or genome not in viral_ids:
            continue
        est[genome] += cov * length_by_contig.get(contig, 0)
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
