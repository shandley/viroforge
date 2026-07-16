"""Shared contig-to-genome alignment for assembly and taxonomy benchmarking.

Aligns contigs to the true genomes with minimap2 (mappy) and returns, per contig,
its primary genome, chimera status, and aligned genomes, plus per-genome coverage
intervals for completeness. Both the assembly module (genome recovery) and the
taxonomy module (contig-based ground truth) consume this.
"""

from __future__ import annotations

import re
from collections import defaultdict

import mappy

ASM_PRESET = "asm5"          # contigs vs their source genome: near-identical
MIN_IDENTITY = 0.90          # drop alignments below this (spurious)
CHIMERA_MIN_SEGMENT = 0.20   # each genome must cover >=20% of the contig
CHIMERA_MAX_SEGMENT_OVERLAP = 0.50  # segments must be largely non-overlapping

_SPADES_COV = re.compile(r"_cov_([0-9]+\.?[0-9]*)")
_MEGAHIT_COV = re.compile(r"multi=([0-9]+\.?[0-9]*)")


def parse_coverage(header: str) -> float | None:
    """Per-contig coverage from a SPAdes or MEGAHIT contig header, else None."""
    m = _SPADES_COV.search(header) or _MEGAHIT_COV.search(header)
    return float(m.group(1)) if m else None


def read_fasta(path):
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


def merge_len(intervals: list[tuple[int, int]]) -> int:
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


def overlap_frac(a: tuple[int, int], b: tuple[int, int]) -> float:
    """Fraction of the shorter interval that overlaps the other."""
    ov = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    shorter = min(a[1] - a[0], b[1] - b[0]) or 1
    return ov / shorter


def align_contigs(contigs_fasta, genome_fasta) -> dict:
    """Align contigs to genomes; return per-contig assignments and genome coverage.

    Returns {
      "contigs": [{name, length, coverage, primary_genome, is_chimera,
                   chimera_segments, aligned_genomes}],
      "genome_intervals": {genome: [(r_st, r_en)]},
      "genome_match_bp": {genome: int}, "genome_block_bp": {genome: int},
    }
    """
    aligner = mappy.Aligner(str(genome_fasta), preset=ASM_PRESET)
    if not aligner:
        raise RuntimeError(f"failed to index genome FASTA: {genome_fasta}")

    genome_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    genome_match_bp: dict[str, int] = defaultdict(int)
    genome_block_bp: dict[str, int] = defaultdict(int)
    contigs: list[dict] = []

    for name, header, seq in read_fasta(contigs_fasta):
        clen = len(seq)
        per_genome_qspans: dict[str, list[tuple[int, int]]] = defaultdict(list)
        best_blen = 0
        primary = None
        for hit in aligner.map(seq):
            if hit.blen == 0 or hit.mlen / hit.blen < MIN_IDENTITY:
                continue
            genome_intervals[hit.ctg].append((hit.r_st, hit.r_en))
            genome_match_bp[hit.ctg] += hit.mlen
            genome_block_bp[hit.ctg] += hit.blen
            per_genome_qspans[hit.ctg].append((hit.q_st, hit.q_en))
            if hit.blen > best_blen:
                best_blen = hit.blen
                primary = hit.ctg

        is_chimera = False
        chimera_segments = None
        if primary is not None:
            spans = {}
            for g, qs in per_genome_qspans.items():
                covered = merge_len(qs)
                if covered / clen >= CHIMERA_MIN_SEGMENT:
                    spans[g] = (min(s for s, _ in qs), max(e for _, e in qs), covered)
            if len(spans) >= 2:
                items = sorted(spans.items(), key=lambda kv: -kv[1][2])[:2]
                (g1, s1), (g2, s2) = items
                if overlap_frac(s1[:2], s2[:2]) < CHIMERA_MAX_SEGMENT_OVERLAP:
                    is_chimera = True
                    chimera_segments = [
                        {"genome": g1, "contig_span": [s1[0], s1[1]]},
                        {"genome": g2, "contig_span": [s2[0], s2[1]]},
                    ]

        contigs.append({
            "name": name, "length": clen, "coverage": parse_coverage(header),
            "primary_genome": primary, "is_chimera": is_chimera,
            "chimera_segments": chimera_segments,
            "aligned_genomes": sorted(per_genome_qspans),
        })

    return {
        "contigs": contigs,
        "genome_intervals": genome_intervals,
        "genome_match_bp": genome_match_bp,
        "genome_block_bp": genome_block_bp,
    }
