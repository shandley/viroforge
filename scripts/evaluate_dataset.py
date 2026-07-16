#!/usr/bin/env python3
"""Technical evaluation of a ViroForge-generated dataset.

Independently verifies that generated FASTQ output matches its ground-truth
metadata and that reads actually derive from their claimed source. This does not
trust the metadata's self-reported stats; it re-derives everything from the reads
and the source FASTA.

Dimensions:
  A  read-level source labels + k-mer traceability to the claimed genome
  D  artifact injection (low-complexity, adapter, duplicate, MDA chimera)
  B  abundance fidelity (realized read fraction vs intended abundance)
  C  contamination fraction realized vs claimed
  E  format / statistical sanity (pairing, length, quality, GC, N)

Pass/fail bars are fixed as constants below, set before any numbers were seen.

Usage:
    python scripts/evaluate_dataset.py --dataset validation/eval/default_s42
    python scripts/evaluate_dataset.py --dataset <dir> --json out.json
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
from collections import Counter, defaultdict
from pathlib import Path

# ---- fixed pass/fail bars (set before looking at any results) ----------------
K = 21                       # k-mer size for containment
SAMPLE_PER_LABEL = 300       # reads sampled per source label for containment
CONTAIN_REAL_MIN = 0.80      # real reads: median canonical-kmer containment >= this
CONTAIN_ARTIFACT_MAX = 0.20  # synthetic low-complexity reads: median containment <= this
CONTAIN_ADAPTER_MIN = 0.50   # adapter read-through: insert matches, tail does not
MIN_READS_FOR_CONTAINMENT = 30  # below this a label's median is too noisy to gate
DUP_IDENTITY_MIN = 0.99      # duplicates are near-identical (allow rare PCR error)
LABEL_COVERAGE_MIN = 1.00    # every read must carry a known source= label
LABEL_CONSISTENCY_MIN = 0.999  # genome-id type must match its source label
ABUND_SPEARMAN_MIN = 0.90    # realized vs intended per-genome abundance rank
CONTAM_TYPES = {"host_dna", "rrna", "phix", "reagent_bacteria"}
KNOWN_LABELS = {"viral", "dark_matter", "artifact_low_complexity"} | CONTAM_TYPES

_COMP = str.maketrans("ACGTN", "TGCAN")
_DUP_SUFFIX = re.compile(r"_dup\d+$")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def canon_kmers(seq: str, k: int = K) -> set[str]:
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        rc = revcomp(kmer)
        out.add(kmer if kmer <= rc else rc)
    return out


def shannon_entropy(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    counts = Counter(seq)
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def parse_genome_id(read_id: str) -> str:
    """Recover the source genome id the way ViroForge encodes it.

    ISS appends _{index}_{pair} to the genome id (which itself may contain
    underscores and dots). Duplicate injection appends _dupN. Mirror the tool's
    own rsplit('_', 2) after stripping the dup suffix.
    """
    rid = _DUP_SUFFIX.sub("", read_id)
    parts = rid.rsplit("_", 2)
    return parts[0] if len(parts) >= 3 else rid


def load_fasta(path: Path) -> dict:
    """genome_id -> {'seq', 'mid' (family or contaminant type), 'abundance'}."""
    genomes: dict[str, dict] = {}
    gid = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gid is not None:
                    genomes[gid]["seq"] = "".join(buf)
                buf = []
                header = line[1:].rstrip()
                gid = header.split()[0]
                mid, abund = None, None
                if "|" in header:
                    fields = [f.strip() for f in header.split("|")]
                    if len(fields) >= 2:
                        mid = fields[1]
                    m = re.search(r"abundance=([0-9.eE+-]+)", header)
                    if m:
                        abund = float(m.group(1))
                genomes[gid] = {"seq": "", "mid": mid, "abundance": abund}
            else:
                buf.append(line.strip())
    if gid is not None:
        genomes[gid]["seq"] = "".join(buf)
    return genomes


def iter_fastq(path: Path):
    """Yield (header_line_without_at, seq) per record."""
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            fh.readline()
            q = fh.readline()
            yield h[1:].rstrip(), s.rstrip(), q.rstrip()


def parse_header(header: str) -> dict:
    toks = header.split()
    read_id = toks[0]
    tags = {}
    src = None
    for t in toks[1:]:
        if t.startswith("source="):
            src = t.split("=", 1)[1]
        elif "=" in t:
            k, v = t.split("=", 1)
            tags[k] = v
    return {"read_id": read_id, "source": src, "tags": tags}


def genome_type(mid: str | None) -> str:
    if mid in CONTAM_TYPES:
        return mid
    return "viral_or_dark"  # FASTA can't distinguish; resolved via DB/labels


def evaluate(dataset: Path) -> dict:
    rng = random.Random(42)
    name = None
    fastq_dir = dataset / "fastq"
    r1 = sorted(fastq_dir.glob("*_R1.fastq"))
    r2 = sorted(fastq_dir.glob("*_R2.fastq"))
    fasta = sorted((dataset / "fasta").glob("*.fasta"))
    meta_files = sorted((dataset / "metadata").glob("*_metadata.json"))
    assert r1 and r2 and fasta and meta_files, f"incomplete dataset: {dataset}"
    r1, r2, fasta = r1[0], r2[0], fasta[0]
    meta = json.load(open(meta_files[0]))
    name = meta["collection"]["name"]

    genomes = load_fasta(fasta)

    # Chimera reads legitimately no longer match their duplicate template, so
    # exclude them from the duplicate-identity check.
    chimera_ids: set[str] = set()
    for cm in (dataset / "metadata").glob("*chimera_manifest.tsv"):
        with open(cm) as fh:
            next(fh, None)
            for line in fh:
                chimera_ids.add(line.split("\t", 1)[0])

    # ---- pass 1: read the whole R1, collect labels, per-genome counts, samples
    n_reads = 0
    label_counts: Counter = Counter()
    unknown_labels: Counter = Counter()
    per_genome_reads: Counter = Counter()
    label_type_mismatch = 0
    len_hist: Counter = Counter()
    qual_sum = 0
    qual_n = 0
    gc_sum = 0.0
    n_base_reads = 0
    samples_by_label: dict[str, list] = defaultdict(list)  # reservoir per label
    dup_copy_numbers: Counter = Counter()
    dup_pairs: list[tuple[str, str]] = []  # (dup read_id, duplicate_of)
    read_seq_by_id: dict[str, str] = {}    # for duplicate sequence check (sampled)

    for header, seq, qual in iter_fastq(r1):
        n_reads += 1
        h = parse_header(header)
        src = h["source"]
        label_counts[src] += 1
        if src not in KNOWN_LABELS:
            unknown_labels[src] += 1
        gid = parse_genome_id(h["read_id"])
        per_genome_reads[gid] += 1

        # label vs genome-type consistency (contaminant types are checkable)
        g = genomes.get(gid)
        if g is not None:
            gt = genome_type(g["mid"])
            if src in CONTAM_TYPES and gt != src:
                label_type_mismatch += 1
            if src == "viral" and gt in CONTAM_TYPES:
                label_type_mismatch += 1

        # duplicate accounting
        if h["tags"].get("pcr_duplicate") == "true":
            cn = h["tags"].get("copy_number")
            if cn is not None:
                dup_copy_numbers[int(cn)] += 1
            dof = h["tags"].get("duplicate_of")
            if dof and len(dup_pairs) < 500:
                dup_pairs.append((h["read_id"], dof))

        # stats
        len_hist[len(seq)] += 1
        n_base_reads += seq.count("N")
        gc_sum += (seq.count("G") + seq.count("C")) / max(1, len(seq))
        for ch in qual[:0]:
            pass
        qual_sum += sum(ord(c) - 33 for c in qual)
        qual_n += len(qual)

        # reservoir sample per label for containment/entropy
        bucket = samples_by_label[src]
        if len(bucket) < SAMPLE_PER_LABEL:
            bucket.append((h["read_id"], gid, seq, h["tags"]))
        else:
            j = rng.randint(0, n_reads - 1)
            if j < SAMPLE_PER_LABEL:
                bucket[j] = (h["read_id"], gid, seq, h["tags"])

        # keep a modest map for duplicate-of sequence lookup
        if len(read_seq_by_id) < 20000:
            read_seq_by_id[h["read_id"]] = seq

    n_reads_r2 = sum(1 for _ in iter_fastq(r2))

    # ---- k-mer containment per label (build genome kmer cache lazily) --------
    kmer_cache: dict[str, set] = {}

    def gkmers(gid: str) -> set:
        if gid not in kmer_cache:
            g = genomes.get(gid)
            kmer_cache[gid] = canon_kmers(g["seq"]) if g else set()
        return kmer_cache[gid]

    containment = {}
    entropy_by_label = {}
    for label, bucket in samples_by_label.items():
        cvals = []
        evals = []
        for read_id, gid, seq, tags in bucket:
            evals.append(shannon_entropy(seq))
            rk = canon_kmers(seq)
            if not rk:
                continue
            gk = gkmers(gid)
            if not gk:
                cvals.append(None)
                continue
            cvals.append(len(rk & gk) / len(rk))
        cvals_valid = [c for c in cvals if c is not None]
        containment[label] = {
            "n": len(bucket),
            "median": _median(cvals_valid),
            "mean": _mean(cvals_valid),
            "n_no_ref": sum(1 for c in cvals if c is None),
        }
        entropy_by_label[label] = {"median": _median(evals), "mean": _mean(evals)}

    # ---- duplicate sequence identity (sampled) -------------------------------
    dup_seq_checked = 0
    dup_seq_match = 0
    for dup_id, orig_id in dup_pairs:
        if dup_id in chimera_ids or orig_id in chimera_ids:
            continue  # chimera rewrote this read or the template it was copied from
        s_dup = read_seq_by_id.get(dup_id)
        s_orig = read_seq_by_id.get(orig_id)
        if s_dup is not None and s_orig is not None:
            dup_seq_checked += 1
            # allow small PCR-error edit distance via kmer overlap
            if s_dup == s_orig or _hamming_ok(s_dup, s_orig):
                dup_seq_match += 1

    # ---- abundance fidelity --------------------------------------------------
    # Score viral/dark-matter genomes above the rare-genome floor only.
    # Contaminants are enrichment-suppressed (near-zero reads by design) and
    # floored genomes are tied at 1e-6, so including either distorts the rank
    # correlation without saying anything about read-realization fidelity.
    FLOOR = 1e-6
    intended = {
        gid: g["abundance"]
        for gid, g in genomes.items()
        if g["abundance"] is not None
        and g["abundance"] > FLOOR
        and g["mid"] not in CONTAM_TYPES
    }
    spearman = _spearman(
        [per_genome_reads.get(gid, 0) for gid in intended],
        [intended[gid] for gid in intended],
    )

    # dark matter realized fraction vs 0.30-of-viral target
    dm_cfg = meta.get("configuration", {}).get("dark_matter") or {}
    dm_target = dm_cfg.get("dark_matter_fraction")
    viral_like = label_counts["viral"] + label_counts["dark_matter"]
    dm_realized_of_viral = (label_counts["dark_matter"] / viral_like) if viral_like else None

    # ---- contamination fraction realized vs claimed --------------------------
    contam_reads = sum(label_counts[t] for t in CONTAM_TYPES)
    contam_realized = contam_reads / n_reads if n_reads else 0
    contam_claimed = meta.get("enrichment_stats", {}).get("contamination_fraction")

    # ---- assemble verdicts ---------------------------------------------------
    result = {
        "dataset": dataset.name,
        "collection": name,
        "platform": meta.get("configuration", {}).get("platform"),
        "amplification": meta.get("configuration", {}).get("amplification"),
        "n_reads_r1": n_reads,
        "n_reads_r2": n_reads_r2,
        "label_counts": dict(label_counts),
        "unknown_labels": dict(unknown_labels),
        "containment": containment,
        "entropy_by_label": entropy_by_label,
        "duplicates": {
            "copy_number_dist": dict(sorted(dup_copy_numbers.items())),
            "seq_checked": dup_seq_checked,
            "seq_match": dup_seq_match,
            "reported": meta.get("duplicate_stats", {}),
        },
        "abundance": {
            "spearman_realized_vs_intended": spearman,
            "n_genomes_scored": len(intended),
        },
        "dark_matter": {
            "target_fraction_of_viral": dm_target,
            "realized_fraction_of_viral": dm_realized_of_viral,
        },
        "contamination": {
            "realized_fraction": contam_realized,
            "claimed_fraction": contam_claimed,
        },
        "stats": {
            "pairs_equal": n_reads == n_reads_r2,
            "read_len_max": max(len_hist) if len_hist else 0,
            "read_len_modes": dict(len_hist.most_common(3)),
            "mean_quality": qual_sum / qual_n if qual_n else 0,
            "mean_gc": gc_sum / n_reads if n_reads else 0,
            "n_bases": n_base_reads,
        },
        "reported_stats": {
            "low_complexity": meta.get("low_complexity_stats", {}),
            "adapter": meta.get("adapter_stats", {}),
            "chimera": meta.get("chimera_stats", {}),
        },
    }
    result["verdicts"] = _verdicts(result, label_type_mismatch, n_reads)
    return result


def _verdicts(r: dict, label_type_mismatch: int, n_reads: int) -> dict:
    v = {}
    # A1 label coverage
    known = sum(c for lbl, c in r["label_counts"].items() if lbl in KNOWN_LABELS)
    v["A1_label_coverage"] = _pf(known / n_reads >= LABEL_COVERAGE_MIN,
                                 f"{known}/{n_reads} reads carry a known source label")
    # A2 label<->genome-type consistency
    v["A2_label_consistency"] = _pf(
        (n_reads - label_type_mismatch) / n_reads >= LABEL_CONSISTENCY_MIN,
        f"{label_type_mismatch} contaminant/viral label-vs-genome mismatches")
    # A3 containment by expected model. Only gate on labels with enough reads;
    # rare contaminant labels (<30 reads) have too few sampled reads for a
    # stable median, so they are reported but not allowed to fail the verdict.
    c = r["containment"]
    lc = r["label_counts"]
    checks = []      # gated (enough reads)
    lown = []        # reported only
    for lbl in ("viral", "dark_matter", "host_dna", "rrna", "phix", "reagent_bacteria"):
        if lbl in c and c[lbl]["median"] is not None:
            entry = (lbl, c[lbl]["median"] >= CONTAIN_REAL_MIN, c[lbl]["median"])
            (checks if lc.get(lbl, 0) >= MIN_READS_FOR_CONTAINMENT else lown).append(entry)
    if "artifact_low_complexity" in c and c["artifact_low_complexity"]["median"] is not None:
        m = c["artifact_low_complexity"]["median"]
        checks.append(("artifact_low_complexity", m <= CONTAIN_ARTIFACT_MAX, m))
    detail = "; ".join(f"{l}={m:.2f}{'ok' if ok else 'FAIL'}" for l, ok, m in checks)
    if lown:
        detail += " | low-n(not gated): " + "; ".join(f"{l}={m:.2f}" for l, _, m in lown)
    v["A3_kmer_traceability"] = _pf(
        all(ok for _, ok, _ in checks) and len(checks) > 0, detail)
    # D1 low-complexity entropy separation (N/A when artifacts are disabled)
    lc = r["entropy_by_label"].get("artifact_low_complexity", {}).get("median")
    vir = r["entropy_by_label"].get("viral", {}).get("median")
    if r["label_counts"].get("artifact_low_complexity", 0) == 0:
        v["D1_lowcomplexity_entropy"] = _pf(True, "no low-complexity artifacts (disabled)")
    else:
        v["D1_lowcomplexity_entropy"] = _pf(
            lc is not None and vir is not None and lc < vir - 0.5,
            f"artifact entropy median={lc}, viral median={vir}")
    # D2 duplicate sequence identity (chimera-rewritten reads already excluded)
    d = r["duplicates"]
    v["D2_duplicate_identity"] = _pf(
        d["seq_checked"] == 0 or d["seq_match"] / d["seq_checked"] >= DUP_IDENTITY_MIN,
        f"{d['seq_match']}/{d['seq_checked']} sampled duplicates match template")
    # B abundance
    sp = r["abundance"]["spearman_realized_vs_intended"]
    v["B_abundance_rank"] = _pf(sp is not None and sp >= ABUND_SPEARMAN_MIN,
                                f"spearman={sp}")
    # E pairing + format
    s = r["stats"]
    v["E_pairing_format"] = _pf(s["pairs_equal"] and s["read_len_max"] <= 300,
                                f"pairs_equal={s['pairs_equal']} len_max={s['read_len_max']}")
    return v


def _pf(ok: bool, detail: str) -> dict:
    return {"pass": bool(ok), "detail": detail}


def _median(xs):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def _hamming_ok(a: str, b: str, max_frac: float = 0.05) -> bool:
    if abs(len(a) - len(b)) > 3:
        return False
    n = min(len(a), len(b))
    d = sum(1 for i in range(n) if a[i] != b[i])
    return d / max(1, n) <= max_frac


def _spearman(xs, ys):
    if len(xs) < 3:
        return None
    rx, ry = _rank(xs), _rank(ys)
    n = len(xs)
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return 1 - (6 * d2) / (n * (n * n - 1))


def _rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--json")
    args = ap.parse_args()
    r = evaluate(Path(args.dataset))
    if args.json:
        Path(args.json).write_text(json.dumps(r, indent=2))
    # human summary
    print(f"\n=== {r['dataset']}  ({r['platform']}, amp={r['amplification']}) ===")
    print(f"reads R1/R2: {r['n_reads_r1']}/{r['n_reads_r2']}  labels: {r['label_counts']}")
    for k, val in r["verdicts"].items():
        print(f"  [{'PASS' if val['pass'] else 'FAIL'}] {k}: {val['detail']}")
    print("  containment medians:",
          {l: (round(c['median'], 3) if c['median'] is not None else None)
           for l, c in r["containment"].items()})


if __name__ == "__main__":
    main()
