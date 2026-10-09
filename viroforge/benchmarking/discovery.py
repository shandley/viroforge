"""Novel virus discovery benchmarking (Module 8).

Scores a discovery tool's ability to detect viral sequences absent from its
reference database.  ViroForge datasets contain two strata:

- **known** (is_known=True): ICTV family assigned, likely in reference DBs.
- **dark matter** (is_known=False): family='Unknown', representing the large
  fraction of real virome reads with no database match.

A good discovery tool flags dark-matter reads as novel/uncharacterized instead of
silently dropping them.  Module 8 measures:

1. **Dark matter detection rate** — fraction of dark-matter reads the tool
   flagged as present (classified, binned, or called as anything).
2. **Per-genome sensitivity** — detection rate per dark-matter genome, since one
   missed genome with 1000 reads differs from 10 genomes missing 1 read each.
3. **Detection by depth** — how many reads a genome needs before the tool finds it.

Input: a classifier/discovery tool's output listing ALL reads with a
classified/unclassified status (same files Module 4 accepts), plus the ViroForge
ground-truth metadata.  The existing taxonomy parsers are reused.
"""

from __future__ import annotations

import re
from collections import Counter

_DUP = re.compile(r"_dup\d+$")


def _parse_genome_id(read_id: str) -> str:
    return _DUP.sub("", read_id).rsplit("_", 2)[0]


def benchmark_discovery(
    assignments: dict[str, int | None],
    taxonomy_gt: dict,
) -> dict:
    """Score novel virus detection against ViroForge ground truth.

    Args:
        assignments: {read_id: assigned_taxid_or_None} from a classifier.
            Same format as taxonomy.benchmark_taxonomy input — every read the
            classifier processed, with None for unclassified reads.
        taxonomy_gt: {genome_id: {ncbi_taxid, is_known, family, ...}} from
            the dataset metadata benchmarking.taxonomy block.

    Returns:
        dict of discovery metrics.
    """
    known_total = 0
    known_detected = 0
    dark_total = 0
    dark_detected = 0

    per_genome_total: Counter = Counter()
    per_genome_detected: Counter = Counter()
    per_genome_known: dict[str, bool] = {}
    per_genome_family: dict[str, str] = {}

    non_viral = 0

    for rid, assigned in assignments.items():
        gid = _parse_genome_id(rid)
        gt = taxonomy_gt.get(gid)
        if gt is None:
            non_viral += 1
            continue

        is_known = gt.get("is_known", True)
        is_detected = assigned is not None

        per_genome_total[gid] += 1
        per_genome_known[gid] = is_known
        per_genome_family[gid] = gt.get("family", "Unknown")

        if is_known:
            known_total += 1
            if is_detected:
                known_detected += 1
        else:
            dark_total += 1
            if is_detected:
                dark_detected += 1
                per_genome_detected[gid] += 1

    dark_genomes = sorted(gid for gid, known in per_genome_known.items() if not known)
    genome_detection = []
    for gid in dark_genomes:
        total = per_genome_total[gid]
        det = per_genome_detected.get(gid, 0)
        genome_detection.append({
            "genome_id": gid,
            "total_reads": total,
            "detected_reads": det,
            "detection_rate": det / total if total else 0.0,
            "family": per_genome_family.get(gid, "Unknown"),
        })

    detected_genomes = sum(1 for g in genome_detection if g["detection_rate"] > 0)
    n_viral = known_total + dark_total

    result = {
        "reliable": n_viral > 0,
        "n_assignments": len(assignments),
        "n_viral_reads": n_viral,
        "n_non_viral_reads": non_viral,
        "n_known_reads": known_total,
        "n_dark_matter_reads": dark_total,
        "n_known_detected": known_detected,
        "n_dark_detected": dark_detected,
        "known_detection_rate": known_detected / known_total if known_total else None,
        "dark_matter_detection_rate": dark_detected / dark_total if dark_total else None,
        "n_dark_genomes": len(dark_genomes),
        "n_dark_genomes_detected": detected_genomes,
        "genome_detection_rate": (detected_genomes / len(dark_genomes)
                                  if dark_genomes else None),
        "per_genome": genome_detection,
    }

    bins = [(1, 10), (11, 50), (51, 200), (201, 1000), (1001, None)]
    bin_stats = []
    for lo, hi in bins:
        in_bin = [g for g in genome_detection
                  if g["total_reads"] >= lo and (hi is None or g["total_reads"] <= hi)]
        if in_bin:
            detected_in_bin = sum(1 for g in in_bin if g["detection_rate"] > 0)
            label = f"{lo}-{hi}" if hi else f"{lo}+"
            bin_stats.append({
                "bin": label,
                "n_genomes": len(in_bin),
                "n_detected": detected_in_bin,
                "detection_rate": detected_in_bin / len(in_bin),
            })
    result["detection_by_depth"] = bin_stats

    return result
