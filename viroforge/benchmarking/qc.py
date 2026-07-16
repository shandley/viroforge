"""QC benchmarking metrics engine (Module 1).

Compares a QC pipeline's cleaned reads against ViroForge's per-read source labels
and reports how well contamination was removed and, critically, how well viral
reads were retained.

Scope (v1): measures read *removal*, R1 only. It does not measure trimming
quality. Adapter read-through reads keep their `source=viral` label and correctly
count as keep-and-trim, so they need no special handling. Paired-end reads are
scored once per fragment via R1; mates share the label.

PCR duplicates inherit their template's source (a duplicated host read is still
host), so removing a viral duplicate is dedup, not over-filtering. Contamination
metrics are therefore computed over non-duplicate reads only, and duplicate
handling is reported separately as a dedup rate.
"""

from __future__ import annotations

from collections import Counter

from .parsers import strip_mate

# Explicit, overridable ground-truth policy: which source labels a good virome
# QC pipeline should keep versus remove.
DEFAULT_KEEP_REMOVE: dict[str, str] = {
    "viral": "keep",
    "dark_matter": "keep",
    "erv_exogenous": "keep",       # exogenous retrovirus: treat as viral signal
    "host_dna": "remove",
    "rrna": "remove",
    "phix": "remove",
    "reagent_bacteria": "remove",  # the spec calls this "bacterial"
    "artifact_low_complexity": "remove",
    "erv_endogenous": "remove",    # endogenous retrovirus: host-derived
}

MATCH_GATE = 0.99  # min fraction of cleaned reads that must map back to raw


def _resolve_kept(raw_names: set[str], kept_names: set[str]) -> tuple[set[str], float, bool]:
    """Map surviving read names to raw read names; return (kept_raw, match_rate, stripped)."""
    if not kept_names:
        return set(), 1.0, False
    exact = kept_names & raw_names
    rate = len(exact) / len(kept_names)
    if rate >= MATCH_GATE:
        return exact, rate, False
    # Retry ignoring a mate suffix on both sides (some tools drop /1,/2).
    raw_by_stripped: dict[str, str] = {}
    for n in raw_names:
        raw_by_stripped.setdefault(strip_mate(n), n)
    stripped_kept = set()
    for n in kept_names:
        rn = raw_by_stripped.get(strip_mate(n))
        if rn is not None:
            stripped_kept.add(rn)
    rate2 = len(stripped_kept) / len(kept_names)
    if rate2 > rate:
        return stripped_kept, rate2, True
    return exact, rate, False


def benchmark_qc(
    raw: dict[str, tuple[str, bool]],
    kept_names: set[str],
    keep_remove: dict[str, str] | None = None,
) -> dict:
    """Compute QC benchmark metrics.

    Args:
        raw: {read_name: (source_label, is_pcr_duplicate)} from the raw reads.
        kept_names: set of read names present in the cleaned reads.
        keep_remove: source -> "keep"|"remove" policy (defaults to DEFAULT_KEEP_REMOVE).

    Returns a dict of metrics; see keys below. Includes `match_rate` and a
    `reliable` flag that is False when the name match rate is below MATCH_GATE.
    """
    policy = dict(DEFAULT_KEEP_REMOVE)
    if keep_remove:
        policy.update(keep_remove)

    raw_names = set(raw)
    kept_raw, match_rate, stripped = _resolve_kept(raw_names, kept_names)

    # Confusion matrix over non-duplicate reads. Positive class = should-remove.
    tp = fp = fn = tn = 0
    per_type_total: Counter = Counter()
    per_type_removed: Counter = Counter()
    unknown_total = unknown_removed = 0
    kept_total = kept_removed = 0  # keep-class reads (viral retention)

    # Dedup metric over duplicate reads (any source).
    dup_total = dup_removed = 0

    for name, (source, is_dup) in raw.items():
        removed = name not in kept_raw
        if is_dup:
            dup_total += 1
            dup_removed += removed
            continue
        decision = policy.get(source)
        if decision is None:
            unknown_total += 1
            unknown_removed += removed
            continue
        if decision == "remove":
            per_type_total[source] += 1
            per_type_removed[source] += removed
            if removed:
                tp += 1
            else:
                fn += 1
        else:  # keep
            kept_total += 1
            kept_removed += removed
            if removed:
                fp += 1
            else:
                tn += 1

    def rate(num, den):
        return num / den if den else None

    per_type = {
        s: {
            "n": per_type_total[s],
            "removal_rate": rate(per_type_removed[s], per_type_total[s]),
        }
        for s in per_type_total
    }

    precision = rate(tp, tp + fp)
    recall = rate(tp, tp + fn)
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision and recall
        else None
    )

    return {
        "match_rate": match_rate,
        "match_rate_after_mate_strip": stripped,
        "reliable": match_rate >= MATCH_GATE,
        "n_raw_reads": len(raw),
        "n_cleaned_reads": len(kept_names),
        # headline virome metrics
        "viral_retention": rate(tn, tn + fp),
        "over_filtering_rate": rate(fp, tn + fp),
        "removal_rate_by_type": per_type,
        # aggregate contamination-removal confusion matrix (non-duplicate reads)
        "contamination": {
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": precision, "recall": recall, "f1": f1,
        },
        # orthogonal QC operations
        "dedup": {
            "n_duplicates": dup_total,
            "dedup_rate": rate(dup_removed, dup_total),
        },
        "unknown_source": {
            "n": unknown_total,
            "removal_rate": rate(unknown_removed, unknown_total),
        },
        "keep_remove_policy": policy,
    }
