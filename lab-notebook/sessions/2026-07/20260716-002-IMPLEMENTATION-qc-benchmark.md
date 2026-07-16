# QC Benchmarking Module (Phase 13B, Module 1)

**Date**: 2026-07-16
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Start the benchmarking framework (Phase 13B) with Module 1: QC benchmarking.
Validate a QC pipeline's contamination removal and viral retention against the
per-read `source=` ground truth, which the data-quality evaluation confirmed is
accurate.

## What was done

- New package `viroforge/benchmarking/`:
  - `parsers.py`: read raw reads to `{read_name: (source, is_duplicate)}` and
    cleaned reads to a set of surviving names; plain and gzipped FASTQ.
  - `qc.py`: the metrics engine. Read-name match-rate gate (flags unreliable
    below 99 percent, tolerates a trailing mate suffix, never strips `.version`
    or `_dupN`); confusion matrix over non-duplicate reads (positive class =
    should-be-removed); per-type removal rate; viral retention and over-filtering;
    separate dedup rate; unknown-source bucket; explicit overridable keep/remove
    policy.
  - `report.py`: JSON plus a markdown summary that leads with viral retention.
- CLI `viroforge benchmark qc` (`viroforge/cli/benchmark.py`, wired into
  `cli/__init__.py`).
- `tests/test_benchmark_qc.py`: an independently-counted oracle (fixed drop
  counts per class, assert the engine recovers the hand-computed rates) plus
  mate-strip, low-match-rate, gzip, and unknown-source cases. 5 tests pass.

## Design decisions (advisor-reviewed)

- Match rate is a hard gate: refuse to present confident numbers when cleaned
  reads do not map back to raw.
- Duplicates inherit their template's source (host and rrna duplicates exist in
  the output), so contamination is scored on non-duplicate reads and dedup is a
  separate metric.
- Lead the report with viral retention and over-filtering, the failure mode for
  low-biomass viromes, not aggregate F1.
- R1 only for v1; gzip accepted; keep/remove policy explicit and overridable.

## Finding

Low-complexity artifact reads carry two `source=` tags
(`source=viral source=artifact_low_complexity`) because the injector appends
rather than replaces. The benchmark uses the last tag (the read's current
nature). A naive consumer that reads the first tag would misclassify these as
viral. Worth cleaning up in the generator later; not blocking.

## Validation

- 5 oracle tests pass; benchmarking-metadata and phase13a tests still pass (16
  total).
- End-to-end on a real dataset (default_s42) with a mock good-QC cleaned file:
  viral retention 99.0 percent, over-filtering 1.0 percent, all contaminant types
  removed at 100 percent, dedup reported separately.
