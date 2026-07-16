# Taxonomy Benchmarking Module (Module 4, read-based v1)

**Date**: 2026-07-16
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Add Module 4 (taxonomy benchmarking), read-based, chosen scope: export per-genome
taxonomy into the metadata, then ship taxid-exact scoring plus an abundance
profile, with the unclassified/dark-matter handling done right.

## What was done

- Metadata export (prerequisite): `scripts/generate_fastq_dataset.py` now exports
  `benchmarking.taxonomy` = {genome_id: {ncbi_taxid, realm..species, is_known}}
  for every viral genome, via a new `_export_taxonomy()` DB lookup. `is_known`
  is False when the ICTV family is Unknown. Makes the taxonomy benchmark
  self-contained (no 500 MB DB needed). `ncbi_taxid` is populated for all 14,423
  genomes, which is what enables taxid-level comparison.
- `viroforge/benchmarking/taxonomy.py`: parsers for Kraken2 per-read output
  (including `--use-names`) and a generic `read_id\ttaxid` TSV; taxid-exact
  scoring stratified into known viruses vs dark matter; abundance profile
  (Bray-Curtis, Pearson, Spearman, MAE) over NCBI taxids. Read names embed the
  genome accession, so ground truth comes from the read name + the metadata
  export; no raw reads needed.
- `report.py`: `taxonomy_to_markdown`. CLI: `benchmark taxonomy` subcommand.
- `tests/test_benchmark_taxonomy.py`: independent oracle (fixed correct /
  unclassified / misclassified mix per stratum). 5 tests; 17 benchmark tests total.

## The key design decision (unclassified done right)

Reads are split into known viruses (ICTV family assigned) and dark matter
(family Unknown). For dark matter, a high unclassified rate is reported as
correct behavior (classifiers should not force a call on novel content), never as
a penalty. Known-virus sensitivity/precision is the real performance signal. This
avoids the trap where every classifier looks terrible on the dark matter that
viromes are largely made of.

## Validation

- 5 oracle tests: per-stratum correct/unclassified/misclassified and
  sensitivity/precision/unclassified-rate recovered exactly.
- End-to-end on real data (tax_test, synthetic Kraken2 output): known viruses
  91.6% sensitivity / 96.2% precision / 4.8% unclassified; dark matter 83.6%
  correctly unclassified; abundance Bray-Curtis 0.37, Pearson 0.74. Numbers track
  the synthetic construction.

## Scope / deferred

Taxid-exact only (higher-rank genus/family metrics need lineage resolution, and
are the next addition). Read-based only. Kraken2 + generic formats; Centrifuge,
MMseqs2, DIAMOND deferred. Docs: `docs/PHASE13C_TAXONOMY_BENCHMARK.md`.
