# Contig-based Taxonomy + Shared Aligner (Module 4 extension)

**Date**: 2026-07-16
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Add contig-based taxonomy benchmarking with a chimera-handling option (report
separately or LCA), and refactor the assembly aligner into a shared helper both
the assembly and taxonomy modules call.

## What was done

- `viroforge/benchmarking/align.py` (new): shared contig-to-genome aligner
  (mappy). Returns per contig its primary genome, chimera status/segments, and
  aligned genomes, plus per-genome coverage intervals. Moved the FASTA reader,
  interval merge, overlap, coverage parser, and thresholds here.
- `assembly.py`: refactored to consume `align_contigs()`; behavior unchanged (all
  7 assembly tests still pass).
- `taxonomy.py`: extracted a unit-agnostic `_score(items, ...)` shared by
  read-based and contig-based scoring. Added `benchmark_taxonomy_contigs()`:
  derives each contig's true taxid from its primary genome (via the shared
  aligner), counts novel (unaligned) and contaminant contigs separately, and
  handles chimeras by `chimera_handling`: "exclude" (reported, not scored) or
  "lca" (true taxid = LCA of the two segment genomes, so a call is credited only
  to the shared rank). Added `_lca` and `NcbiTree.path_to_root`.
- `report.py`: `taxonomy_contig_to_markdown` (shares strata/per-rank rendering).
- CLI: `benchmark taxonomy --mode contig-based` with `--contigs`, `--genomes`,
  `--contig-taxonomy`, `--chimera-handling`.
- tests: contig oracle (mock contigs of known origin + one chimera, one novel,
  one contaminant); exclude vs lca modes; LCA credits the chimera at family.
  3 new tests; 23 benchmark tests total.

## Validation

- Refactor: assembly + read-based taxonomy tests unchanged (15 pass).
- Contig oracle: classification counts (viral/novel/contaminant/chimeric),
  exclude-mode strata, and lca-mode chimera scoring (known misclassified 1->2,
  family per-rank correct=2) all exact.
- End-to-end (tax_test, mock contigs + synthetic taxdump): 25 viral contigs,
  known sensitivity 57.9% / precision 73.3%, dark stratum reported.

## Notes / deferred

Chimera "lca" requires a taxdump. Contig abundance profile is contig-count based
(not length-weighted) in v1. Docs: `docs/PHASE13C_TAXONOMY_BENCHMARK.md`.
