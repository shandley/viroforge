# Assembly Benchmarking Module (Phase 13B, Module 2)

**Date**: 2026-07-16
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Add Module 2 (assembly benchmarking) to the framework: align a user's assembled
contigs to the true ViroForge genomes and score genome recovery, chimeras,
contiguity, and abundance-estimation accuracy. Approach chosen with Scott:
mappy (minimap2) plus a core-first scope (JSON/markdown, no HTML/coverage plots).

## What was done

- `viroforge/benchmarking/assembly.py`: aligns contigs with mappy (asm5),
  computes per-genome completeness from merged alignment intervals, recovery
  categories, identity, chimera detection (two genomes each covering a distinct
  non-overlapping segment of a contig), N50/L50, observed-vs-expected
  completeness (using the metadata's Lander-Waterman `expected_completeness`),
  and abundance accuracy.
- Abundance accuracy: per-contig coverage from SPAdes/MEGAHIT headers is depth
  (proportional to abundance/length), so per-genome abundance is estimated as
  depth times length, normalized, and compared to true relative abundance.
- `report.py`: `assembly_to_markdown`; `benchmark assembly` CLI subcommand.
- mappy added as the `benchmark` extra in setup.py.
- `tests/test_benchmark_assembly.py`: mock contigs of known recovery plus a
  deliberate chimera; independent oracle. 7 tests pass (12 with QC).

## Validation

- Oracle tests: completeness 1.0 / 0.6 / 0.3 / 0 recovered correctly; chimera
  flagged with the right two genomes; N50/classification/abundance correct.
- End-to-end on real data (default_s42, 33 viral genomes assembled as a mock):
  28/33 recovered (19 complete, 9 partial, 5 missing), 0 chimeras, N50 71kb.
  Abundance Spearman 0.76 when contig coverage is realistic depth
  (expected_coverage), which surfaced and confirmed the depth-times-length
  conversion; a first mock that set coverage proportional to abundance gave 0.02,
  the expected wrong answer for the wrong input.

## Docs

`docs/PHASE13B_ASSEMBLY_BENCHMARK.md`.
