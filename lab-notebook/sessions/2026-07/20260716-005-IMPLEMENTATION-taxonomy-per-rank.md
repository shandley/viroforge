# Higher-rank Taxonomy Metrics (Module 4 extension)

**Date**: 2026-07-16
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Extend Module 4 taxonomy benchmarking beyond taxid-exact to genus/family
precision/recall/F1, crediting correct-but-higher-rank classifier calls.

## Approach

Higher-rank scoring must resolve the classifier's arbitrary assigned taxid up the
NCBI tree to learn its family/genus. Rather than add a dependency (taxopy), added
a minimal self-contained NCBI taxdump parser so the logic is fully testable with a
synthetic tree and the user just supplies the taxdump they already have for
Kraken2.

## What was done

- `viroforge/benchmarking/ncbi_tree.py`: parse nodes.dmp (and optional names.dmp);
  `rank_taxid(taxid, rank)` resolves a taxid to its ancestor at a rank.
- `taxonomy.py`: `benchmark_taxonomy(..., ncbi_tree=...)` adds a `per_rank` block.
  For each rank, over the known-virus stratum, resolve true and assigned taxids to
  that rank and classify correct / misclassified / unclassified (a too-shallow
  call counts against recall). Dark matter excluded, so correctly-unclassified
  novel content is still not penalized.
- CLI `--taxdump-dir`; report per-rank table.
- Tests: synthetic mini-taxdump (species->genus->family) with a fixed assignment
  mix; asserts exact per-rank correct/mis/unclassified and precision/recall.
  8 taxonomy tests, 20 benchmark tests total.

## Validation

- Oracle: species precision 0.25/recall 1/6, genus 0.6/0.5, family 0.8/(4/6),
  matching the hand-constructed tree and assignments.
- End-to-end (tax_test, synthetic taxdump + Kraken2): recall rises species 61%
  -> genus 81% -> family 92%, the expected shape (a higher-rank call is correct at
  the higher ranks but not species).

## Note

A taxid absent from the supplied taxdump resolves to None and counts as
unclassified-at-rank (not misclassified). With a real NCBI taxdump every real
classifier taxid resolves, so this only affects synthetic/edge inputs.
