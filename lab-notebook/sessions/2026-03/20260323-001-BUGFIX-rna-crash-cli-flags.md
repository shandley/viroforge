# Bug Fixes from virome-qc Reference Dataset Generation

**Date**: 2026-03-23
**Session Type**: BUGFIX
**Status**: Complete

---

## Context

Generated 12 reference datasets across sample types and contamination levels to calibrate virome-qc expected ranges. 10/12 DNA datasets succeeded; both RNA datasets failed.

## Bugs Fixed

### Bug 1 (blocking): RNA virome template switching crash

**File**: `viroforge/workflows/rna_virome.py:210`
**Error**: `ValueError: setting an array element with a sequence`
**Cause**: `self.rng.choice(sequences)` where sequences is a list of BioPython SeqRecords with different lengths. NumPy tries to create an array from them and fails on inhomogeneous shapes.
**Fix**: `sequences[self.rng.integers(len(sequences))]`

### Bug 2 (cosmetic): Viral fraction 100% with --no-vlp

**Cause**: Fraction computed from `viral_abundances.sum()` before contamination was added. With heavy contamination (27.1%), the log showed "Viral fraction: 100.00%" and "Contamination: 27.10%".
**Fix**: Compute fraction from combined normalized abundances after contamination is mixed in.

### Bug 3 (missing CLI): viroforge generate missing key flags

**Cause**: The `viroforge generate` CLI only exposed --preset, --collection-id, --platform, --coverage, --seed. The underlying script supported many more flags.
**Fix**: Added to both arg parser and generate_with_params: --contamination-level, --vlp-protocol, --no-vlp, --molecule-type, --rna-depletion, --adapter-rate, --low-complexity-rate, --duplicate-rate.

## Also

- Bundled herv_consensus.fasta (55 HERV families, 319 KB) from Dfam
- Updated ERV test to use bundled reference instead of expecting it to be absent
