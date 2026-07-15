# MDA (phi29) Duplicate and Chimera Modeling

**Date**: 2026-07-15
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Land the amplification-aware duplicate model as a standalone feature, extracted
from the deferred PR #39 (collection-specific contamination). PR #39 bundled
three unrelated things; this is the one cleanly separable, mergeable piece.

## What was done

- `viroforge/simulators/duplicates.py`: `add_pcr_duplicates()` gains an
  `amplification_method` parameter (`"pcr"` default, backward compatible).
  `"mda"` uses a power-law copy distribution and phi29 GC-bias priming hotspots
  instead of PCR's geometric distribution. New `add_mda_chimeras()` models phi29
  branch-migration chimeras (reads joining sequence from two genomic regions).
- `scripts/generate_fastq_dataset.py`: maps `--amplification mda`/`mda-long` to
  the MDA duplicate profile, and adds `--mda-chimera-rate` (defaults to 0.15 for
  MDA when unset, 0 otherwise). PCR/linker/RdAB behavior is unchanged.

## Why this was safe to extract

- Main's `duplicates.py` had not diverged from PR #39's base, so the model
  applied cleanly.
- The change has zero references to the contamination feature (bacterial/fungal
  background), which stays blocked on unverified reference accessions.
- `amplification_method` defaults to `"pcr"`, so existing calls and default
  generation are unchanged.

## Validation

- Compiles; the default (pcr) generation path and `--help` are intact.
- Functional test on a mock paired FASTQ: `amplification_method="mda"` uses the
  power-law distribution and is reproducible for a fixed seed; `add_mda_chimeras`
  injects chimeras and writes a manifest.

## Context

The other two parts of PR #39 remain open: the collection-specific contamination
+ bacterial/fungal background (blocked; `/verify-references` found ~54% of its
230 reference accessions resolve to a different organism, so the reference set
needs programmatic regeneration from NCBI and the FASTA moved out of git).
