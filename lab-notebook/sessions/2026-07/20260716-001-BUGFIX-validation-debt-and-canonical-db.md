# Validation Debt Paydown, setup-db Fixes, and Canonical DB Adoption

**Date**: 2026-07-16
**Session Type**: BUGFIX
**Status**: Complete

---

## Objective

Pay down the v0.13.0 validation debt: build a real runtime environment, exercise
the realistic-default generation path and the MDA amplification path end to end,
run the test suite, and confirm setup-db reproduces the 20 collections. Fix
whatever the exercise surfaces.

## Environment

Built a fresh venv (CPython 3.12) with `viroforge[web]` plus InSilicoSeq 2.0.1.
Clean install, no compatibility issues. PBSIM3 not built; pbccs is Linux x86-64
only, so PacBio HiFi end to end cannot run on this macOS arm64 machine. RefSeq
and the ICTV VMR are cached locally, so setup-db skips its large download.

## What was done

- Generated a real dataset with the v0.13.0 defaults (dark matter 0.30 plus
  adapter, low-complexity, and duplicate artifacts). Per-read source labels,
  artifact manifests, and metadata all present and correct.
- Ran the test suite: 292 passed, 9 failed. The 9 failures were stale tests, not
  a product bug. The low-complexity injector gained rare-genome protection (skip
  any genome with under 50 reads) after those tests were written, and the
  fixtures used `read_{i}/1` headers that the injector treats as unique
  single-read genomes, so every read was protected and nothing was injected.
  Fixed the two fixtures to use ISS-style headers. Suite now 301 passing.
- Ran a real MDA run (`--amplification mda --mda-chimera-rate 0.02`). phi29
  duplicates and MDA chimeras behave correctly. Found and fixed a metadata bug:
  `duplicate_stats` was built from the raw CLI args, so MDA runs misreported the
  PCR defaults (max_copies 5, error 0.001) instead of the values the injector
  actually upgraded to (20, 0.0001, power_law).
- Ran setup-db end to end. It builds all 20 collections but surfaced two defects:
  - `curate_vaginal_virome_collection` inserted collection 16 without deleting
    existing rows first, so it hit a UNIQUE constraint and exited non-zero.
    Added the DELETE guard the other curation scripts already use.
  - Stored `n_genomes` overcounted the actual associations on collections 5 and 6
    because a few duplicate genome_ids collapse on insert. Added a reconcile step
    to setup-db that sets `n_genomes` to the real association count for all
    collections.

## Canonical DB decision

The live DB predates the seeding fix (a854bae), so it reflects one unseeded
`ORDER BY RANDOM()` draw and setup-db could not reproduce it: about 10
collections differed by 1 to 3 genomes. The ICTV VMR was identical between the
two (the rebuild reused the cached VMR), ruling out taxonomy drift, so the
difference is the seeding change alone. Adopted a fresh seeded rebuild as
canonical, verified reproducible across two rebuilds. The GCF_000837465.1
taxonomy-update error is benign (Mycobacterium phage TM4, ends family Unknown,
in no collection; a doubled RefSeq entry the parser de-duplicated).

## Provenance

The 500 MB database is not tracked in git, so added
`scripts/export_collection_manifest.py` plus `data/collection_membership.tsv` and
a sha256 summary as the diffable record of exact collection membership.
`--check` verifies the committed manifest against the database. Canonical total
is 2025 associations across the 20 collections.

## Validation

- Full suite green: 301 passed.
- MDA metadata now records amplification_method=mda, max_copies=20,
  error_rate=0.0001, copy_distribution=power_law.
- Vaginal curator runs twice against a DB copy with no error; collection 16 = 26.
- setup-db rebuild builds all 20 collections with zero curation failures after
  the fixes; n_genomes consistent for all 20.
- `.claude/CLAUDE.md` per-collection counts synced to the canonical DB.

## Files

- `tests/test_realistic_contamination.py` - fixture headers (stale test fix)
- `scripts/generate_fastq_dataset.py` - duplicate_stats from injector return
- `scripts/curate_vaginal_virome_collection.py` - DELETE guard
- `viroforge/cli/setup_db.py` - `_reconcile_collection_counts()` after cleanups
- `scripts/export_collection_manifest.py` (new) - manifest generator and checker
- `data/collection_membership.tsv`, `data/collection_membership.sha256.txt` (new)
- `data/body_site_collections/*` - refreshed composition exports
- `.claude/CLAUDE.md` - collection counts
