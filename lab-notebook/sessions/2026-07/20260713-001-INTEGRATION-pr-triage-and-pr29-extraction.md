# PR Triage and Extraction of Code Fixes from PR #29

**Date**: 2026-07-13
**Session Type**: INTEGRATION
**Status**: In progress

---

## Objective

Review Leran10's 30 open issues and 27 open PRs, then begin landing the
accepted fixes without pulling in the data and schema changes that need
separate decisions.

## Approach

Three parallel reviewers read the actual diffs: one for the cumulative
bug-fix chain (PRs 4-32), two for the independent scientific PRs (39-58).
Findings were cross-checked against the live database and the PR branches
before any code was touched.

## Key findings from the review

- The cumulative chain (PRs 4-32) targets real, reproducible bugs. PR #29
  is the correct code tip for the long-read and web fixes, not PR #32:
  PR #32's PacBio path reads a single `.sam` file, but PBSIM3 emits
  per-genome `{prefix}_0001.bam`, so PR #32 silently reintroduces issue
  #33. PR #29 globs the BAM files and merges them with `samtools merge`.
- The database has 20 collections at IDs 9-28 (IDs 1-8 empty). The
  CLAUDE.md claim of 28 collections with VLP-comparison collections at
  9-16 is not backed by the database. Issue #5's gap complaint is real.
- PR #39 hardcodes about 107 bacterial GCF accessions that were never
  verified. At least four are wrong (for example, GCF_000284795.1 labeled
  Megasphaera elsdenii is actually Salmonella enterica). Blocks PR #39
  until the accessions are verified.

## What was done this session

1. Committed three local fixes already on main: #34 (RNA metadata
   rebuild), #35 (hybrid version string and SPAdes flag), #36 (validator
   sequences alias). Commits db189c0, a91888f, d5867a4.

2. Audited PR #29's full diff (40 files) against current main and found
   it is entangled: it bundles the deferred TSV accession swaps, the
   collection-ID renumbering baked into 18 curate scripts, and it would
   revert the #35 fix. So it cannot be merged wholesale.

3. Extracted only the code fixes from PR #29 by applying its authored
   changes (merge-base to branch tip) for the code, CLI, web, and doc
   files on top of current main. The RNA-workflow section conflicted
   because PR #29 fixed #34 independently; resolved in favor of the
   committed #34 version.

## Files changed (this commit)

- scripts/generate_fastq_dataset.py (PacBio BAM merge, ERRHMM names,
  per-genome depth, graceful tool handling)
- viroforge/simulators/longread.py (per-genome depth)
- viroforge/cli/{setup_db,db_utils,__init__,browse,batch,report,compare,generate}.py
- viroforge/web/app.py and two templates
- README.md, lab-notebook/INDEX.md

## Excluded from PR #29 (handled separately)

- data/body_site_collections/*.tsv accession swaps (pending review)
- scripts/curate_*.py collection-ID renumbering (breaking migration)
- generate_hybrid_dataset.py and validate_fastq_dataset.py (would revert
  #35 and re-do #36)
- docs/ANALYSIS_VALIDATION_TEST.md (uses the renumbered IDs)

## Testing

- py_compile passes on all extracted Python files.
- generate_fastq_dataset.py --help parses.
- setup-db subcommand is registered in the CLI parser.
- CLI/web imports depend on rich/flask, which are absent from
  .venv_test; not run there.

## Next Steps

- [ ] Cherry-pick PR #23 (pbccs Linux x86-64 note) on top.
- [ ] Smoke-test the PacBio HiFi path on a Linux x86-64 host.
- [ ] Clean up the orphaned SAM-reading `_run_pbsim3_clr` in longread.py.
- [ ] Decide on the TSV accession swaps and the collection-ID renumbering.
- [ ] Verify PR #39's bacterial accessions before it can merge.

---

**Status**: In progress
**Impact**: Lands the accepted long-read, CLI, and web fixes from PR #29
while keeping the contested data and schema changes out.
