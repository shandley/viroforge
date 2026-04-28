# ERV Integration into Contamination Profile Factory

**Date**: 2026-04-28
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Make ERV injection seamless by integrating it into the standard contamination profile factory instead of having it as a separate bolt-on step.

## Changes

- Added erv_endogenous_pct and erv_exogenous_pct to contamination profile configs
- create_contamination_profile() and create_rna_contamination_profile() now call add_erv_endogenous/exogenous when rates > 0
- Removed ad-hoc ERV injection block from generate_fastq_dataset.py main()
- ERV kwargs flow through prepare_genomes -> create_contamination_profile
- Exposed --erv-endogenous-rate and --erv-exogenous-rate in viroforge generate CLI
- Wired through generate.py build_command
