# Collection Composition Corrections (from the accuracy review)

**Date**: 2026-07-16
**Session Type**: CURATION
**Status**: Complete

---

## Objective

Act on the recommendations from the biological-accuracy review (entry 008): fix the
phage-detection blind spot, then curate the collections that were genuinely off.

## What was done

- **Metric fix** (`scripts/evaluate_composition.py`): added a name-based phage
  heuristic so bacteriophages carried as `family=Unknown` count as phage instead of
  host-unknown. Split the site verdict into a BIOLOGY verdict and a separate
  DATA-QUALITY (taxonomy-coverage) verdict. This alone corrected several false
  "deviations": skin went from 0% to 51% phage (unclassified Propionibacterium
  phages), and CF/vaginal/HIV+/mouse-gut phage fractions became honest.
- **Curation** (`scripts/fix_collection_composition.py`, idempotent): corrected the
  four genuinely-off collections by reweighting existing genomes and adding only
  real DB genomes (no invented accessions):
  - Blood (17): Anelloviridae -> dominant (0.60).
  - Lung (19): Anelloviridae -> dominant (0.45); unclassified-phage share reduced.
  - Nasopharynx (4): added 6 Anelloviridae (none were present); eukaryote-dominant.
  - Wastewater (9): rebalanced to phage-majority (0.85 phage / 0.15 eukaryotic).
- **Folded into setup-db** (`viroforge/cli/setup_db.py`): the composition correction
  runs post-curation, before renormalization; all three data corrections
  (genome_type, composition, renormalize) are idempotent so rebuilds stay correct.
- Updated `docs/BIOLOGICAL_ACCURACY_REVIEW.md` and regenerated
  `docs/figures/composition_review.{png,pdf}`.

## Result

Biology: 19/20 collections OK; only vaginal (16) minor (documented dairy-proxy
Lactobacillus-phage limitation). The four curated collections (4, 9, 17, 19) all
moved to biology OK. Remaining flags are data-quality (unclassified phages =
taxonomy coverage), reported separately and left visible rather than hidden.

## Validation

- Evaluator re-run after each step; all four targets -> biology OK.
- `fix_collection_composition.py` idempotent (re-run skips re-adding, same targets).
- `renormalize_abundances.py` = no-op after correction (sums already 1.0).
- setup_db.py compiles; all corrections accept the args passed.
- DB backed up to `viral_genomes.db.preComposition.bak`.

## Notes

Data-quality (unknown-family) flags reflect the DB-wide 42.9% `family=Unknown`
taxonomy gap, addressable by classification not curation. See [[composition_review]]
memory and entry 008.
