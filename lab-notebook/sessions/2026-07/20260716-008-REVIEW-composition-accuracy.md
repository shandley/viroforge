# Biological Accuracy Review of Per-Site Viral Composition

**Date**: 2026-07-16
**Session Type**: REVIEW
**Status**: Complete (Layers 1-2)

---

## Objective

Evaluate whether each of the 20 collections' viral composition (phage vs
eukaryotic-virus balance, nucleic-acid class, dominant/signature families) matches
current literature expectations, using a verified, reproducible, tunable framework.
Apply the clear-cut data fix and fold it into `setup-db`.

## What was done

- **Verified ICTV property map** (`scripts/build_composition_reference.py` ->
  `data/reference_profiles/family_properties.tsv`): per-family Baltimore/NA type
  (incl. a coarse ds/ss + DNA/RNA type with its own purity) and host type (phage
  if bacteria/archaea, else eukaryotic) from ICTV VMR MSL40. Covers 108/108 named
  collection families. Used INSTEAD of the unreliable `genomes.genome_type`.
- **Evaluator** (`scripts/evaluate_composition.py` observe|evaluate): observed
  composition in 3 currencies (abundance/count/prevalence), scored against a
  tunable cited reference profile (`data/reference_profiles/virome_composition.yaml`)
  of per-site expected BANDS + a global `strictness` dial. Emits a severity-ranked
  scorecard (`validation/composition_scorecard.json`).
- **Verified literature** (`validation/dossiers/*.json` -> `scripts/build_bibliography.py`
  -> `docs/composition_references.bib`): 50 references, every DOI+PMID confirmed via
  `/verify-references` (3 off-by-one years corrected). See `validation/VERIFICATION_SUMMARY.md`.
- **Report + figure**: `docs/BIOLOGICAL_ACCURACY_REVIEW.md`,
  `scripts/plot_composition_review.py` -> `docs/figures/composition_review.{png,pdf}`.
- **Clear-cut fix applied**: `scripts/fix_genome_type.py --apply` relabeled 3516
  genomes whose `genome_type` mislabeled RNA viruses as dsDNA (DB backed up).
- **Folded into setup-db** (`viroforge/cli/setup_db.py`): post-curation now runs the
  genome_type correction and `scripts/renormalize_abundances.py --apply` (the
  abundance-normalization fix from the web-export work), so rebuilds stay correct.

## Findings

11/20 collections biologically sound. Deviations: nasopharynx (4, major - balance
inverted); skin (3), wastewater (9), vaginal (16), blood (17), lung (19) (moderate);
mouse gut (8), HIV+ gut (11), urinary (20) (minor). Confirmed data bug: `genome_type`
silently defaulted RNA viruses to dsDNA (arbovirus collection was 100% dsDNA; now
100% ssRNA). Composition deviations left as ranked curation recommendations (they
need new genome selection + abundance design, not a mechanical fix).

## Validation

- `/verify-references`: all 50 identifiers resolve to the correct papers.
- Evaluator verdicts match the independent per-site dossier verdicts.
- Post genome_type fix: `PRAGMA integrity_check` = ok; collection 14 = 100% ssRNA;
  collection 15 shows rotavirus dsRNA + norovirus ssRNA; DB-wide ssDNA 6 -> 1520.
- `setup_db.py` compiles; both folded scripts accept the args passed.

## Notes

Scorecard derives NA/host type from ICTV, so it is independent of the genome_type
fix. Tune via the `strictness` dial or any band in `virome_composition.yaml`, then
re-run `evaluate_composition.py evaluate`.
