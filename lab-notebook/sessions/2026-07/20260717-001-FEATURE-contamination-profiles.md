# Collection-Specific Contamination Profiles (PR #39 Part 1, rebuilt)

**Date**: 2026-07-17
**Session Type**: FEATURE
**Status**: Complete (core)

---

## Objective

Rebuild the clean, valuable core of PR #39 on current main: sample-type-specific
contamination baselines per collection, with `--contamination-level` acting as a
multiplier. Deliberately excludes the parts of PR #39 that blocked it (4.4M lines
of FASTA in git, unverified accessions, unrelated duplicates/low_complexity
changes) - those are a separate Part 2 (issue #37).

## What was done

- **Schema** (`viroforge/data/database_schema.py`): added `default_host_pct`,
  `default_rrna_pct`, `default_reagent_pct`, `default_phix_pct`, `host_organism`
  to body_site_collections (NULL = fall back to global preset).
- **Logic** (`viroforge/core/contamination.py`): `create_contamination_profile`
  gained `collection_defaults`. With it, the level is a MULTIPLIER
  (`LEVEL_MULTIPLIERS` = clean 0.25 / realistic 1.0 / heavy 2.5 / failed 4.0) on
  the collection baseline, with per-column NULL fallback and sane clamps (host/
  rRNA <=95%, reagent <=50%; PhiX unscaled). Without it, behaviour is unchanged
  (backwards compatible).
- **Cited values** (`data/reference_profiles/contamination_defaults.tsv`):
  literature-informed per-site baselines (blood 40% host cfDNA, marine 0.05%
  host-free, CF 25%, ocular kitome-heavy reagent, RNA collections higher rRNA).
  Anchored to Salter 2014 (kitome, PMID 25387460, verified) and Zolfo 2019
  (ViromeQC, PMID 31748692, verified). Documented as semi-quantitative modeling
  defaults, not per-site empirical measurements.
- **Populate** (`scripts/populate_contamination_defaults.py`): idempotent; adds
  the columns if an older DB predates them, then writes values by collection_id.
  Folded into `setup-db` post-curation.
- **Wiring** (`scripts/generate_fastq_dataset.py`): `prepare_genomes` passes the
  collection's defaults (gated on a non-NULL host default) to both DNA
  contamination calls. RNA path (`create_rna_contamination_profile`) unchanged for
  now. Closes #38 conceptually.

## Validation

- Multiplier/clamp/fallback verified: blood realistic=40% host, heavy/failed clamp
  to 95%, marine ~host-free, NULL columns -> preset, no-defaults -> original preset.
- End-to-end via the real DB load path: collection 17 -> 40% host, 5 -> 0%, 1 -> 5%.
- Populate idempotent (re-apply is a no-op). setup_db.py compiles.
- Tests: `tests/test_contamination_profiles.py` (8, independent-oracle). Full run
  of realistic+rna contamination + composition + new = 100 passed, no regressions.

## Not done (tracked separately)

- Part 2 (issue #37): bacterial/fungal background fragments - keep FASTA out of git
  (build on demand, resolver-discovered), verify every accession (reclassification
  vs genuine error), own PR.
- RNA-path contamination multiplier (create_rna_contamination_profile).

See [[composition_review]] and the PR #39 review.
