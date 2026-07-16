# Harden the Web-Data Export Contract

**Date**: 2026-07-16
**Session Type**: INTEGRATION
**Status**: Complete

---

## Objective

Turn `scripts/export_web_data.py` (the boundary between ViroForge and the HVP web
visualizations) into a durable, tested, provenance-honest contract, and fix a
per-collection abundance-normalization inconsistency in the DB. Applied from a
prepared draft bundle (`~/Code/research/active/HVP/viroforge-pr-draft/`), reconciled
against the current v0.14.0 state (the bundle targeted 0.13.0).

## What was done

- CI (new, first `.github/` in the repo): `.github/workflows/export-web-data.yml`
  runs `tests/test_export_web_data.py` on PRs/pushes that touch the export script
  or the package. The test builds its own minimal SQLite fixture (the tracked
  `test_viral_genomes.db` has no collections), runs the export as a subprocess,
  and asserts exit 0, both JSON files, the provenance block + `schema_version==1`,
  all four `source_sha256` module hashes as 64-hex, source counts, per-collection
  genome fields, per-collection abundance sum == 1.0, and that sequences are never
  exported.
- Version single-source: grafted a regex `_read_version()` onto the current
  `setup.py` so it reads `__version__` from `viroforge/__init__.py` (kept 0.14.0
  and the `benchmark`/mappy extra; did NOT use the bundle's 0.13.0 setup.py, which
  would have regressed both). Bumped the two stale spots: `CITATION.cff`
  (0.4.0 -> 0.14.0, date 2026-07-16) and README BibTeX (0.12.0 -> 0.14.0). The
  export provenance now stamps 0.14.0 (was 0.1.0 before the earlier bump).
- Abundance renormalizer (new): `scripts/renormalize_abundances.py`, idempotent,
  audit-by-default. Applied `--apply` to the canonical DB: collections 12 (CF
  Respiratory 0.9508), 13 (Human Respiratory RNA 0.8527), and 15 (Fecal RNA
  0.8239) were rescaled to sum to 1.0. Cause: post-hoc `fix_*` scripts removed
  genomes without renormalizing the remainder.

## Validation

- `test_export_web_data.py` passes (0.25s).
- `viroforge.__version__` and `setup.py` both resolve 0.14.0; export against the
  real DB stamps `viroforge_version: 0.14.0`.
- DB backed up to `viral_genomes.db.preRenorm.bak` before `--apply`.
  Post-apply: all 20 collections sum to 1.0 (5 and 6 sit at 0.9992/0.9999 from
  float accumulation over hundreds of tiny abundances, within the 1e-3 tolerance,
  left untouched by design). `PRAGMA integrity_check` = ok; export re-run clean.

## Notes / deferred

The DB is not tracked, so the renormalization is a local data change and does not
enter the PR; re-run the renormalizer as the final step of any future `setup-db`
rebuild. A real release tag should still wait until the DB rebuild settles (per
the PR notes) so a clean, non-dirty SHA can be stamped.
