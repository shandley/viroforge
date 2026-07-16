# NCBI Higher-Rank Taxonomy Enrichment (clear data-quality flags)

**Date**: 2026-07-16
**Session Type**: TAXONOMY
**Status**: Complete

---

## Objective

Clear the remaining data-quality flags (unclassified phages in skin/mouse/HIV+/
vaginal) by improving taxonomy coverage - without fabricating families.

## Key finding

The 6,182 `family=Unknown` genomes were Unknown at EVERY rank (realm..genus), not
just family. NCBI classifies them to class/genus, but assigns NO family to most
Caudoviricetes phage genera because ICTV abolished the morphology-based phage
families (Siphoviridae/Myoviridae/Podoviridae) in 2021. So `family=Unknown` is
correct current taxonomy for these phages; the fix is to fill the ranks ICTV/NCBI
actually provide, not to invent families. (Verified from primary source: efetch
shows Propionibacterium phage P100_A as class Caudoviricetes -> genus Pahexavirus,
no family.)

## What was done

- **`scripts/enrich_taxonomy_from_ncbi.py`**: batched NCBI efetch (taxonomy) over
  the unclassified taxids (`--fetch`), cached taxid->lineage to a tracked TSV
  (`data/reference_profiles/ncbi_lineage_cache.tsv`), then `--apply` fills only
  EMPTY ranks in the DB (never overwrites an existing rank or `species`). XML
  parsed with ElementTree (nested LineageEx). Idempotent; DB backed up
  (`viral_genomes.db.preTaxEnrich.bak`).
  - Higher-rank win (large): 4,577 genomes gained a class; DB class coverage ~89%.
  - Family win (smaller): 2,102 genomes gained a real NCBI family (family coverage
    57% -> 72%) - eukaryotic viruses that missed the name-match (e.g. bocavirus ->
    Parvoviridae) + phage genera that do have families. Families not invented.
- **Metric refactor** (`scripts/evaluate_composition.py`): data-quality metric
  changed from unknown-FAMILY to UNCLASSIFIED (no class). unknown_family kept as
  informational. Profile field `max_unknown_family` -> `max_unclassified`.
- **Folded into setup-db**: enrichment `--apply` (from tracked cache, offline) runs
  first in the post-curation corrections block.
- Updated report + figure (Panel B now shows family-unassigned vs genuinely
  unclassified overlay).

## Result

Biology 19/20 OK (only vaginal minor). All data-quality flags cleared. 1,481 genomes
DB-wide remain genuinely unclassified (no class) - the metric still measures this
real residual; curated collections have <0.3% unclassified. Guard against a vacuous
metric satisfied: the residual exists and would flag if a collection drew from it.

## Validation

- Spot-checked: bocavirus -> Parvoviridae; Pahexavirus phage -> Caudoviricetes,
  family blank. Integrity ok. Enrichment idempotent (re-run fills 0).
- setup_db.py compiles; enrichment applies from cache with no network.

See [[composition_review]] and entries 008/009.
