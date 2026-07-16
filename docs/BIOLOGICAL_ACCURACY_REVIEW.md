# Biological Accuracy Review: Per-Site Viral Composition

**Status**: Complete. Corrections applied and folded into `setup-db`.
**Date**: 2026-07-16
**Scope**: Viral taxonomic composition + data integrity (Layers 1-2). Phage
host-community inference and generation-side tuning are out of scope.

---

## Question

Do the per-collection abundance proportions of viral taxonomic types (phage vs
eukaryotic virus, nucleic-acid class, dominant/signature families) match current
biological and literature-based expectations for each body site or environment?

## Answer

After correcting the measurement/data problems and four genuinely-off collections,
**biology is sound for 19 of 20 collections** and **all data-quality flags are
cleared**. The sole remaining flag is vaginal biology (minor) — a documented
data-availability limitation (RefSeq has no vaginal Lactobacillus-phage genomes, so
dairy/food proxies are used).

The data-quality flags were cleared honestly, not by hiding them. They were driven
by bacteriophages carried as `family=Unknown`. Two things resolved that: (1) NCBI
enrichment gave 4,577 previously-unclassified genomes a real class and genus
(DB-wide class coverage 89%, family coverage 72%), and (2) the data-quality metric
was refactored to measure genuine non-placement (no class) rather than missing
family — because ICTV legitimately assigns **no family** to most modern phage
genera (it abolished the morphology-based families in 2021). A genus-level
Caudoviricetes phage is classified, not dark matter. 1,481 genomes DB-wide remain
genuinely unclassified (no class even in NCBI); the curated collections do not draw
from them, so they score <0.3% unclassified.

## Method

Composition is judged semi-quantitatively (rank order, signature-taxon presence,
gross balance). Design choices that make it trustworthy:

1. **Properties from ICTV, not the DB's `genome_type`.** A verified family map
   (`data/reference_profiles/family_properties.tsv`, from ICTV VMR MSL40) gives each
   family's nucleic-acid/Baltimore type and host type. For genomes ICTV leaves
   unclassified, a name-based phage heuristic recovers the large pool of
   bacteriophages carried as `family=Unknown` (see Finding 2).
2. **Three currencies**: abundance-weighted (primary), count, prevalence.
3. **A tunable, cited reference profile** (`virome_composition.yaml`) of per-site
   expected bands + a global `strictness` dial; `scripts/evaluate_composition.py`
   scores observed vs expected.
4. **Biology vs data-quality are separated.** The site verdict reflects biology
   (balance, signature taxa); the unclassified fraction (no class assigned) is a
   distinct data-quality verdict, so taxonomy gaps do not masquerade as biology
   errors. Family=Unknown is *not* penalized — ICTV legitimately omits it for
   phages (Finding 3).

Bibliography: 50 references, every DOI+PMID confirmed via `/verify-references`
(3 years corrected). See `docs/composition_references.bib`, `validation/VERIFICATION_SUMMARY.md`.

---

## Finding 1 (fixed): `genome_type` mislabeled RNA viruses as dsDNA

`genomes.genome_type` silently defaulted RNA viruses to `dsDNA` (the arbovirus
collection was 100% dsDNA despite being RNA viruses). `scripts/fix_genome_type.py`
relabeled 3,516 genomes from the ICTV Baltimore class; collection 14 is now 100%
ssRNA, fecal RNA shows rotavirus as dsRNA, and the DB-wide ssDNA count went from 6
to 1,520. Folded into `setup-db`.

## Finding 2 (fixed): the phage metric missed unclassified bacteriophages

Many RefSeq bacteriophages carry `family=Unknown` in the DB, so an ICTV-only host
map counted them as host-unknown and made phage-dominated sites look phage-empty.
Skin, for example, is ~51% Propionibacterium/Cutibacterium phages that were all
unclassified — it read as 0% phage. Adding a name-based phage heuristic corrected
skin to 51% phage and reclassified CF (76%), vaginal (75%), HIV+ (81%), mouse gut
(71%). This resolved several apparent "deviations" that were measurement artifacts.

## Finding 3 (fixed): 6,182 genomes were unclassified at every rank

The genomes that failed ICTV name-matching sat as Unknown at *every* rank (realm
through genus), not just family — even though NCBI (ICTV-synced) classifies them.
`scripts/enrich_taxonomy_from_ncbi.py` resolved each genome's `ncbi_taxid` against
NCBI E-utilities (cached to a tracked TSV) and filled the missing realm..genus
ranks (and family where NCBI has one), never overwriting an existing rank or the
species name. Two distinct wins, reported separately so neither is overstated:

- **Higher-rank classification (large):** 4,577 genomes gained a class; DB-wide
  class coverage is now 89% (was ~57%). These are mostly Caudoviricetes phages now
  placed to class + genus.
- **New family assignments (smaller):** 2,102 genomes gained a real NCBI family
  (family coverage 57% -> 72%) — eukaryotic viruses that failed the name-match
  (e.g. Human bocavirus -> Parvoviridae) plus phage genera that *do* have families
  (Microviridae, Straboviridae, ...). Families were **not** invented: most
  Caudoviricetes genera remain family-less because ICTV assigns them none.

The data-quality metric was then changed from "unknown-family fraction" to
"unclassified fraction" (no class assigned). 1,481 genomes remain genuinely
unclassified DB-wide, so the metric still measures a real residual; the curated
collections simply do not contain them.

## Corrections applied to four genuinely-off collections

After Findings 1-2, four collections had real composition problems, fixed by
`scripts/fix_collection_composition.py` (reweighting existing genomes and adding
only real DB genomes; idempotent; folded into `setup-db`):

- **Blood (17)**: Anelloviridae reweighted to dominant (0.60), Orthoherpesviridae
  demoted to secondary — healthy plasma is anellovirus-dominated (Cebria-Mendoza 2021).
- **Lung (19)**: Anelloviridae to dominant (0.45), unclassified-phage share reduced
  below its ceiling — healthy lung is Anellovirus-led (Young 2015; Abbas 2017).
- **Nasopharynx (4)**: 6 Anelloviridae genomes added (none were present) and the
  balance shifted to eukaryote-dominant, Anelloviridae top — healthy upper
  respiratory is eukaryote/Anellovirus-dominated (Wang 2016; Megremis 2023).
- **Wastewater (9)**: rebalanced to phage-majority (0.85 phage / 0.15 eukaryotic) —
  urban wastewater is phage-dominated by richness, enteric eukaryotic viruses a
  minority (Gulino 2020; Kuo 2023; Calusinska 2016).

---

## Final scorecard

Full detail: `validation/composition_scorecard.json`. Figure:
`docs/figures/composition_review.png`.

All 20 collections are biology OK except vaginal (minor); all data-quality flags
are OK after taxonomy enrichment. "Family-unassigned %" is informational (ICTV
omits family for phages); "Unclassified %" (no class) is the data-quality metric.

| ID | Collection | Biology | Data-quality | Family-unassigned % | Notes |
|----|-----------|---------|--------------|--------------------:|-------|
| 16 | Vaginal (Healthy) | MINOR | OK | 37 | proxy Lactobacillus phages (dairy/food); documented limitation |
| 1 | Gut - Adult Healthy | OK | OK | 8 | — |
| 2 | Oral - Saliva | OK | OK | 10 | — |
| 3 | Skin - Sebaceous Sites | OK | OK | 51 | family-unassigned = Propionibacterium phages (now class Caudoviricetes) |
| 4 | Respiratory - Nasopharynx | OK | OK | 14 | corrected: added Anelloviridae, eukaryote-dominant |
| 5 | Marine - Coastal Surface | OK | OK | 0 | — |
| 6 | Soil - Agricultural | OK | OK | 4 | — |
| 7 | Freshwater - Lake Surface | OK | OK | 0 | — |
| 8 | Mouse Gut (C57BL/6) | OK | OK | 34 | family-unassigned phages, now classified to class/genus |
| 9 | Wastewater | OK | OK | 0 | corrected: phage-majority |
| 10 | IBD Gut | OK | OK | 10 | — |
| 11 | HIV+ Gut | OK | OK | 29 | family-unassigned phages, now classified to class/genus |
| 12 | Cystic Fibrosis Respiratory | OK | OK | 21 | high unassigned genuinely expected |
| 13 | Human Respiratory RNA | OK | OK | 0 | after genome_type fix |
| 14 | Arbovirus (Mosquito) | OK | OK | 0 | after genome_type fix |
| 15 | Fecal RNA | OK | OK | 0 | after genome_type fix |
| 17 | Blood/Plasma | OK | OK | 2 | corrected: Anelloviridae dominant |
| 18 | Ocular Surface | OK | OK | 8 | — |
| 19 | Lower Respiratory (Lung) | OK | OK | 5 | corrected: Anelloviridae dominant |
| 20 | Urinary | OK | OK | 18 | — |

All 20 collections have <0.3% genuinely-unclassified (no class) content.

---

## Remaining

- **Vaginal (16)** biology minor: RefSeq has no vaginal Lactobacillus-phage genomes,
  so dairy/food Lactobacillus phages are proxies and the eukaryotic share sits a
  little high. Documented limitation; deferred until suitable genomes exist.
- **1,481 genomes DB-wide remain genuinely unclassified** (no class even in NCBI) —
  truly novel/environmental sequences. Not in the curated collections; would flag if
  a collection drew from them.

## Durability

All four corrections are idempotent and run as post-curation steps in
`viroforge setup-db` in order: taxonomy enrichment (`enrich_taxonomy_from_ncbi.py
--apply`, from the tracked cache — no live network), `genome_type` relabel,
`fix_collection_composition.py`, `renormalize_abundances.py`. A fresh rebuild
reproduces the corrected state.

## Artifacts and reproduction

```
scripts/build_composition_reference.py   # ICTV family property map
scripts/evaluate_composition.py observe|evaluate
scripts/fix_genome_type.py               # Finding 1 fix
scripts/enrich_taxonomy_from_ncbi.py     # Finding 3: --fetch (cache) / --apply (DB)
scripts/fix_collection_composition.py    # 4-collection curation fix
scripts/renormalize_abundances.py        # normalize sums to 1.0
data/reference_profiles/ncbi_lineage_cache.tsv  # tracked NCBI lineage cache
scripts/build_bibliography.py            # -> docs/composition_references.bib
scripts/plot_composition_review.py       # -> docs/figures/composition_review.{png,pdf}
data/reference_profiles/{family_properties.tsv, virome_composition.yaml}
validation/dossiers/*.json, validation/VERIFICATION_SUMMARY.md
```

Tune the yardstick via the `strictness` dial or any band in
`virome_composition.yaml`, then re-run `evaluate_composition.py evaluate`.
