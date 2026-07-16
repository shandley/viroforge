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

After correcting two measurement/data problems and four genuinely-off collections,
**biology is sound for 19 of 20 collections**. The sole remaining biology flag is
vaginal (minor) — a documented data-availability limitation (RefSeq has no vaginal
Lactobacillus-phage genomes, so dairy/food proxies are used). Four collections
carry data-quality (taxonomy-coverage) flags: skin, mouse gut, HIV+ gut, and vaginal
have large fractions of bacteriophages that are real but sit unclassified
(`family=Unknown`) in the database. That is a DB-wide taxonomy-coverage issue
(42.9% of genomes are unclassified), not a composition error, and is reported
separately from the biology verdict.

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
   (balance, signature taxa); the unknown-family fraction is reported as a distinct
   data-quality verdict so taxonomy gaps do not masquerade as biology errors.

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

| ID | Collection | Biology | Data-quality | Notes |
|----|-----------|---------|--------------|-------|
| 16 | Vaginal (Healthy) | MINOR | MINOR | proxy Lactobacillus phages (dairy/food); documented limitation |
| 1 | Gut - Adult Healthy | OK | OK | — |
| 2 | Oral - Saliva | OK | OK | — |
| 3 | Skin - Sebaceous Sites | OK | MODERATE | 51% unclassified Propionibacterium phages (taxonomy coverage) |
| 4 | Respiratory - Nasopharynx | OK | OK | corrected: added Anelloviridae, eukaryote-dominant |
| 5 | Marine - Coastal Surface | OK | OK | — |
| 6 | Soil - Agricultural | OK | OK | — |
| 7 | Freshwater - Lake Surface | OK | OK | — |
| 8 | Mouse Gut (C57BL/6) | OK | MINOR | 37% unclassified phages (taxonomy coverage) |
| 9 | Wastewater | OK | OK | corrected: phage-majority |
| 10 | IBD Gut | OK | OK | — |
| 11 | HIV+ Gut | OK | MINOR | 29% unclassified phages (taxonomy coverage) |
| 12 | Cystic Fibrosis Respiratory | OK | OK | high unknown genuinely expected |
| 13 | Human Respiratory RNA | OK | OK | after genome_type fix |
| 14 | Arbovirus (Mosquito) | OK | OK | after genome_type fix |
| 15 | Fecal RNA | OK | OK | after genome_type fix |
| 17 | Blood/Plasma | OK | OK | corrected: Anelloviridae dominant |
| 18 | Ocular Surface | OK | OK | — |
| 19 | Lower Respiratory (Lung) | OK | OK | corrected: Anelloviridae dominant |
| 20 | Urinary | OK | OK | — |

---

## Remaining, and why not "fixed"

- **Data-quality flags (skin, mouse gut, HIV+ gut, vaginal)** are unclassified
  bacteriophages (`family=Unknown`). They are real phages and biologically correct
  to include; the flag reflects the DB's taxonomy-coverage gap, addressable only by
  improving classification, not by curation. Left visible rather than hidden.
- **Vaginal (16)** biology minor: RefSeq has no vaginal Lactobacillus-phage genomes,
  so dairy/food Lactobacillus phages are proxies and the eukaryotic share sits a
  little high. Documented limitation; deferred until suitable genomes exist.

## Durability

All three corrections are idempotent and run as post-curation steps in
`viroforge setup-db` (`genome_type` relabel, `fix_collection_composition.py`,
`renormalize_abundances.py`), so a fresh rebuild reproduces the corrected state.

## Artifacts and reproduction

```
scripts/build_composition_reference.py   # ICTV family property map
scripts/evaluate_composition.py observe|evaluate
scripts/fix_genome_type.py               # Finding 1 fix
scripts/fix_collection_composition.py    # 4-collection curation fix
scripts/renormalize_abundances.py        # normalize sums to 1.0
scripts/build_bibliography.py            # -> docs/composition_references.bib
scripts/plot_composition_review.py       # -> docs/figures/composition_review.{png,pdf}
data/reference_profiles/{family_properties.tsv, virome_composition.yaml}
validation/dossiers/*.json, validation/VERIFICATION_SUMMARY.md
```

Tune the yardstick via the `strictness` dial or any band in
`virome_composition.yaml`, then re-run `evaluate_composition.py evaluate`.
