# Biological Accuracy Review: Per-Site Viral Composition

**Status**: Complete for this pass (Layers 1-2). Bibliography verified against source.
**Date**: 2026-07-16
**Scope**: Viral taxonomic composition + data integrity. Phage host-community
inference (Layer 3) and generation-side tuning are out of scope here.

---

## Question

Do the per-collection abundance proportions of viral taxonomic types (phage vs
eukaryotic virus, nucleic-acid class, dominant/signature families) match current
biological and literature-based expectations for each body site or environment?

## Answer in one paragraph

Eleven of twenty collections are biologically sound. Nine deviate from the
literature: one major (nasopharynx), five moderate (skin, wastewater, vaginal,
blood, lung), and three minor (mouse gut, HIV+ gut, urinary). Separately, a
confirmed data-integrity bug (`genome_type` mislabels every RNA virus as `dsDNA`)
affects all RNA-aware analysis and is the first fix. The recurring biological
problems are: bacteriophages under-represented at sites that should be
phage-dominated (skin, wastewater, vaginal, urinary), a eukaryote-dominant site
modeled as phage-dominant (nasopharynx), a dominant-family rank error (blood:
Orthoherpesviridae over Anelloviridae), and inflated "unknown-family" fractions
that bury signature taxa (mouse gut, lung, HIV+, with lung and CF sharing a
suspicious ~33% placeholder).

## Method

Composition is judged **semi-quantitatively** (rank order, presence/absence of
signature taxa, gross balance), because virome literature rarely reports clean
numeric family percentages. Three choices make it trustworthy:

1. **Properties from ICTV, not the DB's `genome_type`.** A verified family-level
   map (`data/reference_profiles/family_properties.tsv`, built by
   `scripts/build_composition_reference.py` from ICTV VMR MSL40) supplies each
   family's nucleic-acid/Baltimore type and host type (phage if bacteria/archaea,
   else eukaryotic). It covers 108/108 named families across the collections.
2. **Three currencies**: abundance-weighted (primary), genome count, prevalence.
3. **A tunable, cited reference profile** (`virome_composition.yaml`) with per-site
   expected *bands* and a global `strictness` dial;
   `scripts/evaluate_composition.py` scores observed vs expected.

Two "unknown" axes are reported separately and never conflated: **Unknown-family%**
(taxonomically unclassified, a data-quality metric) and **Unknown-host** (named
families ICTV leaves without a host; excluded from the phage:eukaryote ratio).

**Bibliography**: 50 references, every DOI verified via CrossRef and every PMID via
PubMed with `/verify-references`. All 50 resolve to the correct papers; three
publication years were off by one (online-vs-print) and were corrected. See
`docs/composition_references.bib` and `validation/VERIFICATION_SUMMARY.md`.

---

## Finding 1 (confirmed data-integrity bug): `genome_type` mislabels RNA viruses as dsDNA

`genomes.genome_type` marks RNA-virus genomes as `dsDNA`. The arbovirus collection
(14) is the clearest case: its genomes are correct arboviruses (Dengue, Yellow
fever, West Nile, La Crosse, Saint Louis encephalitis, Mayaro; Flaviviridae,
Peribunyaviridae, Togaviridae) but every one is stamped `dsDNA`. RNA collections
13 and 15 are ~98% / 90% `dsDNA` by this column. The field is populated (non-NULL)
for all 14,423 genomes, so the error is silent, and any analysis trusting it is
wrong for RNA viruses. This review derives nucleic-acid type from ICTV instead.
**Remediation**: relabel `genome_type` from the ICTV Baltimore class per family
(clear-cut, backed up, reversible).

---

## Scorecard

Full detail: `validation/composition_scorecard.json`. Figure:
`docs/figures/composition_review.png` (Panel A: observed phage fraction vs
expected band; Panel B: unknown-family fraction vs ceiling).

| ID | Collection | Verdict | Key discrepancies |
|----|-----------|---------|-------------------|
| 4 | Respiratory - Nasopharynx (Healthy) | MAJOR | Anelloviridae not dominant; phage-dominant where eukaryote-dominant expected |
| 3 | Skin - Sebaceous Sites (Healthy) | MODERATE | phage 0.00 vs band 0.40-0.90; unknown-fam 0.51 > 0.25 ceiling |
| 9 | Wastewater - Urban Treatment Plant | MODERATE | phage 0.35 vs band 0.60-0.95; RNA 0.45 above band |
| 16 | Vaginal (Healthy) | MODERATE | phage 0.61 vs band 0.75-0.98 (eukaryotic over-weighted) |
| 17 | Blood/Plasma (Healthy) | MODERATE | Anelloviridae present but not dominant (top is Orthoherpesviridae) |
| 19 | Lower Respiratory (Lung, Healthy) | MODERATE | Anelloviridae not dominant; unknown-fam 0.34 > 0.20 |
| 8 | Mouse Gut - Laboratory (C57BL/6) | MINOR | phage slightly below band; high unknown-host |
| 11 | HIV+ Gut | MINOR | phage slightly below band; unknown buries signature taxa |
| 20 | Urinary (Healthy) | MINOR | phage below band (bacteriophages under-represented) |
| 1 | Gut - Adult Healthy (Western Diet) | OK | — |
| 2 | Oral - Saliva (Healthy) | OK | (Microviridae top-rank nuance, see below) |
| 5 | Marine - Coastal Surface Water | OK | — |
| 6 | Soil - Agricultural | OK | — |
| 7 | Freshwater - Lake Surface Water | OK | — |
| 10 | IBD Gut | OK | — |
| 12 | Cystic Fibrosis Respiratory | OK | high unknown is genuinely expected here |
| 13 | Human Respiratory RNA | OK | (after genome_type fix) |
| 14 | Arbovirus (Mosquito) | OK | (after genome_type fix) |
| 15 | Fecal RNA | OK | (after genome_type fix) |
| 18 | Ocular Surface (Healthy) | OK | — |

---

## Flagged sites (with cited expectations)

**Nasopharynx (4) - MAJOR.** Observed phage-dominant (phage 50% > euk 35%), top
family a phage (Aliceevansviridae). The healthy upper-respiratory virome is
eukaryote-dominated and led by Anelloviridae, with bacteriophages comparatively
low in health (Wang 2016; Megremis 2023; Xu 2017). The balance is inverted and the
dominant family is wrong. Note the strongest sources are pediatric/asthma cohorts;
adult healthy-nasopharynx data are sparser, but Anellovirus dominance is robust.

**Skin (3) - MODERATE.** Observed 0% phage, 51% unknown-family. The healthy
sebaceous skin virome is bacteriophage-dominated (Cutibacterium/Staphylococcus/
Corynebacterium phages) with Papillomaviridae/Polyomaviridae a minority (Hannigan
2015; Graham 2023). Zero phage is not credible; skin phages appear missing from the
collection, and the dark matter that is present is itself overwhelmingly phage.

**Wastewater (9) - MODERATE.** Observed 65% eukaryotic / 45% RNA. Urban wastewater
viromes are strongly bacteriophage-dominated by richness and abundance (~60-95% of
classifiable reads), with enteric eukaryotic viruses typically ~1-2% (Gulino 2020;
Kuo 2023; Calusinska 2016; Martinez-Puchol 2020; Cantalupo 2011). ViroForge
correctly includes the right enteric eukaryotic families but inverts their
proportion. Internal tell: the collection's own top family is Microviridae (an
ssDNA phage), inconsistent with phage totalling only 35%.

**Vaginal (16) - MODERATE.** Observed phage fraction 0.61 of classified; eukaryotic
viruses ~25% where literature reports ~4% (Jakobsen 2020; Kaelin 2022). Top family
Herelleviridae (large Firmicutes phages) is the documented dairy/food-Lactobacillus
proxy artifact, not the temperate Lactobacillus siphophages expected in a healthy
Lactobacillus-CST vaginal virome. This is a known limitation (RefSeq lacks vaginal
Lactobacillus-phage genomes); documenting it precisely is the action here.

**Blood/plasma (17) - MODERATE.** Structure correct (eukaryote-dominant, DNA-based,
phage ~0), but the dominant family is wrong: healthy plasma is Anelloviridae-
dominated (~97% of viral reads), herpesviruses secondary (Cebria-Mendoza 2021,
2023; Segura-Wang 2018). Observed top is Orthoherpesviridae; Anelloviridae is
present but not the most abundant.

**Lung (19) - MODERATE.** Eukaryote-dominant direction is right, but the top family
should be Anelloviridae/Torque teno virus (Young 2015; Abbas 2017), not Unknown.
Its 33.6% unknown-family is implausibly high for the comparatively well-defined
healthy lung and is near-identical to CF's 33.4% - a tell of a shared placeholder
default rather than a site-specific value.

**Mouse gut (8), HIV+ gut (11), urinary (20) - MINOR.** Mouse gut is phage-dominated
in reality (Caudovirales + Microviridae; Bao 2022; Moltzau Anderson 2023), but the
observed 37% phage / 59% unknown-host suggests phages sitting unlabeled in the dark
pool. HIV+ gut correctly shows eukaryotic-virus expansion (Monaco 2016; Boukadida
2024), but 29% unknown-family buries the diagnostic Adenoviridae/Anelloviridae/
Papillomaviridae signature. Urinary correctly tops with Papillomaviridae (Santiago-
Rodriguez 2015; Garretto 2018; Maqsood 2024) but under-represents bacteriophages.

## Sites that pass

Gut (1: Pargin 2023; Van Espen 2021), oral (2: Pride 2012; Ly 2014 - with a
top-rank nuance that saliva is led by tailed Caudovirales, not Microviridae), marine
(5), soil (6), freshwater (7), IBD (10: Caudovirales-expansion signature, Norman
2015; Zuo 2019), CF (12: high unknown genuinely expected, Willner 2009), the three
RNA collections (13/14/15) once `genome_type` is fixed, and ocular (18: Anellovirus-
led, Doan 2016; Siegal 2021).

## Cross-cutting themes

1. **Bacteriophages under-represented at phage-dominated sites** (skin, wastewater,
   vaginal, urinary) - the single most common failure.
2. **Inflated unknown-family fractions bury signature taxa** (mouse gut, lung,
   HIV+), and at least one shared placeholder (~33%) is reused across CF and lung.
3. **A few dominant-family rank errors** (blood, lung, nasopharynx, oral).

---

## Remediation

Applied this pass (clear-cut, backed up, reversible; see the accompanying
lab-notebook entry):
- **`genome_type` relabel** from ICTV Baltimore class (fixes Finding 1).

Recommended follow-ups (larger or judgment-dependent, tracked separately):
- Re-query **skin (3)** to add Cutibacterium/Staphylococcus/Corynebacterium phages.
- Rebalance **wastewater (9)** toward phage-majority / DNA-majority.
- Replace the shared ~33% unknown default in **lung (19)** with a site-specific
  value and add Anelloviridae as the dominant family; same for **nasopharynx (4)**.
- Trim/annotate the dark pool in **mouse gut (8)** and **HIV+ gut (11)** so phages
  and diagnostic eukaryotic viruses are labeled.
- Re-rank **blood (17)** and **oral (2)** dominant families.
- Keep **vaginal (16)** documented as a proxy limitation until vaginal Lactobacillus
  phage genomes are available.

---

## Artifacts and reproduction

```
scripts/build_composition_reference.py   # ICTV family property map
scripts/evaluate_composition.py observe  # -> validation/observed_composition.json
scripts/evaluate_composition.py evaluate # -> validation/composition_scorecard.json
scripts/build_bibliography.py            # -> docs/composition_references.bib
scripts/plot_composition_review.py       # -> docs/figures/composition_review.{png,pdf}
data/reference_profiles/family_properties.tsv     # verified property map
data/reference_profiles/virome_composition.yaml   # tunable cited profile (strictness dial)
docs/composition_references.bib                    # verified bibliography
validation/dossiers/*.json                         # per-site literature dossiers
validation/VERIFICATION_SUMMARY.md                 # /verify-references outcome
```

Adjust the yardstick with the `strictness` dial or any per-site band in
`virome_composition.yaml`, then re-run `evaluate_composition.py evaluate`.
