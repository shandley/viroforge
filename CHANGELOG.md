# Changelog

Notable changes to ViroForge. This file starts at 0.15.0; for earlier history,
see the git log.

Versions follow [semantic versioning](https://semver.org/), with the caveat that
ViroForge generates data. A change that alters generator output for a fixed seed
is treated as breaking even when the API is untouched, and is called out under
"Reproducibility" below.

## [0.20.0] - 2026-07-30

### Reproducibility

**Every dataset using `linker` or `rdab` amplification changes, substantially.**
`_apply_linker_bias` and `_apply_rdab_bias` raised GC efficiency to the power of
the cycle count (20 and 40). That function returns a total relative efficiency,
not a per-cycle one, so the bias compounded absurdly: a 36% GC genome was
suppressed **134x** against a 50% GC one, and a 29% GC one by **61,000x**. 36%
is the median GC of the gut collection, so this distorted the relative abundance
of every default-amplification dataset. `linker` is the default. MDA was already
correct and is unaffected.

Nothing measured is remotely that strong. Parras-Moltó et al. 2018 report 6-7%
of contigs past 10x for MDA, the harshest method available.

### Added

- **Fungal and archaeal background**, completing issue #37. `--fungal-fraction`
  and `--archaeal-fraction` with per-collection baselines, plus bundled
  reference sets (15 fragments, 150 KB each, all 10 taxa resolved).

  Measured on collection 1, `--no-vlp --contamination-level heavy`, no flags:

  | source | % |
  |---|---|
  | bacterial | 71.8 |
  | host | 12.7 |
  | rRNA | 7.8 |
  | viral | 4.1 |
  | fungal | 1.5 |
  | archaeal | 1.1 |

  Fungal plus archaeal at 2.5% against a real bulk-stool 1-5%.

  The three domains are separate contaminant types, not one bucket: a classifier
  can handle one and miss another, and merging them would hide that in a QC
  benchmark.

### Fixed

- Background fragment count now scales with abundance. A fixed 200 gave a 2%
  domain 0.0001 per fragment, which rounds to zero reads at any realistic depth,
  so the mycobiome and archaeome never appeared in output at all.
- The curation script no longer selects organelle genomes. The *Trichoderma
  reesei* search returned its 42 kb mitochondrion as the top hit, which has
  atypical base composition and would supply only four fragments.

## [0.19.0] - 2026-07-30

### Fixed

- **Every bundled preset pointed at the wrong collection.** The 9-28 to 1-20
  renumbering in PR #6 never updated them, so `--preset gut-standard` and
  `--preset gut-bulk` generated a **wastewater** virome, `marine-standard`
  generated respiratory RNA, and `respiratory-rna` pointed at collection 21,
  which does not exist. All eight are corrected (`new = old - 8`).

  The wastewater case is the dangerous one: collection 9 exists, so generation
  succeeded and silently produced the wrong dataset. Anyone who benchmarked
  against a gut preset since 2026-07-14 was benchmarking against wastewater.

- **`hybrid-standard` could never run through `viroforge generate`.** It sets
  `short_platform` and `long_platform`, which belong to
  `scripts/generate_hybrid_dataset.py`. Presets can now declare
  `metadata.requires_script`, and the CLI prints what to run instead of failing
  on unrecognised parameters.

### Added

- **Sample-type presets** completing issue #37: `gut-vlp` (viral-dominated,
  96.7% viral-origin) and `stool-clinical` (host-dominated), alongside the
  corrected `gut-bulk`. `gut-bulk` and `gut-vlp` share a collection so the pair
  isolates the effect of enrichment.

  Measured against the compositions proposed in the issue:

  | preset | host | bacterial | viral |
  |---|---|---|---|
  | `gut-bulk` | 16.4% | 67.5% | 7.9% |
  | `stool-clinical` | 40.0% | 36.5% | 16.0% |
  | `gut-vlp` | 0.0% | 1.8% | 96.7% |

- `--host-fraction` pins host DNA absolutely, bypassing the
  `--contamination-level` multiplier. Needed for host-dominated sample types:
  40% host is unreachable from the gut collection's 5% baseline even at the
  highest level. Mirrors `--bacterial-fraction`.

- `tests/test_presets.py` checks every preset loads, resolves to a real
  collection, points at a collection matching its own name, and converts to a
  valid generator namespace. The stale IDs survived two years because nothing
  loaded a preset and checked where it pointed.

## [0.18.0] - 2026-07-30

### Reproducibility

**Bacterial background now applies by default, so every collection's output
changes.** Each collection carries a `default_bacterial_pct` baseline and it is
used unless `--bacterial-fraction` says otherwise. Pass `--bacterial-fraction 0`
to restore 0.17.0 behaviour.

### Added

- **Per-collection bacterial baselines.** `default_bacterial_pct` and
  `bacterial_community` columns on `body_site_collections`, populated from
  `data/reference_profiles/contamination_defaults.tsv` by the existing
  `populate_contamination_defaults.py` (already folded into `setup-db`, so
  rebuilds carry them).

  The baseline is the sample's microbiome **before** any VLP step, so it tracks
  microbial biomass per site rather than prep quality: soil 80%, wastewater and
  marine 75%, gut 70%, vaginal 60%, down to lung 20%, ocular 10% and blood 1%.
  Anchored to Sender et al. 2016, PLoS Biol (PMID 27541692) and the Human
  Microbiome Project Consortium 2012, Nature (PMID 22699609), both verified.

  One number serves both workflows because VLP does the physical work. Measured
  on collection 1 with no flags at all:

  |  | default (VLP) | `--no-vlp` |
  |---|---|---|
  | bacterial | 1.8% | 74.3% |
  | viral-origin | 96.7% | 15.8% |

  A default VLP run keeps a realistic few-percent bacterial residual, which real
  VLP preps do carry; `--no-vlp` keeps the full bulk-metagenome background.

  RNA collections (13, 14, 15) are set to 0.0: bacterial background is wired
  into the DNA contamination path only, and an RNA library's bacterial signal is
  ribosomal, already covered by `rrna_pct`.

### Changed

- `--bacterial-fraction` now defaults to `None` rather than `0.0`, so an
  explicit `0` can switch the feature off. With a `0.0` default there was no way
  to distinguish "not given" from "off", and a nonzero baseline could not be
  disabled.
- `populate_contamination_defaults.py` builds its UPDATE statement from its
  `COLUMNS` map instead of a hand-written column list. The old list silently
  left both new columns NULL on the first run.

## [0.17.0] - 2026-07-29

### Added

- **Bacterial background**, so `--no-vlp` produces a realistic bulk metagenome.
  `--bacterial-fraction` adds the sample's own microbiome and
  `--bacterial-community` selects which one (ten profiles; soil is GC-rich,
  marine is not). Off by default. Measured on the IBD gut collection with
  `--bacterial-fraction 0.70`, same seed:

  |            | `--no-vlp` | `--vlp tangential` | real bulk stool |
  |---|---|---|---|
  | bacterial  | 68.9% | 1.2%  | 60-80% |
  | viral      | 8.0%  | 96.9% | 1-5%   |
  | host       | 9.8%  | 0.2%  | 10-30% |
  | rRNA       | 10.4% | 1.2%  | 5-15%  |

  The same run was 81.6% viral and 1.7% bacterial in 0.16.0. VLP enrichment now
  removes 98.2% of the background, so it finally models what VLP prep is for.
  Closes the core of issue #37.

  Bacterial background is a separate contaminant type from reagent bacteria, and
  is not scaled by `--contamination-level`. That dial models how well the prep
  went; this models how much microbiome the sample held.

- `scripts/curate_bacterial_background.py` builds the reference set from RefSeq
  at build time, resolving genomes by taxon **name** so no accession is stored in
  the repository, and fetching only 10 kb slices so whole genomes are never
  downloaded. 9 gut taxa produce 18 fragments at 180 KB.

### Reproducibility

**Contamination fractions are now delivered as requested, which shifts existing
output slightly.** Previously the viral and contaminant blocks were concatenated
and the total normalised, so a requested fraction `c` arrived as `c/(1+c)`: a
requested 8.6% landed at 7.9%. The viral community is now scaled to fill
`1 - contamination` instead. Existing datasets shift by well under a percentage
point; a requested 70% would have landed at 41%.

### Changed

- `--dark-matter-fraction` is documented as a fraction **of the viral portion**,
  not of all reads. This was always true of the implementation but only mattered
  once contamination could be large: at 70% bacterial, a 0.30 dark-matter
  fraction is 9% of reads.

- **The reference set now ships bundled**:
  `viroforge/data/references/bacterial_fragments.fasta`, 120 fragments of 10 kb
  across all 10 communities, 1.2 MB, 0.017% ambiguous bases. All 40 taxa resolved
  against RefSeq. `--bacterial-fraction` works out of the box; the synthetic
  fallback remains for anyone replacing the set via
  `VIROFORGE_BACTERIAL_FRAGMENTS`.

  The recommended bulk stool recipe is `--no-vlp --contamination-level heavy
  --bacterial-fraction 0.75`, measured at 75.0% bacterial, 3.6% viral, 11.8%
  rRNA, 9.5% host. Bacterial fraction alone does not push viral into the 1-5%
  band; host and rRNA have to occupy the remainder, which the heavy level
  provides.

### Changed

- `BACTERIAL_COMMUNITY_PROFILES` GC values are now **measured** from the bundled
  reference rather than estimated, so the synthetic fallback matches the real
  set. Several were well off: wastewater 55 to 44.4, urinary 45 to 36.9, gut 45
  to 48.6. A test remeasures the reference and fails if the two drift apart.

### Known issues

- Per-collection bacterial baselines (a `default_bacterial_pct` column) are not
  implemented, so the fraction must be given explicitly per run.
- Marine uses *Synechococcus elongatus* rather than a marine strain such as
  WH 8102, because strain-qualified names do not resolve through an `[Organism]`
  search. The genus and GC profile are representative; the ecology is not exact.

## [0.16.0] - 2026-07-29

### Reproducibility

**MDA-amplified datasets change.** Runs using `--amplification mda` or
`mda-long`, and PCR duplicate injection under MDA, produce different output than
0.15.0 for the same seed. Other amplification methods are unaffected.

### Fixed

- **MDA GC bias peaked at the wrong GC content.** The φ29 model centred its
  efficiency curve at 40% GC, so it penalised genomes at 45-60% GC by up to 6x.
  Parras-Moltó et al. 2018, Microbiome 6:119 (PMID 29954453) measured the
  opposite in saliva DNA viromes: MDA *over*-amplifies contigs at 45-60% GC and
  under-represents both extremes. The optimum moves to 50% GC in
  `viroforge/amplification.py` and in the MDA branch of
  `viroforge/simulators/duplicates.py`. RdAB and Linker already used 50%.

  A 60% GC genome went from a 6.05x penalty to 1.57x; a 30% GC genome now takes
  the 6.05x penalty it previously escaped. Low-GC-dominated collections (gut,
  median 36% GC) shift more than mid-GC ones.

### Added

- `scripts/benchmark_amplification_bias.py` calibrates the MDA GC model against
  the same paper, which reports 6.2-7.6% of contigs exceeding 10x or 0.1x fold
  change under MDA. It reports, for real collections, the fraction of genomes
  past that threshold, the fold-bias across the published over-amplified band,
  and the shift in contaminant efficiency. Defaults are read from the shipped
  model rather than restated, so the script cannot drift from what it measures.

### Notes

`gc_bias_strength` stays at 3.0. The benchmark does not support lowering it: at
3.0 the model puts 3.3% of genomes past 10x against a published 4.5-7.6%, while
1.5 puts none there at all, producing less bias than measured. The parameter
that was wrong was the peak position, not the strength.

### Known issues

- MDA and RdAB now produce nearly identical GC efficiency curves at their
  defaults (within 2% at 30% and 70% GC). `MDAAmplification` is documented as
  having 2-5x stronger bias than PCR, and most of that difference previously
  came from the mis-centred optimum rather than from `gc_bias_strength`.
  Whether φ29 should fall off more sharply than RdAB needs its own measurement;
  `test_mda_and_rdab_curves_are_currently_similar` pins the current state so the
  question is not lost.

## [0.15.0] - 2026-07-28

### Reproducibility

**Datasets generated with 0.14.0 and earlier will not reproduce byte-for-byte
under 0.15.0, even with the same seed.** The ground truth is still exact and
still deterministic; the specific reads differ. Three changes are responsible:

- Contamination generation moved off the global `random.seed()` onto local
  `random.Random(seed)` instances, which draw a different stream.
- Duplicate template selection moved to `numpy.random.Generator`. The weighted
  sampling distribution is unchanged, verified against the previous sequential
  draw over 40,000 trials, but the specific templates chosen for a given seed
  differ.
- Three Gaussian weights in the MDA and PCR duplicate models used the literal
  `2.718` in place of `math.exp()`, which shifts values in the low decimals.

If you have published results tied to a generated dataset, pin the version that
produced it rather than regenerating.

### Fixed

- **VLP generation crashed.** Any run combining `--vlp-protocol` with
  contamination raised `NameError: collection_defaults`. This affected every VLP
  protocol and was present throughout 0.14.0.
- **Filtration size bias ran backwards.** The sigmoid, step and linear retention
  curves gave high recovery to particles *larger* than the pore size. ViroForge
  models dead-end filtration, where the VLP fraction is collected in the
  filtrate, so recovery must fall as virion diameter rises. Losing giant viruses
  and jumbo phages to a 0.22 um membrane is the documented cost of VLP prep. The
  measured size/enrichment correlation moves from roughly +0.97 to -0.86, so
  size-stratified results from earlier versions are inverted.
- Synthetic contaminant sequences (reagent bacteria, and PhiX when no reference
  resolves) are seeded again. `add_phix_control` gained a `random_seed`
  parameter.
- `--seed` had no effect on `scripts/viroforge_subset.py`, which used SQLite's
  unseeded `ORDER BY RANDOM()`. Note that `random_subset` now returns rows in
  database order when the requested count meets or exceeds the number of
  matching rows.
- `PARSERS["auto"]` held `None`, so any caller indexing the registry with
  `auto` outside the CLI raised `TypeError`. Format choices moved to
  `FORMAT_CHOICES`.

### Added

- **Collection-specific contamination baselines.** Each collection carries its
  own contamination profile (blood 40% host cfDNA, marine effectively host-free,
  CF 25%, ocular reagent-heavy). `--contamination-level` is now a multiplier on
  that baseline (clean 0.25x, realistic 1.0x, heavy 2.5x) rather than a fixed
  global fraction. Collections without defaults fall back to the previous preset
  behaviour. Baselines and their citations are in
  `data/reference_profiles/contamination_defaults.tsv`.
- **Classifier format auto-detection** for `viroforge benchmark taxonomy`.
  `--format` defaults to `auto` and recognises Kraken2, Kaiju, Centrifuge,
  DIAMOND and MMseqs2. Unrecognised files are rejected with instructions rather
  than guessed at; use `--format generic` with `--read-id-column` and
  `--taxid-column`.
- **Composition review tooling**: `scripts/evaluate_composition.py` scores each
  collection against per-site literature bands in
  `data/reference_profiles/virome_composition.yaml`, with a strictness dial.

### Changed

- **Taxonomy coverage improved to 71.7% family and 88.9% class** (from 57.1%)
  after enriching ranks from NCBI. No families were invented: ICTV assigns no
  family to most Caudoviricetes genera, so those genomes correctly remain
  `family='Unknown'`.
- **`genome_type` was wrong for 3,516 genomes**, which defaulted to `dsDNA`. The
  entire arbovirus collection was labelled dsDNA. Values are now derived from
  the ICTV Baltimore class. This field is used for display and statistics, not
  for molecule-type logic, which reads `--molecule-type`.
- Four collections were recurated against the literature: blood, lung and
  nasopharynx are now Anelloviridae-dominant, and wastewater is phage-majority.

### Performance

- Duplicate template selection went from O(n*k) to O(n) by replacing a
  `list.index()` scan per draw with a single vectorised call. The effect grows
  with read count and is large at the hundreds-of-thousands scale.
- `pandas` is imported lazily in `viroforge/core/contamination.py`, since only
  `get_contamination_table()` uses it. Worth roughly 90-110 ms per CLI
  invocation on an M-series Mac. That is the marginal cost of pandas on top of
  the numpy and Biopython imports the module already needs, not the ~120 ms
  pandas costs on its own.

### Known issues

- `LEVEL_MULTIPLIERS` defines a `failed` contamination level at 4.0x and the
  docstrings describe it, but both CLI parsers accept only `clean`, `realistic`
  and `heavy`, so it cannot be selected.
- The pore-less and column branches of VLP enrichment still favour larger
  particles, the opposite of the filtration curves. Pelleting and column capture
  are not sieving, so this may be correct, but it has not been checked against
  the literature.
