# Changelog

Notable changes to ViroForge. This file starts at 0.15.0; for earlier history,
see the git log.

Versions follow [semantic versioning](https://semver.org/), with the caveat that
ViroForge generates data. A change that alters generator output for a fixed seed
is treated as breaking even when the API is untouched, and is called out under
"Reproducibility" below.

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
