# ViroForge - Development Context

**Last Updated**: 2026-07-29
**Current Version**: v0.16.0
**Status**: Main green at 371 tests; PR backlog cleared; CHANGELOG started at 0.15.0

---

## Session Handoff (2026-07-28 / 07-29)

origin/main at `ef877c7`, working tree clean, CI green. Earlier handoffs
(2026-07-14, 07-16, 07-17) are in git history and in the memory files under
`~/.claude/projects/-Users-shandley-Code-software-viroforge/memory/`.

### Main was broken, and is not any more

`677dab7` shipped a `NameError: collection_defaults` in `_apply_vlp_enrichment`,
so every run combining `--vlp-protocol` with contamination crashed. Eleven
integration tests were failing. PR #66 fixed it. The suite went 341 -> 371
passing over the two days.

### PR backlog cleared

Ten of Leran10's PRs merged: #61, #62, #63, #64, #66, #67, #68, #69, #70, #71.
Two were closed: #39 (superseded, see below) and #65 (bad citation, see below).
Several needed fixing forward on main rather than in the PR:

- **#67** was a reproducibility regression, not just an incomplete fix. It added
  an `rng` parameter to `_generate_sequence_with_gc` that no call site passed,
  so synthetic sequences became unseeded where the old global `random.seed()`
  had made them deterministic. It fires on every default run, since
  `add_reagent_contamination` is always called with no database path.
- **#71** reintroduced the `collection_defaults` NameError in the extracted
  `generator.py`, and its `_params_to_namespace` restated all 41 parser defaults,
  having already drifted to 6 missing attributes and 11 wrong values.
  `viroforge generate` raised AttributeError. The parser now lives in
  `viroforge/generator.py` as `build_parser()` and the CLI takes defaults from
  `parse_args([])`, so drift is impossible. The script is down to 100 lines.
- **#61** inverted the VLP filtration curves so small virions are recovered,
  which is correct for dead-end filtration. Nothing covered `FiltrationCurve`,
  so the direction had never been tested.

### Two modelling corrections

- **VLP size bias runs the other way now.** Recovery falls with virion diameter;
  losing giant viruses and jumbo phages to a 0.22 um membrane is the documented
  cost of VLP prep. The size/enrichment correlation moved from about +0.97 to
  -0.86, so size-stratified results from before 0.15.0 are inverted.
- **MDA GC bias peaked at the wrong GC content.** The curve was centred at 40%
  GC and penalised 45-60% GC genomes up to 6x, where Parras-Molto et al. 2018
  (PMID 29954453) measured MDA *over*-amplifying that band. Optimum moved to
  50%. `scripts/benchmark_amplification_bias.py` calibrates this against the
  paper's measured 6.2-7.6% of contigs past 10x fold change.

### Releases

0.15.0 and 0.16.0, both output-changing. `CHANGELOG.md` now exists and starts at
0.15.0. For a generator whose value is reproducible ground truth, any change to
output for a fixed seed is treated as breaking and recorded there. No git tags
and nothing published to PyPI, so versions are repo state only.

### Retired

The lab-notebook workflow is gone: `.git/hooks/pre-commit` archived to
`pre-commit.disabled`, and the four unregistered `.claude/hooks/` scripts
deleted. The global `~/.githooks/pre-commit` (gitleaks) still runs and still
chains to any repo-local hook.

### Watch out for

Leran10's citations have not been reliable. Every PMID she has cited (2 of 2)
resolved to an unrelated paper, and in #65 the number attributed to the paper
was invented as well. `/verify-references` catches the first problem, not the
second. Check that the cited paper actually contains the number being used.

### Still open

- **PR #51** host_associations, deferred.
- **Issue #37** bacterial/fungal background. PR #39 was closed rather than
  rebased: its Part 1 was superseded by `677dab7`, and the branch carried 4.4M
  lines of FASTA in git, unrelated `duplicates.py` hunks that now conflict, and
  at least one wrong accession. Build the FASTA on demand and verify accessions.
- **Issue #30** web visual form builder, not started.
- MDA and RdAB now produce nearly identical GC curves at their defaults. The
  docstring claim that MDA bias is 2-5x stronger than PCR no longer holds and
  needs its own measurement.
- `--contamination-level failed` is defined in `LEVEL_MULTIPLIERS` and
  documented, but both CLI parsers accept only clean/realistic/heavy, so it
  cannot be selected.

### Environment

A full `.venv` (CPython 3.12) with `viroforge[web]` + InSilicoSeq 2.0.1 + mappy.
Run tests with `.venv/bin/python -m pytest`, and `export PATH="$PWD/.venv/bin:$PATH"`
so `iss` resolves. The full suite takes about 5 minutes. PBSIM3 is not built and
`pbccs` is Linux-x86-64 only, so long-read generation cannot run on this machine.

---

## Project Overview

ViroForge is a comprehensive mock metavirome data generator for benchmarking virome analysis pipelines. It generates synthetic FASTQ datasets with complete ground truth for validation.

**Core Capabilities**:
- 20 curated virome collections (host-associated, environmental, disease states)
- 14,423 RefSeq viral genomes with ICTV taxonomy (71.7% family, 88.9% class coverage after 2026-07-16 NCBI rank enrichment)
- 5 sequencing platforms (NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore)
- DNA and RNA virome workflows (RT, rRNA depletion, degradation)
- VLP enrichment modeling (5 protocols)
- Real reference contamination (rRNA from NCBI, host DNA from T2T, PhiX NC_001422.1, Illumina adapters)
- Per-read source labels (source=viral/host_dna/rrna/phix/reagent/erv_endogenous/erv_exogenous)
- Sequencing artifact injection (adapter read-through, low-complexity, PCR duplicates, ERVs)
- Complete ground truth metadata for validation (v1.1)

---

## Current Status

Phases 1 through 13 (metadata, QC, assembly and taxonomy benchmarking) are
delivered. What follows is the feature history; see the handoff above for where
things stand today and `CHANGELOG.md` for what changed in each release.

**Phase 13A: Benchmarking Metadata - COMPLETE**

Metadata enhancements (2025-11-11):
- Metadata schema v1.1 with benchmarking support
- Contamination manifest export (all contaminant sources tracked)
- Expected coverage calculation per genome (Lander-Waterman formula)
- Coverage categorization (complete/high/partial/fragmented/missing)
- HTML report format implemented
- HTML comparison format implemented

**QC Validation Toolkit - COMPLETE (2026-03-22)**

Real reference contamination, sequencing artifact injection, and per-read ground truth labeling for QC pipeline validation.

Contamination realism:
- rRNA: 23 real sequences from NCBI RefSeq (E. coli, human, gut bacteria 16S/18S/23S/28S)
- Host DNA: 48 fragments from T2T-CHM13v2.0 human genome (all chromosomes, 10 kb each)
- PhiX174: Real NC_001422.1 genome (5,386 bp)
- Adapters: TruSeq and Nextera Illumina adapter sequences (11 sequences)
- PhiX curation fix: GCF_000819615.1 removed from 8 collections (lab phage, not natural)

Sequencing artifact modules:
- Adapter read-through (`--adapter-rate`, `--adapter-type`): Normal insert size distribution, header tags, manifest TSV
- Low-complexity artifacts (`--low-complexity-rate`): Homopolymers, dinucleotide repeats, simple repeats, low-entropy
- Controlled entropy (`--entropy-range 0.3-0.7`): Gray-zone reads for complexity filter threshold testing
- PCR duplicates (`--duplicate-rate`): Geometric copy distribution, PCR error modeling, paired-end consistency
- Retroviral reads (`--erv-endogenous-rate`, `--erv-exogenous-rate`): Endogenous HERV + exogenous Retroviridae

Per-read ground truth:
- Source labels on all FASTQ headers (`source=viral/host_dna/rrna/phix/reagent/erv_*`)
- Adapter, low-complexity, and duplicate manifest TSV files
- All artifact stats saved to metadata JSON

Key files:
- `viroforge/data/references/` - Bundled reference FASTA files (~528 KB total)
- `viroforge/data/references/resolver.py` - Auto-discovers references (user > env var > bundled > synthetic)
- `viroforge/simulators/adapters.py` - Adapter read-through post-processor
- `viroforge/simulators/low_complexity.py` - Low-complexity artifact injector with entropy control
- `viroforge/simulators/duplicates.py` - PCR duplicate injector
- `viroforge/core/contamination.py` - ERV contamination types (endogenous + exogenous)
- `scripts/curate_reference_sequences.py` - Provenance script for regenerating references

**Phase 12: CLI & Web Enhancements - COMPLETE**

All subphases delivered:

**Phase 12.1**: Full generate command with progress reporting
- Preset-based generation system (8 built-in presets)
- Real-time progress bars via rich.progress
- Parameter override system
- Verbose mode for debugging

**Phase 12.2**: Batch generation and result reporting
- YAML-based batch configuration
- Parameter sweep support (cartesian product expansion)
- Sequential and parallel execution modes
- Dataset quality reporting (viroforge report)
- Intelligent dataset comparison (viroforge compare)
- 5 example batch configurations

**Phase 12.3**: Web interface
- Flask-based web application (viroforge web)
- Bootstrap 5 responsive UI with 9 pages
- RESTful API (10+ endpoints)
- Visual collection browser with search/filter
- Interactive dataset generation with progress monitoring
- Batch configuration builder
- Dataset reporting and comparison dashboards

**Commands Available**:
```bash
viroforge browse              # Interactive collection browser (TUI)
viroforge generate            # Generate datasets (CLI)
viroforge report              # View dataset reports
viroforge compare             # Compare multiple datasets
viroforge batch               # Batch generation from YAML
viroforge presets             # Manage configuration presets
viroforge web                 # Launch web interface
```

---

## Database

**Location**: `viroforge/data/viral_genomes.db` (SQLite, NOT in git)
**Size**: ~500 MB

**Key Tables**:
- `genomes` - 14,423 RefSeq viral genome sequences
- `taxonomy` - ICTV/NCBI taxonomy mappings (71.7% family, 88.9% class matched)
- `body_site_collections` - 20 collection metadata records
- `collection_genomes` - Collection-genome associations with abundances

**Regeneration**: Database can be regenerated from RefSeq if needed (scripts in `scripts/` directory)

---

## Key File Locations

### Core Generation
- `scripts/generate_fastq_dataset.py` - Main dataset generation script (all platforms)
- `viroforge/simulators/longread.py` - PacBio HiFi and Nanopore simulator (850+ lines)
- `viroforge/simulators/illumina.py` - Short-read simulator wrapper
- `viroforge/simulators/adapters.py` - Adapter read-through post-processor
- `viroforge/simulators/low_complexity.py` - Low-complexity artifact injector with entropy control
- `viroforge/simulators/duplicates.py` - PCR duplicate injector
- `viroforge/workflows/rna_virome.py` - RNA workflow (RT, rRNA depletion, degradation)
- `viroforge/vlp_enrichment.py` - VLP protocol modeling

### Contamination References
- `viroforge/data/references/resolver.py` - Reference file auto-discovery
- `viroforge/data/references/phix174.fasta` - Real PhiX174 genome (NC_001422.1)
- `viroforge/data/references/rrna_representatives.fasta` - 23 real rRNA sequences
- `viroforge/data/references/host_fragments.fasta` - 48 human genome fragments (T2T-CHM13v2.0)
- `viroforge/data/references/adapters.fasta` - Illumina TruSeq/Nextera adapters
- `scripts/curate_reference_sequences.py` - Downloads and curates reference sequences from NCBI

### CLI Application
- `viroforge/cli/__init__.py` - Main CLI entry point and argument parsing
- `viroforge/cli/browse.py` - Interactive TUI collection browser
- `viroforge/cli/generate.py` - Dataset generation command (340 lines)
- `viroforge/cli/batch.py` - Batch generation (358 lines)
- `viroforge/cli/report.py` - Dataset reporting (319 lines)
- `viroforge/cli/compare.py` - Dataset comparison (261 lines)
- `viroforge/cli/presets.py` - Preset management
- `viroforge/cli/web.py` - Web server launcher (80 lines)

### Web Interface
- `viroforge/web/app.py` - Flask application with API routes (280 lines)
- `viroforge/web/templates/` - 9 HTML templates with Bootstrap 5 (~1,250 lines)

### Database Utilities
- `viroforge/database/utils.py` - Database query utilities for CLI/web
- `viroforge/database/manager.py` - Collection management

### Collection Curation
- `scripts/curate_*_virome_collection.py` - 20 collection curation scripts
- Follow proactive taxonomy rescanning protocol (see below)

### Documentation
- `docs/PHASE12.1_SUMMARY.md` - Generate command documentation
- `docs/PHASE12.2_SUMMARY.md` - Batch/report/compare documentation
- `docs/PHASE12.3_SUMMARY.md` - Web interface documentation
- `docs/LONGREAD_TUTORIAL.md` - PacBio HiFi and Nanopore guide
- `docs/HYBRID_ASSEMBLY_TUTORIAL.md` - Hybrid assembly workflows
- `docs/COLLECTION_IMPLEMENTATION_GUIDE.md` - All 20 collections documented
- `docs/TAXONOMY_BUG_FIX.md` - Critical taxonomy fix documentation

### Example Configurations
- `examples/presets/` - 8 built-in configuration presets
- `examples/batch_configs/` - 5 batch YAML examples

---

## Collection Curation Protocol

**When adding new collections, always follow this protocol to avoid taxonomy issues:**

### 1. Research Phase
- Literature review for virome composition
- Identify dominant viruses and prevalence
- Document clinical/research applications

### 2. Script Creation
Create `scripts/curate_[collection]_virome_collection.py` with:
- Separate methods for each viral category
- SQL queries targeting specific viruses
- Abundance assignment based on literature

### 3. Proactive Taxonomy Rescanning (CRITICAL)
**Always verify expected viruses are present after initial run**

**Common Query Issues**:
- Non-human matches: Use specific patterns like `"Human bocavirus%"` not `"%bocavirus%"`
- Substring matches: Avoid broad patterns (e.g., `"%HIV%"` matches "Phives")
- Family name changes: ICTV taxonomy evolves (e.g., `Orthoherpesviridae` not `Herpesviridae`)
- Missing families: Some viruses have `family="Unknown"` (e.g., HCV)
- Word boundaries: Use `"Influenza A virus%"` to exclude parainfluenza

**Query Fix Patterns**:
```sql
-- Too broad
WHERE g.genome_name LIKE '%bocavirus%'
-- Better (human-specific)
WHERE g.genome_name LIKE 'Human bocavirus%'

-- Substring issues
WHERE g.genome_name LIKE '%HIV%'
-- Better (word boundary)
WHERE g.genome_name LIKE 'Human immunodeficiency virus%'

-- Explicit exclusions when needed
WHERE g.genome_name LIKE 'Respiratory syncytial virus%'
  AND g.genome_name NOT LIKE '%Bovine%'
  AND g.genome_name NOT LIKE '%Potato%'
```

### 4. Verification
- Check all CRITICAL biomarker viruses are present
- Verify genome counts match targets
- Confirm database insertion

### 5. Documentation
- Update README.md with collection details
- Include applications, literature basis, special notes

---

## Testing Workflow

**Unit Tests**: `pytest tests/ -v`
**Integration Tests**: Require full environment (InSilicoSeq, PBSIM3, etc.)

**Test Coverage**:
- 80+ tests passing
- RNA workflow tests (40+ tests)
- Contamination tests (30+ tests)
- Long-read simulator tests
- VLP enrichment tests

---

## Git Workflow

**Database**: `viral_genomes.db` is in `.gitignore` (too large)
**Commits**: Use detailed commit messages with documentation
**Branches**: Main branch for stable releases

**Recent Commits**:
- d57e58c: Phase 12.3 - Web Interface
- 8f7169e: Phase 12.2 - Batch Generation and Result Reporting
- cf1fd64: Phase 12.1 - Complete Core Features

---

## Next Steps / Future Work

### Immediate (No Additional Phases Required)
- User testing and feedback collection
- Bug fixes and polish based on usage
- Performance optimization if needed
- Documentation improvements

### Optional Future Phases (Not Planned)

**Phase 13: Advanced Visualizations** (if requested)
- Interactive composition charts (plotext or matplotlib)
- Assembly quality dashboards
- Benchmark result visualization
- Coverage depth heatmaps

**Phase 14: Animal Model Collections** (if requested)
- Zebrafish virome
- Pig virome (agricultural/swine flu)
- Chicken virome (poultry health)
- Primate virome (research models)

**Phase 15: Environmental Diversity** (if requested)
- Hot spring virome (extremophiles)
- Hypersaline environment
- Hospital environment (nosocomial)
- Plant virome

**Web Interface Enhancements** (if needed for production)
- User authentication (Flask-Login)
- Job queue system (Redis/Celery)
- Real-time progress streaming (WebSockets)
- Result file downloads via UI
- Dataset browser (browse generated datasets)
- Preset editor (create/edit in UI)

---

## Important Notes

### Security Considerations
**Web interface is for LOCAL USE or TRUSTED NETWORKS only**:
- No authentication implemented
- No HTTPS
- In-memory job tracking (not persistent)
- Localhost binding by default (127.0.0.1)

For production deployment: Add authentication, HTTPS, job queue, rate limiting.

### Taxonomy Limitations
- 88.9% of genomes are classified to class, 71.7% to family (after the 2026-07-16
  NCBI rank enrichment; see docs/BIOLOGICAL_ACCURACY_REVIEW.md)
- 28.3% (4,080) have `family="Unknown"` - largely modern bacteriophages that ICTV
  assigns to a genus but no family (families abolished 2021), so this is correct
  taxonomy, not a gap. 1,481 genomes are genuinely unclassified (no class).
- Enhanced fuzzy matching fixed 469 genomes; NCBI enrichment then classified 4,577
  more to class and 2,102 to family
- See `docs/TAXONOMY_BUG_FIX.md` for details

### Platform Requirements
**Short-read (Illumina)**: InSilicoSeq required
**Long-read (PacBio HiFi)**: PBSIM3 + pbccs (samtools) required
**Long-read (Nanopore)**: PBSIM3 required
**RNA workflows**: All dependencies above
**Web interface**: Flask >= 2.0.0

---

## References

**Installation**: `pip install -e .` (core) or `pip install -e ".[web]"` (with Flask)
**Repository**: https://github.com/shandley/viroforge
**Principal Investigator**: Scott Handley, Washington University in St. Louis
**Lab Website**: https://www.handleylab.org
**License**: MIT

---

## Quick Command Reference

```bash
# CLI Commands
viroforge browse                              # Browse collections (TUI)
viroforge generate --preset gut-standard      # Generate with preset
viroforge generate --collection-id 1 --output data/gut --platform novaseq --coverage 30
viroforge batch config.yaml --parallel 4      # Batch generation
viroforge report data/gut-standard            # View report
viroforge compare data/gut_* --format json    # Compare datasets
viroforge presets list                        # List presets
viroforge web                                 # Launch web UI

# Generation Script (Direct)
python scripts/generate_fastq_dataset.py \
    --collection-id 1 \
    --output data/gut \
    --platform novaseq \
    --coverage 30 \
    --vlp-protocol tangential_flow

# Hybrid Assembly
python scripts/generate_hybrid_dataset.py \
    --collection-id 1 \
    --output data/gut_hybrid \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 \
    --depth 15 \
    --seed 42

# Database Query
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT collection_id, collection_name, n_genomes
   FROM body_site_collections
   ORDER BY collection_id;"
```

---

## Known Issues

**Interactive Preset Creation**
- Status: Not implemented (low priority)
- Workaround: Use `viroforge presets create <name> --from-dataset <path>`

**ISS Read Count Deviation for Small Genomes**
- InSilicoSeq `--mode basic` can produce fewer reads than expected for very small
  genomes due to rounding and minimum-read effects. Example: PhiX174 at 0.06%
  abundance with 537K total reads should produce ~322 reads but ISS generated 276.
- Impact: Small (< 15% deviation), only affects genomes with very low abundance.
- Workaround: Use actual output read counts rather than abundance-predicted counts
  when computing benchmarking metrics.

**PhiX174 Curation Fix (2026-03-22)**
- GCF_000819615.1 (Escherichia phage phiX174) was removed from 8 collections
  (9, 10, 12, 13, 14, 16, 17, 18) because it is the same lab phage used as
  Illumina spike-in control. Including it as a natural virus causes QC tools
  that remove PhiX to also remove "natural" community reads.

**Vaginal Virome: Proxy Phages from Dairy/Food Lactobacillus (2026-06-02)**
- The vaginal virome collection (Collection 16) uses Lactobacillus phages from
  RefSeq as proxies for vaginal Lactobacillus phages. However, nearly all 60
  Lactobacillus phages in RefSeq infect dairy/food species (L. plantarum,
  L. delbrueckii, L. casei), not vaginal species (L. crispatus, L. iners,
  L. gasseri, L. jensenii).
- Only ~1 phage (Lv-1, infecting L. jensenii) is from a vaginal isolate.
- No Gardnerella, Prevotella, Atopobium, or Sneathia phages exist in RefSeq,
  preventing accurate modeling of CST IV (bacterial vaginosis) communities.
- Impact: The vaginal virome is biologically plausible at the family level
  (Lactobacillus phages dominate, HPV present) but not species-specific.
  Benchmarking results for vaginal datasets should note this limitation.
- CST (Community State Type) modeling is deferred until vaginal-specific
  phage genomes become available in public databases (RefSeq, INPHARED, GPD).
- Reference: Lv-1 vaginal phage (Lactobacillus jensenii) described in
  Miller-Ensminger et al. 2020 (PLOS One, PMC7289420).

**Dark Matter Sequences: Random Selection, Not Body-Site-Specific (2026-06-05)**
- The `--dark-matter-fraction` flag adds unclassified viral genomes (family='Unknown')
  to simulate the large fraction of reads in real viromes that do not match known
  references. As of v0.13.0 it is on by default (0.30).
- Dark matter genomes are selected randomly from the Unknown pool, NOT filtered by
  body site. A gut virome and a vaginal virome draw dark matter from the same pool.
- In reality dark matter is body-site-specific, but most Unknown-family genomes lack
  host/habitat metadata, making body-site filtering impractical.
- Exclusions applied: animal viruses, insect viruses, and known human viruses are
  excluded from the pool. Plant viruses are intentionally kept (dietary plant viruses
  are a real component of human gut viromes).
- Impact on QC benchmarking: none (QC tools process reads the same regardless of dark
  matter origin). Body-site specificity only matters for taxonomy benchmarking.

---

## Recent Bug Fixes (2025-11-16)

**HTML Report/Comparison Empty (FIXED)**
- Issue: HTML reports generated mostly empty content
- Root cause: Metadata field name mismatch between v1.1 schema and HTML templates
  - Templates accessed `metadata['platform']` - not present in v1.1
  - Templates accessed `metadata['vlp_enrichment']` - renamed to `enrichment_stats` in v1.1
- Fix: Updated both `report.py` and `compare.py` to use correct v1.1 field names
  - Platform info: `configuration['platform']`
  - VLP/enrichment info: `enrichment_stats['viral_enrichment']`, etc.
- Status: RESOLVED
- Files: `viroforge/cli/report.py:355-432`, `viroforge/cli/compare.py:294-337`

---
## Benchmarking Framework (Phase 13)

Modules 1, 2 and 4 are built and shipped in `viroforge/benchmarking/`, each with
a `viroforge benchmark <module>` subcommand, JSON and markdown reports, and an
independent-oracle test suite.

- **Module 1 QC** (`benchmark qc`): contamination removal and viral retention
  scored against per-read `source=` labels, with a read-name match-rate gate.
  Duplicates are scored separately. See `docs/PHASE13B_QC_BENCHMARK.md`.
- **Module 2 Assembly** (`benchmark assembly`): genome recovery, completeness
  and identity, chimeras, N50/L50, observed-vs-expected completeness, abundance
  accuracy. minimap2 via mappy, shared in `benchmarking/align.py` (the
  `benchmark` extra). See `docs/PHASE13B_ASSEMBLY_BENCHMARK.md`.
- **Module 4 Taxonomy** (`benchmark taxonomy`): read- and contig-based,
  taxid-exact plus genus/family via an NCBI taxdump, abundance profile,
  known/dark-matter stratification. Kraken2, Kaiju, Centrifuge, DIAMOND, MMseqs2
  and generic formats, auto-detected by default. See
  `docs/PHASE13C_TAXONOMY_BENCHMARK.md`.

Ground truth comes from metadata `benchmarking.{expected_coverage,taxonomy,
contamination_manifest}`, the per-read `source=` labels, and the source genome
FASTA.

**Not built**: Module 5 (completeness across coverage), HTML reports and
visualizations, length-weighted contig taxonomy metrics, and Modules 3 and 6-9
(binning, annotation, host prediction, novel discovery, end-to-end). See
`ROADMAP.md` and `docs/PHASE13_BENCHMARKING_FRAMEWORK.md`.

**Why it matters**: the field has no CAMI-equivalent benchmark for viromes.
Completing generate -> benchmark -> publish is what ViroForge is for.
