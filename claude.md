# ViroForge - Development Context

**Last Updated**: 2026-07-16
**Current Version**: v0.14.0 (benchmarking framework)
**Status**: Canonical DB adopted; data quality evaluated; benchmarking framework Modules 1/2/4 built

---

## Session Handoff (2026-07-16)

origin/main is at `6edef85`. Working tree clean except the pre-existing untracked
`data/web/` + `scripts/export_web_data.py` (a prior static web export) and the
gitignored `validation/` scratch dir.

### Environment

A full `.venv` (CPython 3.12) now exists with `viroforge[web]` + InSilicoSeq 2.0.1
+ mappy (the `benchmark` extra). Real generation, the test suite, and the
benchmarks all run in it. PBSIM3 is not built and `pbccs` is Linux-x86-64 only, so
long-read generation is not exercisable on this macOS-arm64 machine.

### Canonical database

Adopted a fresh seeded `setup-db` rebuild as the canonical DB (the prior DB
predated the seeding fix and was not regenerable). It is reproducible via
`setup-db`; provenance is pinned by `data/collection_membership.tsv` +
`scripts/export_collection_manifest.py` (the 500 MB DB is gitignored). Also fixed:
vaginal-curation idempotency and an `n_genomes` reconcile step.

### Data quality evaluation

Ran a full pre-QC technical evaluation (`docs/DATA_QUALITY_EVALUATION.md`,
`scripts/evaluate_dataset.py`). Verdict: per-read `source=` labels and artifact
tags are accurate, traceable, and reproducible. Fixes from it: the default 0.30
dark-matter fraction is now actually delivered (was 0.11-0.48; amplification
reweighting); the three Illumina platforms are documented as interchangeable in
ISS basic mode (identical 125 bp reads); MDA duplicate metadata records real values.

### Benchmarking framework (Phase 13B/C) - `viroforge/benchmarking/`

Three modules built, each with `viroforge benchmark <module>`, JSON + markdown
reports, and an independent-oracle test suite (26 tests):
- **Module 1 QC** (`benchmark qc`): contamination removal + viral retention vs
  per-read labels; match-rate gate; dedup scored separately.
- **Module 2 Assembly** (`benchmark assembly`): genome recovery/completeness,
  chimeras, N50/L50, observed-vs-expected completeness, abundance accuracy. Uses
  minimap2 via mappy (shared `benchmarking/align.py`).
- **Module 4 Taxonomy** (`benchmark taxonomy`): read- and contig-based;
  taxid-exact + genus/family (via NCBI taxdump); abundance profile;
  known/dark-matter stratification; Kraken2/Centrifuge/DIAMOND/MMseqs2/generic
  parsers. Metadata now exports `benchmarking.taxonomy` (ncbi_taxid + lineage).

Deferred: Module 5 (completeness across coverage), HTML reports + visualizations,
length-weighted contig taxonomy metrics.

### Still open (unchanged)

- PR #39 (collection-specific contamination, blocked - bad accessions), PR #51
  (host_associations, deferred). Issues #37, #38 (tied to #39), #30 (web form builder).

---

## Session Handoff (2026-07-14)

Reviewed Leran10's 30 issues and 27 PRs (she is a developer on the team, using
Claude Code). Full per-PR triage and the verified facts are in the memory file
`session_handoff.md`.

### Landed on main (pushed, origin up to 3482003)

- Three local fixes committed: #34 (RNA metadata rebuild), #35 (hybrid version
  and SPAdes flag), #36 (validator sequences alias).
- Code from PR #29 extracted (commit 13eced4): PacBio HiFi BAM merge (#17, #33),
  Nanopore ERRHMM names (#15), per-genome depth (#19, #31), setup-db command
  (#2, #3), browse checks (#7, #13), batch time estimate (#26), report/compare
  metadata (#24), web fixes (#28). Applied on top of main so #34/#35/#36 are
  preserved; the TSV accession swaps, the curate-script renumbering, and two
  regressed scripts were deliberately excluded.
- PR #23 pbccs Linux x86-64 note (commit c525e33, #22).
- Merged scientific PRs: #48 (vaginal docs), #47 (VLP pore size), #58 (poly-G
  low-complexity), #53 (rare-genome floor).
- Closed 14 PRs and 19 issues. Rework feedback posted on PRs #41, #50, #51, #55,
  #56.

### Second session (2026-07-14): PR backlog resolved

Scott decided to adopt the 1-20 renumbering and move to realistic default
output with a version bump. Worked through the open PRs one at a time:

- **Merged/landed** (credited to Leran10): #56 (Okabe-Ito palette, `fdd5682`),
  #6 (renumber to 1-20 + migration, `1e3af14`), #41 (remove animal/plant
  viruses, `9eae691`), #50 (replace non-site phages, `68439a5`), and the
  v0.13.0 realistic-defaults cluster `74e33b4` (#43 dark matter, #45 artifact
  defaults, #55 host filter folded in).
- **Deferred with feedback (still open)**: #51 (host_associations - consumer-less
  infra, needs a real body_site column) and #39 (collection-specific
  contamination - commits ~4.4M lines of FASTA into git and hardcodes unverified
  bacterial accessions; needs the FASTA out of git and /verify-references).
- **Issues**: closed #5, #40, #42, #44, #49, #54 (resolved by the merges).
  Still open: #37, #38 (both tied to blocked #39), #30 (web form builder).

### Verified facts

- The live database has **20 collections at contiguous IDs 1-20** (renumbered
  from 9-28). There are no separate "VLP comparison" collections; VLP is applied
  via `--vlp-protocol`. The animal-virus and phage-host cleanups run as a
  post-curation step in `viroforge setup-db`.
- Test interpreter is `.venv_test/bin/python` (the older `.venv` path is gone).
  `.venv_test` lacks rich and flask, so CLI and web runtime imports fail there.
  It has numpy + viroforge, enough for dry-run generation.

### Open PRs

39, 51 (both intentionally deferred with feedback).

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

**Phase 13A: Benchmarking Metadata - COMPLETE**

Metadata enhancements (2025-11-11):
- Metadata schema v1.1 with benchmarking support
- Contamination manifest export (all contaminant sources tracked)
- Expected coverage calculation per genome (Lander-Waterman formula)
- Coverage categorization (complete/high/partial/fragmented/missing)
- CLI version updated to 0.11.0
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

## Next Phase: Benchmarking Framework (Phase 13)

**Status**: Phase 13A Complete | All CLI Features Working | Ready for Phase 13B
**Timeline**: 6-8 weeks for MVP (2 weeks elapsed)
**Priority**: VERY HIGH

### Phase 13A Complete (v0.11.0)

**Delivered**:
- Metadata version 1.1 with benchmarking support
- Contamination manifest export
- Expected coverage calculation per genome
- CLI version updated to 0.11.0
- HTML report format (implemented and tested)
- HTML comparison format (implemented and tested)

See `docs/PHASE13A_IMPLEMENTATION_SUMMARY.md` for details.

### Overview

Transform ViroForge from a data generator to a **complete validation platform** by adding comprehensive benchmarking tools.

**Phase 13A Achievement**: ViroForge now exports enhanced ground truth metadata that benchmarking tools (Phase 13B-D) will use to validate pipeline performance.

**Solution**: Modular benchmarking framework matching virome analysis workflow:
- Module 1: QC Benchmarking (contamination removal validation)
- Module 2: Assembly Benchmarking (genome recovery, completeness)
- Module 3: Binning Benchmarking (vMAG reconstruction)
- Module 4: Taxonomy Benchmarking (read + contig classification)
- Module 5: Completeness Analysis (genome recovery across coverage)
- Module 6: Annotation Benchmarking (gene calling + function)
- Module 7: Host Prediction Benchmarking
- Module 8: Novel Discovery Benchmarking
- Module 9: End-to-End Pipeline Benchmarking

### Architecture

```
viroforge/
├── benchmarking/              # NEW (optional: pip install viroforge[benchmark])
│   ├── parsers/              # Kraken2, Centrifuge, DIAMOND, etc.
│   ├── metrics/              # Precision, recall, F1, completeness, etc.
│   ├── visualizations/       # Plots for HTML reports
│   └── reports/              # HTML + JSON report generation
```

### CLI Commands (Planned)

```bash
viroforge benchmark qc          # Test contamination removal
viroforge benchmark assembly    # Test assembly quality
viroforge benchmark taxonomy    # Test classification accuracy
viroforge benchmark pipeline    # Comprehensive end-to-end report
```

### ViroForge Enhancements Needed

**Metadata Enhancements** (v0.11.0):
1. Add `contamination_manifest` - Track all contaminant sources
2. Add `expected_coverage` - Per-genome expected coverage/completeness
3. Add `read_manifest` - Which genome each read came from (optional)
4. Add `gene_annotations` - Export RefSeq CDS annotations (future)

### Why This Matters

- **For Pipeline Developers**: Rigorous validation before publication
- **For Bioinformaticians**: Choose best tool for their data type
- **For ViroForge**: Completes "generate → benchmark → publish" workflow
- **For the Field**: Standardized virome pipeline benchmarking (CAMI-equivalent for viromes)

### Implementation Phases

**Phase 13A** (Weeks 1-2): Foundation + metadata enhancements
**Phase 13B** (Weeks 3-5): QC + Assembly benchmarking
**Phase 13C** (Weeks 6-8): Taxonomy + Completeness benchmarking
**Phase 13D** (Future): Advanced modules (binning, annotation, host, discovery)

See `ROADMAP.md` for complete Phase 13 specifications.
