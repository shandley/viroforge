# ViroForge - Development Context

**Last Updated**: 2025-11-10
**Current Version**: v0.10.0
**Status**: Phase 12 Complete (Production Ready)

---

## Project Overview

ViroForge is a comprehensive mock metavirome data generator for benchmarking virome analysis pipelines. It generates synthetic FASTQ datasets with complete ground truth for validation.

**Core Capabilities**:
- 28 curated virome collections (host-associated, environmental, disease states)
- 14,423 RefSeq viral genomes with ICTV taxonomy (57.1% coverage)
- 5 sequencing platforms (NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore)
- DNA and RNA virome workflows (RT, rRNA depletion, degradation)
- VLP enrichment modeling (5 protocols)
- Complete ground truth metadata for validation

---

## Current Status

**Phase 12: CLI & Web Enhancements - COMPLETE**

All three subphases delivered (3 days total):

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
- `taxonomy` - ICTV taxonomy mappings (57.1% matched)
- `body_site_collections` - 28 collection metadata records
- `collection_genomes` - Collection-genome associations with abundances

**Regeneration**: Database can be regenerated from RefSeq if needed (scripts in `scripts/` directory)

---

## Key File Locations

### Core Generation
- `scripts/generate_fastq_dataset.py` - Main dataset generation script (all platforms)
- `viroforge/simulators/longread.py` - PacBio HiFi and Nanopore simulator (850+ lines)
- `viroforge/simulators/illumina.py` - Short-read simulator wrapper
- `viroforge/workflows/rna_virome.py` - RNA workflow (RT, rRNA depletion, degradation)
- `viroforge/vlp_enrichment.py` - VLP protocol modeling

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
- `scripts/curate_*_virome_collection.py` - 28 collection curation scripts
- Follow proactive taxonomy rescanning protocol (see below)

### Documentation
- `docs/PHASE12.1_SUMMARY.md` - Generate command documentation
- `docs/PHASE12.2_SUMMARY.md` - Batch/report/compare documentation
- `docs/PHASE12.3_SUMMARY.md` - Web interface documentation
- `docs/LONGREAD_TUTORIAL.md` - PacBio HiFi and Nanopore guide
- `docs/HYBRID_ASSEMBLY_TUTORIAL.md` - Hybrid assembly workflows
- `docs/COLLECTION_IMPLEMENTATION_GUIDE.md` - All 28 collections documented
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
- 57.1% of genomes have ICTV taxonomy matches
- 42.9% have `family="Unknown"` (strain-specific nomenclature issues)
- Enhanced fuzzy matching fixed 469 genomes (7.1%)
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
viroforge generate --collection-id 9 --output data/gut --platform novaseq --coverage 30
viroforge batch config.yaml --parallel 4      # Batch generation
viroforge report data/gut-standard            # View report
viroforge compare data/gut_* --format json    # Compare datasets
viroforge presets list                        # List presets
viroforge web                                 # Launch web UI

# Generation Script (Direct)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut \
    --platform novaseq \
    --coverage 30 \
    --vlp-protocol tangential_flow

# Hybrid Assembly
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
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
