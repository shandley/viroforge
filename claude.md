# ViroForge - Claude Development Context

**Last Updated**: 2025-11-10
**Current Version**: v0.10.0 (Phase 12 complete - all 3 subphases)
**Current Phase**: Phase 12 - CLI & Web Enhancements 🎨🌐 (Complete)

## Project Overview

ViroForge is a comprehensive mock metavirome data generator for benchmarking virome analysis pipelines. It generates synthetic FASTQ datasets with complete ground truth for validation.

**Key Features**:
- 28 curated virome collections (1,559 genomes total)
- DNA + RNA virome workflows
- **5 sequencing platforms**: NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore
- VLP enrichment modeling (adapted for long reads)
- Amplification bias (RdAB, MDA, Linker)
- Platform-specific error models (all 5 platforms)
- Literature-validated compositions

## Current Status

**Phase 9**: ✅ COMPLETE (November 2025)
- Collections 24-28: Vaginal, Blood, Ocular, Lung, Urinary viromes
- All collections use proactive taxonomy rescanning protocol
- Total: 28 collections across diverse environments

**Phase 10**: ✅ COMPLETE (November 2025) - Long-Read Sequencing Support
- ✅ PBSIM3 research and architecture design
- ✅ Core simulator module (`viroforge/simulators/longread.py` - 850+ lines)
- ✅ PacBio HiFi (two-step CCS workflow) and Nanopore support
- ✅ VLP modeling updates for long reads (60% size bias reduction)
- ✅ Integration with generate_fastq_dataset.py (full routing logic)
- ✅ Comprehensive testing (`tests/test_longread_simulator.py`)
- ✅ Complete documentation (`docs/LONGREAD_TUTORIAL.md`)
- **Timeline**: Completed in 3 weeks (originally planned 3-4 weeks)

**Phase 11**: ✅ COMPLETE (November 2025) - Hybrid Assembly Support
- ✅ Convenience script for matched short+long datasets
- ✅ Composition validation utility
- ✅ Complete hybrid assembly tutorial
- ✅ Support for Unicycler, SPAdes hybrid mode, MaSuRCA
- **Timeline**: Completed in 1 day

**Phase 12**: ✅ COMPLETE (November 2025) - CLI Enhancements
- ✅ Unified `viroforge` CLI command
- ✅ Interactive collection browser (`viroforge browse`)
- ✅ Configuration presets system (8 built-in presets)
- ✅ Beautiful terminal UI with `rich` library
- ✅ Database utilities for CLI
- ✅ **Phase 12.1 COMPLETE**: Full generate command with progress reporting
  - ✅ Preset-based generation (`viroforge generate --preset`)
  - ✅ Real-time progress bars and status updates
  - ✅ Parameter override system
  - ✅ Verbose mode for detailed output
- ✅ **Phase 12.2 COMPLETE**: Batch generation and result reporting
  - ✅ Batch generation from YAML (`viroforge batch`)
  - ✅ Parameter sweep support (itertools.product)
  - ✅ Sequential and parallel execution
  - ✅ Result reporting (`viroforge report`)
  - ✅ Dataset comparison (`viroforge compare`)
  - ✅ 5 example batch configurations
  - ✅ Intelligent recommendations (hybrid assembly, platform comparison)
- ✅ **Phase 12.3 COMPLETE**: Web interface
  - ✅ Flask-based web application (`viroforge web`)
  - ✅ Modern Bootstrap 5 responsive UI
  - ✅ Visual collection browser with search/filter
  - ✅ Interactive dataset generation with progress monitoring
  - ✅ Batch configuration builder with YAML editor
  - ✅ Dataset reporting and comparison dashboards
  - ✅ RESTful API (10+ endpoints)
  - ✅ 9 HTML templates with JavaScript interactivity
- **Timeline**: Phase 12.1 (1 day), Phase 12.2 (1 day), Phase 12.3 (1 day) - Total 3 days

## Database

**Location**: `viroforge/data/viral_genomes.db`
**Type**: SQLite
**Contents**:
- 14,423 RefSeq viral genomes
- ICTV taxonomy (57.1% matched via enhanced fuzzy matching)
- 28 body site collections with metadata

**Key Tables**:
- `genomes`: Genome sequences and metadata
- `taxonomy`: ICTV taxonomy mappings
- `body_site_collections`: Collection metadata
- `collection_genomes`: Collection-genome associations with abundances

## Collection Curation Workflow

All new collections follow this established protocol:

### 1. Research Phase
- Literature review for virome composition
- Identify dominant viruses and prevalence data
- Document clinical/research applications

### 2. Script Creation
Create `scripts/curate_[collection]_virome_collection.py`:
- Separate methods for each viral category
- SQL queries targeting specific viruses
- Abundance assignment based on literature

### 3. Proactive Taxonomy Rescanning
**CRITICAL**: Always verify expected viruses are present after initial run.

**Common Issues & Fixes**:
- **Non-human matches**: Use specific patterns like `"Human bocavirus%"` not `"%bocavirus%"`
- **Substring matches**: Avoid broad patterns (e.g., `"%HIV%"` matches "Phives")
- **Family name changes**: ICTV taxonomy evolves (e.g., `Orthoherpesviridae` not `Herpesviridae`)
- **Missing families**: Some viruses have `family="Unknown"` (e.g., HCV)
- **Word boundaries**: Use `"Influenza A virus%"` to exclude parainfluenza

**Query Fix Patterns**:
```sql
-- ❌ Too broad
WHERE g.genome_name LIKE '%bocavirus%'
-- ✅ Specific
WHERE g.genome_name LIKE 'Human bocavirus%'

-- ❌ Substring issues
WHERE g.genome_name LIKE '%HIV%'
-- ✅ Word boundary
WHERE g.genome_name LIKE 'Human immunodeficiency virus%'

-- ❌ May include non-human
WHERE t.family = 'Anelloviridae'
-- ✅ Human-specific
WHERE g.genome_name LIKE 'Torque teno virus%'
```

**Explicit Exclusions** when needed:
```sql
WHERE g.genome_name LIKE 'Respiratory syncytial virus%'
  AND g.genome_name NOT LIKE '%Bovine%'
  AND g.genome_name NOT LIKE '%Potato%'
```

### 4. Verification
- Check all CRITICAL biomarker viruses are present (e.g., CMV for transplant, BK virus for nephritis)
- Verify genome counts match targets
- Confirm database insertion

### 5. Documentation
- Update `README.md` with collection details
- Include applications, literature basis, special notes
- Update badges (collection count, phase status)

### 6. Commit
Use detailed commit messages documenting:
- Collection composition
- Taxonomy fixes applied
- Literature basis
- Database updates

## Recent Collections (Phase 9)

### Collection 24: Vaginal Virome (26 genomes)
- **Key Fix**: Lactobacillus phages taxonomy (`Siphoviridae`/`Myoviridae` → modern genera)
- **Issue**: Herpesviruses used old family name (`Herpesviridae` → `Orthoherpesviridae`)

### Collection 25: Blood/Plasma Virome (21 genomes)
- **Key Fix**: Human HBV not in RefSeq → Used HCV/HIV instead
- **Issue**: Broad "HIV" pattern matched "Phives" and "Crohivirus"

### Collection 26: Ocular Surface Virome (17 genomes)
- **Key Fix**: Targeted clinically relevant herpesviruses (HSV-1, VZV, EBV, CMV)
- **Issue**: Random selection missed HSV-1 (most important for keratitis)

### Collection 27: Lower Respiratory/Lung Virome (31 genomes)
- **Key Fixes**:
  1. RSV matched bovine RSV and potato virus
  2. Influenza matched parainfluenza (wrong family!)
  3. Bocavirus matched bat bocavirus
  4. Random selection missed CMV (critical for transplant)

### Collection 28: Urinary Virome (20 genomes)
- **Key Fixes**:
  1. Anelloviruses matched seal/simian/rodent TTVs → `"Torque teno virus%"`
  2. Herpesviruses: Split queries to ensure CMV + EBV (not random)
- **Critical**: BK polyomavirus (60-100% in transplant complications)

## Key Development Patterns

### Testing Changes
```bash
# Run curation script
python scripts/curate_[collection]_virome_collection.py

# Verify database insertion
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT collection_id, collection_name, n_genomes, version
   FROM body_site_collections WHERE collection_id = XX;"
```

### Investigating Taxonomy
```bash
# Check virus names in database
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT genome_name, family, genus
   FROM genomes g
   LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
   WHERE genome_name LIKE '%virus_pattern%'
   LIMIT 20;"
```

### Git Workflow
```bash
git add scripts/curate_*.py README.md
git commit -m "feat: Phase X - Add Collection YY"
git push
```

**Note**: Database file (`viral_genomes.db`) is in `.gitignore`

## Phase 10 Summary - COMPLETE ✅

### All Components Completed
1. ✅ Research long-read simulators → **PBSIM3 selected**
2. ✅ Design integration architecture → `docs/PHASE10_ARCHITECTURE.md`
3. ✅ Implement `viroforge/simulators/longread.py` (850+ lines)
4. ✅ PacBio HiFi support (two-step: PBSIM3 CLR → ccs consensus)
5. ✅ Nanopore support (single-step: PBSIM3 with ONT error model)
6. ✅ Configuration classes (PacBioHiFiConfig, NanoporeConfig)
7. ✅ Ground truth tracking extended for long reads
8. ✅ VLP modeling updates for long-read size bias (60% reduction)
9. ✅ Full integration with `generate_fastq_dataset.py`
10. ✅ Comprehensive testing (`tests/test_longread_simulator.py`)
11. ✅ Complete user tutorial (`docs/LONGREAD_TUTORIAL.md`)
12. ✅ Documentation updates (README, ROADMAP, claude.md)

### Deliverables
- **Code**: 4 files modified, 2 new docs, 1 new test file
- **Tests**: 80+ unit tests covering all configurations
- **Documentation**: 20+ page tutorial with benchmarking workflows
- **Timeline**: 3 weeks (on schedule)

### Long-Read Simulator API

```python
from viroforge.simulators import (
    generate_long_reads,
    LongReadPlatform,
    PacBioHiFiConfig,
    NanoporeConfig
)

# PacBio HiFi with custom configuration
config = PacBioHiFiConfig(passes=15, read_length_mean=20000)
output = generate_long_reads(
    composition=composition,
    output_prefix='data/hifi_reads',
    platform=LongReadPlatform.PACBIO_HIFI,
    depth=10.0,
    platform_config=config,
    random_seed=42
)
# Returns: {'reads': Path('hifi.fastq.gz'), 'ground_truth': Path('ground_truth.tsv')}

# Nanopore with defaults
output = generate_long_reads(
    composition=composition,
    output_prefix='data/nanopore_reads',
    platform=LongReadPlatform.NANOPORE,
    depth=15.0,
    random_seed=42
)
# Returns: {'reads': Path('nanopore.fastq'), 'ground_truth': Path('ground_truth.tsv')}
```

## Important Files

- `ROADMAP.md`: Development roadmap and timelines
- `README.md`: User-facing documentation
- `docs/COLLECTION_IMPLEMENTATION_GUIDE.md`: Collection curation details
- `docs/PHASE10_LONGREAD_RESEARCH.md`: Long-read simulator evaluation
- `docs/PHASE10_ARCHITECTURE.md`: Long-read integration architecture
- `scripts/curate_*_virome_collection.py`: Collection curation scripts
- `scripts/generate_fastq_dataset.py`: Main data generation script
- `viroforge/simulators/longread.py`: PacBio HiFi and Nanopore simulator
- `viroforge/simulators/illumina.py`: Short-read simulator
- `viroforge/data/viral_genomes.db`: SQLite database (NOT in git)
- `lab-notebook/sessions/2025-11/20251110-001-IMPLEMENTATION-phase10-longread-simulator.md`: Phase 10 lab notebook

## References

See individual collection scripts and `docs/COLLECTION_IMPLEMENTATION_GUIDE.md` for literature citations.
