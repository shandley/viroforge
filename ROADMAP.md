# ViroForge Development Roadmap

**Version**: 0.4.0 → 1.0.0
**Timeline**: 3-6 months (complete)
**Goal**: Comprehensive virome simulation covering DNA/RNA, short/long-read, host/environmental, healthy/disease states

---

## Current Status (v0.10.0)

**Phase 12 Complete** - CLI Enhancements & Production Ready 🎉

✅ 14,423 RefSeq viral genomes with ICTV taxonomy (57.1% ICTV coverage)
✅ 28 curated collections (host-associated, environmental, disease states)
✅ **5 sequencing platforms** - NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore
✅ **Long-read simulation** - PBSIM3-based PacBio HiFi and Nanopore
✅ **Hybrid assembly support** - Matched short + long read datasets
✅ VLP enrichment (5 protocols) adapted for long reads
✅ RNA virome workflow with reverse transcription and rRNA depletion
✅ Amplification bias modeling (RdAB, MDA, Linker amplification)
✅ **Complete CLI** - browse, generate, batch, report, compare, presets
✅ **Web interface** - Bootstrap 5 UI with RESTful API
✅ **Configuration presets** - 8 built-in presets for common use cases
✅ **Batch generation** - YAML-based parameter sweeps
✅ Complete ground truth tracking for all platforms
✅ Comprehensive documentation with literature validation
✅ 80+ comprehensive tests

**Next Phase: Benchmarking Framework** - Transform ViroForge from data generator to complete validation platform

---

## Development Path B: Comprehensive Enhancement

### **PHASE 6: Enhanced Library Preparation**

**Timeline**: 1-2 weeks
**Status**: ✅ Complete

#### Objectives
- Integrate existing amplification bias code into user workflow
- Enable comparison of library prep methods
- Model realistic PCR artifacts

#### Tasks
- [x] Create development roadmap
- [ ] Add `--amplification` flag to `generate_fastq_dataset.py`
- [ ] Integrate RdAB, MDA, Linker amplification into workflow
- [ ] Add amplification to batch generation presets
- [ ] Create amplification comparison preset
- [ ] Update metadata to include amplification parameters
- [ ] Add integration tests for amplification
- [ ] Update documentation (QUICKSTART, PHASE4, tutorials)
- [ ] Create amplification comparison tutorial

#### Deliverables
- `--amplification {none,rdab,mda,linker}` flag in all scripts
- Amplification statistics in metadata.json
- Tutorial: "Comparing Library Preparation Methods"
- Tests: Integration tests for all 4 amplification methods

#### Impact
- **HIGH**: Affects all virome studies
- **Effort**: LOW (code exists, needs integration)
- **Users**: Pipeline developers, method comparison studies

---

### **PHASE 7: Critical Collections**

**Timeline**: 3-4 weeks
**Status**: ✅ Complete

#### 7.1 Wastewater Virome (Week 1)

**Rationale**: Highly timely for epidemiological surveillance

**Collection Design**:
- **Name**: Wastewater Virome - Urban Treatment Plant
- **Size**: 300-500 genomes
- **Composition**:
  - Human enteric viruses (40%): Norovirus, rotavirus, adenovirus, astrovirus, sapovirus
  - Bacteriophages (35%): Gut bacterial phages
  - Environmental viruses (15%): Plant viruses, insect viruses
  - Emerging pathogens (10%): SARS-CoV-2, mpox, polio, etc.
- **Literature**: Crits-Christoph et al. 2021, Crank et al. 2022

**Tasks**:
- [x] Literature review of wastewater virome composition
- [x] Curate genome list from RefSeq
- [x] Assign literature-validated abundances
- [x] Add to database as collection ID 17
- [x] Validate collection
- [ ] Create wastewater-specific contamination profile (future enhancement)
- [x] Document in Collection Implementation Guide

#### 7.2 Disease State Collections (Weeks 2-4)

**Collection 18: IBD Gut Virome** (Week 2)
- **Comparison**: vs healthy gut (collection 9)
- **Key differences**: Lower diversity, altered Caudovirales, increased temperate phages
- **Size**: 80-100 genomes
- **Literature**: Norman et al. 2015, Zuo et al. 2019

**Collection 19: HIV+ Gut Virome** (Week 3)
- **Comparison**: vs healthy gut (collection 9)
- **Key differences**: Dramatically reduced diversity, altered composition, increased eukaryotic viruses
- **Size**: 40-60 genomes
- **Literature**: Handley et al. 2012, Monaco et al. 2016

**Collection 20: Cystic Fibrosis Respiratory Virome** (Week 4)
- **Comparison**: vs healthy respiratory (collection 12)
- **Key differences**: Chronic viral infections, altered bacteriophage composition
- **Size**: 25-35 genomes
- **Literature**: Lim et al. 2014, Willner et al. 2009

#### Deliverables
- 4 new collections (1 wastewater, 3 disease states)
- Disease vs healthy comparison tutorials
- Dysbiosis modeling examples
- Updated Collection Implementation Guide

#### Impact
- **VERY HIGH**: Clinical validation, public health surveillance
- **Effort**: MEDIUM (1 week per collection)
- **Users**: Clinical researchers, epidemiologists, public health

---

### **PHASE 8: RNA Virome Workflow**

**Timeline**: 3-4 weeks
**Status**: ✅ Complete

#### Objectives
- Support RNA virus sequencing workflows
- Model reverse transcription step
- Implement rRNA depletion for RNA viromes

#### 8.1 RNA Virus Collections (Weeks 1-2)

**Collection 21: Human Respiratory RNA Virome**
- **Viruses**: Influenza, RSV, HMPV, coronaviruses, rhinoviruses, enteroviruses
- **Size**: 50-70 genomes
- **Type**: ssRNA (negative and positive sense)

**Collection 22: Arbovirus Environmental (Mosquito)**
- **Viruses**: Flaviviruses, alphaviruses, bunyaviruses
- **Size**: 30-50 genomes
- **Type**: Mixed ssRNA, dsRNA

**Collection 23: Fecal RNA Virome**
- **Viruses**: Norovirus, rotavirus, astrovirus, picornaviruses
- **Size**: 40-60 genomes
- **Type**: Mixed ssRNA, dsRNA

#### 8.2 RNA Workflow Implementation (Weeks 3-4)

**Tasks**:
- [ ] Create RNAViromeWorkflow class
- [ ] Implement reverse transcription modeling
  - RT efficiency by virus type
  - Random vs oligo-dT vs specific primers
  - RT artifacts (template switching, truncation)
- [ ] Implement rRNA depletion
  - Different from DNA rRNA removal
  - Ribo-Zero or Ribominus modeling
  - Variable efficiency
- [ ] Add RNA degradation modeling
  - More fragile than DNA
  - RNase contamination
- [ ] Create `--molecule-type {dna,rna}` flag
- [ ] Update contamination profiles for RNA
- [ ] Add RNA-specific tests
- [ ] Create RNA virome tutorial

#### Deliverables
- RNA virome workflow
- 3 RNA virus collections
- RNA-specific documentation
- Tutorial: "Generating RNA Virome Datasets"

#### Impact
- **HIGH**: Opens entirely new user base
- **Effort**: MEDIUM-HIGH (new workflow components)
- **Users**: Respiratory researchers, clinical diagnostics, environmental surveillance

---

### **PHASE 9: Additional Host Niches**

**Timeline**: 2-3 weeks
**Status**: ✅ Complete

#### Collections to Add

**Collection 24: Vaginal Virome** (Week 1)
- **Size**: 20-30 genomes
- **Composition**: Papillomaviruses, herpesviruses, bacteriophages
- **Literature**: Wylie et al. 2014, Dols et al. 2016
- **Impact**: Women's health, pregnancy outcomes

**Collection 25: Blood/Plasma Virome** (Week 1)
- **Size**: 15-25 genomes
- **Composition**: Anelloviruses, herpesviruses, parvoviruses
- **Literature**: Moustafa et al. 2017, Young et al. 2016
- **Impact**: Viremia, transplant monitoring, diagnostics

**Collection 26: Ocular Surface Virome** (Week 2)
- **Size**: 10-20 genomes
- **Composition**: Adenoviruses, herpesviruses, bacteriophages
- **Literature**: Doan et al. 2016
- **Impact**: Ophthalmology, infectious keratitis

**Collection 27: Lower Respiratory (Lung) Virome** (Week 2)
- **Size**: 25-40 genomes
- **Composition**: Different from nasopharynx, includes lung-resident viruses
- **Literature**: Kitsios et al. 2018
- **Impact**: Pneumonia, COPD, lung transplant

**Collection 28: Urinary Virome** (Week 3)
- **Size**: 15-25 genomes
- **Composition**: BK virus, JC virus, adenoviruses
- **Literature**: Santiago-Rodriguez et al. 2015
- **Impact**: Transplant monitoring, UTI

#### Deliverables
- 5 new host-associated collections
- Coverage of major human body sites
- Updated collection comparison tools

#### Impact
- **MEDIUM**: Expands applicable research niches
- **Effort**: LOW-MEDIUM (similar to existing collections)
- **Users**: Clinical researchers across specialties

---

### **PHASE 10: Long-Read Sequencing Support**

**Timeline**: 3 weeks
**Status**: ✅ Complete

#### Objectives
- ✅ Support PacBio HiFi and Nanopore platforms
- ✅ Enable complete genome assembly benchmarking
- ✅ Model long-read specific artifacts

#### Tasks
- [x] Research long-read simulators (pbsim3, NanoSim, PBSIM2)
  - Selected PBSIM3 (supports both PacBio and Nanopore)
  - Documented in `docs/PHASE10_LONGREAD_RESEARCH.md`
- [x] Integrate PacBio HiFi simulator
  - High accuracy (>99.9%, QV20+)
  - Two-step workflow: PBSIM3 CLR → ccs consensus
  - Configurable passes (3-20), read lengths (10-30kb)
- [x] Integrate Nanopore simulator
  - Homopolymer errors (hp_del_bias)
  - Ultra-long reads (10kb-2Mb)
  - R9.4 and R10.4 chemistry support
- [x] Add `--platform {novaseq,miseq,hiseq,pacbio-hifi,nanopore}` options
- [x] Update VLP modeling for long reads (60% size bias reduction)
  - Long reads span entire genomes → less sequencing bias
- [x] Create long-read specific tests (`tests/test_longread_simulator.py`)
- [x] Create comprehensive long-read tutorial (`docs/LONGREAD_TUTORIAL.md`)

#### Deliverables
- ✅ PacBio HiFi support with realistic CCS workflow
- ✅ Nanopore support with characteristic homopolymer errors
- ✅ Long-read assembly benchmarking capability
- ✅ Tutorial: "ViroForge Long-Read Sequencing Tutorial"
- ✅ 80+ unit tests including long-read configurations
- ✅ Architecture documentation (`docs/PHASE10_ARCHITECTURE.md`)

#### Key Files
- `viroforge/simulators/longread.py` (850+ lines)
- `scripts/generate_fastq_dataset.py` (updated with long-read routing)
- `viroforge/enrichment/vlp.py` (updated with read_type parameter)
- `tests/test_longread_simulator.py` (comprehensive test suite)
- `docs/LONGREAD_TUTORIAL.md` (complete user guide)

#### Impact
- **HIGH**: Future-proofing, better assemblies
- **Effort**: MEDIUM-HIGH (new external tools)
- **Users**: Assembly developers, complete genome studies

---

### **PHASE 11: Hybrid Assembly Support**

**Timeline**: 1 week
**Status**: ✅ Complete

#### Objectives
- ✅ Enable matched short + long read dataset generation
- ✅ Support hybrid assemblers (Unicycler, SPAdes hybrid, MaSuRCA)
- ✅ Provide composition validation utilities

#### Tasks
- [x] Create `generate_hybrid_dataset.py` convenience script
- [x] Ensure identical composition across short and long reads
- [x] Add `validate_hybrid_composition.py` validation utility
- [x] Create comprehensive tutorial (`docs/HYBRID_ASSEMBLY_TUTORIAL.md`)
- [x] Add hybrid preset (`hybrid-standard`)
- [x] Document supported hybrid assemblers

#### Deliverables
- ✅ Hybrid dataset generation with matched compositions
- ✅ Validation utilities for composition verification
- ✅ Tutorial: "Hybrid Assembly with ViroForge"
- ✅ Example workflows for Unicycler, SPAdes, MaSuRCA

#### Key Files
- `scripts/generate_hybrid_dataset.py` - Convenience wrapper
- `scripts/validate_hybrid_composition.py` - Composition validator
- `docs/HYBRID_ASSEMBLY_TUTORIAL.md` - Complete user guide
- `viroforge/presets/hybrid-standard.yaml` - Hybrid preset

#### Impact
- **HIGH**: Superior assemblies from hybrid approaches
- **Effort**: LOW (leverages existing infrastructure)
- **Users**: Assembly researchers, complete genome reconstruction

---

### **PHASE 12: CLI Enhancements & Web Interface**

**Timeline**: 3 days (distributed across 3 subphases)
**Status**: ✅ Complete

#### Objectives
- ✅ Create unified CLI interface with subcommands
- ✅ Add interactive collection browser
- ✅ Implement configuration preset system
- ✅ Build batch generation framework
- ✅ Add dataset reporting and comparison tools
- ✅ Create modern web interface

#### Subphases

**Phase 12.1: Full Generate Command** (Day 1)
- [x] Preset-based generation system (8 built-in presets)
- [x] Real-time progress bars via rich.progress
- [x] Parameter override system
- [x] Verbose mode for debugging

**Phase 12.2: Batch Generation & Reporting** (Day 2)
- [x] YAML-based batch configuration
- [x] Parameter sweep support (cartesian product expansion)
- [x] Sequential and parallel execution modes
- [x] Dataset quality reporting (`viroforge report`)
- [x] Intelligent dataset comparison (`viroforge compare`)
- [x] 5 example batch configurations

**Phase 12.3: Web Interface** (Day 3)
- [x] Flask-based web application
- [x] Bootstrap 5 responsive UI (9 pages)
- [x] RESTful API (10+ endpoints)
- [x] Visual collection browser with search/filter
- [x] Interactive dataset generation with progress monitoring
- [x] Batch configuration builder
- [x] Dataset reporting and comparison dashboards

#### Deliverables
- ✅ Unified CLI: `viroforge {browse,generate,batch,report,compare,presets,web}`
- ✅ Interactive TUI collection browser
- ✅ 8 configuration presets
- ✅ Batch generation framework with YAML configs
- ✅ Dataset reporting and comparison tools
- ✅ Web interface with Bootstrap 5 UI
- ✅ Complete documentation for all features

#### Key Files
- `viroforge/cli/__init__.py` - Main CLI entry point
- `viroforge/cli/browse.py` - Interactive TUI browser (340 lines)
- `viroforge/cli/generate.py` - Generation command (340 lines)
- `viroforge/cli/batch.py` - Batch generation (358 lines)
- `viroforge/cli/report.py` - Dataset reporting (319 lines)
- `viroforge/cli/compare.py` - Dataset comparison (261 lines)
- `viroforge/cli/presets.py` - Preset management
- `viroforge/cli/web.py` - Web server launcher
- `viroforge/web/app.py` - Flask application (280 lines)
- `viroforge/web/templates/` - 9 HTML templates (~1,250 lines)
- `examples/batch_configs/` - 5 example configurations
- `docs/PHASE12.1_SUMMARY.md` - Generate command docs
- `docs/PHASE12.2_SUMMARY.md` - Batch/report/compare docs
- `docs/PHASE12.3_SUMMARY.md` - Web interface docs

#### Impact
- **VERY HIGH**: Production-ready UX
- **Effort**: LOW-MEDIUM (3 days total)
- **Users**: All users benefit from improved usability

---

### **PHASE 13: Benchmarking & Validation Framework**

**Timeline**: 6-8 weeks
**Status**: 📋 In Planning

#### Objectives
- Transform ViroForge from data generator to **complete validation platform**
- Enable comprehensive pipeline benchmarking across virome analysis workflow
- Provide modular benchmarking tools for each workflow stage
- Support multiple pipeline tools (Kraken2, Centrifuge, DIAMOND, etc.)

#### Strategic Rationale

**Current Gap**: ViroForge generates world-class synthetic data with perfect ground truth, but provides **zero tools** to help users benchmark their pipelines against that ground truth.

**Problem**: Users must manually:
1. Parse pipeline outputs (every tool has different format)
2. Match pipeline taxa to ground truth genomes
3. Calculate accuracy metrics
4. Create comparison plots
5. Write benchmark reports

**Solution**: Integrated benchmarking framework that takes pipeline output + ground truth → comprehensive HTML report with metrics, plots, error analysis.

#### Architecture

**Modular Design** - Match virome analysis workflow stages:

```
viroforge/
├── benchmarking/              # NEW MODULE (optional install)
│   ├── __init__.py
│   ├── parsers/              # Parse pipeline outputs
│   │   ├── kraken2.py        # Kraken2/Bracken format
│   │   ├── centrifuge.py     # Centrifuge format
│   │   ├── diamond.py        # DIAMOND BLAST output
│   │   ├── metaphlan.py      # MetaPhlAn profile
│   │   ├── kaiju.py          # Kaiju output
│   │   └── generic.py        # Generic TSV
│   ├── metrics/              # Calculate performance metrics
│   │   ├── taxonomic.py      # Classification accuracy
│   │   ├── abundance.py      # Quantification accuracy
│   │   ├── assembly.py       # Assembly quality
│   │   ├── completeness.py   # Genome recovery
│   │   └── contamination.py  # QC validation
│   ├── visualizations/       # Diagnostic plots
│   │   ├── composition.py    # Stacked bars, pie charts
│   │   ├── scatter.py        # Abundance correlations
│   │   ├── confusion.py      # Confusion matrices
│   │   ├── assembly.py       # Coverage, completeness
│   │   └── sankey.py         # Read flow diagrams
│   ├── reports/              # Report generation
│   │   ├── html_report.py    # Publication-quality HTML
│   │   ├── json_export.py    # Machine-readable JSON
│   │   └── templates/        # Jinja2 HTML templates
│   └── utils.py              # Taxonomy matching utilities
```

**Installation**: `pip install viroforge[benchmark]`

#### Benchmarking Modules

**Module 1: QC Benchmarking** (CRITICAL for viromes)
- Validate contamination removal (host DNA, bacterial, rRNA, PhiX)
- Measure false positive rate (viral reads wrongly removed)
- Measure viral retention rate
- **Metrics**: Precision, recall, F1 for each contamination type

**Module 2: Assembly Benchmarking** (MOST CRITICAL for viromes)
- Align assembled contigs to true viral genomes
- Calculate genome recovery (complete, high-quality, partial, missing)
- Detect chimeric contigs (segments from multiple genomes)
- Measure coverage uniformity and bias
- **Metrics**: Completeness, identity, N50, chimera rate, coverage bias

**Module 3: Binning Benchmarking**
- Validate vMAG (viral MAG) reconstruction
- Measure bin purity (one genome per bin)
- Detect mixed bins (multiple genomes wrongly binned together)
- **Metrics**: Bin purity, precision, recall, bin quality (CheckV-style)

**Module 4: Taxonomy Benchmarking** (read + contig-based)
- Compare pipeline classifications to ground truth
- Calculate accuracy at multiple taxonomic levels (species, genus, family)
- Identify false positives/negatives by taxonomy
- **Metrics**: Precision, recall, F1, accuracy at each level

**Module 5: Completeness Benchmarking**
- Measure genome recovery across coverage ranges
- Analyze completeness by genome length, GC content
- Identify systematic biases in recovery
- **Metrics**: Mean completeness, completeness by coverage/length

**Module 6: Annotation Benchmarking**
- Validate gene calling accuracy (ORF prediction)
- Validate functional annotation accuracy
- **Metrics**: Gene calling precision/recall, function precision/recall

**Module 7: Host Prediction Benchmarking**
- Validate virus-host linkage predictions
- **Metrics**: Accuracy by host type (bacteriophage, human virus, etc.)

**Module 8: Novel Discovery Benchmarking** (Advanced)
- Test ability to find viruses NOT in reference databases
- Holdout testing: Withhold 20% of genomes from DB
- **Metrics**: Discovery rate, false discovery rate

**Module 9: End-to-End Pipeline Benchmarking**
- Comprehensive benchmarking of entire pipeline
- Combines all modules into single comprehensive report
- **Output**: Multi-page HTML report with all 8 modules

#### CLI Interface

```bash
# Modular approach - test individual components
viroforge benchmark qc          # Contamination removal
viroforge benchmark assembly    # Assembly quality
viroforge benchmark binning     # MAG reconstruction
viroforge benchmark taxonomy    # Classification
viroforge benchmark completeness # Genome recovery
viroforge benchmark annotation  # Gene calling + function
viroforge benchmark host        # Host prediction
viroforge benchmark discovery   # Novel virus detection

# Comprehensive end-to-end
viroforge benchmark pipeline    # All modules together

# Multi-pipeline comparison
viroforge benchmark \
  --ground-truth data/gut/metadata/metadata.json \
  --pipelines kraken2:results/k2.txt centrifuge:results/cf.txt \
  --output reports/comparison.html

# Batch benchmarking across collections
viroforge benchmark-batch \
  --datasets data/*_standard/ \
  --pipeline-results results/ \
  --format kraken2
```

#### ViroForge Metadata Enhancements

**Phase 13A: Essential Enhancements**

Add to metadata.json:
```json
{
  "benchmarking": {
    "version": "1.0",
    "capabilities": ["qc", "assembly", "taxonomy", "completeness"],

    "contamination_manifest": {
      "host_dna": {"n_sequences": 10, "total_abundance": 0.035, "sequences": [...]},
      "rrna": {...},
      "bacterial": {...},
      "phix": {...}
    },

    "expected_coverage": {
      "mean_coverage": 30.0,
      "by_genome": {
        "GCF_015160975.1": {
          "expected_coverage": 68.7,
          "expected_completeness": 0.999,
          "expected_n_reads": 950000
        }
      }
    },

    "read_manifest": {
      "enabled": true,
      "path": "metadata/read_manifest.tsv.gz",
      "total_reads": 10000000
    }
  }
}
```

**Phase 13B: Advanced Enhancements**
- Gene annotations export (from RefSeq CDS)
- Expected assembly output (optional, expensive)
- Strain variant tracking (future)

#### Implementation Roadmap

**Phase 13A: Foundation + Minimal Enhancements** (Weeks 1-2) - ✅ COMPLETE (2025-11-11)
- [x] Add `enable_benchmarking` flag to metadata generation
- [x] Implement contamination manifest export
- [x] Calculate expected coverage per genome
- [ ] Create `viroforge/benchmarking/` module structure
- [ ] Implement metadata loader (parses enhanced metadata)
- [ ] Build core utilities (taxonomy matching, alignment)

**Phase 13B: QC + Assembly Benchmarking** (Weeks 3-5)
- [ ] Implement read manifest tracking (optional)
- [ ] Module 1: QC benchmarking (contamination removal validation)
- [ ] Module 2: Assembly benchmarking (contig-to-genome alignment)
- [ ] Basic HTML report templates
- [ ] Unit tests for QC and assembly modules

**Phase 13C: Taxonomy + Completeness** (Weeks 6-8)
- [ ] Module 4: Taxonomy benchmarking (Kraken2, Centrifuge, DIAMOND parsers)
- [ ] Module 5: Completeness analysis (genome recovery metrics)
- [ ] Comprehensive HTML reports with visualizations
- [ ] Integration tests with example datasets
- [ ] Documentation and tutorials

**Phase 13D: Advanced Modules** (Future)
- [ ] Module 3: Binning benchmarking
- [ ] Module 6: Annotation benchmarking (requires gene annotation export)
- [ ] Module 7: Host prediction benchmarking
- [ ] Module 8: Novel discovery benchmarking
- [ ] Module 9: End-to-end pipeline benchmarking

#### Deliverables

**MVP (Phase 13A-C)**:
- ✅ Enhanced metadata with contamination manifest, expected coverage
- ✅ Module 1: QC benchmarking
- ✅ Module 2: Assembly benchmarking
- ✅ Module 4: Taxonomy benchmarking
- ✅ Module 5: Completeness analysis
- ✅ HTML reports with metrics + visualizations
- ✅ CLI: `viroforge benchmark {qc,assembly,taxonomy,completeness}`
- ✅ Documentation: Benchmarking tutorial

**Full Implementation (Phase 13D)**:
- All 9 benchmarking modules
- Multi-pipeline comparison
- Batch benchmarking
- Comprehensive end-to-end reports
- Publication-ready figures

#### Dependencies

New optional dependencies:
```python
extras_require={
    "benchmark": [
        "matplotlib>=3.5.0",      # Plotting
        "seaborn>=0.11.0",        # Statistical visualizations
        "scikit-learn>=1.0.0",    # Metrics calculations
        "plotly>=5.0.0",          # Interactive plots
        "jinja2>=3.0.0",          # HTML report templates
    ]
}
```

#### Impact

**For Pipeline Developers**:
- ✅ Rigorous validation before publication
- ✅ Identify weak spots (low-abundance taxa, marine viruses, etc.)
- ✅ Optimize parameters (database, thresholds, filters)
- ✅ Publication-ready benchmark figures

**For Bioinformaticians**:
- ✅ Choose best tool for their data type
- ✅ Understand trade-offs (speed vs accuracy)
- ✅ Confidence in results

**For ViroForge Project**:
- ✅ Completes "generate → benchmark → publish" workflow
- ✅ Positions as CAMI-equivalent for viromes
- ✅ Dramatically increases utility and adoption
- ✅ Enables community-wide standardized benchmarking

**For the Field**:
- ✅ Standardized virome pipeline benchmarking
- ✅ Reproducible method comparisons
- ✅ Better tools through rigorous testing
- ✅ Increased confidence in virome studies

**Comparison to Bacterial Tools**:
- CAMI/OPAL focuses on bacterial metagenomes
- ViroForge benchmarking will be virome-specific:
  - Assembly-centric (not just taxonomy)
  - Completeness-focused (viral genomes are small)
  - QC validation (contamination removal critical)
  - Novel discovery (most viruses unknown)

**Priority**: **VERY HIGH** - This is arguably more important than additional collections because it unlocks the value of all existing data generation capabilities.

**Effort**: MEDIUM-HIGH (6-8 weeks for MVP, more for full implementation)

**Users**: All ViroForge users (pipeline developers, bioinformaticians, researchers)

---

## Version Milestones

### **v0.5.0 - Enhanced Library Prep** (2 weeks from now)
- Amplification bias integration
- Amplification comparison tools
- Updated documentation

### **v0.6.0 - Critical Collections** (6 weeks from now)
- Wastewater virome
- 3 disease state collections
- Disease vs healthy comparison framework

### **v0.7.0 - RNA Virome Support** (10 weeks from now)
- RNA workflow implementation
- 3 RNA virus collections
- RNA-specific documentation

### **v0.8.0 - Expanded Host Niches** (13 weeks from now)
- 5 additional host-associated collections
- Comprehensive human body site coverage

### **v1.0.0 - Feature Complete** (16-20 weeks from now)
- Long-read support (PacBio HiFi minimum)
- Temporal dynamics framework
- 30+ curated collections
- Comprehensive benchmarking platform

---

## Publication Strategy

### **Manuscript Timeline**

**Phase 1: Initial Submission** (After v0.7.0, ~10 weeks)
- Title: "ViroForge: A Comprehensive Synthetic Virome Data Generator for Pipeline Benchmarking"
- Content: DNA viromes, VLP modeling, amplification bias, disease states, RNA support
- Target: Bioinformatics, BMC Bioinformatics, or Microbiome

**Phase 2: Major Update** (After v1.0.0, ~20 weeks)
- Update with long-read support, temporal dynamics
- Expanded collection catalog
- Comprehensive validation studies

### **Preprint Strategy**

- Post to bioRxiv after v0.6.0 (wastewater + disease collections)
- Update preprint with v1.0.0 features
- Submit to peer review after v0.7.0 (RNA support)

---

## Success Metrics

### **Technical Metrics**
- [ ] 30+ curated virome collections
- [ ] DNA + RNA workflow support
- [ ] 5+ sequencing platforms (Illumina + PacBio + Nanopore)
- [ ] 4+ library prep methods
- [ ] Temporal dynamics framework
- [ ] 100% test coverage on critical paths

### **Community Metrics**
- [ ] 100+ GitHub stars
- [ ] 10+ citations in peer-reviewed literature
- [ ] Used in published benchmarking studies
- [ ] 50+ downloads per month
- [ ] Active community contributions

### **Scientific Impact**
- [ ] Enables standardized virome pipeline benchmarking
- [ ] Cited in clinical virome studies
- [ ] Used for method comparison papers
- [ ] Adopted by pipeline developers (Hecatomb, others)

---

## Resource Requirements

### **Computational**
- Database storage: ~500 MB (current) → ~2 GB (v1.0.0)
- Test datasets: ~10 GB
- CI/CD: GitHub Actions (free tier sufficient)

### **Time Commitment**
- **Phase 6**: 1-2 weeks (amplification)
- **Phase 7**: 3-4 weeks (critical collections)
- **Phase 8**: 3-4 weeks (RNA workflow)
- **Phase 9-13**: 8-10 weeks (additional collections, long-read, temporal)
- **Total**: 16-20 weeks for v1.0.0

### **Expertise Needed**
- Virome biology (curating collections)
- Bioinformatics (workflow implementation)
- Literature review (validation)
- Technical writing (documentation)

---

## Risk Assessment

### **Technical Risks**

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Long-read simulator integration issues | Medium | Medium | Start with PacBio HiFi only, well-documented tools |
| RNA workflow complexity | Medium | High | Leverage existing RT-PCR knowledge, literature |
| Database size growth | Low | Low | SQLite handles 2GB easily, can migrate to PostgreSQL |
| Test suite maintenance | Medium | Medium | Incremental testing, dry-run modes for CI/CD |

### **Scientific Risks**

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Literature gaps for collections | Low | Medium | Start with well-studied conditions, expert consultation |
| Parameter validation challenges | Low | High | Conservative estimates, sensitivity analysis |
| User adoption | Medium | High | Excellent documentation, tutorials, examples |

---

## Community Engagement

### **Documentation**
- Comprehensive tutorials for each phase
- Video tutorials (YouTube)
- Workshop materials
- Example workflows

### **Outreach**
- Present at conferences (ISME, ASM, Microbiome meetings)
- Twitter/X announcements for new features
- Engage with virome research community
- Collaborate with pipeline developers

### **Support**
- GitHub Discussions for Q&A
- Responsive issue tracking
- Monthly release notes
- Changelog documentation

---

## Conclusion

ViroForge v1.0.0 will be the **most comprehensive virome simulation platform** available, covering:

- **DNA + RNA viromes**
- **Host-associated + environmental samples**
- **Healthy + disease states**
- **Short-read + long-read sequencing**
- **Multiple library prep methods**
- **Temporal dynamics**
- **30+ curated collections**
- **Complete ground truth for rigorous benchmarking**

**Timeline**: 16-20 weeks to v1.0.0
**Impact**: Essential tool for virome research community
**Publication**: Strong methods paper in top-tier bioinformatics journal

---

**Let's build the future of virome benchmarking! 🦠🔨**
