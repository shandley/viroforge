# Phase 13: Implementation Checklist

**Author**: ViroForge Development Team
**Date**: 2025-11-11
**Status**: Ready to implement

---

## Quick Reference

**Timeline**: 6-8 weeks for MVP (Phases 13A-C)
**Goal**: Transform ViroForge into complete validation platform
**Priority**: VERY HIGH

**Phases**:
- Phase 13A (Weeks 1-2): Foundation + metadata enhancements
- Phase 13B (Weeks 3-5): QC + Assembly benchmarking
- Phase 13C (Weeks 6-8): Taxonomy + Completeness benchmarking
- Phase 13D (Future): Advanced modules

---

## Phase 13A: Foundation + Minimal Enhancements (Weeks 1-2)

### ViroForge Metadata Enhancements

#### 1. Add Benchmarking Flag
- [ ] `scripts/generate_fastq_dataset.py`: Add `--enable-benchmarking` flag
- [ ] `scripts/generate_fastq_dataset.py`: Add `enable_benchmarking` parameter to `save_metadata()`
- [ ] Default: `True` (enabled by default in v0.11.0)

#### 2. Contamination Manifest Export
- [ ] `scripts/generate_fastq_dataset.py`: Create `_export_contamination_manifest()` method
- [ ] Export `ContaminationProfile` data to metadata JSON
- [ ] Include: host_dna, rrna, bacterial, phix sequences and abundances
- [ ] Test: Verify manifest matches actual contamination added

#### 3. Expected Coverage Calculation
- [ ] `scripts/generate_fastq_dataset.py`: Create `_calculate_expected_coverage()` method
- [ ] Calculate per-genome expected coverage from abundances + coverage parameter
- [ ] Calculate expected completeness (based on coverage)
- [ ] Calculate expected number of reads per genome
- [ ] Test: Verify calculations match observed coverage

#### 4. Enhanced Metadata Schema
- [ ] Update metadata version to `1.1`
- [ ] Add `benchmarking` section to metadata.json:
  ```json
  {
    "benchmarking": {
      "version": "1.0",
      "capabilities": ["qc", "assembly", "taxonomy", "completeness"],
      "contamination_manifest": {...},
      "expected_coverage": {...}
    }
  }
  ```
- [ ] Maintain backward compatibility (old format still works)

### Benchmarking Module Infrastructure

#### 5. Module Structure
- [ ] Create `viroforge/benchmarking/` directory
- [ ] Create `viroforge/benchmarking/__init__.py`
- [ ] Create `viroforge/benchmarking/parsers/__init__.py`
- [ ] Create `viroforge/benchmarking/metrics/__init__.py`
- [ ] Create `viroforge/benchmarking/visualizations/__init__.py`
- [ ] Create `viroforge/benchmarking/reports/__init__.py`
- [ ] Create `viroforge/benchmarking/reports/templates/` directory

#### 6. Core Utilities
- [ ] `viroforge/benchmarking/utils.py`: Taxonomy matching utilities
  - [ ] `match_taxonomy()`: Match pipeline taxonomy to ground truth
  - [ ] `normalize_taxon_name()`: Handle name variations
  - [ ] `taxonomy_distance()`: Calculate taxonomic distance
- [ ] `viroforge/benchmarking/utils.py`: Sequence alignment utilities
  - [ ] `align_contig_to_genome()`: Wrapper for minimap2/BLAST
  - [ ] `calculate_identity()`: Sequence identity percentage
  - [ ] `detect_chimera()`: Identify chimeric contigs
- [ ] `viroforge/benchmarking/utils.py`: Metadata loader
  - [ ] `load_ground_truth()`: Parse ViroForge metadata.json
  - [ ] `validate_metadata()`: Check for required fields
  - [ ] `get_contamination_manifest()`: Extract contamination data

#### 7. Base Classes
- [ ] `viroforge/benchmarking/metrics/__init__.py`: `BaseMetric` class
- [ ] `viroforge/benchmarking/visualizations/__init__.py`: `BaseVisualization` class
- [ ] `viroforge/benchmarking/reports/__init__.py`: `BaseReport` class

#### 8. Unit Tests
- [ ] `tests/test_benchmarking_utils.py`: Test taxonomy matching
- [ ] `tests/test_benchmarking_utils.py`: Test alignment utilities
- [ ] `tests/test_benchmarking_utils.py`: Test metadata loader
- [ ] Create test fixtures: `tests/fixtures/test_metadata.json`

#### 9. setup.py Updates
- [ ] Add `[benchmark]` extras_require:
  ```python
  "benchmark": [
      "matplotlib>=3.5.0",
      "seaborn>=0.11.0",
      "scikit-learn>=1.0.0",
      "plotly>=5.0.0",
      "jinja2>=3.0.0",
  ]
  ```

### Deliverables (Phase 13A)
- ✅ Enhanced metadata format (v1.1) with contamination manifest + expected coverage
- ✅ Core benchmarking infrastructure (~200 lines)
- ✅ Taxonomy matching and alignment utilities
- ✅ Unit tests for core functionality

---

## Phase 13B: QC + Assembly Benchmarking (Weeks 3-5)

### Module 1: QC Benchmarking

#### 10. QC Metrics
- [ ] `viroforge/benchmarking/metrics/contamination.py`: Create module
- [ ] `calculate_removal_rate()`: Contamination removal metrics
- [ ] `calculate_viral_retention()`: Viral read retention
- [ ] `calculate_false_positive_rate()`: Over-filtering detection
- [ ] `calculate_precision_recall()`: Per-contaminant-type metrics
- [ ] Support contamination types: host_dna, rrna, bacterial, phix

#### 11. Read Manifest (Optional)
- [ ] `scripts/generate_fastq_dataset.py`: Add `--read-manifest` flag
- [ ] Track which genome each read comes from
- [ ] Compress manifest: gzip TSV format
- [ ] `read_manifest.tsv.gz`: columns = [read_id, genome_id, sequence_type]

#### 12. QC Visualizations
- [ ] `viroforge/benchmarking/visualizations/qc_plots.py`: Create module
- [ ] `plot_removal_rates()`: Bar chart of removal rates by contaminant type
- [ ] `plot_viral_retention()`: Scatter plot of viral retention
- [ ] `plot_precision_recall()`: Precision-recall bar chart

#### 13. QC Report
- [ ] `viroforge/benchmarking/reports/templates/qc_report.html`: Jinja2 template
- [ ] `viroforge/benchmarking/reports/html_report.py`: `generate_qc_report()`
- [ ] Include: Metrics tables, plots, recommendations

### Module 2: Assembly Benchmarking

#### 14. Assembly Metrics
- [ ] `viroforge/benchmarking/metrics/assembly.py`: Create module
- [ ] `align_contigs_to_genomes()`: Align contigs to ground truth genomes
- [ ] `calculate_genome_recovery()`: Complete, high-quality, partial, missing
- [ ] `calculate_completeness()`: Per-genome completeness percentage
- [ ] `detect_chimeras()`: Identify contigs matching multiple genomes
- [ ] `calculate_coverage_bias()`: GC bias, length bias
- [ ] `calculate_assembly_stats()`: N50, L50, longest contig

#### 15. Assembly Visualizations
- [ ] `viroforge/benchmarking/visualizations/assembly_plots.py`: Create module
- [ ] `plot_recovery_heatmap()`: Genome recovery quality heatmap
- [ ] `plot_completeness_distribution()`: Histogram of completeness values
- [ ] `plot_coverage_uniformity()`: Coverage along genome position
- [ ] `plot_gc_bias()`: Scatter plot of GC% vs coverage

#### 16. Assembly Report
- [ ] `viroforge/benchmarking/reports/templates/assembly_report.html`: Jinja2 template
- [ ] `viroforge/benchmarking/reports/html_report.py`: `generate_assembly_report()`
- [ ] Include: Genome recovery table, completeness heatmap, coverage plots

### CLI Commands

#### 17. Benchmark CLI
- [ ] `viroforge/cli/benchmark.py`: Create module
- [ ] Add subparser for `viroforge benchmark` command
- [ ] `viroforge benchmark qc`: QC benchmarking subcommand
- [ ] `viroforge benchmark assembly`: Assembly benchmarking subcommand
- [ ] Common args: `--ground-truth`, `--output`, `--format`

### HTML Report Infrastructure

#### 18. Base Template
- [ ] `viroforge/benchmarking/reports/templates/base.html`: Base Jinja2 template
- [ ] Include Bootstrap 5 CSS
- [ ] Include Plotly/Matplotlib integration
- [ ] Header/footer navigation

#### 19. Report Generator
- [ ] `viroforge/benchmarking/reports/html_report.py`: `HTMLReportGenerator` class
- [ ] `render_template()`: Jinja2 template rendering
- [ ] `embed_plot()`: Embed matplotlib/plotly figures
- [ ] `generate_toc()`: Table of contents generation

### Testing

#### 20. Unit Tests
- [ ] `tests/test_qc_metrics.py`: Test contamination removal metrics
- [ ] `tests/test_assembly_metrics.py`: Test genome recovery calculations
- [ ] `tests/test_assembly_metrics.py`: Test chimera detection
- [ ] `tests/test_visualizations.py`: Test plot generation (no display)

#### 21. Integration Tests
- [ ] `tests/test_qc_benchmark_integration.py`: End-to-end QC benchmark
- [ ] `tests/test_assembly_benchmark_integration.py`: End-to-end assembly benchmark
- [ ] Create test fixtures: fake pipeline outputs

### Deliverables (Phase 13B)
- ✅ Module 1: QC Benchmarking (functional)
- ✅ Module 2: Assembly Benchmarking (functional)
- ✅ CLI: `viroforge benchmark {qc,assembly}`
- ✅ HTML reports with plots
- ✅ Unit + integration tests (~800 lines of code)

---

## Phase 13C: Taxonomy + Completeness (Weeks 6-8)

### Module 4: Taxonomy Benchmarking

#### 22. Pipeline Parsers
- [ ] `viroforge/benchmarking/parsers/kraken2.py`: Kraken2 report parser
  - [ ] Parse Kraken2 report format
  - [ ] Extract taxonomy + abundances
  - [ ] Return standardized format
- [ ] `viroforge/benchmarking/parsers/centrifuge.py`: Centrifuge output parser
- [ ] `viroforge/benchmarking/parsers/diamond.py`: DIAMOND BLAST parser
  - [ ] Parse BLAST tabular format
  - [ ] Extract best hit per query
- [ ] `viroforge/benchmarking/parsers/generic.py`: Generic TSV parser
  - [ ] Columns: genome_id, abundance, taxonomy
  - [ ] Flexible format support

#### 23. Taxonomic Metrics
- [ ] `viroforge/benchmarking/metrics/taxonomic.py`: Create module
- [ ] `calculate_classification_metrics()`: TP/FP/TN/FN at species/genus/family levels
- [ ] `calculate_precision_recall_f1()`: Classification accuracy metrics
- [ ] `identify_false_positives()`: List of wrongly called taxa
- [ ] `identify_false_negatives()`: List of missed taxa
- [ ] `calculate_confusion_matrix()`: Confusion matrix by taxonomy

#### 24. Abundance Metrics
- [ ] `viroforge/benchmarking/metrics/abundance.py`: Create module
- [ ] `calculate_correlation()`: Pearson, Spearman correlation
- [ ] `calculate_errors()`: MAE, RMSE, MAPE
- [ ] `calculate_bray_curtis()`: Beta diversity dissimilarity
- [ ] `calculate_residuals()`: Observed - expected abundances

#### 25. Taxonomy Visualizations
- [ ] `viroforge/benchmarking/visualizations/taxonomy_plots.py`: Create module
- [ ] `plot_composition_comparison()`: Stacked bar charts (ground truth vs observed)
- [ ] `plot_abundance_scatter()`: Scatter plot with regression line
- [ ] `plot_confusion_matrix()`: Heatmap of confusion matrix
- [ ] `plot_residuals()`: Residual plot for systematic bias

#### 26. Taxonomy Report
- [ ] `viroforge/benchmarking/reports/templates/taxonomy_report.html`: Jinja2 template
- [ ] `viroforge/benchmarking/reports/html_report.py`: `generate_taxonomy_report()`
- [ ] Include: Classification metrics table, abundance scatter, confusion matrix

### Module 5: Completeness Analysis

#### 27. Completeness Metrics
- [ ] `viroforge/benchmarking/metrics/completeness.py`: Create module
- [ ] `calculate_completeness_by_coverage()`: Group by coverage ranges
- [ ] `calculate_completeness_by_length()`: Group by genome length
- [ ] `calculate_completeness_by_gc()`: Group by GC content
- [ ] `identify_biases()`: Detect systematic recovery biases

#### 28. Completeness Visualizations
- [ ] `viroforge/benchmarking/visualizations/completeness_plots.py`: Create module
- [ ] `plot_completeness_vs_coverage()`: Scatter plot with trend line
- [ ] `plot_completeness_distribution()`: Histogram by category
- [ ] `plot_recovery_by_category()`: Bar charts by length/GC

#### 29. Completeness Report
- [ ] `viroforge/benchmarking/reports/templates/completeness_report.html`: Jinja2 template
- [ ] `viroforge/benchmarking/reports/html_report.py`: `generate_completeness_report()`

### CLI Updates

#### 30. Additional CLI Commands
- [ ] `viroforge benchmark taxonomy`: Taxonomy benchmarking subcommand
  - [ ] `--format {kraken2,centrifuge,diamond,generic}`: Pipeline format
  - [ ] `--pipeline-output`: Path to pipeline results
  - [ ] `--mode {read-based,contig-based}`: Classification mode
- [ ] `viroforge benchmark completeness`: Completeness analysis subcommand
  - [ ] `--contigs`: Path to assembled contigs
  - [ ] `--coverage`: Path to coverage file

### Testing

#### 31. Parser Tests
- [ ] `tests/test_kraken2_parser.py`: Test Kraken2 report parsing
- [ ] `tests/test_centrifuge_parser.py`: Test Centrifuge output parsing
- [ ] `tests/test_diamond_parser.py`: Test DIAMOND BLAST parsing
- [ ] Create test fixtures: example pipeline outputs

#### 32. Metrics Tests
- [ ] `tests/test_taxonomic_metrics.py`: Test classification metrics
- [ ] `tests/test_abundance_metrics.py`: Test abundance correlations
- [ ] `tests/test_completeness_metrics.py`: Test completeness calculations

#### 33. Integration Tests
- [ ] `tests/test_taxonomy_benchmark_integration.py`: End-to-end taxonomy benchmark
- [ ] `tests/test_completeness_benchmark_integration.py`: End-to-end completeness analysis

### Documentation

#### 34. Tutorial
- [ ] `docs/BENCHMARKING_TUTORIAL.md`: Complete user tutorial
  - [ ] Quick start example
  - [ ] QC benchmarking walkthrough
  - [ ] Assembly benchmarking walkthrough
  - [ ] Taxonomy benchmarking walkthrough
  - [ ] Multi-pipeline comparison
  - [ ] Interpreting results

#### 35. API Documentation
- [ ] `docs/BENCHMARKING_API.md`: Programmatic usage reference
  - [ ] Importing modules
  - [ ] Using parsers programmatically
  - [ ] Calculating metrics in Python
  - [ ] Generating custom reports

#### 36. README Updates
- [ ] Update main `README.md` with benchmarking examples
- [ ] Update `docs/QUICKSTART.md` with benchmarking workflow
- [ ] Add benchmarking section to table of contents

### Deliverables (Phase 13C - MVP COMPLETE)
- ✅ Module 4: Taxonomy Benchmarking (functional)
- ✅ Module 5: Completeness Analysis (functional)
- ✅ Parsers: Kraken2, Centrifuge, DIAMOND, Generic TSV
- ✅ CLI: `viroforge benchmark {qc,assembly,taxonomy,completeness}`
- ✅ Comprehensive HTML reports with plots
- ✅ Tutorial documentation
- ✅ 10+ unit tests, 5+ integration tests
- ✅ ~1000 lines of code

**🎉 MVP COMPLETE at end of Phase 13C 🎉**

---

## Phase 13D: Advanced Modules (Future)

### Module 3: Binning Benchmarking
- [ ] `viroforge/benchmarking/metrics/binning.py`
- [ ] Bin purity metrics (one genome per bin)
- [ ] Mixed bin detection
- [ ] Bin quality metrics (CheckV-style)

### Module 6: Annotation Benchmarking
- [ ] Requires: Gene annotation export from RefSeq
- [ ] `viroforge/benchmarking/metrics/annotation.py`
- [ ] Gene calling accuracy
- [ ] Functional annotation accuracy

### Module 7: Host Prediction Benchmarking
- [ ] `viroforge/benchmarking/metrics/host_prediction.py`
- [ ] Virus-host linkage accuracy
- [ ] Accuracy by host type

### Module 8: Novel Discovery Benchmarking
- [ ] Holdout testing implementation
- [ ] Discovery rate metrics
- [ ] False discovery detection

### Module 9: End-to-End Pipeline Benchmarking
- [ ] `viroforge benchmark pipeline`: Comprehensive subcommand
- [ ] Combines all 8 modules
- [ ] Multi-page HTML report

### Additional Features
- [ ] Multi-pipeline comparison
- [ ] Batch benchmarking (`viroforge benchmark-batch`)
- [ ] Additional parsers (MMseqs2, Kaiju, MetaPhlAn)
- [ ] Interactive Plotly visualizations
- [ ] PDF report generation

---

## Testing Checklist

### Unit Tests (Target: 10+)
- [ ] `tests/test_benchmarking_utils.py` (3 tests)
- [ ] `tests/test_qc_metrics.py` (2 tests)
- [ ] `tests/test_assembly_metrics.py` (3 tests)
- [ ] `tests/test_taxonomic_metrics.py` (3 tests)
- [ ] `tests/test_abundance_metrics.py` (2 tests)
- [ ] `tests/test_completeness_metrics.py` (2 tests)
- [ ] `tests/test_parsers.py` (4 tests - one per parser)
- [ ] `tests/test_visualizations.py` (2 tests)

### Integration Tests (Target: 5+)
- [ ] `tests/test_qc_benchmark_integration.py`
- [ ] `tests/test_assembly_benchmark_integration.py`
- [ ] `tests/test_taxonomy_benchmark_integration.py`
- [ ] `tests/test_completeness_benchmark_integration.py`
- [ ] `tests/test_end_to_end.py`

### Test Fixtures
- [ ] `tests/fixtures/test_metadata.json`: Small test dataset metadata
- [ ] `tests/fixtures/test_kraken2_report.txt`: Example Kraken2 output
- [ ] `tests/fixtures/test_centrifuge_output.txt`: Example Centrifuge output
- [ ] `tests/fixtures/test_contigs.fasta`: Example assembled contigs
- [ ] `tests/fixtures/test_coverage.txt`: Example coverage file

---

## Documentation Checklist

### User Documentation
- [ ] `docs/BENCHMARKING_TUTORIAL.md` (Complete tutorial)
- [ ] `docs/BENCHMARKING_API.md` (API reference)
- [ ] Update `README.md` (Add benchmarking examples)
- [ ] Update `docs/QUICKSTART.md` (Add benchmarking workflow)

### Developer Documentation
- [ ] `docs/PHASE13_BENCHMARKING_FRAMEWORK.md` ✅ (Already created)
- [ ] `docs/PHASE13_IMPLEMENTATION_CHECKLIST.md` ✅ (This document)
- [ ] Inline code documentation (docstrings for all functions)
- [ ] Type hints for all public functions

### Project Documentation Updates
- [ ] `ROADMAP.md` ✅ (Already updated with Phase 13)
- [ ] `claude.md` ✅ (Already updated with Phase 13 summary)
- [ ] `.claude/claude.md` ✅ (Already updated with current status)

---

## Success Metrics

### Code Metrics
- [ ] ≥ 1500 lines of new code (modules + CLI)
- [ ] ≥ 500 lines of tests
- [ ] ≥ 80% test coverage on critical paths
- [ ] All type hints present (mypy compliance)

### Functionality Metrics
- [ ] 4+ modules functional (QC, Assembly, Taxonomy, Completeness)
- [ ] 4+ pipeline formats supported (Kraken2, Centrifuge, DIAMOND, Generic)
- [ ] HTML reports with embedded plots
- [ ] JSON export for programmatic access

### Documentation Metrics
- [ ] Tutorial with 5+ examples
- [ ] API documentation for all public functions
- [ ] README updated with benchmarking section

### User Experience Metrics
- [ ] Single command to run each benchmark module
- [ ] Clear error messages when inputs invalid
- [ ] Progress indicators for long operations
- [ ] Beautiful HTML reports (Bootstrap 5)

---

## Dependencies Installation

```bash
# Install ViroForge with benchmarking support
pip install -e ".[benchmark]"

# Or manually install dependencies
pip install matplotlib>=3.5.0 seaborn>=0.11.0 scikit-learn>=1.0.0 plotly>=5.0.0 jinja2>=3.0.0

# External tools (optional, for alignment)
conda install -c bioconda minimap2
```

---

## Quick Start (Post-Implementation)

```bash
# Generate dataset with benchmarking enabled
viroforge generate --preset gut-standard --output benchmarks/gut --enable-benchmarking

# Run your pipeline
your_pipeline benchmarks/gut/fastq/*.fastq --output results/

# Benchmark QC
viroforge benchmark qc \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --raw-reads benchmarks/gut/fastq/*.fastq \
  --cleaned-reads results/qc/*.fastq \
  --output reports/qc.html

# Benchmark assembly
viroforge benchmark assembly \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --contigs results/assembly/contigs.fasta \
  --output reports/assembly.html

# Benchmark taxonomy
viroforge benchmark taxonomy \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --pipeline-output results/kraken2_report.txt \
  --format kraken2 \
  --output reports/taxonomy.html
```

---

## Timeline Summary

| Phase | Duration | Deliverables | Code (est.) |
|-------|----------|--------------|-------------|
| 13A | 2 weeks | Foundation + metadata | ~200 lines |
| 13B | 3 weeks | QC + Assembly | ~800 lines |
| 13C | 3 weeks | Taxonomy + Completeness | ~1000 lines |
| **Total MVP** | **8 weeks** | **4 modules functional** | **~2000 lines** |
| 13D | Future | Advanced modules | ~1500 lines |

---

## Next Steps

1. ✅ Review and approve Phase 13 specifications
2. ✅ Update planning documents (ROADMAP, claude.md) - **DONE**
3. ✅ Create detailed documentation (PHASE13_BENCHMARKING_FRAMEWORK.md) - **DONE**
4. ✅ Create implementation checklist (this document) - **DONE**
5. ⏭️ Begin Phase 13A implementation
6. ⏭️ Implement metadata enhancements
7. ⏭️ Build core benchmarking infrastructure
8. ⏭️ Continue with Phase 13B

---

**Document Version**: 1.0
**Last Updated**: 2025-11-11
**Status**: Ready for implementation
