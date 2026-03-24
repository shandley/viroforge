# Phase 13: Benchmarking & Validation Framework

**Author**: ViroForge Development Team
**Date**: 2025-11-11
**Status**: In Planning
**Timeline**: 6-8 weeks for MVP
**Priority**: VERY HIGH

---

## Executive Summary

**Transform ViroForge from a data generator to a complete validation platform.**

### The Problem

ViroForge generates world-class synthetic virome data with perfect ground truth, but provides **zero tools** to help users benchmark their pipelines against that ground truth.

**Current user workflow** (broken):
1. ✅ Generate ViroForge dataset with ground truth
2. ✅ Run pipeline → get taxonomic classifications
3. ❌ Manually parse pipeline output (every tool has different format)
4. ❌ Manually match pipeline taxa to ground truth genomes
5. ❌ Manually calculate accuracy metrics
6. ❌ Manually create comparison plots
7. ❌ Manually write benchmark report

**Ideal workflow** (with benchmarking framework):
1. ✅ Generate ViroForge dataset
2. ✅ Run pipeline
3. ✅ `viroforge benchmark --ground-truth metadata.json --pipeline-output results.txt`
4. ✅ Get comprehensive HTML report with metrics, plots, error analysis

### The Solution

**Modular benchmarking framework** matching the virome analysis workflow:

```
viroforge benchmark qc          # Contamination removal validation
viroforge benchmark assembly    # Genome recovery & completeness
viroforge benchmark taxonomy    # Classification accuracy
viroforge benchmark completeness # Recovery across coverage ranges
viroforge benchmark binning     # vMAG reconstruction
viroforge benchmark annotation  # Gene calling + function
viroforge benchmark host        # Virus-host prediction
viroforge benchmark discovery   # Novel virus detection
viroforge benchmark pipeline    # Comprehensive end-to-end
```

### Why This Matters

#### Virome Analysis ≠ Bacterial Metagenome Analysis

**Bacterial microbiome tools** (CAMI/OPAL):
- Read-based taxonomy is often sufficient
- Assembly is nice-to-have
- Contamination less critical (high biomass)

**Virome analysis** (ViroForge benchmarking):
- **Assembly is essential** (most viruses not in databases)
- **Completeness matters** (viral genomes are small, 5-100kb)
- **QC is make-or-break** (low biomass, 90%+ contamination)
- **Novel discovery is primary goal** (not just classification)

### Impact

**For Pipeline Developers**:
- ✅ Rigorous validation before publication
- ✅ Identify weak spots (low-abundance taxa, specific environments)
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

**For the Viromics Field**:
- ✅ Standardized pipeline benchmarking
- ✅ Reproducible method comparisons
- ✅ Better tools through rigorous testing
- ✅ Increased confidence in virome studies

---

## Architecture

### Integration with ViroForge

**Recommendation**: Integrated with optional dependencies

```python
# setup.py
extras_require={
    "benchmark": [
        "matplotlib>=3.5.0",      # Plotting
        "seaborn>=0.11.0",        # Statistical visualizations
        "scikit-learn>=1.0.0",    # Metrics calculations
        "plotly>=5.0.0",          # Interactive plots
        "jinja2>=3.0.0",          # HTML report templates
    ],
    "web": ["flask>=2.0.0"],
    "all": ["matplotlib", "seaborn", "scikit-learn", "plotly", "jinja2", "flask"]
}
```

**Installation**:
```bash
# Core only
pip install viroforge

# With benchmarking
pip install viroforge[benchmark]

# Everything
pip install viroforge[all]
```

### Module Structure

```
viroforge/
├── benchmarking/              # NEW MODULE (optional install)
│   ├── __init__.py
│   ├── parsers/              # Parse pipeline outputs
│   │   ├── __init__.py
│   │   ├── kraken2.py        # Kraken2/Bracken format
│   │   ├── centrifuge.py     # Centrifuge format
│   │   ├── diamond.py        # DIAMOND BLAST output
│   │   ├── metaphlan.py      # MetaPhlAn profile
│   │   ├── kaiju.py          # Kaiju output
│   │   ├── mmseqs2.py        # MMseqs2 output
│   │   └── generic.py        # Generic TSV (genome_id, abundance, taxonomy)
│   ├── metrics/              # Calculate performance metrics
│   │   ├── __init__.py
│   │   ├── taxonomic.py      # Classification accuracy (TP/FP/TN/FN, precision, recall, F1)
│   │   ├── abundance.py      # Quantification accuracy (Pearson, Spearman, RMSE, Bray-Curtis)
│   │   ├── assembly.py       # Assembly quality (completeness, identity, chimeras)
│   │   ├── completeness.py   # Genome recovery (by coverage, length, GC)
│   │   ├── contamination.py  # QC validation (contamination removal rates)
│   │   └── diversity.py      # Ecological metrics (Shannon, Simpson, richness)
│   ├── visualizations/       # Generate diagnostic plots
│   │   ├── __init__.py
│   │   ├── composition.py    # Stacked bars, pie charts
│   │   ├── scatter.py        # Abundance correlations, residual plots
│   │   ├── confusion.py      # Confusion matrix heatmaps
│   │   ├── assembly.py       # Coverage plots, completeness heatmaps
│   │   ├── sankey.py         # Read flow diagrams
│   │   └── barplots.py       # Error breakdowns, metric comparisons
│   ├── reports/              # Report generation
│   │   ├── __init__.py
│   │   ├── html_report.py    # Publication-quality HTML reports
│   │   ├── json_export.py    # Machine-readable JSON export
│   │   └── templates/        # Jinja2 HTML templates
│   │       ├── base.html
│   │       ├── qc_report.html
│   │       ├── assembly_report.html
│   │       ├── taxonomy_report.html
│   │       └── pipeline_report.html
│   └── utils.py              # Taxonomy matching, alignment utilities
├── cli/
│   └── benchmark.py          # NEW: viroforge benchmark commands
```

---

## Benchmarking Modules

### Module 1: QC Benchmarking (CRITICAL for viromes)

**Purpose**: Validate contamination removal and quality filtering

#### Why Critical for Viromes

Low biomass virome samples are 90%+ contamination:
- Host DNA: 2-30%
- Bacterial DNA: 5-20%
- rRNA: 5-90% (RNA viromes)
- PhiX spike-in: 0.5-2%

**Bad QC = failed virome study**. Must validate:
1. Contamination removed correctly
2. Viral reads retained (not over-filtered)

#### ViroForge Ground Truth

```json
{
  "contamination_manifest": {
    "host_dna": {
      "sequences": [
        {"id": "HOST_DNA_001", "name": "Human_chr1_fragment", "length": 5000, "abundance": 0.02}
      ],
      "total_abundance": 0.035
    },
    "rrna": {
      "sequences": [
        {"id": "RRNA_16S_001", "name": "E.coli_16S_rRNA", "abundance": 0.008}
      ],
      "total_abundance": 0.013
    },
    "bacterial": {...},
    "phix": {...}
  }
}
```

#### Metrics

```python
qc_metrics = {
    'phix_removal_rate': 0.98,      # 98% PhiX removed
    'host_removal_rate': 0.92,      # 92% host DNA removed
    'bacterial_removal_rate': 0.85, # 85% bacterial DNA removed
    'rrna_removal_rate': 0.95,      # 95% rRNA removed (RNA viromes)

    'viral_retention': 0.98,        # 98% viral reads kept (CRITICAL)
    'false_positive_rate': 0.02,    # 2% viral reads wrongly removed

    'by_contaminant_type': {
        'host_dna': {'precision': 0.95, 'recall': 0.92, 'f1': 0.935},
        'rrna': {'precision': 0.97, 'recall': 0.95, 'f1': 0.960},
        'bacterial': {'precision': 0.90, 'recall': 0.85, 'f1': 0.875},
        'phix': {'precision': 0.99, 'recall': 0.98, 'f1': 0.985}
    }
}
```

#### CLI

```bash
viroforge benchmark qc \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --raw-reads data/gut_standard/fastq/raw_*.fastq \
  --cleaned-reads results/qc/cleaned_*.fastq \
  --qc-stats results/qc/fastp.json \
  --output reports/qc_benchmark.html
```

#### Visualizations

1. **Contamination Removal Bar Chart**
```
Removal Rate                Host    Bacterial  rRNA    PhiX
(%)         100 |            ▓▓▓▓    ▓▓▓▓      ▓▓▓▓▓   ▓▓▓▓▓
             75 |            ▓▓▓▓    ▓▓▓▓      ▓▓▓▓▓   ▓▓▓▓▓
             50 |            ▓▓▓▓    ▓▓▓       ▓▓▓▓▓   ▓▓▓▓▓
             25 |            ▓▓▓▓    ▓▓▓       ▓▓▓▓▓   ▓▓▓▓▓
              0 |____________▓▓▓▓____▓▓▓_______▓▓▓▓▓___▓▓▓▓▓
                             92%     85%       95%     98%
                           Target: 90%+
```

2. **Viral Retention Scatter**
```
Observed                     Perfect retention line
Viral     ●                 /
Reads   ●   ●             /
      ●       ●         /
    ●           ●     /
  ●               ● /
                  /●
                /
              /
            /
          Expected Viral Reads

Retention Rate: 98% ✓
False Positive Rate: 2%
```

---

### Module 2: Assembly Benchmarking (MOST CRITICAL for viromes)

**Purpose**: Compare assembled contigs to true viral genomes

#### Why Most Critical for Viromes

- Viral genomes are small (5-100kb) → complete recovery realistic
- Most viruses not in databases → assembly-first approach essential
- Novel virus discovery relies on assembly quality
- Chimeras are problematic → need detection

#### ViroForge Ground Truth

```json
{
  "sequences": [
    {
      "genome_id": "GCF_015160975.1",
      "genome_name": "CrAssphage cr118_1",
      "length": 92600,
      "relative_abundance": 0.229,
      "expected_coverage": 68.7,
      "expected_completeness": 0.999
    }
  ]
}
```

Plus: Full genome sequences in FASTA for alignment.

#### Metrics

**A. Genome Recovery**
```python
assembly_metrics = {
    'contigs_total': 450,
    'contigs_viral': 380,           # Actually viral (vs contaminant)
    'contigs_matched': 320,         # Match known genomes

    'genomes_recovered': {
        'complete': 45,              # ≥95% genome recovered
        'high_quality': 78,          # ≥75% genome recovered
        'partial': 105,              # ≥50% genome recovered
        'fragmented': 42,            # <50% genome recovered
        'missing': 64                # Not recovered at all
    },

    'genome_recovery_rate': 0.76,   # 320/420 expected genomes
    'mean_completeness': 0.68,      # Average % genome recovered
    'mean_identity': 0.998,         # Sequence identity

    'assembly_contiguity': {
        'n50': 15400,
        'l50': 12,
        'longest_contig': 94500
    }
}
```

**B. Chimera Detection**
```python
chimera_analysis = {
    'chimeric_contigs': 8,          # Contigs matching multiple genomes
    'chimera_rate': 0.021,          # 2.1% of contigs
    'example_chimeras': [
        {
            'contig_id': 'NODE_45',
            'length': 12500,
            'segments': [
                {'genome': 'Crass_phage_1', 'pos': '1-8000', 'identity': 0.999},
                {'genome': 'Microvirus_2', 'pos': '8001-12500', 'identity': 0.998}
            ],
            'verdict': 'CHIMERA'
        }
    ]
}
```

**C. Coverage Bias**
```python
coverage_bias = {
    'mean_coverage': 45.2,
    'coverage_std': 12.3,
    'coverage_cv': 0.27,            # Coefficient of variation
    'gc_bias_correlation': -0.15,   # Negative = low GC undercovered
    'length_bias_correlation': 0.08,
    'regions_zero_coverage': 245
}
```

#### CLI

```bash
viroforge benchmark assembly \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --contigs results/assembly/contigs.fasta \
  --coverage results/assembly/coverage.txt \
  --output reports/assembly_benchmark.html
```

#### Visualizations

1. **Assembly Quality Heatmap**
```
Genome         Completeness  Identity  Coverage  Verdict
─────────────────────────────────────────────────────────
Crass_1        ████████████  99.9%     45x      ✓ Complete
Microvirus_2   ████████      99.8%     38x      ✓ High-quality
Coliphage_3    ████          99.7%     22x      ~ Partial
Bocavirus_1    █             98.2%     8x       ✗ Fragmented
Rotavirus_A    <empty>       N/A       2x       ✗ Missing

Summary: 45 complete, 78 high-quality, 105 partial, 64 missing
```

2. **Coverage Uniformity Plot**
```
Coverage                   Expected (horizontal line)
   ^                ──────────────────────────────────
   |        ╱╲    ╱  ╲         ╱╲
   |       ╱  ╲  ╱    ╲       ╱  ╲
   |      ╱    ╲╱      ╲     ╱    ╲
   |─────╱              ╲───╱      ╲─────────────────
   |    ╱                ╲ ╱        ╲
   |   ╱                  ╱          ╲
   +───────────────────────────────────────────────→
       Genome position

GC bias: -0.15 (low GC regions undercovered)
```

---

### Module 3: Binning Benchmarking

**Purpose**: Validate vMAG (viral MAG) reconstruction

#### Metrics

```python
binning_metrics = {
    'bins_total': 85,
    'bins_high_quality': 52,        # ≥90% complete, <5% contamination
    'bins_medium_quality': 28,      # ≥50% complete, <10% contamination

    'bin_purity': {
        'pure': 75,                  # One genome per bin
        'mixed': 8,                  # Multiple genomes wrongly binned
        'chimeric': 2                # Fragments from different genomes
    },

    'genome_binning_rate': 0.82,    # 82% of genomes binned
    'precision': 0.88,               # Correct bins / total bins
    'recall': 0.82                   # Binned genomes / total genomes
}
```

#### CLI

```bash
viroforge benchmark binning \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --bins results/binning/bins/*.fasta \
  --bin-stats results/binning/checkv_quality.txt \
  --output reports/binning_benchmark.html
```

---

### Module 4: Taxonomy Benchmarking

**Purpose**: Compare pipeline classifications to ground truth

**Two sub-modes**:
- **Read-based**: Benchmark Kraken2, Centrifuge, MMseqs2, DIAMOND
- **Contig-based**: Benchmark contig taxonomy assignment

#### Metrics

**At Multiple Taxonomic Levels** (species, genus, family):

```python
taxonomic_metrics = {
    'species_level': {
        'true_positives': 120,
        'false_positives': 15,
        'false_negatives': 8,
        'true_negatives': 14000,
        'precision': 0.889,      # TP / (TP + FP)
        'recall': 0.9375,        # TP / (TP + FN)
        'f1_score': 0.9123,      # Harmonic mean
        'accuracy': 0.9984       # (TP + TN) / Total
    },
    'genus_level': {...},
    'family_level': {...}
}
```

**Abundance Accuracy**:
```python
abundance_metrics = {
    'pearson_r': 0.95,         # Linear correlation
    'spearman_rho': 0.93,      # Rank correlation
    'r_squared': 0.90,         # Coefficient of determination
    'mae': 0.012,              # Mean absolute error
    'rmse': 0.018,             # Root mean squared error
    'mape': 8.5,               # Mean absolute percentage error
    'bray_curtis': 0.12        # Dissimilarity (0=identical)
}
```

#### Supported Tools

| Tool | Format | Support |
|------|--------|---------|
| Kraken2/Bracken | kraken2 report | ✅ Phase 13C |
| Centrifuge | centrifuge output | ✅ Phase 13C |
| DIAMOND | BLAST tabular | ✅ Phase 13C |
| MMseqs2 | taxonomy TSV | ✅ Phase 13C |
| MetaPhlAn | profile | ⏸ Phase 13D |
| Kaiju | output table | ⏸ Phase 13D |
| Generic | TSV (genome_id, abundance, taxonomy) | ✅ Phase 13C |

#### CLI

```bash
# Read-based classification
viroforge benchmark taxonomy \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --pipeline-output results/kraken2_report.txt \
  --format kraken2 \
  --mode read-based \
  --output reports/taxonomy_benchmark.html

# Contig-based classification
viroforge benchmark taxonomy \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --contig-taxonomy results/taxonomy/contig_tax.txt \
  --format generic \
  --mode contig-based \
  --output reports/taxonomy_benchmark.html
```

#### Visualizations

1. **Confusion Matrix** (Family-level)
```
                    Predicted Family
                Micro  Sino  Inoviridae  Unknown
Actual   Micro   120    2      0          8
Family   Sino     1    95      0          4
         Ino      0     0     82          3
         Other    3     2      1         45
```

2. **Abundance Scatter Plot**
```
Observed                  ●  Perfect correlation
Abundance              ●     ●
(log10)           ●  ●     ●
                ●     ●  ●
              ●    ●  ●
            ●   ●
          ●  ●
         ●●
        ●
     Ground Truth Abundance (log10)

Pearson r = 0.95, p < 0.001
R² = 0.90
```

---

### Module 5: Completeness Benchmarking

**Purpose**: Measure genome recovery across coverage ranges

#### Metrics

```python
completeness_by_coverage = {
    '< 5x': {'mean_completeness': 0.15, 'n_genomes': 45},
    '5-10x': {'mean_completeness': 0.45, 'n_genomes': 38},
    '10-20x': {'mean_completeness': 0.75, 'n_genomes': 52},
    '20-50x': {'mean_completeness': 0.92, 'n_genomes': 68},
    '> 50x': {'mean_completeness': 0.98, 'n_genomes': 85}
}

completeness_by_length = {
    '< 10kb': {'mean_completeness': 0.85},
    '10-50kb': {'mean_completeness': 0.78},
    '50-100kb': {'mean_completeness': 0.65},
    '> 100kb': {'mean_completeness': 0.52}
}

completeness_by_gc = {
    '< 30% GC': {'mean_completeness': 0.62},
    '30-40% GC': {'mean_completeness': 0.75},
    '40-50% GC': {'mean_completeness': 0.80},
    '> 50% GC': {'mean_completeness': 0.68}
}
```

#### CLI

```bash
viroforge benchmark completeness \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --contigs results/assembly/contigs.fasta \
  --coverage results/assembly/coverage.txt \
  --output reports/completeness_benchmark.html
```

---

### Module 6: Annotation Benchmarking (Phase 13D)

**Purpose**: Validate gene calling and functional annotation

**Requires**: Gene annotation export from RefSeq CDS

#### Metrics

```python
gene_calling = {
    'genes_true': 2450,
    'genes_predicted': 2380,
    'genes_correct': 2100,
    'sensitivity': 0.857,
    'precision': 0.882,
    'f1_score': 0.869
}

functional_annotation = {
    'functions_true': 1850,
    'functions_predicted': 1620,
    'functions_correct': 1420,
    'precision': 0.877,
    'recall': 0.768,
    'f1_score': 0.819
}
```

---

### Module 7: Host Prediction Benchmarking (Phase 13D)

**Purpose**: Validate virus-host linkage predictions

#### Metrics

```python
host_prediction = {
    'predictions_total': 320,
    'predictions_correct': 245,
    'overall_accuracy': 0.766,
    'by_host_type': {
        'bacteriophage': {'accuracy': 0.88, 'n': 180},
        'human_virus': {'accuracy': 0.92, 'n': 85},
        'environmental': {'accuracy': 0.65, 'n': 55}
    }
}
```

---

### Module 8: Novel Discovery Benchmarking (Phase 13D - Advanced)

**Purpose**: Test ability to find viruses NOT in reference databases

#### How It Works

1. ViroForge withholds 20% of genomes from reference DB
2. User builds pipeline with 80% DB
3. Pipeline should still find the 20% "novel" viruses via assembly
4. Compare recovered genomes to withheld truth

#### Metrics

```python
discovery_metrics = {
    'novel_viruses_true': 28,       # Withheld from DB
    'novel_viruses_found': 22,      # Found via assembly
    'discovery_rate': 0.786,        # 22 / 28
    'false_discoveries': 5,         # Claimed novel but known
    'precision': 0.815,
    'recall': 0.786
}
```

#### CLI

```bash
# Step 1: Generate with holdout
viroforge generate --preset gut-standard --holdout 0.2 --holdout-db exclude.fasta

# Step 2: User builds DB excluding holdouts
# Step 3: User runs pipeline
# Step 4: Benchmark discovery
viroforge benchmark discovery \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --holdout-list data/gut_standard/metadata/holdout_genomes.txt \
  --discovered-genomes results/novel/*.fasta \
  --output reports/discovery_benchmark.html
```

---

### Module 9: End-to-End Pipeline Benchmarking

**Purpose**: Comprehensive benchmarking of entire pipeline

**Combines all modules** into one comprehensive report.

#### CLI

```bash
viroforge benchmark pipeline \
  --ground-truth data/gut_standard/metadata/metadata.json \
  --raw-reads data/gut_standard/fastq/raw_*.fastq \
  --cleaned-reads results/qc/*.fastq \
  --contigs results/assembly/contigs.fasta \
  --bins results/binning/bins/*.fasta \
  --taxonomy results/taxonomy/classifications.txt \
  --annotations results/annotation/functions.txt \
  --output reports/comprehensive_benchmark.html
```

#### Report Structure

```
ViroForge Pipeline Benchmark Report
═══════════════════════════════════════════════════════════

Executive Summary
─────────────────────────────────────────────────────────
Overall Score: B+ (83/100)
✓ PASS: QC (92%), Assembly (85%), Taxonomy (88%)
⚠ WARN: Binning (72%), Annotation (68%)
✗ FAIL: Host prediction (45%)

Detailed Metrics
─────────────────────────────────────────────────────────
1. QC Benchmarking ........................... 92% ✓
2. Assembly Benchmarking ..................... 85% ✓
3. Binning Benchmarking ...................... 72% ⚠
4. Taxonomy Benchmarking ..................... 88% ✓
5. Completeness Analysis ..................... 78% ✓
6. Annotation Benchmarking ................... 68% ⚠
7. Host Prediction Benchmarking .............. 45% ✗

Recommendations
─────────────────────────────────────────────────────────
• Improve binning: Consider metaBAT2 instead of CONCOCT
• Host prediction weak: Update host prediction database
• Overall strong performance on core tasks
```

---

## ViroForge Metadata Enhancements

### Current Metadata (v0.10.0)

```json
{
  "generation_info": {...},
  "collection": {...},
  "configuration": {...},
  "enrichment_stats": {...},
  "sequences": [
    {
      "genome_id": "GCF_015160975.1",
      "genome_name": "CrAssphage cr118_1",
      "length": 92600,
      "relative_abundance": 0.229,
      "family": "Suoliviridae",
      "genus": "Besingivirus",
      "species": "Besingivirus coli"
    }
  ]
}
```

### Enhanced Metadata (v0.11.0)

#### Phase 13A: Essential Enhancements

**1. Contamination Manifest** (HIGH priority, easy)

```json
{
  "benchmarking": {
    "version": "1.0",
    "contamination_manifest": {
      "host_dna": {
        "n_sequences": 10,
        "total_abundance": 0.035,
        "sequences": [
          {"id": "HOST_DNA_001", "name": "Human_chr1_fragment", "length": 5000, "abundance": 0.02}
        ]
      },
      "rrna": {
        "n_sequences": 15,
        "total_abundance": 0.013,
        "sequences": [
          {"id": "RRNA_16S_001", "name": "E.coli_16S_rRNA", "length": 1542, "abundance": 0.008}
        ]
      },
      "bacterial": {...},
      "phix": {...}
    }
  }
}
```

**Implementation**: Export `ContaminationProfile` data to metadata.

**2. Expected Coverage** (HIGH priority, trivial)

```json
{
  "benchmarking": {
    "expected_coverage": {
      "mean_coverage": 30.0,
      "by_genome": {
        "GCF_015160975.1": {
          "relative_abundance": 0.229,
          "expected_coverage": 68.7,
          "expected_completeness": 0.999,
          "expected_n_reads": 950000
        }
      }
    }
  }
}
```

**Implementation**: Calculate from abundances + coverage parameter.

**3. Read Manifest** (MEDIUM priority, moderate effort - optional)

```json
{
  "benchmarking": {
    "read_manifest": {
      "enabled": true,
      "path": "metadata/read_manifest.tsv.gz",
      "format": "tsv_gzip",
      "total_reads": 10000000,
      "viral_reads": 9500000,
      "contamination_reads": 500000
    }
  }
}
```

**read_manifest.tsv.gz**:
```tsv
read_id                  genome_id           sequence_type
M00001:1:000:1          GCF_015160975.1     viral
M00001:1:000:2          GCF_015160975.1     viral
M00001:1:000:3          HOST_DNA_001        host_contamination
M00001:1:000:4          RRNA_16S_001        rrna_contamination
```

**Implementation**: Modify InSilicoSeq/PBSIM3 wrapper to track read origins.

#### Phase 13B: Advanced Enhancements

**4. Gene Annotations** (LOW priority for MVP, harder)

```json
{
  "benchmarking": {
    "gene_annotations": {
      "enabled": true,
      "path": "metadata/gene_annotations.gff3",
      "total_genes": 12450,
      "by_genome": {
        "GCF_015160975.1": {
          "n_genes": 85,
          "n_proteins": 82
        }
      }
    }
  }
}
```

**Implementation**: Parse RefSeq GenBank CDS annotations.

**5. Expected Assembly** (LOW priority, expensive - optional)

```json
{
  "benchmarking": {
    "expected_assembly": {
      "enabled": false,
      "path": "metadata/expected_contigs.fasta",
      "n_contigs": 145,
      "n50": 85000
    }
  }
}
```

**Implementation**: Run megahit/metaSPAdes internally, save expected output.

---

## Implementation Roadmap

### Phase 13A: Foundation + Minimal Enhancements (Weeks 1-2)

**ViroForge Enhancements**:
- [ ] Add `--enable-benchmarking` flag to `save_metadata()`
- [ ] Implement `_export_contamination_manifest()` method
- [ ] Implement `_calculate_expected_coverage()` method
- [ ] Add `benchmarking` section to metadata.json schema
- [ ] Update metadata version to 1.1

**Benchmarking Module**:
- [ ] Create `viroforge/benchmarking/` module structure
- [ ] Implement `benchmarking/__init__.py` with core imports
- [ ] Implement `benchmarking/utils.py` (taxonomy matching, sequence alignment)
- [ ] Implement metadata loader (parses enhanced metadata)
- [ ] Create base classes for modules, metrics, visualizations
- [ ] Add unit tests for core utilities

**Deliverables**:
- Enhanced metadata format (v1.1)
- Core benchmarking infrastructure
- Taxonomy matching utilities
- ~200 lines of code

### Phase 13B: QC + Assembly Benchmarking (Weeks 3-5)

**ViroForge Enhancements**:
- [ ] Implement read manifest tracking (optional)
- [ ] Compress and store read manifest efficiently (gzip)

**Benchmarking Module**:
- [ ] Implement Module 1: QC Benchmarking
  - [ ] `benchmarking/metrics/contamination.py`
  - [ ] Read manifest parser (if enabled)
  - [ ] Contamination removal metrics
  - [ ] Viral retention metrics
  - [ ] `benchmarking/visualizations/qc_plots.py`
- [ ] Implement Module 2: Assembly Benchmarking
  - [ ] `benchmarking/metrics/assembly.py`
  - [ ] Contig-to-genome aligner (minimap2 wrapper?)
  - [ ] Genome recovery metrics
  - [ ] Chimera detection
  - [ ] Coverage bias analysis
  - [ ] `benchmarking/visualizations/assembly_plots.py`
- [ ] Create basic HTML report templates
  - [ ] `benchmarking/reports/templates/base.html`
  - [ ] `benchmarking/reports/templates/qc_report.html`
  - [ ] `benchmarking/reports/templates/assembly_report.html`
  - [ ] `benchmarking/reports/html_report.py`
- [ ] Add CLI commands
  - [ ] `viroforge/cli/benchmark.py`
  - [ ] `viroforge benchmark qc`
  - [ ] `viroforge benchmark assembly`
- [ ] Unit tests for QC and assembly modules

**Deliverables**:
- Module 1: QC Benchmarking (functional)
- Module 2: Assembly Benchmarking (functional)
- Basic HTML reports
- CLI: `viroforge benchmark {qc,assembly}`
- ~800 lines of code

### Phase 13C: Taxonomy + Completeness (Weeks 6-8)

**Benchmarking Module**:
- [ ] Implement Module 4: Taxonomy Benchmarking
  - [ ] `benchmarking/parsers/kraken2.py`
  - [ ] `benchmarking/parsers/centrifuge.py`
  - [ ] `benchmarking/parsers/diamond.py`
  - [ ] `benchmarking/parsers/generic.py`
  - [ ] `benchmarking/metrics/taxonomic.py`
  - [ ] `benchmarking/metrics/abundance.py`
  - [ ] Classification metrics (TP/FP/TN/FN, precision, recall, F1)
  - [ ] Abundance correlation metrics (Pearson, Spearman, RMSE)
  - [ ] Error analysis (false positives/negatives by taxon)
  - [ ] `benchmarking/visualizations/taxonomy_plots.py`
- [ ] Implement Module 5: Completeness Analysis
  - [ ] `benchmarking/metrics/completeness.py`
  - [ ] Completeness by coverage range
  - [ ] Completeness by genome length, GC content
  - [ ] Systematic bias detection
  - [ ] `benchmarking/visualizations/completeness_plots.py`
- [ ] Enhanced HTML reports with visualizations
  - [ ] `benchmarking/reports/templates/taxonomy_report.html`
  - [ ] `benchmarking/reports/templates/completeness_report.html`
  - [ ] Matplotlib/Plotly plots embedded in HTML
- [ ] Add CLI commands
  - [ ] `viroforge benchmark taxonomy`
  - [ ] `viroforge benchmark completeness`
- [ ] Integration tests with example datasets
- [ ] Documentation
  - [ ] Benchmarking tutorial (`docs/BENCHMARKING_TUTORIAL.md`)
  - [ ] API documentation

**Deliverables**:
- Module 4: Taxonomy Benchmarking (functional)
- Module 5: Completeness Analysis (functional)
- Comprehensive HTML reports with plots
- CLI: `viroforge benchmark {qc,assembly,taxonomy,completeness}`
- Tutorial documentation
- ~1000 lines of code

**MVP COMPLETE** at end of Phase 13C

### Phase 13D: Advanced Modules (Future - Weeks 9+)

- [ ] Module 3: Binning Benchmarking
- [ ] Module 6: Annotation Benchmarking
- [ ] Module 7: Host Prediction Benchmarking
- [ ] Module 8: Novel Discovery Benchmarking
- [ ] Module 9: End-to-End Pipeline Benchmarking
- [ ] Multi-pipeline comparison
- [ ] Batch benchmarking
- [ ] Additional parsers (MMseqs2, Kaiju, MetaPhlAn)
- [ ] Interactive plots (Plotly)
- [ ] PDF report generation

---

## Example Use Cases

### Use Case 1: Validate Hecatomb Pipeline

```bash
# Generate benchmark dataset
viroforge generate --preset gut-standard --output benchmarks/gut --enable-benchmarking

# Run Hecatomb
hecatomb run \
  --reads benchmarks/gut/fastq/*.fastq \
  --output hecatomb_results/

# Benchmark each component
viroforge benchmark qc \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --raw-reads benchmarks/gut/fastq/*.fastq \
  --cleaned-reads hecatomb_results/QC/*.fastq \
  --output reports/hecatomb_qc.html

viroforge benchmark assembly \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --contigs hecatomb_results/assembly/contigs.fasta \
  --output reports/hecatomb_assembly.html

viroforge benchmark taxonomy \
  --ground-truth benchmarks/gut/metadata/metadata.json \
  --contig-taxonomy hecatomb_results/taxonomy/contig_tax.txt \
  --format generic \
  --output reports/hecatomb_taxonomy.html
```

### Use Case 2: Compare 3 Taxonomic Classifiers

```bash
# Generate benchmark dataset
viroforge generate --preset marine-standard --enable-benchmarking

# Run all 3 pipelines
kraken2 --db viral_db --paired reads_R1.fastq reads_R2.fastq --report kraken2.txt
centrifuge -x viral_db -1 reads_R1.fastq -2 reads_R2.fastq --report-file centrifuge.txt
diamond blastx -d viral_db -q reads.fastq -o diamond.txt

# Compare all 3
viroforge benchmark taxonomy \
  --ground-truth data/marine-standard/metadata/metadata.json \
  --pipelines \
    kraken2:kraken2.txt \
    centrifuge:centrifuge.txt \
    diamond:diamond.txt \
  --output reports/classifier_comparison.html
```

### Use Case 3: Benchmark Across Disease States

```bash
# Generate datasets
viroforge generate --preset gut-healthy --output data/gut_healthy --enable-benchmarking
viroforge generate --collection-id 18 --output data/gut_ibd --enable-benchmarking
viroforge generate --collection-id 19 --output data/gut_hiv --enable-benchmarking

# Run pipeline on all
run_pipeline.sh data/gut_*/fastq/

# Batch benchmark
viroforge benchmark-batch \
  --datasets data/gut_* \
  --pipeline-results results/ \
  --format kraken2 \
  --output reports/disease_state_benchmark/
```

---

## Testing Strategy

### Unit Tests

- [ ] Test taxonomy matching utilities
- [ ] Test contamination metrics calculations
- [ ] Test assembly metrics calculations
- [ ] Test classification metrics (TP/FP/TN/FN)
- [ ] Test abundance correlation metrics
- [ ] Test parsers for each tool format
- [ ] Test visualization generation (no display, just file creation)

### Integration Tests

- [ ] Generate small test dataset with benchmarking enabled
- [ ] Run mock pipeline (create fake results)
- [ ] Run each benchmark module
- [ ] Verify HTML reports generated correctly
- [ ] Verify JSON exports match expected format

### Example Datasets

Create small reference datasets for testing:
- `tests/data/gut_test/` - 10 genome test dataset with ground truth
- `tests/data/mock_results/` - Fake pipeline outputs for testing parsers

---

## Documentation

### User Documentation

- [ ] `docs/BENCHMARKING_TUTORIAL.md` - Complete user tutorial
- [ ] `docs/BENCHMARKING_API.md` - API reference for programmatic use
- [ ] Update `README.md` with benchmarking examples
- [ ] Update `QUICKSTART.md` with benchmarking workflow

### Developer Documentation

- [ ] `docs/PHASE13_IMPLEMENTATION_CHECKLIST.md` - Detailed implementation tasks
- [ ] `docs/BENCHMARKING_MODULE_DESIGN.md` - Architecture and design decisions
- [ ] Inline code documentation (docstrings)
- [ ] Type hints for all functions

---

## Success Criteria

### MVP (Phase 13A-C Complete)

- ✅ Enhanced metadata with contamination manifest, expected coverage
- ✅ Module 1: QC Benchmarking functional
- ✅ Module 2: Assembly Benchmarking functional
- ✅ Module 4: Taxonomy Benchmarking functional (Kraken2, Centrifuge, DIAMOND)
- ✅ Module 5: Completeness Analysis functional
- ✅ HTML reports with metrics + visualizations
- ✅ CLI: `viroforge benchmark {qc,assembly,taxonomy,completeness}`
- ✅ Tutorial documentation
- ✅ 10+ unit tests, 5+ integration tests

### Full Implementation (Phase 13D)

- All 9 benchmarking modules functional
- 10+ pipeline format parsers
- Multi-pipeline comparison
- Batch benchmarking
- Interactive plots
- Publication-ready figures
- Comprehensive documentation

---

## Risks & Mitigation

### Risk 1: Complexity

**Risk**: Benchmarking framework is complex, could take longer than 6-8 weeks

**Mitigation**:
- MVP approach (4 modules first)
- Modular design (one module at a time)
- Skip advanced modules initially (Phase 13D)

### Risk 2: Pipeline Format Variability

**Risk**: Every pipeline has different output format

**Mitigation**:
- Start with most common (Kraken2, Centrifuge, DIAMOND)
- Provide generic TSV format for custom pipelines
- Community can contribute parsers

### Risk 3: Alignment Performance

**Risk**: Aligning contigs to 14,423 genomes might be slow

**Mitigation**:
- Use minimap2 (very fast)
- Only align contigs to genomes in current collection (~100-400)
- Cache alignments for reuse

### Risk 4: User Adoption

**Risk**: Users might not adopt benchmarking tools

**Mitigation**:
- Excellent documentation and tutorials
- Hecatomb integration (dogfood internally)
- Publication with benchmark examples
- Community outreach

---

## Future Enhancements (Beyond Phase 13)

### Phase 14: Advanced Visualizations

- Interactive Plotly dashboards
- Real-time progress monitoring
- Comparison matrices (N pipelines × M collections)

### Phase 15: Benchmarking Presets

- Predefined benchmarking workflows
- "Standard virome pipeline benchmark"
- "Assembly quality benchmark"
- "Novel discovery benchmark"

### Phase 16: Community Benchmark Repository

- Public repository of benchmark results
- Compare your pipeline to published baselines
- Leaderboards (optional, community-driven)

---

## Conclusion

Phase 13 transforms ViroForge from a **data generator** to a **complete validation platform**, enabling the viromics community to rigorously benchmark pipelines with standardized methods.

**This is the missing piece** that completes the "generate → benchmark → publish" workflow.

**Timeline**: 6-8 weeks for MVP (4 core modules)
**Priority**: VERY HIGH
**Impact**: Positions ViroForge as the gold standard for virome pipeline validation

---

## References

### Related Tools

**CAMI/OPAL**: Assessment framework for bacterial metagenomes
- Focus: Bacterial metagenomes, assembly, binning
- Gap: Not virus-specific, no VLP/RNA modeling
- URL: https://github.com/CAMI-challenge

**AMBER**: Assessment of Metagenome BinnERs
- Focus: Genome binning quality
- Gap: Not virome-specific
- URL: https://github.com/CAMI-challenge/AMBER

**ViromeQC**: Virome enrichment assessment
- Focus: QC metric for real samples (no ground truth)
- Gap: Cannot benchmark classification accuracy
- URL: https://github.com/SegataLab/viromeqc

**ViroForge Benchmarking**: Combines all above + virome-specific features

### Literature

**Virome Analysis Workflows**:
- Roux et al. 2016 - Virome analysis workflow best practices
- Emerson et al. 2018 - Host-associated viral ecology framework
- Gregory et al. 2019 - Marine virome analysis

**Pipeline Benchmarking**:
- Sczyrba et al. 2017 - CAMI challenge results
- McIntyre et al. 2017 - Benchmarking taxonomic classifiers
- Meyer et al. 2019 - CAMI2 challenge

**VLP Enrichment & QC**:
- Lim et al. 2020 - VLP protocol comparison
- Thurber et al. 2009 - Viral metagenomics methods
- Reyes et al. 2012 - DNase treatment efficiency

---

**Document Version**: 1.0
**Last Updated**: 2025-11-11
**Next Review**: Start of Phase 13A implementation
