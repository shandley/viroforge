# ViroForge Long-Read Sequencing Tutorial

<!-- collection-id-scheme-note -->
> **Collection IDs use the 1-20 scheme.** ViroForge renumbered its collections to a
> contiguous 1-20 layout (healthy gut = 1). Example commands in this document may still
> show legacy IDs from the older 9-28 numbering. For the current collection-to-ID map,
> run `viroforge browse` or query the `body_site_collections` table.


**Date**: 2026-09-22
**Phase**: 10 - Long-Read Sequencing Support

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Platform Details](#platform-details)
5. [Configuration Guide](#configuration-guide)
6. [Usage Examples](#usage-examples)
7. [Output Files](#output-files)
8. [Benchmarking Workflows](#benchmarking-workflows)
9. [Troubleshooting](#troubleshooting)
10. [FAQ](#faq)

---

## Introduction

ViroForge supports long-read sequencing simulation for **Oxford Nanopore** using PBSIM3, enabling:

- **Complete Viral Genome Assembly**: Long reads (10-200kb) can span entire viral genomes
- **Structural Variant Detection**: Identify insertions, deletions, rearrangements
- **Technology Comparison**: Compare short-read vs long-read performance
- **Pipeline Validation**: Benchmark assembly and variant calling tools with complete ground truth

### Why Long Reads for Viromics?

| Feature | Short Reads (150-300bp) | Long Reads (10-200kb) |
|---------|------------------------|----------------------|
| Genome Assembly | Fragmented, gaps | Complete, circular |
| Structural Variants | Difficult to detect | Easy to detect |
| Repeat Regions | Collapsed | Resolved |
| Strain Diversity | Mixed, ambiguous | Phased, separated |
| Cost per Gb | Low | Medium |
| Accuracy | >99.9% | ~95% (R10.4) |

---

## Installation

### Prerequisites

ViroForge long-read sequencing requires PBSIM3:

```bash
# Create conda environment with all dependencies
conda create -n viroforge-longread \
    python=3.9 \
    biopython \
    numpy \
    pandas \
    pbsim3

# Activate environment
conda activate viroforge-longread

# Install ViroForge
pip install -e .
```

### Dependency Details

| Tool | Purpose | Installation |
|------|---------|--------------|
| **pbsim3** | Nanopore read simulator | `conda install -c bioconda pbsim3` |

### Verify Installation

```bash
# Check PBSIM3
pbsim --version
# Expected: PBSIM3 v3.0.0 or later
```

---

## Quick Start

### Nanopore: 5-Minute Example

```bash
# Generate soil virome with Nanopore
viroforge generate \
    --collection-id 6 \
    --output data/quickstart_nanopore \
    --depth 15 \
    --platform nanopore

# Output:
#   data/quickstart_nanopore/fastq/soil_virome.fastq
#   data/quickstart_nanopore/metadata/soil_virome_ground_truth.tsv
```

---

## Platform Details

### Oxford Nanopore

**Overview**: Single-molecule sequencing with characteristic homopolymer errors.

**Characteristics**:
- **Accuracy**: ~95% (R10.4 chemistry), improving to 99%+
- **Read Length**: 10kb-2Mb (ultra-long possible)
- **Error Profile**: Homopolymer indels (deletions > insertions)
- **Throughput**: ~10-50 Gb per flow cell (MinION)
- **Cost**: $$ (lower per Gb)

**Best For**:
- Ultra-long reads for spanning repeats
- Structural variant detection
- Real-time sequencing (portable MinION)
- Complete viral genome assembly
- Large viral genome assembly (herpesviruses, poxviruses)
- Phage genomics (large genomes with repeats)
- Rapid outbreak response

---

## Configuration Guide

### Nanopore Parameters

```bash
viroforge generate \
    --platform nanopore \
    --depth 15 \                     # Sequencing depth (15x typical)
    --ont-chemistry R10.4 \          # Chemistry version (R9.4 or R10.4)
    --ont-read-length 20000          # Mean read length in bp
```

#### Parameter Details

| Parameter | Default | Options | Effect |
|-----------|---------|---------|--------|
| `--depth` | 15 | 10-100 | Higher depth = better coverage, larger files |
| `--ont-chemistry` | R10.4 | R9.4, R10.4 | R10.4 has lower error rate (~5% vs ~10%) |
| `--ont-read-length` | 20000 | 5000-200000 | Longer reads = better assembly, more chimeras |

#### Common Configurations

**Standard Sequencing** (R10.4, ~95% accuracy):
```bash
--ont-chemistry R10.4 --ont-read-length 20000
```

**Long Reads** (R10.4, 30-50kb):
```bash
--ont-chemistry R10.4 --ont-read-length 40000
```

**Ultra-Long Reads** (R10.4, >100kb):
```bash
--ont-chemistry R10.4 --ont-read-length 100000
```

**Older Chemistry** (R9.4, ~90% accuracy):
```bash
--ont-chemistry R9.4 --ont-read-length 15000
```

---

## Usage Examples

### Example 1: Gut Virome Assembly Benchmark

**Goal**: Generate Nanopore dataset for evaluating viral genome assemblers.

```bash
viroforge generate \
    --collection-id 1 \
    --output data/gut_nanopore_benchmark \
    --depth 20 \
    --platform nanopore \
    --ont-chemistry R10.4 \
    --ont-read-length 25000 \
    --vlp-protocol tangential_flow \
    --contamination-level realistic \
    --seed 42

# Expected output:
# - 20x coverage of gut virome genomes
# - ~25kb mean read length
# - ~95% accuracy (R10.4)
# - Realistic contamination (bacterial/host DNA)
```

**Use Case**: Benchmark Flye, Canu, and other long-read assemblers.

### Example 2: Structural Variant Detection

**Goal**: Generate ultra-long Nanopore reads for detecting large insertions/deletions.

```bash
viroforge generate \
    --collection-id 6 \
    --output data/soil_nanopore_sv \
    --depth 30 \
    --platform nanopore \
    --ont-chemistry R10.4 \
    --ont-read-length 50000 \
    --vlp-protocol tangential_flow \
    --seed 42

# Expected output:
# - 30x coverage
# - ~50kb mean read length (ultra-long)
# - ~95% accuracy
# - Can span large viral genomes entirely
```

**Use Case**: Benchmark Sniffles, cuteSV, SVIM structural variant callers.

### Example 3: Technology Comparison (Short vs Long)

**Goal**: Compare short-read (Illumina) vs long-read (Nanopore) assembly quality.

```bash
# Generate Illumina NovaSeq dataset
viroforge generate \
    --collection-id 5 \
    --output data/marine_novaseq \
    --coverage 50 \
    --platform novaseq \
    --vlp-protocol tangential_flow \
    --seed 42

# Generate Nanopore dataset (SAME collection, SAME seed)
viroforge generate \
    --collection-id 5 \
    --output data/marine_nanopore \
    --depth 20 \
    --platform nanopore \
    --ont-chemistry R10.4 \
    --vlp-protocol tangential_flow \
    --seed 42
```

**Analysis**:
- Compare N50, L50, genome completeness
- Compare misassembly rates (using ground truth)
- Compare computational time and memory

### Example 4: VLP Enrichment vs Bulk Metagenome

**Goal**: Compare VLP-enriched vs bulk metagenome sequencing with long reads.

```bash
# VLP-enriched (high viral fraction)
viroforge generate \
    --collection-id 1 \
    --output data/gut_vlp_nanopore \
    --depth 15 \
    --platform nanopore \
    --vlp-protocol tangential_flow \
    --contamination-level realistic \
    --seed 42

# Bulk metagenome (low viral fraction)
viroforge generate \
    --collection-id 1 \
    --output data/gut_bulk_nanopore \
    --depth 15 \
    --platform nanopore \
    --no-vlp \
    --contamination-level heavy \
    --seed 42
```

**Expected Difference**:
- VLP: ~85-95% viral reads
- Bulk: ~5-15% viral reads

### Example 5: Multi-Platform Dataset Generation

**Goal**: Generate datasets for all 4 platforms (NovaSeq, MiSeq, HiSeq, Nanopore).

```bash
#!/bin/bash
# Generate multi-platform benchmark suite

COLLECTION=1  # Gut virome
OUTPUT_BASE=data/multiplatform_benchmark
SEED=42

# Short-read platforms
for PLATFORM in novaseq miseq hiseq; do
    viroforge generate \
        --collection-id $COLLECTION \
        --output ${OUTPUT_BASE}/${PLATFORM} \
        --coverage 30 \
        --platform $PLATFORM \
        --vlp-protocol tangential_flow \
        --seed $SEED
done

# Nanopore
viroforge generate \
    --collection-id $COLLECTION \
    --output ${OUTPUT_BASE}/nanopore \
    --depth 20 \
    --platform nanopore \
    --ont-chemistry R10.4 \
    --vlp-protocol tangential_flow \
    --seed $SEED
```

---

## Output Files

### Directory Structure

```
data/my_longread_dataset/
├── fasta/
│   └── collection_name.fasta              # Input genomes (with abundances)
├── fastq/
│   └── collection_name.fastq              # Nanopore reads (uncompressed)
└── metadata/
    ├── collection_name_metadata.json      # Complete dataset metadata
    ├── collection_name_composition.tsv    # Genome composition table
    └── collection_name_ground_truth.tsv   # Read-to-genome ground truth
```

### Ground Truth Format

**File**: `*_ground_truth.tsv`

```tsv
genome_id       genome_type     length  relative_abundance      platform        read_type
NC_001416       viral           5386    0.25                    nanopore        long
NC_001422       viral           5386    0.20                    nanopore        long
NC_007605       viral           48502   0.15                    nanopore        long
contam_001      bacterial       2500000 0.05                    nanopore        long
```

**Columns**:
- `genome_id`: Unique genome identifier
- `genome_type`: `viral` or `contaminant` (bacterial, host, reagent)
- `length`: Genome length in bp
- `relative_abundance`: Fraction of total reads (0-1, sums to 1.0)
- `platform`: Sequencing platform used (`nanopore`)
- `read_type`: `long` (distinguishes from short-read ground truth)

### Metadata JSON

**File**: `*_metadata.json`

```json
{
  "generation_info": {
    "timestamp": "2026-09-22T10:30:00",
    "viroforge_version": "0.20.0",
    "random_seed": 42
  },
  "collection": {
    "id": 1,
    "name": "Human Gut Virome",
    "n_viral_genomes": 25,
    "n_contaminants": 8
  },
  "configuration": {
    "platform": "nanopore",
    "depth": 15.0,
    "vlp_protocol": "tangential_flow",
    "contamination_level": "realistic"
  },
  "enrichment_stats": {
    "viral_fraction": 0.89,
    "contamination_fraction": 0.11
  }
}
```

---

## Benchmarking Workflows

### Workflow 1: Viral Genome Assembly

**Goal**: Evaluate assembler performance on long-read data.

```bash
# 1. Generate Nanopore dataset
viroforge generate \
    --collection-id 1 \
    --output data/assembly_benchmark \
    --depth 20 \
    --platform nanopore \
    --seed 42

# 2. Run assemblers
READS=data/assembly_benchmark/fastq/gut_virome.fastq

# Flye
flye --nano-hq $READS --out-dir results/flye --threads 8

# 3. Evaluate assemblies against ground truth
GROUND_TRUTH=data/assembly_benchmark/metadata/gut_virome_ground_truth.tsv

# Calculate metrics:
# - N50, L50
# - Genome completeness (% of ground truth genomes assembled)
# - Misassembly rate
# - Chimeric contig rate
```

**Metrics to Report**:
- Assembly contiguity: N50, L50, largest contig
- Completeness: % genomes >90% assembled
- Accuracy: Misassemblies per 100kb
- Computational: Time, memory, CPU

### Workflow 2: Structural Variant Detection

**Goal**: Benchmark SV callers on long-read data.

```bash
# 1. Generate Nanopore ultra-long reads
viroforge generate \
    --collection-id 6 \
    --output data/sv_benchmark \
    --depth 30 \
    --platform nanopore \
    --ont-read-length 50000 \
    --seed 42

# 2. Align reads to reference
READS=data/sv_benchmark/fastq/soil_virome.fastq
REF=data/sv_benchmark/fasta/soil_virome.fasta

minimap2 -ax map-ont -t 8 $REF $READS | samtools sort > aligned.bam
samtools index aligned.bam

# 3. Call structural variants
# Sniffles
sniffles -i aligned.bam -v sniffles.vcf

# cuteSV
cuteSV aligned.bam $REF cutesv.vcf ./cutesv_temp

# SVIM
svim alignment ./svim aligned.bam $REF

# 4. Evaluate SV calls against ground truth
# (Since genomes are known, any large indels are false positives)
```

### Workflow 3: Technology Comparison

**Goal**: Compare short-read vs long-read assembly quality.

```bash
# Generate datasets (from Example 3)
# ...

# Assemble short reads (SPAdes)
spades.py --meta \
    -1 data/marine_novaseq/fastq/marine_virome_R1.fastq \
    -2 data/marine_novaseq/fastq/marine_virome_R2.fastq \
    -o results/spades

# Assemble long reads (Flye)
flye --nano-hq \
    data/marine_nanopore/fastq/marine_virome.fastq \
    --out-dir results/flye

# Compare:
# - Short reads: fragmented, many small contigs
# - Long reads: complete genomes, circular
```

---

## Troubleshooting

### Problem: "PBSIM3 (pbsim) not found in PATH"

**Solution**:
```bash
# Install PBSIM3
conda install -c bioconda pbsim3

# Verify installation
pbsim --version
```

### Problem: Very slow generation (hours for high depth)

**Solutions**:
1. **Reduce depth**: Use `--depth 5` for quick tests
2. **Shorter reads**: Use `--ont-read-length 15000` instead of 50000
3. **Smaller collection**: Use a collection with fewer genomes

### Problem: Very large output files (>10 GB)

**Solutions**:
1. **Reduce depth**: Lower `--depth` parameter
2. **Smaller collection**: Use collection with fewer genomes
3. **Compress**: Nanopore output is uncompressed by default
   ```bash
   gzip data/output/fastq/*.fastq
   ```

### Problem: "Abundances do not sum to 1.0"

**Cause**: Numerical precision issue in abundance normalization.

**Solution**: This warning is harmless (automatic renormalization occurs).

---

## FAQ

### Q: Can I combine short-read and long-read datasets?

**A**: Yes! Use the **same collection and seed** to ensure identical genome composition:

```bash
# Short reads
viroforge generate --collection-id 1 --platform novaseq --seed 42 --output data/short

# Long reads (same collection, same seed)
viroforge generate --collection-id 1 --platform nanopore --seed 42 --output data/long
```

Then perform **hybrid assembly** with tools like Unicycler, MaSuRCA, or SPAdes hybrid mode.

### Q: What depth should I use?

**A**:
- **Quick tests**: 5-10x
- **Standard benchmarking**: 15-20x
- **High-quality assembly**: 30-50x
- **Rare variant detection**: 100x+

### Q: How long does generation take?

**A**: Nanopore generation typically takes ~5-20 minutes depending on depth and number of genomes.

### Q: Can I use ViroForge long-read data for training ML models?

**A**: Yes! The complete ground truth enables supervised learning:
- Read classification (viral vs contamination)
- Assembly graph resolution
- Error correction models
- Variant calling models

### Q: Are the error profiles realistic?

**A**: Yes, PBSIM3 uses empirically-derived error models from real Nanopore data. Error rates, homopolymer biases, and quality scores match real sequencing.

### Q: Can I simulate mixed viral strains?

**A**: Yes! ViroForge collections include multiple strains of the same virus (e.g., multiple influenza strains in respiratory virome). The long reads enable **strain phasing** benchmarking.

### Q: What if I need even longer reads (>100kb)?

**A**: For ultra-long Nanopore reads:
```bash
viroforge generate \
    --platform nanopore \
    --ont-read-length 150000 \
    --ont-chemistry R10.4
```

---

## Next Steps

1. **Try the Quick Start example** to familiarize yourself with the workflow
2. **Generate a multi-platform dataset** for your favorite collection
3. **Benchmark your assembly/variant calling pipeline** using ground truth
4. **Publish your findings** using ViroForge as the simulation framework

## References

1. **PBSIM3**: Ono Y, et al. PBSIM3: a simulator for all types of PacBio and ONT long reads. NAR Genomics Bioinformatics 2022;4(4):lqac092.

2. **Nanopore R10.4**: Oxford Nanopore Technologies. R10.4.1 chemistry technical note. 2023.

3. **ViroForge**: [Add publication when available]

---

## Support

- **Issues**: https://github.com/shandley/viroforge/issues
- **Documentation**: `docs/` directory
- **Examples**: `scripts/` directory
