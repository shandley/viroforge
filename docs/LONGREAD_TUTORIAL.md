# ViroForge Long-Read Sequencing Tutorial

**Version**: 0.9.0
**Date**: 2025-11-10
**Phase**: 10 - Long-Read Sequencing Support

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Platform Comparison](#platform-comparison)
5. [Configuration Guide](#configuration-guide)
6. [Usage Examples](#usage-examples)
7. [Output Files](#output-files)
8. [Benchmarking Workflows](#benchmarking-workflows)
9. [Troubleshooting](#troubleshooting)
10. [FAQ](#faq)

---

## Introduction

ViroForge now supports long-read sequencing simulation for **PacBio HiFi** and **Oxford Nanopore** platforms, enabling:

- **Complete Viral Genome Assembly**: Long reads (10-30kb) can span entire small viral genomes
- **Structural Variant Detection**: Identify insertions, deletions, rearrangements
- **Technology Comparison**: Compare short-read vs long-read performance
- **Pipeline Validation**: Benchmark assembly and variant calling tools with complete ground truth

### Why Long Reads for Viromics?

| Feature | Short Reads (150-300bp) | Long Reads (10-30kb) |
|---------|------------------------|----------------------|
| Genome Assembly | Fragmented, gaps | Complete, circular |
| Structural Variants | Difficult to detect | Easy to detect |
| Repeat Regions | Collapsed | Resolved |
| Strain Diversity | Mixed, ambiguous | Phased, separated |
| Cost per Gb | Low | Medium |
| Accuracy | >99.9% | 95-99.9% |

---

## Installation

### Prerequisites

ViroForge long-read sequencing requires additional dependencies:

```bash
# Create conda environment with all dependencies
conda create -n viroforge-longread \
    python=3.9 \
    biopython \
    numpy \
    pandas \
    pbsim3 \
    pbccs \
    samtools

# Activate environment
conda activate viroforge-longread

# Install ViroForge
pip install -e .
```

### Dependency Details

| Tool | Purpose | Platform | Installation |
|------|---------|----------|--------------|
| **pbsim3** | Long-read simulator | Both | `conda install -c bioconda pbsim3` |
| **pbccs** | HiFi consensus calling | PacBio only | `conda install -c bioconda pbccs` |
| **samtools** | SAM/BAM manipulation | PacBio only | `conda install -c bioconda samtools` |

### Verify Installation

```bash
# Check PBSIM3
pbsim --version
# Expected: PBSIM3 v3.0.0 or later

# Check PacBio ccs (for HiFi)
ccs --version
# Expected: ccs 6.4.0 or later

# Check samtools
samtools --version
# Expected: samtools 1.15 or later
```

---

## Quick Start

### PacBio HiFi: 5-Minute Example

```bash
# Generate gut virome with PacBio HiFi
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/quickstart_hifi \
    --depth 10 \
    --platform pacbio-hifi

# Output:
#   data/quickstart_hifi/fastq/gut_virome_hifi.fastq.gz
#   data/quickstart_hifi/metadata/gut_virome_ground_truth.tsv
```

### Nanopore: 5-Minute Example

```bash
# Generate soil virome with Nanopore
python scripts/generate_fastq_dataset.py \
    --collection-id 1 \
    --output data/quickstart_nanopore \
    --depth 15 \
    --platform nanopore

# Output:
#   data/quickstart_nanopore/fastq/soil_virome.fastq
#   data/quickstart_nanopore/metadata/soil_virome_ground_truth.tsv
```

---

## Platform Comparison

### PacBio HiFi (Circular Consensus Sequencing)

**Overview**: Multi-pass sequencing of the same molecule generates >99.9% accuracy reads.

**Characteristics**:
- **Accuracy**: >99.9% (QV20+)
- **Read Length**: 10-25kb typical, 30kb possible
- **Error Profile**: Random errors, no systematic bias
- **Throughput**: ~100-200 Gb per SMRT Cell 8M
- **Cost**: $$$ (higher per Gb)

**Best For**:
- Complete viral genome assembly
- High-accuracy variant calling
- SNP detection
- Projects requiring perfect consensus

**Example Use Cases**:
- Clinical diagnostics (need accuracy)
- Reference genome generation
- Mixed viral populations (need precise SNPs)

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
- Projects needing read length over accuracy

**Example Use Cases**:
- Field sequencing (portability)
- Large viral genome assembly (herpesviruses, poxviruses)
- Phage genomics (large genomes with repeats)
- Rapid outbreak response

---

## Configuration Guide

### PacBio HiFi Parameters

```bash
python scripts/generate_fastq_dataset.py \
    --platform pacbio-hifi \
    --depth 10 \                      # Sequencing depth (10x typical)
    --pacbio-passes 10 \              # CCS passes (10 default, 15-20 for higher accuracy)
    --pacbio-read-length 15000        # Mean read length in bp (10-25kb range)
```

#### Parameter Details

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `--depth` | 10 | 5-50 | Higher depth = more coverage, larger files |
| `--pacbio-passes` | 10 | 3-20 | More passes = higher accuracy, shorter reads |
| `--pacbio-read-length` | 15000 | 5000-30000 | Longer reads = better assembly, lower throughput |

**Tradeoff**: More CCS passes (higher accuracy) → Shorter read lengths (lower throughput)

#### Common Configurations

**Standard Accuracy** (QV20, >99.9%):
```bash
--pacbio-passes 10 --pacbio-read-length 15000
```

**High Accuracy** (QV30, >99.99%):
```bash
--pacbio-passes 15 --pacbio-read-length 12000
```

**Ultra-High Accuracy** (QV40, >99.999%):
```bash
--pacbio-passes 20 --pacbio-read-length 10000
```

**Long Reads** (QV20, but longer):
```bash
--pacbio-passes 8 --pacbio-read-length 20000
```

### Nanopore Parameters

```bash
python scripts/generate_fastq_dataset.py \
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

**Goal**: Generate PacBio HiFi dataset for evaluating viral genome assemblers.

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_hifi_benchmark \
    --depth 20 \
    --platform pacbio-hifi \
    --pacbio-passes 15 \
    --pacbio-read-length 18000 \
    --vlp-protocol tangential_flow \
    --contamination-level realistic \
    --seed 42

# Expected output:
# - 20x coverage of gut virome genomes
# - ~18kb mean read length
# - >99.9% accuracy (15 CCS passes)
# - Realistic contamination (bacterial/host DNA)
```

**Use Case**: Benchmark Flye, Canu, HiCanu, hifiasm-meta assemblers.

### Example 2: Nanopore Structural Variant Detection

**Goal**: Generate ultra-long Nanopore reads for detecting large insertions/deletions.

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 1 \
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

**Goal**: Compare short-read (Illumina) vs long-read (PacBio HiFi) assembly quality.

```bash
# Generate Illumina NovaSeq dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output data/marine_novaseq \
    --coverage 50 \
    --platform novaseq \
    --vlp-protocol tangential_flow \
    --seed 42

# Generate PacBio HiFi dataset (SAME collection, SAME seed)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output data/marine_hifi \
    --depth 20 \
    --platform pacbio-hifi \
    --pacbio-passes 15 \
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
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_vlp_hifi \
    --depth 15 \
    --platform pacbio-hifi \
    --vlp-protocol tangential_flow \
    --contamination-level realistic \
    --seed 42

# Bulk metagenome (low viral fraction)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_bulk_hifi \
    --depth 15 \
    --platform pacbio-hifi \
    --no-vlp \
    --contamination-level heavy \
    --seed 42
```

**Expected Difference**:
- VLP: ~85-95% viral reads
- Bulk: ~5-15% viral reads

### Example 5: Multi-Platform Dataset Generation

**Goal**: Generate datasets for all 5 platforms (NovaSeq, MiSeq, HiSeq, PacBio, Nanopore).

```bash
#!/bin/bash
# Generate multi-platform benchmark suite

COLLECTION=9  # Gut virome
OUTPUT_BASE=data/multiplatform_benchmark
SEED=42

# Short-read platforms
for PLATFORM in novaseq miseq hiseq; do
    python scripts/generate_fastq_dataset.py \
        --collection-id $COLLECTION \
        --output ${OUTPUT_BASE}/${PLATFORM} \
        --coverage 30 \
        --platform $PLATFORM \
        --vlp-protocol tangential_flow \
        --seed $SEED
done

# PacBio HiFi
python scripts/generate_fastq_dataset.py \
    --collection-id $COLLECTION \
    --output ${OUTPUT_BASE}/pacbio_hifi \
    --depth 15 \
    --platform pacbio-hifi \
    --pacbio-passes 15 \
    --vlp-protocol tangential_flow \
    --seed $SEED

# Nanopore
python scripts/generate_fastq_dataset.py \
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
│   ├── collection_name_hifi.fastq.gz      # PacBio HiFi reads (gzipped)
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
NC_001416       viral           5386    0.25                    pacbio-hifi     long
NC_001422       viral           5386    0.20                    pacbio-hifi     long
NC_007605       viral           48502   0.15                    pacbio-hifi     long
contam_001      bacterial       2500000 0.05                    pacbio-hifi     long
```

**Columns**:
- `genome_id`: Unique genome identifier
- `genome_type`: `viral` or `contaminant` (bacterial, host, reagent)
- `length`: Genome length in bp
- `relative_abundance`: Fraction of total reads (0-1, sums to 1.0)
- `platform`: Sequencing platform used (`pacbio-hifi`, `nanopore`)
- `read_type`: `long` (distinguishes from short-read ground truth)

### Metadata JSON

**File**: `*_metadata.json`

```json
{
  "generation_info": {
    "timestamp": "2025-11-10T10:30:00",
    "viroforge_version": "0.9.0",
    "random_seed": 42
  },
  "collection": {
    "id": 9,
    "name": "Human Gut Virome",
    "n_viral_genomes": 25,
    "n_contaminants": 8
  },
  "configuration": {
    "platform": "pacbio-hifi",
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
# 1. Generate PacBio HiFi dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/assembly_benchmark \
    --depth 20 \
    --platform pacbio-hifi \
    --seed 42

# 2. Run assemblers
READS=data/assembly_benchmark/fastq/gut_virome_hifi.fastq.gz

# Flye
flye --pacbio-hifi $READS --out-dir results/flye --threads 8

# HiCanu
canu -p gut_hifi -d results/hicanu \
    genomeSize=50k \
    -pacbio-hifi $READS

# hifiasm-meta
hifiasm_meta -o results/hifiasm/gut_hifi $READS

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
python scripts/generate_fastq_dataset.py \
    --collection-id 1 \
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
flye --pacbio-hifi \
    data/marine_hifi/fastq/marine_virome_hifi.fastq.gz \
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

### Problem: "PacBio ccs not found in PATH"

**Solution**:
```bash
# Install pbccs (PacBio consensus caller)
conda install -c bioconda pbccs

# Verify installation
ccs --version
```

### Problem: "SAMtools not found"

**Solution**:
```bash
# Install samtools
conda install -c bioconda samtools

# Verify installation
samtools --version
```

### Problem: Very slow generation (hours for 10x depth)

**Cause**: PBSIM3 and ccs are computationally intensive, especially for PacBio HiFi.

**Solutions**:
1. **Reduce depth**: Use `--depth 5` for quick tests
2. **Reduce passes**: Use `--pacbio-passes 8` instead of 15
3. **Shorter reads**: Use `--pacbio-read-length 12000` instead of 20000
4. **Use Nanopore**: Faster than PacBio HiFi (no CCS step)

### Problem: Very large output files (>10 GB)

**Cause**: High depth with many genomes produces many reads.

**Solutions**:
1. **Reduce depth**: Lower `--depth` parameter
2. **Smaller collection**: Use collection with fewer genomes
3. **Compress**: PacBio HiFi automatically gzips, Nanopore doesn't
   ```bash
   gzip data/output/fastq/*.fastq
   ```

### Problem: "Abundances do not sum to 1.0"

**Cause**: Numerical precision issue in abundance normalization.

**Solution**: This warning is harmless (automatic renormalization occurs).

### Problem: CCS produces no reads

**Cause**: Insufficient CCS passes (min_passes too high).

**Solution**: Ensure `--pacbio-passes >= 3` and check CLR generation succeeded.

---

## FAQ

### Q: Which platform should I use for my benchmarking study?

**A**: It depends on your goal:

- **Variant calling accuracy**: PacBio HiFi (>99.9% accuracy)
- **Complete genome assembly**: PacBio HiFi or Nanopore (both good)
- **Structural variants**: Nanopore ultra-long reads (>50kb)
- **Cost-effective**: Nanopore (lower cost per Gb)
- **Maximum accuracy**: PacBio HiFi with 15-20 passes

### Q: Can I combine short-read and long-read datasets?

**A**: Yes! Use the **same collection and seed** to ensure identical genome composition:

```bash
# Short reads
python scripts/generate_fastq_dataset.py --collection-id 9 --platform novaseq --seed 42 --output data/short

# Long reads (same collection, same seed)
python scripts/generate_fastq_dataset.py --collection-id 9 --platform pacbio-hifi --seed 42 --output data/long
```

Then perform **hybrid assembly** with tools like Unicycler, MaSuRCA, or SPAdes hybrid mode.

### Q: What depth should I use?

**A**:
- **Quick tests**: 5-10x
- **Standard benchmarking**: 15-20x
- **High-quality assembly**: 30-50x
- **Rare variant detection**: 100x+

### Q: How long does generation take?

**A**:
- **PacBio HiFi**: ~10-60 minutes (depends on depth, passes, genomes)
  - Slower: Higher depth, more passes, more genomes
- **Nanopore**: ~5-20 minutes (faster than PacBio, no CCS step)

### Q: Can I use ViroForge long-read data for training ML models?

**A**: Yes! The complete ground truth enables supervised learning:
- Read classification (viral vs contamination)
- Assembly graph resolution
- Error correction models
- Variant calling models

### Q: Are the error profiles realistic?

**A**: Yes, PBSIM3 uses empirically-derived error models from real PacBio and Nanopore data. Error rates, homopolymer biases, and quality scores match real sequencing.

### Q: Can I simulate mixed viral strains?

**A**: Yes! ViroForge collections include multiple strains of the same virus (e.g., multiple influenza strains in respiratory virome). The long reads enable **strain phasing** benchmarking.

### Q: What if I need even longer reads (>100kb)?

**A**: For ultra-long Nanopore reads:
```bash
python scripts/generate_fastq_dataset.py \
    --platform nanopore \
    --ont-read-length 150000 \
    --ont-chemistry R10.4
```

---

## Next Steps

1. **Try the Quick Start examples** to familiarize yourself with the workflow
2. **Generate a multi-platform dataset** for your favorite collection
3. **Benchmark your assembly/variant calling pipeline** using ground truth
4. **Publish your findings** using ViroForge as the simulation framework

## References

1. **PBSIM3**: Ono Y, et al. PBSIM3: a simulator for all types of PacBio and ONT long reads. NAR Genomics Bioinformatics 2022;4(4):lqac092.

2. **PacBio HiFi**: Wenger AM, et al. Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. Nat Biotechnol 2019;37:1155-1162.

3. **Nanopore R10.4**: Oxford Nanopore Technologies. R10.4.1 chemistry technical note. 2023.

4. **ViroForge**: [Add publication when available]

---

## Support

- **Issues**: https://github.com/hecatomb/viroforge/issues
- **Documentation**: `docs/` directory
- **Examples**: `scripts/` directory

**Happy benchmarking!** 🧬🦠
