# ViroForge Hybrid Assembly Tutorial

**Version**: 0.9.0
**Date**: 2025-11-10
**Phase**: 11 - Hybrid Assembly Support

---

## Table of Contents

1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Generation Methods](#generation-methods)
4. [Supported Assemblers](#supported-assemblers)
5. [Benchmarking Workflows](#benchmarking-workflows)
6. [Validation](#validation)
7. [Troubleshooting](#troubleshooting)
8. [FAQ](#faq)

---

## Introduction

### What is Hybrid Assembly?

Hybrid assembly combines **short-read accuracy** with **long-read length** to produce superior genome assemblies. It addresses the weaknesses of each technology:

| Technology | Strength | Weakness | Hybrid Benefit |
|------------|----------|----------|----------------|
| Short reads | High accuracy (>99.9%) | Fragmented assemblies | Provides accuracy |
| Long reads | Complete genomes | Higher error rate (5-10%) | Provides contiguity |
| **Hybrid** | **Best of both** | **Requires more data** | **Complete + accurate** |

### Why Hybrid Assembly for Viromes?

**Viral Genomes are Ideal for Hybrid Assembly**:
- Small size (5-300kb) → manageable data requirements
- Often circular → long reads can span entire genomes
- Strain diversity → short reads resolve SNPs accurately
- Structural variants → long reads detect large indels/rearrangements

**ViroForge Enables**:
- Matched short + long datasets with identical compositions
- Benchmarking hybrid assemblers with complete ground truth
- Technology comparison studies
- Cost-effectiveness analysis

---

## Quick Start

### Method 1: Convenience Script (Recommended)

Generate matched datasets in one command:

```bash
# NovaSeq + PacBio HiFi hybrid dataset
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
    --output data/gut_hybrid \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 \
    --depth 15 \
    --seed 42

# Output structure:
# data/gut_hybrid/
# ├── short_reads/       # NovaSeq reads
# ├── long_reads/        # PacBio HiFi reads
# └── hybrid_metadata.json
```

### Method 2: Separate Runs

Generate datasets separately with matching parameters:

```bash
# Step 1: Generate short reads
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_novaseq \
    --platform novaseq \
    --coverage 30 \
    --seed 42

# Step 2: Generate long reads (SAME collection, SAME seed)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_hifi \
    --platform pacbio-hifi \
    --depth 15 \
    --seed 42
```

**Critical**: Use the **same** `--collection-id` and `--seed` for matching compositions!

---

## Generation Methods

### Convenience Script Parameters

```bash
python scripts/generate_hybrid_dataset.py \
    --collection-id ID \              # Collection to use
    --output OUTPUT_DIR \             # Output directory
    --short-platform {novaseq,miseq,hiseq} \
    --long-platform {pacbio-hifi,nanopore} \
    --coverage DEPTH \                # Short-read coverage
    --depth DEPTH \                   # Long-read depth
    --seed SEED                       # Random seed (critical!)
```

### Recommended Hybrid Combinations

| Short | Long | Coverage | Depth | Use Case |
|-------|------|----------|-------|----------|
| NovaSeq | PacBio HiFi | 30x | 15x | High accuracy, complete genomes |
| NovaSeq | Nanopore | 50x | 20x | Ultra-long scaffolds, cost-effective |
| MiSeq | PacBio HiFi | 50x | 15x | Moderate throughput, high accuracy |
| HiSeq | Nanopore | 30x | 25x | Legacy + long reads |

### Platform-Specific Parameters

#### PacBio HiFi Options
```bash
python scripts/generate_hybrid_dataset.py \
    --long-platform pacbio-hifi \
    --pacbio-passes 15 \              # More passes = higher accuracy
    --pacbio-read-length 18000        # Mean read length (10-30kb)
```

#### Nanopore Options
```bash
python scripts/generate_hybrid_dataset.py \
    --long-platform nanopore \
    --ont-chemistry R10.4 \           # R10.4 > R9.4 (lower error)
    --ont-read-length 30000           # Mean read length (10-200kb)
```

---

## Supported Assemblers

### 1. Unicycler (Recommended)

**Best For**: Bacterial and viral genomes, automated hybrid assembly

**Installation**:
```bash
conda install -c bioconda unicycler
```

**Usage**:
```bash
unicycler \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    -l data/gut_hybrid/long_reads/fastq/*.fastq* \
    -o results/unicycler \
    --threads 8
```

**Output**: `results/unicycler/assembly.fasta`

**Pros**:
- Fully automated
- Conservative approach (few misassemblies)
- Polishes with Pilon using short reads
- Good for circular genomes

**Cons**:
- Slower than alternatives
- Conservative = may miss complex regions

### 2. SPAdes (Hybrid Mode)

**Best For**: Metagenomes, mixed communities

**Installation**:
```bash
conda install -c bioconda spades
```

**Usage**:
```bash
# For PacBio HiFi
spades.py --meta \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    --pacbio data/gut_hybrid/long_reads/fastq/*.fastq* \
    -o results/spades \
    --threads 8

# For Nanopore
spades.py --meta \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    --nanopore data/gut_hybrid/long_reads/fastq/*.fastq \
    -o results/spades \
    --threads 8
```

**Output**: `results/spades/scaffolds.fasta`

**Pros**:
- Fast
- Good for metagenomes (`--meta` mode)
- Handles mixed communities well

**Cons**:
- More aggressive = more misassemblies
- Less polishing than Unicycler

### 3. MaSuRCA

**Best For**: Large genomes, complex repeats

**Installation**:
```bash
conda install -c bioconda masurca
```

**Usage**:
```bash
# Create configuration file
masurca -g config.txt

# Edit config.txt:
# PE= pe 180 20  short_R1.fastq short_R2.fastq
# PACBIO=long_reads.fastq
# END

# Run assembly
./assemble.sh
```

**Pros**:
- Excellent repeat resolution
- Handles ultra-long reads well

**Cons**:
- Complex configuration
- Slower
- High memory usage

### 4. HybridSPAdes

**Best For**: Metagenomic hybrid assembly

**Installation**:
```bash
conda install -c bioconda spades
```

**Usage**:
```bash
spades.py --meta \
    -1 short_R1.fastq \
    -2 short_R2.fastq \
    --pacbio long_reads.fastq \
    -o results/hybridspades
```

**Note**: Same as SPAdes `--meta` mode with long reads.

---

## Benchmarking Workflows

### Workflow 1: Assembler Comparison

**Goal**: Compare Unicycler vs SPAdes vs long-read-only assemblies

```bash
# 1. Generate hybrid dataset
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
    --output data/gut_hybrid \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 --depth 15 --seed 42

# 2. Unicycler hybrid
unicycler \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    -l data/gut_hybrid/long_reads/fastq/*.fastq.gz \
    -o results/unicycler

# 3. SPAdes hybrid
spades.py --meta \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    --pacbio data/gut_hybrid/long_reads/fastq/*.fastq.gz \
    -o results/spades

# 4. Long-read only (Flye for comparison)
flye --pacbio-hifi data/gut_hybrid/long_reads/fastq/*.fastq.gz \
    --out-dir results/flye_longonly

# 5. Short-read only (SPAdes for comparison)
spades.py --meta \
    -1 data/gut_hybrid/short_reads/fastq/*_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/*_R2.fastq \
    -o results/spades_shortonly

# 6. Compare all assemblies against ground truth
GROUND_TRUTH=data/gut_hybrid/short_reads/metadata/*_composition.tsv

# Calculate metrics:
# - N50, L50
# - Genome completeness (% of ground truth genomes assembled)
# - Misassembly rate
# - Computational cost (time, memory)
```

**Expected Results**:
- **Hybrid** > Long-only > Short-only (completeness)
- **Short-only** > Hybrid > Long-only (accuracy)
- **Hybrid** = best overall balance

### Workflow 2: Coverage Optimization

**Goal**: Determine optimal short/long read ratios

```bash
# Test different coverage combinations
for SHORT_COV in 20 30 50 100; do
  for LONG_DEPTH in 10 15 20 30; do
    echo "Testing ${SHORT_COV}x short + ${LONG_DEPTH}x long"

    python scripts/generate_hybrid_dataset.py \
        --collection-id 9 \
        --output data/cov_${SHORT_COV}_${LONG_DEPTH} \
        --short-platform novaseq \
        --long-platform pacbio-hifi \
        --coverage $SHORT_COV \
        --depth $LONG_DEPTH \
        --seed 42

    unicycler \
        -1 data/cov_${SHORT_COV}_${LONG_DEPTH}/short_reads/fastq/*_R1.fastq \
        -2 data/cov_${SHORT_COV}_${LONG_DEPTH}/short_reads/fastq/*_R2.fastq \
        -l data/cov_${SHORT_COV}_${LONG_DEPTH}/long_reads/fastq/*.fastq.gz \
        -o results/cov_${SHORT_COV}_${LONG_DEPTH}
  done
done

# Analyze results to find optimal coverage combination
```

### Workflow 3: Technology Comparison

**Goal**: Compare NovaSeq+HiFi vs NovaSeq+Nanopore

```bash
# HiFi hybrid
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
    --output data/gut_novaseq_hifi \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 --depth 15 --seed 42

# Nanopore hybrid
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
    --output data/gut_novaseq_nanopore \
    --short-platform novaseq \
    --long-platform nanopore \
    --coverage 30 --depth 20 --seed 42

# Assemble both with Unicycler
# Compare cost, time, and assembly quality
```

---

## Validation

### Validate Composition Matching

After generating datasets, verify they have matching compositions:

```bash
# Using hybrid directory
python scripts/validate_hybrid_composition.py \
    --hybrid-dir data/gut_hybrid

# Using explicit files
python scripts/validate_hybrid_composition.py \
    --short data/gut_novaseq/metadata/*_composition.tsv \
    --long data/gut_hifi/metadata/*_composition.tsv
```

**Expected Output**:
```
ViroForge Hybrid Composition Validator
================================================================================

Loading compositions...
  Short-read: 33 genomes
  Long-read:  33 genomes

Validating compositions...

================================================================================
✓ VALIDATION PASSED
================================================================================

Genome compositions match!
  • 33 genomes in both datasets
  • Abundances match within tolerance (1e-06)

These datasets are suitable for hybrid assembly.
```

---

## Troubleshooting

### Problem: Compositions don't match

**Cause**: Different `--seed` or `--collection-id` used

**Solution**:
```bash
# Make sure BOTH use same seed and collection
python scripts/generate_fastq_dataset.py --collection-id 9 --seed 42 ...
python scripts/generate_fastq_dataset.py --collection-id 9 --seed 42 ...
```

### Problem: Assembler fails with "not enough data"

**Cause**: Insufficient coverage/depth

**Solution**: Increase coverage:
```bash
--coverage 50 --depth 25  # Higher values
```

### Problem: Assembly is fragmented despite hybrid

**Cause**: Long reads may be too short to span repeats

**Solution**: Use longer reads:
```bash
--pacbio-read-length 20000  # Increase from 15000
--ont-read-length 40000     # Increase from 20000
```

### Problem: Very slow assembly

**Cause**: Hybrid assembly is computationally intensive

**Solutions**:
1. Reduce coverage: `--coverage 20 --depth 10`
2. Use more threads: `--threads 16`
3. Use faster assembler (SPAdes instead of Unicycler)

---

## FAQ

### Q: What coverage/depth should I use for hybrid assembly?

**A**: Recommended starting points:
- **Short reads**: 30-50x
- **Long reads**: 15-25x
- **Minimum**: 20x short + 10x long
- **High quality**: 50x short + 25x long

### Q: Which assembler should I use?

**A**:
- **Unicycler**: Best for single genomes, circular genomes
- **SPAdes**: Best for metagenomes, faster
- **MaSuRCA**: Best for complex repeats, large genomes

### Q: Can I mix different short-read platforms?

**A**: Not recommended. Each platform has different error profiles. Stick to one short-read platform per hybrid dataset.

### Q: Can I use RNA-seq short reads with DNA long reads?

**A**: No. Use `--molecule-type dna` for both, or generate separate RNA+DNA hybrids.

### Q: How do I know if my hybrid assembly is better than short/long-only?

**A**: Compare against ground truth:
```bash
# Ground truth is in metadata/*_composition.tsv
# Compare:
# - N50 (contiguity)
# - # contigs (fewer = better)
# - Genome completeness (% genomes >90% assembled)
# - Misassemblies (check against known genome structures)
```

### Q: What if I already have separate short/long datasets?

**A**: As long as they're from the **same collection** with the **same seed**, they work for hybrid assembly. Validate with:
```bash
python scripts/validate_hybrid_composition.py \
    --short short_reads/metadata/*_composition.tsv \
    --long long_reads/metadata/*_composition.tsv
```

### Q: Can I generate multiple hybrid datasets in parallel?

**A**: Yes! Each with different seeds:
```bash
for SEED in 42 123 456 789; do
  python scripts/generate_hybrid_dataset.py \
      --collection-id 9 \
      --output data/gut_hybrid_seed${SEED} \
      --seed $SEED &
done
wait
```

---

## Next Steps

1. **Try the Quick Start** to generate your first hybrid dataset
2. **Run Unicycler** to see hybrid assembly in action
3. **Compare assemblers** using benchmarking workflows
4. **Optimize coverage** for your specific use case

## References

1. **Unicycler**: Wick RR, et al. Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 2017;13(6):e1005595.

2. **SPAdes**: Bankevich A, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol 2012;19(5):455-477.

3. **MaSuRCA**: Zimin AV, et al. The MaSuRCA genome assembler. Bioinformatics 2013;29(21):2669-2677.

4. **Hybrid Assembly Review**: Wick RR, Holt KE. Benchmarking of long-read assemblers for prokaryote whole genome sequencing. F1000Research 2019;8:2138.

---

## Support

- **Issues**: https://github.com/hecatomb/viroforge/issues
- **Long-Read Tutorial**: `docs/LONGREAD_TUTORIAL.md`
- **Examples**: `scripts/` directory

**Happy hybrid assembling!** 🧬+🧬=🎉
