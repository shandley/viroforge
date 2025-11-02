# FASTQ Generation Scripts

**Purpose**: Generate realistic synthetic virome FASTQ datasets from curated body site collections

---

## Overview

These scripts provide the Phase 4 database-driven FASTQ generation workflow. They load viral genomes from curated collections in the SQLite database and generate realistic paired-end FASTQ files using InSilicoSeq with platform-specific error models.

**Key Scripts**:
1. `generate_fastq_dataset.py` - Generate FASTQ from a single collection
2. `batch_generate_fastq.py` - Generate multiple datasets with presets

---

## Prerequisites

### Dependencies

```bash
# Python packages
pip install biopython numpy pandas

# InSilicoSeq (for FASTQ generation)
conda install -c bioconda insilicoseq
# or
pip install InSilicoSeq
```

### Database

Requires ViroForge viral genome database populated from Phase 3:
- Location: `viroforge/data/viral_genomes.db`
- 14,423 RefSeq viral genomes
- 8 curated body site collections
- ICTV taxonomy annotations

---

## Script 1: generate_fastq_dataset.py

### Purpose

Generate FASTQ files from a single body site collection with configurable parameters.

### Basic Usage

```bash
# List available collections
python scripts/generate_fastq_dataset.py --list-collections

# Generate gut virome at 10x coverage
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --coverage 10 \
    --platform novaseq
```

### Command-Line Options

#### Required

- `--collection-id INT` - Collection ID to generate (see `--list-collections`)
- `--output PATH` - Output directory for generated files

#### Coverage Options

- `--coverage FLOAT` - Mean coverage depth (default: 10.0)
- `--n-reads INT` - Number of reads (overrides `--coverage`)

#### Read Parameters

- `--read-length INT` - Read length in bp (default: 150)
- `--insert-size INT` - Insert size for paired-end (default: 350)

#### Platform Selection

- `--platform {novaseq,miseq,hiseq}` - Sequencing platform (default: novaseq)

Platform-specific error models:
- `novaseq` - NovaSeq 6000 (higher quality, patterned flow cell)
- `miseq` - MiSeq (cluster-based, no polyG tails)
- `hiseq` - HiSeq 2500 (older platform)

#### VLP Enrichment

- `--no-vlp` - Skip VLP enrichment simulation (generate bulk metagenome)
- `--vlp-efficiency FLOAT` - VLP nuclease efficiency (default: 0.95)

#### Other Options

- `--seed INT` - Random seed for reproducibility (default: 42)
- `--dry-run` - Show what would be generated without running ISS
- `--database PATH` - Path to database (default: viroforge/data/viral_genomes.db)

### Examples

**1. Basic Generation**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_10x \
    --coverage 10
```

**2. High Coverage with MiSeq**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_50x_miseq \
    --coverage 50 \
    --platform miseq
```

**3. Specific Read Count**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 10 \
    --output oral_100k_reads \
    --n-reads 100000
```

**4. Bulk Metagenome (No VLP)**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_bulk \
    --coverage 10 \
    --no-vlp
```

**5. Custom Read Parameters**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 14 \
    --output soil_custom \
    --coverage 5 \
    --read-length 250 \
    --insert-size 500
```

**6. Dry Run (Preview)**

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 16 \
    --output test \
    --coverage 1 \
    --dry-run
```

### Output Structure

```
output_dir/
├── fasta/
│   └── collection_name.fasta          # Reference genomes with abundances
├── fastq/
│   ├── collection_name_R1.fastq       # Forward reads
│   └── collection_name_R2.fastq       # Reverse reads
└── metadata/
    ├── collection_name_metadata.json  # Complete ground truth
    ├── collection_name_composition.tsv # Abundance table
    └── collection_name_abundances.txt  # ISS abundance file
```

### Ground Truth Metadata

**metadata.json** contains:
- Generation timestamp and version
- Collection information
- Configuration parameters
- Complete genome list with:
  - Genome ID, name, length, GC content
  - Taxonomy (family, genus, species)
  - Relative abundance (VLP-adjusted)
  - Original abundance (pre-VLP)

**composition.tsv** provides tab-separated table for easy analysis

### Available Collections

| ID | Name | Genomes | Environment |
|----|------|---------|-------------|
| 9 | Gut Virome - Adult Healthy (Western Diet) | 134 | Human gut |
| 10 | Oral Virome - Saliva (Healthy) | 47 | Human oral cavity |
| 11 | Skin Virome - Sebaceous Sites (Healthy) | 15 | Human skin |
| 12 | Respiratory Virome - Nasopharynx (Healthy) | 41 | Human respiratory |
| 13 | Marine Virome - Coastal Surface Water | 448 | Marine |
| 14 | Soil Virome - Agricultural | 291 | Soil |
| 15 | Freshwater Virome - Lake Surface Water | 200 | Freshwater |
| 16 | Mouse Gut Virome - Laboratory (C57BL/6) | 22 | Mouse gut |

See [Collection Implementation Guide](../docs/COLLECTION_IMPLEMENTATION_GUIDE.md) for detailed composition information.

---

## Script 2: batch_generate_fastq.py

### Purpose

Automate generation of multiple datasets with predefined configurations (presets).

### Basic Usage

```bash
# Quick test with small collections
python scripts/batch_generate_fastq.py \
    --preset quick-test \
    --output data/test_datasets

# Full benchmark suite
python scripts/batch_generate_fastq.py \
    --preset benchmark-standard \
    --output data/benchmark_datasets
```

### Command-Line Options

- `--preset {quick-test,benchmark-standard,vlp-comparison,platform-comparison,coverage-series}` - Use preset configuration
- `--config PATH` - Path to custom YAML configuration file
- `--output PATH` - Output directory for all datasets (required)
- `--dry-run` - Show what would be generated without running

### Available Presets

#### 1. quick-test

**Purpose**: Rapid testing with small datasets

**Configuration**:
- Collections: Mouse gut (22), Skin (15), Oral (47)
- Coverage: 1x
- Platform: MiSeq
- VLP: Enabled

**Use Case**: Quick validation, testing new features

**Runtime**: ~5-10 minutes

**Example**:
```bash
python scripts/batch_generate_fastq.py \
    --preset quick-test \
    --output test_data
```

#### 2. benchmark-standard

**Purpose**: Comprehensive benchmark suite for all body sites

**Configuration**:
- Collections: All 8 collections
- Coverage: 10x
- Platform: NovaSeq
- VLP: Enabled

**Use Case**: Complete pipeline benchmarking across all environments

**Runtime**: Several hours

**Example**:
```bash
python scripts/batch_generate_fastq.py \
    --preset benchmark-standard \
    --output benchmarks
```

#### 3. vlp-comparison

**Purpose**: Compare VLP enrichment vs bulk metagenome

**Configuration**:
- Collections: Gut (9), Marine (13)
- Coverage: 10x
- Platform: NovaSeq
- Generates: Both VLP and non-VLP versions

**Use Case**: Study VLP enrichment effects on viral recovery

**Runtime**: ~2-4 hours

**Example**:
```bash
python scripts/batch_generate_fastq.py \
    --preset vlp-comparison \
    --output vlp_study
```

#### 4. platform-comparison

**Purpose**: Cross-platform reproducibility testing

**Configuration**:
- Collections: Gut (9)
- Coverage: 10x
- Platforms: NovaSeq, MiSeq, HiSeq
- VLP: Enabled

**Use Case**: Assess platform-specific artifacts and assembly quality

**Runtime**: ~1-2 hours

**Example**:
```bash
python scripts/batch_generate_fastq.py \
    --preset platform-comparison \
    --output platform_study
```

#### 5. coverage-series

**Purpose**: Coverage depth sensitivity analysis

**Configuration**:
- Collections: Gut (9)
- Coverages: 1x, 5x, 10x, 20x, 50x
- Platform: NovaSeq
- VLP: Enabled

**Use Case**: Study coverage requirements for assembly/detection

**Runtime**: ~2-3 hours

**Example**:
```bash
python scripts/batch_generate_fastq.py \
    --preset coverage-series \
    --output coverage_study
```

### Custom Configuration Files

Create YAML file with custom dataset specifications:

```yaml
datasets:
  - collection_id: 9
    name: gut_high_coverage
    coverage: 50
    platform: novaseq
    vlp: true
    seed: 42

  - collection_id: 13
    name: marine_low_coverage
    coverage: 1
    platform: miseq
    vlp: false
    seed: 123
```

Use with:
```bash
python scripts/batch_generate_fastq.py \
    --config my_config.yaml \
    --output custom_datasets
```

### Batch Output

**Directory Structure**:
```
output_dir/
├── dataset_1_name/
│   ├── fasta/
│   ├── fastq/
│   └── metadata/
├── dataset_2_name/
│   ├── fasta/
│   ├── fastq/
│   └── metadata/
└── batch_generation_summary.json  # Summary of all generations
```

**Summary File** contains:
- Batch start/end times
- List of all datasets generated
- Success/failure status for each
- Commands used
- Summary statistics

---

## Workflow Integration

### With Hecatomb

```bash
# Generate benchmark dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmarks/gut \
    --coverage 10

# Run Hecatomb
hecatomb run \
    --reads benchmarks/gut/fastq/*_R{1,2}.fastq \
    --outdir results/gut_benchmark

# Compare results to ground truth
# Use benchmarks/gut/metadata/gut_virome_metadata.json
```

### With Custom Pipelines

```python
import json
import pandas as pd

# Load ground truth
with open('output/metadata/collection_metadata.json') as f:
    truth = json.load(f)

# Extract composition
composition = pd.DataFrame(truth['genomes'])

# Compare to pipeline results
# ... your analysis code
```

---

## Troubleshooting

### InSilicoSeq Not Found

**Error**: `iss: command not found`

**Solution**:
```bash
conda install -c bioconda insilicoseq
# or
pip install InSilicoSeq
```

### Database Not Found

**Error**: `FileNotFoundError: Database not found`

**Solution**: Ensure Phase 3 database was populated:
```bash
ls -lh viroforge/data/viral_genomes.db
# Should show ~100+ MB file
```

### Collection Not Found

**Error**: `Collection X not found`

**Solution**: List available collections:
```bash
python scripts/generate_fastq_dataset.py --list-collections
```

### Memory Issues

**Error**: Out of memory during large collection generation

**Solution**: Reduce coverage or split into smaller batches:
```bash
# Use lower coverage
--coverage 5

# Or generate subsets
--n-reads 100000
```

### ISS Generation Slow

**Issue**: Large collections (marine, soil) take hours

**Explanation**: Normal - generating millions of realistic reads with error models is computationally intensive

**Optimization**: Use `--mode basic` in ISS (already default in script)

---

## Technical Details

### Coverage Calculation

Reads needed for target coverage:
```
n_reads = (total_genome_length * coverage) / (2 * read_length)
```

Factor of 2 accounts for paired-end reads covering same fragments.

### VLP Enrichment Simulation

Abundances are adjusted to simulate VLP enrichment:
```python
if apply_vlp:
    abundance = abundance * (1.0 + np.random.normal(0.2, 0.05))
```

Then renormalized to sum to 1.0. This increases viral abundance ~20% with variation, simulating contamination removal.

### Abundance Assignment

Collections use structured random tiered distributions (NOT metagenomic-derived):

| Tier | Range | Purpose |
|------|-------|---------|
| Dominant | 10-30% | Major community members |
| Common | 1-10% | Prevalent species |
| Moderate | 0.1-1% | Detectable species |
| Rare | 0.01-0.1% | Low abundance |
| Very Rare | <0.01% | Marginal detection |

See [Collection Implementation Guide](../docs/COLLECTION_IMPLEMENTATION_GUIDE.md) for detailed methodology.

---

## Related Documentation

- [FASTQ Generation Guide](../docs/PHASE4_FASTQ_GENERATION.md) - Complete Phase 4 documentation
- [Collection Implementation Guide](../docs/COLLECTION_IMPLEMENTATION_GUIDE.md) - How collections were curated
- [Database Design](../docs/GENOME_DATABASE_DESIGN.md) - Database schema
- [Exploration Tools](README_EXPLORATION_TOOLS.md) - Database exploration utilities

---

## Support

**Issues**: Open GitHub issue at https://github.com/shandley/viroforge/issues

**Documentation**: See `docs/` directory

**Questions**: scott.handley@wustl.edu
