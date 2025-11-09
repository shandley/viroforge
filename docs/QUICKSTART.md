# ViroForge Quick Start Guide

**Goal**: Generate synthetic virome FASTQ files in 5 minutes

---

## Installation

```bash
# Clone repository
git clone https://github.com/shandley/viroforge.git
cd viroforge

# Install dependencies
pip install biopython numpy pandas insilicoseq
```

---

## Generate Your First Dataset

**Step 1: See what's available**

```bash
python scripts/generate_fastq_dataset.py --list-collections
```

Output:
```
Available Collections:
ID: 9  - Gut Virome - Adult Healthy (Western Diet) - 134 genomes
ID: 10 - Oral Virome - Saliva (Healthy) - 47 genomes
ID: 11 - Skin Virome - Sebaceous Sites (Healthy) - 15 genomes
ID: 12 - Respiratory Virome - Nasopharynx (Healthy) - 41 genomes
ID: 13 - Marine Virome - Coastal Surface Water - 448 genomes
ID: 14 - Soil Virome - Agricultural - 291 genomes
ID: 15 - Freshwater Virome - Lake Surface Water - 200 genomes
ID: 16 - Mouse Gut Virome - Laboratory (C57BL/6) - 22 genomes
```

**Step 2: Generate FASTQ**

```bash
# Generate gut virome dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output my_first_virome \
    --coverage 10
```

Wait 5-15 minutes depending on collection size.

**Step 3: Check output**

```bash
ls -lh my_first_virome/fastq/
```

You'll see:
- `*_R1.fastq` - Forward reads
- `*_R2.fastq` - Reverse reads

**Step 4: View ground truth**

```bash
# View composition table
cat my_first_virome/metadata/*_composition.tsv
```

Shows exact viral composition with abundances.

---

## Common Use Cases

### Quick Test Dataset

```bash
# Small, fast dataset for testing
python scripts/generate_fastq_dataset.py \
    --collection-id 16 \
    --output test_data \
    --coverage 1
```

Runtime: 2-3 minutes

### Full Benchmark Dataset

```bash
# Comprehensive gut virome
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmarks/gut \
    --coverage 10
```

Runtime: 10-15 minutes

### Different VLP Protocols

```bash
# Tangential flow filtration (default - highest purity)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_tff \
    --coverage 10 \
    --vlp-protocol tangential_flow

# Ultracentrifugation (highest recovery)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_ultra \
    --coverage 10 \
    --vlp-protocol ultracentrifugation

# Syringe filter (field-friendly)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_syringe \
    --coverage 10 \
    --vlp-protocol syringe

# Norgen kit (convenient)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_norgen \
    --coverage 10 \
    --vlp-protocol norgen
```

VLP Protocols: `tangential_flow`, `syringe`, `ultracentrifugation`, `norgen`

### Bulk Metagenome (No VLP)

```bash
# Without VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_bulk \
    --coverage 10 \
    --no-vlp
```

### Contamination Levels

```bash
# Clean VLP enrichment (minimal contamination)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_clean \
    --coverage 10 \
    --contamination-level clean

# Realistic VLP enrichment (typical contamination)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_realistic \
    --coverage 10 \
    --contamination-level realistic

# Heavy contamination (sub-optimal VLP)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_heavy \
    --coverage 10 \
    --contamination-level heavy
```

Contamination Levels: `clean` (~0.7%), `realistic` (~7.4%), `heavy` (~27%)

### Library Preparation Amplification ⭐ NEW (Phase 6)

```bash
# No amplification (default - high biomass samples)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_no_amp \
    --coverage 10 \
    --amplification none

# RdAB with 40 cycles (standard virome protocol)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_rdab \
    --coverage 10 \
    --amplification rdab

# RdAB with 30 cycles (moderate bias)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_rdab30 \
    --coverage 10 \
    --amplification rdab-30

# MDA 4 hours (low biomass samples)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_mda \
    --coverage 10 \
    --amplification mda

# MDA overnight (very low biomass)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_mda_long \
    --coverage 10 \
    --amplification mda-long

# Linker-based (minimal bias, modern kits)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_linker \
    --coverage 10 \
    --amplification linker
```

Amplification Methods:
- `none` - No amplification (control, high biomass)
- `rdab` - RdAB 40 cycles (standard, length + GC bias)
- `rdab-30` - RdAB 30 cycles (moderate bias)
- `mda` - MDA 4h (extreme GC bias + stochasticity)
- `mda-long` - MDA 16h (very extreme bias)
- `linker` - Linker-based (minimal GC bias only)

### Different Sequencing Platform

```bash
# MiSeq instead of NovaSeq
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_miseq \
    --coverage 10 \
    --platform miseq
```

Platforms: `novaseq`, `miseq`, `hiseq`

---

## Generate Multiple Datasets

```bash
# All 8 collections at 10x coverage
python scripts/batch_generate_fastq.py \
    --preset benchmark-standard \
    --output benchmarks

# Compare all VLP protocols (Phase 5)
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output vlp_comparison

# Compare amplification methods (Phase 6) ⭐ NEW
python scripts/batch_generate_fastq.py \
    --preset amplification-comparison \
    --output amplification_comparison
```

Runtime: Several hours (benchmark-standard), 20-30 minutes (vlp-protocol-comparison)

Other presets:
- `quick-test` - Fast test (3 collections, 1x coverage)
- `vlp-protocol-comparison` - All 5 VLP protocols (tangential_flow, syringe, ultra, norgen, bulk)
- `vlp-comparison` - VLP vs bulk for 2 collections
- `platform-comparison` - NovaSeq vs MiSeq vs HiSeq
- `coverage-series` - 1x, 5x, 10x, 20x, 50x

---

## Use with Hecatomb

```bash
# Generate benchmark
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmarks/gut \
    --coverage 10

# Run Hecatomb
hecatomb run \
    --reads benchmarks/gut/fastq/*_R{1,2}.fastq \
    --outdir results/gut

# Compare to ground truth
# See: benchmarks/gut/metadata/*_metadata.json
```

---

## Output Files

Every dataset includes:

```
output_dir/
├── fasta/
│   └── collection.fasta              # Reference genomes
├── fastq/
│   ├── collection_R1.fastq          # Forward reads
│   └── collection_R2.fastq          # Reverse reads
└── metadata/
    ├── collection_metadata.json     # Complete ground truth
    └── collection_composition.tsv   # Abundance table
```

**Ground truth metadata** contains:
- Exact viral composition (viral genomes + contaminants)
- Genome IDs, names, taxonomy
- Relative abundances
- Coverage and platform parameters
- **VLP enrichment statistics** (Phase 5):
  - Viral fraction before/after VLP
  - Contamination reduction by type (host DNA, rRNA, bacteria, PhiX)
  - Overall contamination removal percentage
  - VLP protocol used

Use this to validate your pipeline results and understand VLP effects.

---

## Customization Options

### Coverage

```bash
--coverage 10    # 10x mean coverage (default)
--coverage 50    # High coverage
--coverage 1     # Low coverage (fast)
```

### Read Parameters

```bash
--read-length 150    # bp (default)
--insert-size 350    # bp (default)
```

### VLP Enrichment (Phase 5 Enhanced)

```bash
# With VLP (default: tangential flow filtration)
(no flag needed)

# Specify VLP protocol
--vlp-protocol tangential_flow    # Highest purity (91% contamination reduction)
--vlp-protocol ultracentrifugation # Highest recovery (90% viral recovery)
--vlp-protocol syringe             # Field-friendly (60% recovery, 86% reduction)
--vlp-protocol norgen              # Column-based kit (70% recovery, 87% reduction)

# Without VLP (bulk metagenome)
--no-vlp

# Specify contamination level
--contamination-level clean        # ~0.7% initial contamination
--contamination-level realistic    # ~7.4% initial (default)
--contamination-level heavy        # ~27% initial
```

See [VLP Protocol Comparison Tutorial](VLP_PROTOCOL_COMPARISON_TUTORIAL.md) for detailed protocol comparison.

### Reproducibility

```bash
--seed 42    # Use same random seed for reproducible results
--seed 123   # Different seed = different read sampling
```

---

## Troubleshooting

**InSilicoSeq not found**

```bash
pip install insilicoseq
# or
conda install -c bioconda insilicoseq
```

**Database not found**

Check database exists:
```bash
ls -lh viroforge/data/viral_genomes.db
```

Should show ~100MB file. If missing, Phase 3 database construction wasn't completed.

**Out of memory**

Use lower coverage:
```bash
--coverage 5
```

Or limit read count:
```bash
--n-reads 100000
```

**Slow generation**

Normal for large collections. Marine (448 genomes) and Soil (291 genomes) take longest.

Faster alternatives:
- Use smaller collections (Skin: 15, Mouse Gut: 22)
- Lower coverage (`--coverage 1`)
- Test with dry run first (`--dry-run`)

---

## Collection Details

| ID | Environment | Genomes | Best For |
|----|-------------|---------|----------|
| 9 | Human Gut | 134 | General benchmarking |
| 10 | Human Oral | 47 | Oral microbiome studies |
| 11 | Human Skin | 15 | Quick tests (small) |
| 12 | Human Respiratory | 41 | Respiratory virome |
| 13 | Marine | 448 | Environmental viromics |
| 14 | Soil | 291 | Soil ecology |
| 15 | Freshwater | 200 | Aquatic viromics |
| 16 | Mouse Gut | 22 | Quick tests (smallest) |

See [Collection Implementation Guide](COLLECTION_IMPLEMENTATION_GUIDE.md) for detailed composition information.

---

## Next Steps

### Explore Your Data

```bash
# Count reads
wc -l output/fastq/*_R1.fastq
# Divide by 4 for read count

# Check composition
less output/metadata/*_composition.tsv
```

### Run Through Your Pipeline

```bash
your_pipeline output/fastq/*_R{1,2}.fastq
```

### Compare Results to Ground Truth

```python
import json
with open('output/metadata/*_metadata.json') as f:
    truth = json.load(f)
# Compare to pipeline output
```

---

## Advanced Usage: Python API

For users who need more control, ViroForge also provides a Python API for custom workflows:

```python
from viroforge.core.community import create_body_site_profile
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles

# Create custom composition
community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
# Apply VLP enrichment
vlp = standard_vlp()
vlp.apply(composition)
# Apply amplification
amplification = rdab_40_cycles()
amplification.apply(composition)
```

See [USER_GUIDE.md](USER_GUIDE.md) for complete Python API documentation.

---

## Documentation

**Getting Started**:
- **This guide** - Quick start in 5 minutes
- [VLP Protocol Comparison Tutorial](VLP_PROTOCOL_COMPARISON_TUTORIAL.md) - Compare all VLP protocols (Phase 5)
- [FASTQ Generation Guide](PHASE4_FASTQ_GENERATION.md) - Complete technical documentation
- [Scripts README](../scripts/README_FASTQ_GENERATION.md) - All command-line options

**Advanced Topics**:
- [User Guide](USER_GUIDE.md) - Comprehensive Python API documentation
- [Collection Implementation](COLLECTION_IMPLEMENTATION_GUIDE.md) - How collections were curated
- [VLP Integration Guide](VLP_CONTAMINATION_INTEGRATION.md) - VLP enrichment technical details
- [Phase 5 Validation Report](PHASE5_TASK3_VALIDATION_REPORT.md) - Literature validation results

**Database**:
- 14,423 RefSeq viral genomes
- ICTV taxonomy integration (53.9% coverage)
- 8 curated body site collections
- Literature-validated compositions

---

## Support

**Questions**: scott.handley@wustl.edu

**Issues**: https://github.com/shandley/viroforge/issues

**Full Documentation**: See `docs/` directory

---

**Happy Forging!**

*Generated datasets are ready for immediate use in benchmarking studies.*
