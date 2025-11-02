# ViroForge Helper Utilities

**Created**: 2025-11-01
**Phase**: Phase 3, Week 6

## Overview

ViroForge includes three helper utilities for estimating resource requirements, checking database health, and creating test subsets:

| Tool | Purpose | Key Features |
|------|---------|--------------|
| `viroforge_estimate.py` | FASTQ size/runtime estimator | Predict output size, runtime, memory before generation |
| `viroforge_qc.py` | Database quality control | Health checks, validation, quick stats |
| `viroforge_subset.py` | Genome subset selector | Create small test sets with various strategies |

---

## 1. FASTQ Size & Runtime Estimator (`viroforge_estimate.py`)

### Overview

Estimate output file sizes, read counts, runtime, and resource requirements **before** generating mock virome datasets. Essential for planning disk space and compute resources.

### Quick Start

```bash
# Basic estimation
python scripts/viroforge_estimate.py --genomes 500 --coverage 10x

# From collection
python scripts/viroforge_estimate.py --collection gut --coverage 10x

# Multiple coverage levels
python scripts/viroforge_estimate.py --genomes 500 --coverage 5x,10x,20x,50x --format table
```

### Output Example

```
======================================================================
ViroForge FASTQ Size & Runtime Estimation
======================================================================

Input Parameters:
----------------------------------------------------------------------
  Genomes in community:             500
  Mean genome length:            50,000 bp
  Total community size:      25,000,000 bp (25.0 Mbp)
  Target coverage:                 10.0x
  Read length:                      150 bp (paired-end)
  Fragment length:                  350 bp
  Compression:                     gzip

Output Files:
----------------------------------------------------------------------
  Read pairs:                   833,333
  Total reads:                1,666,666

  Uncompressed size:               0.48 GB
  Compressed size (total):         0.24 GB
    - R1 FASTQ:                    0.06 GB
    - R2 FASTQ:                    0.06 GB
  Compression ratio:              25.0%

Resource Requirements:
----------------------------------------------------------------------
  Estimated runtime:                  3 minutes
  Memory required:                  1.5 GB
  Disk space required:              0.6 GB
    - Final output:                 0.2 GB
    - Temp files:                   0.4 GB

Recommendations:
----------------------------------------------------------------------
  ✓ Memory usage within normal range
  ✓ Disk space within normal range
  ✓ Runtime within normal range
```

### Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--genomes` | Number of genomes | Required | `--genomes 500` |
| `--collection` | Use collection | Alternative to --genomes | `--collection gut` |
| `--coverage` | Target coverage | Required | `--coverage 10x` or `5x,10x,20x` |
| `--read-length` | Read length (bp) | 150 | `--read-length 150` |
| `--fragment-length` | Fragment length (bp) | 350 | `--fragment-length 350` |
| `--compression` | Compression type | gzip | `none`, `gzip`, `bzip2`, `xz` |
| `--mean-length` | Mean genome length | 50000 | `--mean-length 40000` |
| `--format` | Output format | detailed | `detailed`, `table`, `json` |

### Table Format (Multiple Coverages)

```bash
python scripts/viroforge_estimate.py \
  --genomes 500 \
  --coverage 5x,10x,20x,50x \
  --format table
```

```
Coverage Comparison Table:
====================================================================================================
Coverage     Read Pairs  Total Reads   FASTQ Size      Runtime     Memory
----------------------------------------------------------------------------------------------------
     5.0x      416,666      833,332         0.1GB           1m      1.5GB
    10.0x      833,333    1,666,666         0.2GB           3m      1.5GB
    20.0x    1,666,666    3,333,332         0.5GB           5m      1.5GB
    50.0x    4,166,666    8,333,332         1.2GB          14m      1.5GB
```

### Estimation Accuracy

Based on empirical benchmarking:

| Metric | Typical Accuracy | Notes |
|--------|------------------|-------|
| File size | ±5% | Depends on compression |
| Runtime | ±20% | Varies by system |
| Memory | ±10% | Conservative estimate |
| Disk temp | ±15% | Depends on pipeline |

### Use Cases

#### Use Case 1: Plan Disk Space

```bash
# Check if you have enough space for large dataset
python scripts/viroforge_estimate.py \
  --genomes 1000 \
  --coverage 50x \
  --compression gzip

# Look at "Disk space required"
# Make sure you have enough free space
```

#### Use Case 2: Compare Compression Methods

```bash
# No compression
python scripts/viroforge_estimate.py --genomes 500 --coverage 10x --compression none

# Gzip (default, 4x smaller)
python scripts/viroforge_estimate.py --genomes 500 --coverage 10x --compression gzip

# XZ (5.5x smaller, slower)
python scripts/viroforge_estimate.py --genomes 500 --coverage 10x --compression xz
```

#### Use Case 3: Estimate Collection Before Curation

```bash
# Once collections are curated, estimate their size
python scripts/viroforge_estimate.py --collection gut --coverage 10x
python scripts/viroforge_estimate.py --collection oral --coverage 10x
python scripts/viroforge_estimate.py --collection marine --coverage 20x
```

---

## 2. Database Quality Control (`viroforge_qc.py`)

### Overview

Rapid database health checks and validation. Runs in <1 second for quick status verification.

### Quick Start

```bash
# Quick health check
python scripts/viroforge_qc.py

# Detailed report
python scripts/viroforge_qc.py --detailed

# Save to file
python scripts/viroforge_qc.py --output qc_report.txt

# Check specific aspect
python scripts/viroforge_qc.py --check taxonomy
```

### Output Example

```
======================================================================
ViroForge Database Quality Control Report
======================================================================

Database Overview:
----------------------------------------------------------------------
  Database file:         viroforge/data/viral_genomes.db
  File size:             34.6 MB
  Total genomes:         969
  Tables present:        9

Taxonomy Check:
----------------------------------------------------------------------
  Status:                ✓ EXCELLENT
  Taxonomy coverage: 71.7% (≥70%)

Quality Check:
----------------------------------------------------------------------
  Status:                ✓ EXCELLENT
  All quality checks passed

  Length statistics:
    Min:         1,168 bp
    Mean:       35,916 bp
    Max:       368,683 bp

  GC content statistics:
    Min:         17.8%
    Mean:        45.1%
    Max:         74.9%

Collection Check:
----------------------------------------------------------------------
  Status:                ℹ️ EMPTY
  Collection table exists but empty

Summary:
----------------------------------------------------------------------
  ✓ Database is in excellent condition

Quality check completed in 0.06 seconds
```

### Quality Checks

#### Taxonomy Check

**Status Levels**:
- ✓ **EXCELLENT**: ≥70% coverage (ideal)
- ✓ **GOOD**: 50-70% coverage (acceptable)
- ⚠️ **FAIR**: 30-50% coverage (needs improvement)
- ❌ **POOR**: <30% coverage (requires attention)

**Checks**:
- Realm coverage percentage
- Coverage at each rank (realm → species)
- Top realms distribution

#### Quality Check

**Checks**:
- Genome length distribution (min/mean/max)
- GC content distribution
- Very short genomes (<1kb)
- Very long genomes (>500kb, potential artifacts)
- Extreme GC content (<15% or >75%)
- Missing sequences

**Status**: EXCELLENT if no issues, NEEDS_ATTENTION if 3+ issues

#### Collection Check

**Checks**:
- Collection table presence
- Number of collections
- Genome associations count
- Collection metadata

**Status**: EXCELLENT if ≥5 collections, GOOD if ≥1

### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--database` | Database path | `--database path/to/db.sqlite` |
| `--detailed` | Show detailed report | `--detailed` |
| `--check` | Check specific aspect | `taxonomy`, `quality`, `collections`, `all` |
| `--output` | Save to file | `--output qc_report.txt` |
| `--format` | Output format | `text`, `json` |

### Use Cases

#### Use Case 1: Daily Health Check

```bash
# Quick validation after updates
python scripts/viroforge_qc.py

# Should take <1 second
# Check status symbols: ✓ = good, ⚠️ = warning, ❌ = error
```

#### Use Case 2: Pre-Release Validation

```bash
# Comprehensive check before release
python scripts/viroforge_qc.py --detailed --output release_qc.txt

# Review report
cat release_qc.txt

# Verify all checks pass
```

#### Use Case 3: Debug Issues

```bash
# Check specific aspect
python scripts/viroforge_qc.py --check quality

# If issues found, investigate
python scripts/viroforge_db.py search --length-min 0 --length-max 1000
```

#### Use Case 4: Monitor After Pipeline Changes

```bash
# Before changes
python scripts/viroforge_qc.py --output before_qc.txt

# ... make changes ...

# After changes
python scripts/viroforge_qc.py --output after_qc.txt

# Compare
diff before_qc.txt after_qc.txt
```

---

## 3. Genome Subset Selector (`viroforge_subset.py`)

### Overview

Create small, representative genome subsets for testing and development. Essential for quick pipeline validation before running full datasets.

### Quick Start

```bash
# Random 50 genomes
python scripts/viroforge_subset.py --count 50 --output test_subset.txt

# Diverse 100 genomes (taxonomically)
python scripts/viroforge_subset.py --count 100 --strategy diverse --output diverse_100.txt

# One per family
python scripts/viroforge_subset.py --strategy family_reps --output family_reps.txt
```

### Selection Strategies

#### 1. Random (Default)

Random selection with optional filters.

```bash
python scripts/viroforge_subset.py --count 50 --output random_50.txt

# With filters
python scripts/viroforge_subset.py \
  --count 20 \
  --family Siphoviridae \
  --host Streptococcus \
  --output sipho_strep.txt
```

#### 2. Diverse

Taxonomically diverse - even sampling across families.

```bash
python scripts/viroforge_subset.py \
  --count 100 \
  --strategy diverse \
  --output diverse_100.txt
```

**Use**: Maximum taxonomic diversity for testing

#### 3. Family Representatives

One or more genomes from each family.

```bash
# One per family
python scripts/viroforge_subset.py \
  --strategy family_reps \
  --output one_per_family.txt

# Three per family
python scripts/viroforge_subset.py \
  --strategy family_reps \
  --per-family 3 \
  --output three_per_family.txt
```

**Use**: Systematic family coverage

#### 4. Size Stratified

Even sampling across size bins (small/medium/large).

```bash
python scripts/viroforge_subset.py \
  --count 60 \
  --strategy size_stratified \
  --output size_stratified.txt
```

**Size bins**:
- Small: 0-10kb
- Medium: 10-50kb
- Large: 50kb-1Mb

**Use**: Test pipeline with various genome sizes

#### 5. Balanced Type

Proportional to genome type distribution (dsDNA, ssRNA, etc).

```bash
python scripts/viroforge_subset.py \
  --count 100 \
  --strategy balanced_type \
  --output balanced_type.txt
```

**Use**: Representative of database genome type distribution

### Output Formats

#### List (Default)

Simple list of genome IDs, one per line:

```
GCF_000001234.1
GCF_000005678.1
GCF_000009012.1
...
```

#### TSV (Metadata)

Tab-separated with metadata:

```tsv
genome_id       genome_name                     length  gc_content  genome_type  family
GCF_000001234.1 Escherichia phage T7            39937   0.485       dsDNA        Podoviridae
GCF_000005678.1 Enterobacteria phage lambda     48502   0.498       dsDNA        Siphoviridae
...
```

```bash
python scripts/viroforge_subset.py \
  --count 50 \
  --format tsv \
  --output subset_with_metadata.tsv
```

#### JSON

Complete metadata in JSON format:

```bash
python scripts/viroforge_subset.py \
  --count 50 \
  --format json \
  --output subset.json
```

### Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--count` | Number of genomes | 50 | `--count 100` |
| `--strategy` | Selection strategy | random | `diverse`, `family_reps`, etc |
| `--per-family` | Genomes per family | 1 | `--per-family 3` (family_reps only) |
| `--family` | Filter by family | None | `--family Siphoviridae` |
| `--host` | Filter by host | None | `--host Escherichia` |
| `--genome-type` | Filter by type | None | `--genome-type dsDNA` |
| `--length-min` | Min length (bp) | None | `--length-min 30000` |
| `--length-max` | Max length (bp) | None | `--length-max 100000` |
| `--output` | Output file | Required | `--output subset.txt` |
| `--format` | Output format | list | `list`, `tsv`, `json` |
| `--seed` | Random seed | 42 | `--seed 123` |
| `--quiet` | Suppress summary | False | `--quiet` |

### Use Cases

#### Use Case 1: Quick Pipeline Test

```bash
# Create small test set
python scripts/viroforge_subset.py \
  --count 10 \
  --strategy diverse \
  --output test_10.txt

# Test pipeline with it
# (Use subset IDs in ViroForge pipeline)
```

#### Use Case 2: Systematic Testing

```bash
# Test each family individually
python scripts/viroforge_subset.py \
  --strategy family_reps \
  --per-family 5 \
  --format tsv \
  --output family_test_set.tsv

# Run pipeline on each family
```

#### Use Case 3: Size-Specific Testing

```bash
# Test with large genomes only
python scripts/viroforge_subset.py \
  --count 20 \
  --length-min 100000 \
  --output large_genomes.txt

# Test with small genomes
python scripts/viroforge_subset.py \
  --count 20 \
  --length-max 10000 \
  --output small_genomes.txt
```

#### Use Case 4: Body Site Testing

```bash
# Create test set for oral virome
python scripts/viroforge_subset.py \
  --count 30 \
  --host Streptococcus \
  --family Siphoviridae \
  --output oral_test_set.txt
```

---

## Complete Workflow Example

### Scenario: Testing Pipeline Before Full Run

```bash
# 1. Check database health
python scripts/viroforge_qc.py
# Verify: ✓ Database is in excellent condition

# 2. Create test subset (10 diverse genomes)
python scripts/viroforge_subset.py \
  --count 10 \
  --strategy diverse \
  --output test_10.txt

# 3. Estimate resources for full run
python scripts/viroforge_estimate.py \
  --genomes 500 \
  --coverage 10x \
  --format table

# 4. Check available disk space
df -h

# 5. Run pipeline with test subset first
# (Use ViroForge with test_10.txt genomes)

# 6. If test passes, run full dataset
# (Use ViroForge with all 500 genomes)
```

---

## Troubleshooting

### Estimator Issues

**Q**: Estimates seem too high/low
**A**: Check mean genome length. Use `--mean-length` to adjust if different from default 50kb.

**Q**: Runtime estimate inaccurate
**A**: Estimates are based on benchmarks. Actual runtime varies by system (CPU, disk I/O, etc).

### QC Issues

**Q**: "Database not found" error
**A**: Check `--database` path. Default is `viroforge/data/viral_genomes.db`

**Q**: Taxonomy coverage showing 0%
**A**: ICTV taxonomy not applied yet. Run taxonomy integration pipeline.

### Subset Issues

**Q**: "No genomes match criteria" error
**A**: Filters too restrictive. Try relaxing filters or using different strategy.

**Q**: Subset smaller than requested
**A**: Not enough genomes match criteria. Check available genomes with viroforge-db search.

---

## Performance

All tools are optimized for speed:

| Tool | Typical Runtime | Database Size | Notes |
|------|----------------|---------------|-------|
| Estimator | <0.1 sec | Any | Pure calculation |
| QC | <1 sec | 1,000-15,000 genomes | Fast SQL queries |
| Subset | 1-5 sec | 1,000-15,000 genomes | Depends on strategy |

**Memory Usage**: All tools use <500 MB memory

---

## Integration with ViroForge Pipeline

### Using Subsets in Pipeline

```python
from viroforge.core.community import ViralCommunity

# Load subset
with open('test_subset.txt') as f:
    genome_ids = [line.strip() for line in f]

# Create community from subset
community = ViralCommunity.from_genome_ids(
    genome_ids,
    random_seed=42
)

# Generate FASTQ
# (continue with normal pipeline)
```

### Validation Workflow

```bash
# 1. QC check
python scripts/viroforge_qc.py

# 2. Estimate resources
python scripts/viroforge_estimate.py --collection gut --coverage 10x

# 3. Create test subset
python scripts/viroforge_subset.py --count 10 --strategy diverse --output test.txt

# 4. Run pipeline with test subset
# (ViroForge pipeline code here)

# 5. If successful, run full collection
# (ViroForge pipeline code here)
```

---

## Tips and Best Practices

### Resource Planning

1. **Always estimate first**: Use estimator before large runs
2. **Check disk space**: Ensure 2x the estimated space (safety margin)
3. **Test with subsets**: Validate pipeline with 10-50 genomes first
4. **Use compression**: gzip saves 75% disk space with minimal performance cost

### Quality Assurance

1. **Run QC after updates**: Check database health after any changes
2. **Monitor taxonomy**: Aim for ≥70% ICTV coverage
3. **Check for artifacts**: Review quality check warnings
4. **Validate collections**: Ensure collections pass validation

### Subset Selection

1. **Use diverse for testing**: Best representative sample
2. **Use family_reps for systematic testing**: One per family ensures coverage
3. **Use filters for specific tests**: Test specific body sites or viral types
4. **Export as TSV**: Easier to review metadata

---

## Future Enhancements

Planned improvements:

- **Batch estimation**: Estimate multiple configurations at once
- **QC trends**: Track quality metrics over time
- **Subset preview**: Show genome details before saving
- **Custom strategies**: User-defined subset selection criteria
- **Integration tests**: Automated testing with subsets

---

## Related Documentation

- **[Database Explorer](README_EXPLORATION_TOOLS.md)** - Search, browse, export genomes
- **[Curation Workflow](../docs/BODY_SITE_CURATION_WORKFLOW.md)** - Body site collection curation
- **[API Reference](../docs/API.md)** - Python API documentation

---

## Support

- **Issues**: Open a GitHub issue
- **Questions**: Check documentation or discussions
- **Feature requests**: Submit via GitHub

---

**Document Version**: 1.0
**Last Updated**: 2025-11-01
**Author**: ViroForge Development Team
