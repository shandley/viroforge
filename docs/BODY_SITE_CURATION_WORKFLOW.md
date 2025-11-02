# Body Site Virome Collection Curation Workflow

**Date**: 2025-11-01
**Phase**: Phase 3, Week 6-8
**Status**: Ready for implementation

## Overview

This document describes the complete workflow for curating body site-specific virome collections in ViroForge. The workflow takes literature-validated composition specifications and generates realistic, curated viral genome collections for 8 different body sites.

## Collection Specifications

ViroForge includes specifications for 8 body site collections totaling 2,200 genomes:

| Collection | Size | Body Site | Phage % | Dominant Family | Key Feature |
|------------|------|-----------|---------|-----------------|-------------|
| Gut (adult, healthy, Western) | 500 | Gut | 96% | Crassvirales | crAssphage 30% |
| Oral (saliva, healthy) | 200 | Oral cavity | 92.5% | Siphoviridae | Streptococcus phages |
| Skin (sebaceous, healthy) | 150 | Skin | 96% | Cutibacterium | P. acnes phages 47% |
| Respiratory (nasopharynx) | 200 | Respiratory | 87.5% | Mixed | Cutibacterium prevalent |
| Marine (coastal surface) | 500 | Marine | 97% | Caudoviricetes | Pelagiphages |
| Soil (agricultural) | 300 | Soil | 94% | Caudoviricetes | Actinophages |
| Freshwater (lake surface) | 200 | Freshwater | 95% | Caudoviricetes | Cyanophages |
| Mouse gut (C57BL/6) | 150 | Gut | 90% | Siphoviridae | Lactobacillus phages |

See `docs/BODY_SITE_COLLECTIONS.md` for complete specifications.

## Workflow Overview

```
RefSeq Database (14,568 genomes)
         ↓
   [1. Collection Specification]
         ↓
   [2. Genome Selection]
    (Taxonomic + Host filters)
         ↓
   [3. Abundance Modeling]
    (Log-normal + Power-law)
         ↓
   [4. Database Population]
         ↓
   [5. Validation]
         ↓
   [6. Export]
    (TSV + Summary files)
```

## Prerequisites

1. **Complete RefSeq database** with ICTV taxonomy
   - 14,568 viral genomes downloaded
   - Quality filtered (95-97% pass rate)
   - ICTV taxonomy applied (70-75% coverage)

2. **Database schema** with body site collection tables
   - `body_site_collections` - Collection metadata
   - `collection_genomes` - Genome-collection associations

3. **Curation scripts** (created in Phase 3, Week 6)
   - `scripts/curate_body_site_collections.py` - Main curation script
   - `scripts/validate_body_site_collections.py` - Validation script

## Step-by-Step Workflow

### Step 1: Prepare RefSeq Database

Ensure the RefSeq database is complete with ICTV taxonomy:

```bash
# Check database status
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT COUNT(*) as total_genomes FROM genomes"

# Check ICTV taxonomy coverage
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT COUNT(*) as with_taxonomy FROM taxonomy WHERE realm IS NOT NULL"
```

**Expected Results**:
- Total genomes: ~13,800 (after quality filtering)
- ICTV taxonomy: ~10,200 genomes (74% coverage)

### Step 2: Curate Individual Collection

Curate a single body site collection:

```bash
# Curate gut virome (500 genomes)
python scripts/curate_body_site_collections.py \
  --collection gut \
  --output data/body_site_collections \
  --seed 42

# Output:
# - data/body_site_collections/gut_virome_adult_healthy_western_composition.tsv
# - data/body_site_collections/gut_virome_adult_healthy_western_summary.txt
# - Database tables populated
```

**Script Actions**:
1. Loads collection specification from `COLLECTION_SPECS`
2. For each taxonomic target:
   - Queries database with taxonomic/host filters
   - Randomly selects required number of genomes
3. Generates abundance distribution using tiered model
4. Populates database tables
5. Exports composition TSV and summary

### Step 3: Validate Collection

Run validation checks on curated collection:

```bash
# Validate gut virome
python scripts/validate_body_site_collections.py \
  --collection gut_virome_adult_healthy_western \
  --output data/body_site_collections/validation

# Output:
# - data/body_site_collections/validation/gut_virome_adult_healthy_western_validation.txt
```

**Validation Checks**:
- ✓ Genome count matches specification
- ✓ Abundance sum equals 1.0 (±1%)
- ✓ Taxonomy coverage ≥70%
- ✓ Shannon diversity ≥3.0
- ✓ Realistic tier distribution

### Step 4: Curate All Collections

Curate all 8 collections in batch:

```bash
# Curate all collections
python scripts/curate_body_site_collections.py \
  --all \
  --output data/body_site_collections \
  --seed 42

# Validate all collections
python scripts/validate_body_site_collections.py \
  --all \
  --output data/body_site_collections/validation
```

**Expected Execution Time**:
- Curation: ~10-15 minutes for all 8 collections
- Validation: ~2-3 minutes for all 8 collections

### Step 5: Review Results

Check curation results:

```bash
# Count curated collections
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT collection_id, name, genome_count FROM body_site_collections"

# Total genomes in collections
sqlite3 viroforge/data/viral_genomes.db \
  "SELECT COUNT(*) as total_collection_genomes FROM collection_genomes"

# Check validation summaries
cat data/body_site_collections/validation/*_validation.txt | grep "Quality Checks"
```

## Abundance Model Details

### Tiered Abundance Distribution

Collections use a **5-tier abundance model** to simulate realistic viral community structure:

| Tier | Abundance Range | Example Count (n=500) | % of Total Abundance |
|------|----------------|---------------------|---------------------|
| **Tier 1: Dominant** | 10-30% | 8 genomes | ~150% (1.5 dominant) |
| **Tier 2: Common** | 1-10% | 75 genomes | ~30-40% |
| **Tier 3: Moderate** | 0.1-1% | 175 genomes | ~15-20% |
| **Tier 4: Rare** | 0.01-0.1% | 217 genomes | ~2-3% |
| **Tier 5: Very Rare** | <0.01% | 25 genomes | ~0.1% |

**Key Features**:
- Realistic "long tail" distribution
- Matches metagenomic rarefaction curves
- Allows dominant viruses (e.g., crAssphage in gut)
- High diversity in rare fraction

### Mathematical Model

```python
# Gut virome example (n=500)
abundances = []

# Tier 1: Dominant (8 genomes, 10-30% each)
for _ in range(8):
    abundances.append(np.random.uniform(0.10, 0.30))

# Tier 2: Common (75 genomes, 1-10% each)
for _ in range(75):
    abundances.append(np.random.uniform(0.01, 0.10))

# Tier 3: Moderate (175 genomes, 0.1-1% each)
for _ in range(175):
    abundances.append(np.random.uniform(0.001, 0.01))

# Tier 4: Rare (217 genomes, 0.01-0.1% each)
for _ in range(217):
    abundances.append(np.random.uniform(0.0001, 0.001))

# Tier 5: Very rare (25 genomes, <0.01% each)
for _ in range(25):
    abundances.append(np.random.uniform(0.00001, 0.0001))

# Normalize to sum to 1.0
abundances = np.array(abundances)
abundances = abundances / abundances.sum()
```

## Genome Selection Criteria

### Taxonomic Filtering

Genomes are selected based on ICTV taxonomy hierarchy:

```python
# Example: Select Siphoviridae genomes
SELECT g.genome_id, g.species_name, t.family
FROM genomes g
LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.family LIKE '%Siphoviridae%'
ORDER BY RANDOM()
```

### Host Filtering

Some targets require host-specific phages:

```python
# Example: Select Cutibacterium phages (skin virome)
SELECT g.genome_id
FROM genomes g
LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
LEFT JOIN host_associations h ON g.genome_id = h.genome_id
WHERE t.genus LIKE '%Cutibacterium%'
  AND h.host_name LIKE '%Cutibacterium%'
ORDER BY RANDOM()
```

### Genome Type Filtering

Some targets specify genome type (DNA vs RNA):

```python
# Example: Select dsDNA eukaryotic viruses
SELECT g.genome_id
FROM genomes g
LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.family LIKE '%Adenoviridae%'
  AND g.genome_type = 'dsDNA'
ORDER BY RANDOM()
```

## Output Files

### 1. Composition TSV

**File**: `{collection_id}_composition.tsv`

```tsv
genome_id       abundance       rank
NC_024711       0.23451234      1
NC_019724       0.15672341      2
NC_031940       0.08234123      3
...
```

**Columns**:
- `genome_id`: RefSeq accession
- `abundance`: Relative abundance (0-1, sum=1.0)
- `rank`: Abundance rank (1 = most abundant)

### 2. Summary Statistics

**File**: `{collection_id}_summary.txt`

```
Body Site Virome Collection Summary
======================================================================

Collection: Gut Virome - Adult Healthy (Western Diet)
ID: gut_virome_adult_healthy_western
Body Site: gut
Sample Type: stool
Health Status: healthy

Composition:
  Total genomes: 500
  Phage fraction: 96.0%

Taxonomic Targets:
  - Crassvirales: 150 genomes (30.0%)
  - Siphoviridae: 140 genomes (28.0%)
  - Myoviridae: 70 genomes (14.0%)
  ...

Abundance Distribution:
  Model: mixed
  Dominant tier (10-30%): 8 genomes
  Common tier (1-10%): 75 genomes
  ...

Statistics:
  Max abundance: 28.45%
  Mean abundance: 0.20%
  Median abundance: 0.05%
  Min abundance: 0.001%
  Shannon diversity: 5.23
```

### 3. Validation Report

**File**: `{collection_id}_validation.txt`

```
================================================================================
BODY SITE VIROME COLLECTION VALIDATION REPORT
================================================================================

Collection Information:
  ID: gut_virome_adult_healthy_western
  Name: Gut Virome - Adult Healthy (Western Diet)
  Body Site: gut
  ...

Quality Checks:
--------------------------------------------------------------------------------
  ✓ genome_count: PASS
      Expected: 500
      Actual: 500

  ✓ abundance_sum: PASS
      Expected: 1.0
      Actual: 1.0000

  ✓ taxonomy_coverage: PASS
      Threshold: 0.70
      Actual: 0.7420

  ✓ shannon_diversity: PASS
      Threshold: 3.0
      Actual: 5.23

  ...

Abundance Distribution:
--------------------------------------------------------------------------------
  Total Genomes: 500
  Shannon Diversity: 5.23
  Simpson Diversity: 0.9856

  Tier Distribution (genome counts):
    Dominant (10-30%):      8 genomes (143.2% of total abundance)
    Common (1-10%):        75 genomes (35.7% of total abundance)
    Moderate (0.1-1%):    175 genomes (18.4% of total abundance)
    ...

Taxonomic Composition:
--------------------------------------------------------------------------------
  Taxonomy Coverage: 74.2%

  Top 10 Families:
    Crassvirales                              150 (30.0%)
    Siphoviridae                              140 (28.0%)
    Myoviridae                                 70 (14.0%)
    ...
```

## Database Schema

### body_site_collections Table

```sql
CREATE TABLE body_site_collections (
    collection_id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    body_site TEXT NOT NULL,
    sample_type TEXT,
    health_status TEXT,
    genome_count INTEGER,
    phage_fraction REAL,
    curation_date TEXT,
    notes TEXT
);
```

### collection_genomes Table

```sql
CREATE TABLE collection_genomes (
    collection_id TEXT NOT NULL,
    genome_id TEXT NOT NULL,
    abundance REAL NOT NULL,
    inclusion_reason TEXT,
    PRIMARY KEY (collection_id, genome_id),
    FOREIGN KEY (collection_id) REFERENCES body_site_collections(collection_id),
    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id)
);
```

## Using Curated Collections

### Load Collection in Python

```python
from viroforge.core.community import create_body_site_profile

# Load curated gut virome
community = create_body_site_profile(
    body_site='gut',
    collection_name='gut_virome_adult_healthy_western',
    random_seed=42
)

# Genomes loaded with literature-validated abundances
print(f"Loaded {len(community.genomes)} genomes")
print(f"Max abundance: {max(community.abundances):.2%}")
```

### Query Collection in SQL

```sql
-- Get top 10 most abundant genomes in gut virome
SELECT cg.genome_id, g.species_name, cg.abundance
FROM collection_genomes cg
JOIN genomes g ON cg.genome_id = g.genome_id
WHERE cg.collection_id = 'gut_virome_adult_healthy_western'
ORDER BY cg.abundance DESC
LIMIT 10;
```

### Generate FASTQ with Collection

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.reads import FASTQGenerator

# Load collection
community = create_body_site_profile('gut', 'gut_virome_adult_healthy_western')

# Generate reads
generator = FASTQGenerator(
    random_seed=42,
    read_length=150,
    fragment_length_mean=350
)

reads = generator.generate_reads(
    community=community,
    coverage=10.0,
    output_prefix='gut_virome_mock'
)

# Output: gut_virome_mock_R1.fastq.gz, gut_virome_mock_R2.fastq.gz
```

## Quality Metrics

### Expected Validation Results

All collections should pass these quality checks:

| Check | Threshold | Expected | Notes |
|-------|-----------|----------|-------|
| Genome count | Exact match | 100% | Must match specification |
| Abundance sum | 0.99-1.01 | 1.000 | Normalized abundances |
| Taxonomy coverage | ≥70% | 74-75% | ICTV taxonomy |
| Shannon diversity | ≥3.0 | 4.5-6.0 | High diversity |
| Dominant tier | ≥1% | 1-2% | At least 1% dominant |
| Mid-tier distribution | ≥50% | 60-75% | Most in moderate/rare |

### Diversity Indices

**Shannon Diversity**: Measures evenness of abundance distribution
```
H' = -Σ(p_i * ln(p_i))
```
- Range: 0 to ln(N)
- Higher = more even distribution
- Expected: 4.5-6.0 for curated collections

**Simpson Diversity**: Probability two random reads are from different genomes
```
D = 1 - Σ(p_i²)
```
- Range: 0 to 1
- Higher = more diverse
- Expected: 0.98-0.99 for curated collections

## Troubleshooting

### Issue: Low Taxonomy Coverage (<70%)

**Cause**: Insufficient genomes with ICTV taxonomy for specific taxon

**Solution**:
1. Check ICTV taxonomy coverage in database
2. Try broader taxonomic rank (e.g., family instead of genus)
3. Add fallback selection without taxonomy filter

### Issue: Genome Count Mismatch

**Cause**: Not enough genomes match filtering criteria

**Solution**:
1. Check database has sufficient genomes for taxon
2. Relax host filter (try without host constraint)
3. Adjust target counts in specification

### Issue: Failed Shannon Diversity Check

**Cause**: Abundance distribution too uneven

**Solution**:
1. Check tier counts sum to target size
2. Verify abundance model parameters
3. Increase genomes in mid-tiers (moderate/rare)

## Timeline

**Phase 3, Week 6-8**: Body Site Collection Curation

| Day | Collection | Size | Status |
|-----|------------|------|--------|
| Week 6, Day 7 | Gut virome | 500 | Pending |
| Week 7, Day 1 | Oral virome | 200 | Pending |
| Week 7, Day 2 | Skin virome | 150 | Pending |
| Week 7, Day 3 | Respiratory virome | 200 | Pending |
| Week 7, Day 4-5 | Marine virome | 500 | Pending |
| Week 7, Day 6 | Soil virome | 300 | Pending |
| Week 7, Day 7 | Freshwater virome | 200 | Pending |
| Week 8, Day 1 | Mouse gut virome | 150 | Pending |

**Total**: 8 collections, 2,200 genomes

## References

### Gut Virome
- Dutilh et al. (2014) *Nat Commun* - crAssphage discovery
- Shkoporov et al. (2019) *Nat Microbiol* - Gut virome structure
- Guerin et al. (2018) *Cell Host Microbe* - Crassvirales diversity

### Oral Virome
- Edlund et al. (2015) *Cell Host Microbe* - Oral virome landscape
- Baker et al. (2021) *Microbiome* - Saliva virome composition

### Skin Virome
- Hannigan et al. (2015) *mBio* - Skin virome diversity
- Lobb et al. (2020) *Front Cell Infect Microbiol* - P. acnes phages

### General
- Zolfo et al. (2019) *Microbiome* - ViromeQC metrics
- Gregory et al. (2020) *Nat Rev Microbiol* - Viral ecology

---

**Document Version**: 1.0
**Last Updated**: 2025-11-01
**Author**: ViroForge Development Team
