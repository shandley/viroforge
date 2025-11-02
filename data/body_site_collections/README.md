# Body Site Virome Collections

This directory contains curated body site-specific virome collections for ViroForge.

## Overview

ViroForge includes 8 pre-defined body site collections based on literature-validated compositions:

| Collection | Size | Body Site | Key Feature |
|------------|------|-----------|-------------|
| `gut_virome_adult_healthy_western` | 500 | Gut | crAssphage-dominated (30%) |
| `oral_virome_saliva_healthy` | 200 | Oral cavity | Streptococcus phages |
| `skin_virome_sebaceous_healthy` | 150 | Skin | P. acnes phages (47%) |
| `respiratory_virome_nasopharynx_healthy` | 200 | Respiratory | Mixed phage/eukaryotic |
| `marine_virome_coastal_surface` | 500 | Marine | Pelagiphages |
| `soil_virome_agricultural` | 300 | Soil | Actinophages |
| `freshwater_virome_lake_surface` | 200 | Freshwater | Cyanophages |
| `mouse_gut_virome_lab` | 150 | Mouse gut | Lactobacillus phages |

**Total**: 2,200 curated genomes across 8 body sites

## Quick Start

### Curate a Collection

```bash
# Curate gut virome (500 genomes)
python scripts/curate_body_site_collections.py \
  --collection gut \
  --output data/body_site_collections

# Curate all collections
python scripts/curate_body_site_collections.py \
  --all \
  --output data/body_site_collections
```

### Validate a Collection

```bash
# Validate gut virome
python scripts/validate_body_site_collections.py \
  --collection gut_virome_adult_healthy_western \
  --output data/body_site_collections/validation

# Validate all collections
python scripts/validate_body_site_collections.py \
  --all \
  --output data/body_site_collections/validation
```

## Output Files

Each curated collection generates three files:

### 1. Composition TSV
**File**: `{collection_id}_composition.tsv`

Contains genome IDs, abundances, and ranks:
```tsv
genome_id       abundance       rank
NC_024711       0.23451234      1
NC_019724       0.15672341      2
```

### 2. Summary Statistics
**File**: `{collection_id}_summary.txt`

Contains:
- Collection metadata (body site, sample type, health status)
- Taxonomic target breakdown
- Abundance distribution statistics
- Diversity indices (Shannon, Simpson)

### 3. Validation Report
**File**: `validation/{collection_id}_validation.txt`

Contains:
- Quality check results (PASS/WARN/FAIL)
- Taxonomic composition analysis
- Abundance tier distribution
- Genome characteristics (length, GC, types)

## Directory Structure

```
data/body_site_collections/
├── README.md
├── gut_virome_adult_healthy_western_composition.tsv
├── gut_virome_adult_healthy_western_summary.txt
├── oral_virome_saliva_healthy_composition.tsv
├── oral_virome_saliva_healthy_summary.txt
├── ...
└── validation/
    ├── gut_virome_adult_healthy_western_validation.txt
    ├── oral_virome_saliva_healthy_validation.txt
    └── ...
```

## Using Collections in ViroForge

### Load Collection in Python

```python
from viroforge.core.community import create_body_site_profile

# Load curated gut virome
community = create_body_site_profile(
    body_site='gut',
    collection_name='gut_virome_adult_healthy_western',
    random_seed=42
)

print(f"Loaded {len(community.genomes)} genomes")
```

### Generate Mock Dataset

```python
from viroforge.core.community import create_body_site_profile
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000

# Load collection
community = create_body_site_profile('gut', 'gut_virome_adult_healthy_western')

# Apply VLP enrichment
vlp = standard_vlp()
vlp.apply(community)

# Apply amplification
amp = rdab_40_cycles()
amp.apply(community)

# Generate reads (see examples/ for complete workflow)
```

## Abundance Model

Collections use a **5-tier abundance model** for realistic viral community structure:

| Tier | Abundance Range | Typical Count (n=500) | % of Total |
|------|----------------|---------------------|-----------|
| **Dominant** | 10-30% | 8 genomes | ~150% |
| **Common** | 1-10% | 75 genomes | ~35% |
| **Moderate** | 0.1-1% | 175 genomes | ~18% |
| **Rare** | 0.01-0.1% | 217 genomes | ~2% |
| **Very Rare** | <0.01% | 25 genomes | ~0.1% |

This model:
- ✓ Creates realistic "long tail" distribution
- ✓ Allows dominant viruses (e.g., crAssphage)
- ✓ Matches metagenomic rarefaction curves
- ✓ High diversity in rare fraction

## Quality Metrics

All curated collections are validated against these criteria:

| Metric | Threshold | Expected |
|--------|-----------|----------|
| Genome count | Exact | 100% match |
| Abundance sum | 0.99-1.01 | 1.000 |
| ICTV taxonomy | ≥70% | 74-75% |
| Shannon diversity | ≥3.0 | 4.5-6.0 |
| Simpson diversity | - | 0.98-0.99 |

## Collection Specifications

### Gut Virome (500 genomes)
**Target**: Adult healthy gut, Western diet
**Composition**:
- Crassvirales: 150 (30%) - crAssphage dominant
- Siphoviridae: 140 (28%) - Temperate phages
- Myoviridae: 70 (14%) - Lytic phages
- Podoviridae: 40 (8%)
- Microviridae: 50 (10%)
- Others: 50 (10%)

**Key References**: Dutilh 2014, Shkoporov 2019, Guerin 2018

### Oral Virome (200 genomes)
**Target**: Healthy saliva
**Composition**:
- Siphoviridae: 80 (40%) - Streptococcus phages
- Myoviridae: 40 (20%)
- Podoviridae: 30 (15%)
- Eukaryotic: 15 (7.5%) - Herpesviruses, HPV

**Key References**: Edlund 2015, Baker 2021

### Skin Virome (150 genomes)
**Target**: Sebaceous sites (forehead, back)
**Composition**:
- Cutibacterium phages: 70 (47%) - P. acnes dominant
- Staphylococcus phages: 40 (27%)
- Corynebacterium phages: 15 (10%)
- Others: 25 (16%)

**Key References**: Hannigan 2015, Lobb 2020

### Marine Virome (500 genomes)
**Target**: Coastal surface water
**Composition**:
- Caudoviricetes: 400 (80%) - Pelagiphages, cyanophages
- Microviridae: 50 (10%)
- Phycodnaviridae: 30 (6%)
- Others: 20 (4%)

**Key References**: Brum 2015, Gregory 2019

See `docs/BODY_SITE_COLLECTIONS.md` for complete specifications.

## Troubleshooting

### No collections found?
Run the curation script first:
```bash
python scripts/curate_body_site_collections.py --all
```

### Validation failures?
Check:
1. RefSeq database is complete (14,568 genomes)
2. ICTV taxonomy applied (70%+ coverage)
3. Database schema includes collection tables

### Low diversity?
- Increase genomes in moderate/rare tiers
- Check tier counts sum to target size
- Verify abundance model parameters

## Documentation

- **Detailed Workflow**: `docs/BODY_SITE_CURATION_WORKFLOW.md`
- **Collection Specs**: `docs/BODY_SITE_COLLECTIONS.md`
- **Curation Plans**: `docs/GUT_VIROME_CURATION.md`

## Scripts

- **Curation**: `scripts/curate_body_site_collections.py`
- **Validation**: `scripts/validate_body_site_collections.py`

---

**Last Updated**: 2025-11-01
**Phase**: Phase 3, Week 6-8
