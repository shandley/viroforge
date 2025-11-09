# Phase 4: FASTQ Generation Workflows

**Status**: ✅ Completed (Enhanced in Phase 5 & 6)
**Date**: 2025-11-01 (Updated: 2025-11-09)
**ViroForge Version**: 0.5.0-dev

## Overview

Phase 4 implements a complete FASTQ generation pipeline that converts curated body site collections from the database into realistic synthetic sequencing datasets. This infrastructure enables benchmarking of virome analysis pipelines with ground truth data.

## Components

### 1. Core Generation Script

**File**: `scripts/generate_fastq_dataset.py`

Generates FASTQ files from a single body site collection with configurable parameters.

**Key Features**:
- Loads genomes directly from SQLite database
- Enhanced VLP enrichment with 5 protocols (Phase 5)
- Size-based filtration and contamination reduction
- Integrates with InSilicoSeq for realistic error profiles
- Exports complete ground truth metadata (JSON + TSV)
- Supports NovaSeq, MiSeq, and HiSeq platforms

**Usage**:
```bash
# Generate gut virome with VLP enrichment (default: tangential flow)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --coverage 10 \
    --platform novaseq

# Specify VLP protocol
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_ultracentrifuge \
    --coverage 10 \
    --vlp-protocol ultracentrifugation

# List available collections
python scripts/generate_fastq_dataset.py --list-collections

# Generate without VLP (bulk metagenome)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output data/fastq/marine_bulk \
    --no-vlp \
    --coverage 10

# Specify contamination level
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_heavy_contam \
    --coverage 10 \
    --contamination-level heavy

# Add amplification bias (Phase 6) ⭐ NEW
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_rdab \
    --coverage 10 \
    --amplification rdab

# MDA amplification for low biomass samples
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_mda \
    --coverage 10 \
    --amplification mda
```

**Key Classes**:

- **`CollectionLoader`**: Database integration
  - Loads collections and genomes with abundances
  - Joins genome sequences with taxonomy
  - Lists available collections

- **`FASTQGenerator`**: FASTQ generation workflow
  - Prepares genomes and adjusts abundances (VLP simulation)
  - Writes FASTA with abundance annotations
  - Calls InSilicoSeq for read generation
  - Exports metadata

### 2. Batch Generation Script

**File**: `scripts/batch_generate_fastq.py`

Automates generation of multiple datasets with different configurations.

**Preset Configurations**:

1. **`quick-test`**: Small datasets for rapid testing
   - Collections: Mouse gut, Skin, Oral (smallest 3)
   - Coverage: 1x
   - Platform: MiSeq

2. **`benchmark-standard`**: Full benchmark suite
   - Collections: All 8 collections
   - Coverage: 10x
   - Platform: NovaSeq

3. **`vlp-protocol-comparison`**: Compare all VLP protocols (Phase 5)
   - Collection: Gut
   - Protocols: tangential_flow, syringe, ultracentrifugation, norgen, bulk
   - Coverage: 10x
   - Generates 5 datasets for direct comparison

4. **`amplification-comparison`**: Compare amplification methods (Phase 6) ⭐ NEW
   - Collection: Gut
   - Methods: none, rdab, rdab-30, mda, mda-long, linker
   - Coverage: 10x
   - Generates 6 datasets for method comparison

5. **`vlp-comparison`**: VLP vs bulk comparison
   - Collections: Gut and Marine
   - Generates both VLP and non-VLP versions
   - Coverage: 10x

6. **`platform-comparison`**: Cross-platform comparison
   - Collection: Gut
   - Platforms: NovaSeq, MiSeq, HiSeq
   - Coverage: 10x

7. **`coverage-series`**: Coverage depth series
   - Collection: Gut
   - Coverages: 1x, 5x, 10x, 20x, 50x
   - Platform: NovaSeq

**Usage**:
```bash
# Run quick test
python scripts/batch_generate_fastq.py \
    --preset quick-test \
    --output data/test_datasets

# Run full benchmark suite
python scripts/batch_generate_fastq.py \
    --preset benchmark-standard \
    --output data/benchmark_datasets

# Compare VLP protocols (Phase 5)
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output data/vlp_protocols

# Use custom config file
python scripts/batch_generate_fastq.py \
    --config my_config.yaml \
    --output data/custom_datasets
```

### 3. Output Structure

Each dataset generates:

```
output_dir/
├── fasta/
│   └── collection_name.fasta          # Reference genomes with abundances
├── fastq/
│   ├── collection_name_R1.fastq       # Forward reads
│   └── collection_name_R2.fastq       # Reverse reads
└── metadata/
    ├── collection_name_metadata.json  # Complete ground truth
    ├── collection_name_composition.tsv # Composition table
    └── collection_name_abundances.txt  # ISS abundance file
```

## Technical Details

### Enhanced VLP Enrichment Modeling (Phase 5)

VLP (Virus-Like Particle) enrichment now uses a comprehensive biological model with size-based filtration and type-specific contamination reduction.

#### Available VLP Protocols

| Protocol | Filtration Method | Pore Size | Nuclease Efficiency | Viral Recovery | Contamination Reduction |
|----------|------------------|-----------|---------------------|----------------|------------------------|
| **tangential_flow** | 0.2 μm TFF | 0.2 μm | 98% | 85% | 91.2% |
| **syringe** | 0.22 μm syringe | 0.22 μm | 90% | 60% | 85.7% |
| **ultracentrifugation** | Density gradient | N/A | 95% | 90% | 88.4% |
| **norgen** | Column-based | N/A | 92% | 70% | 87.1% |
| **none** (--no-vlp) | Bulk metagenome | N/A | 0% | 100% | 0% |

#### VLP Enrichment Features

**1. Virion Size Estimation**
- Genome length and type (dsDNA, ssDNA, ssRNA, dsRNA) used to estimate virion diameter
- Based on empirical relationships from literature (Cui et al. 2014, Nasir et al. 2017)

**2. Size-Based Filtration**
- Protocol-specific retention curves (sigmoid for TFF, step for syringe)
- Smaller viruses partially lost during filtration
- Larger contaminants (bacteria) highly removed

**3. Type-Specific Contamination Reduction**
- **Host DNA**: 90-98% removal (DNase treatment)
- **rRNA**: 85-95% removal (RNase + size-based)
- **Bacteria**: 94-99% removal (filtration)
- **PhiX**: 10-40% retention (treated as small virus)

**4. Stochastic Variation**
- Realistic biological and technical variability
- 5-15% coefficient of variation depending on mechanism

See [VLP Integration Guide](VLP_CONTAMINATION_INTEGRATION.md) for detailed implementation.

### Coverage Calculation

For paired-end reads, coverage is calculated as:

```python
n_reads = (total_genome_length * coverage) / (2 * read_length)
```

Where the factor of 2 accounts for both mates covering the same fragment.

### InSilicoSeq Integration

ISS parameters used:
- `--model`: Platform-specific error model (novaseq, miseq, hiseq)
- `--mode basic`: Faster generation with simplified error model
- `--abundance_file`: Genome-specific relative abundances
- `--n_reads`: Total number of read pairs
- `--seed`: Reproducible generation

## Test Results

### Test Dataset: Mouse Gut Virome (Collection 16)

**Configuration**:
- Genomes: 22 viral genomes
- Total length: 522,734 bp
- Coverage: 1x
- Platform: MiSeq
- VLP: Enabled

**Output**:
```
✓ FASTQ generation complete!
  Collection: Mouse Gut Virome - Laboratory (C57BL/6)
  Genomes: 22
  Coverage: 1.0x
  Platform: miseq
  VLP enrichment: True

  Output files:
    - R1: data/test_fastq/fastq/mouse_gut_virome_-_laboratory_c57bl6_R1.fastq (237 KB)
    - R2: data/test_fastq/fastq/mouse_gut_virome_-_laboratory_c57bl6_R2.fastq (237 KB)
    - FASTA: data/test_fastq/fasta/mouse_gut_virome_-_laboratory_c57bl6.fasta
    - Metadata: data/test_fastq/metadata/
```

**Verification**:
- Read pairs generated: 871
- Proper FASTQ format: ✓
- Genome IDs in headers: ✓
- Quality scores present: ✓
- Metadata files complete: ✓

### Ground Truth Metadata Example

**JSON** (`metadata.json`) - Enhanced with Phase 5 VLP Statistics:
```json
{
  "generation_info": {
    "timestamp": "2025-11-08T18:22:09.895652",
    "viroforge_version": "0.4.0",
    "random_seed": 42
  },
  "collection": {
    "id": 9,
    "name": "Gut Virome - Adult Healthy (Western Diet)",
    "description": "Gut virome from healthy adults...",
    "n_genomes": 134
  },
  "configuration": {
    "coverage": 10.0,
    "read_length": 150,
    "insert_size": 350,
    "platform": "novaseq",
    "vlp_protocol": "tangential_flow",
    "contamination_level": "realistic"
  },
  "enrichment_statistics": {
    "vlp_protocol": "tangential_flow",
    "viral_fraction_before": 0.926,
    "viral_fraction_after": 0.993,
    "viral_enrichment_fold_change": 1.07,
    "contamination_reduction": {
      "overall_reduction": 0.912,
      "host_dna_removal": 0.964,
      "rrna_removal": 0.891,
      "bacteria_removal": 0.990,
      "phix_retention": 0.600
    }
  },
  "sequences": [
    {
      "genome_id": "GCF_000001405.40",
      "genome_name": "Enterobacteria phage T7",
      "type": "viral",
      "length": 39937,
      "gc_content": 0.485,
      "relative_abundance": 0.0234,
      "taxonomy": {
        "family": "Podoviridae",
        "genus": "T7virus",
        "species": "Enterobacteria phage T7"
      }
    }
  ]
}
```

**TSV** (`composition.tsv`):
| genome_id | genome_name | length | gc_content | genome_type | family | genus | species | relative_abundance | original_abundance |
|-----------|-------------|--------|------------|-------------|--------|-------|---------|-------------------|-------------------|
| GCF_004146805.1 | Salmonella phage SeSz-2 | 45049 | 0.460 | dsDNA | Unknown | | Salmonella phage SeSz-2 | 0.3127 | 0.3107 |
| GCF_018595485.1 | Chrysothrix chrysovirus 1 | 12382 | 0.358 | ssRNA | Chrysoviridae | Alphachrysovirus | ... | 0.2754 | 0.2809 |

## Available Body Site Collections

From validation in Phase 3:

1. **ID 9**: Gut Virome - Adult Healthy (Western Diet) - 134 genomes
2. **ID 10**: Oral Virome - Saliva (Healthy) - 47 genomes
3. **ID 11**: Skin Virome - Sebaceous Sites (Healthy) - 15 genomes
4. **ID 12**: Respiratory Virome - Nasopharynx (Healthy) - 41 genomes
5. **ID 13**: Marine Virome - Coastal Surface Water - 448 genomes
6. **ID 14**: Soil Virome - Agricultural - 291 genomes
7. **ID 15**: Freshwater Virome - Lake Surface Water - 200 genomes
8. **ID 16**: Mouse Gut Virome - Laboratory (C57BL/6) - 22 genomes

## Dependencies

- Python 3.9+
- BioPython
- InSilicoSeq (iss)
- NumPy
- Pandas
- SQLite3

Install InSilicoSeq:
```bash
conda install -c bioconda insilicoseq
# or
pip install InSilicoSeq
```

## Troubleshooting

### InSilicoSeq Not Found
```bash
conda install -c bioconda insilicoseq
```

### Database Not Found
Ensure Phase 3 database population completed:
```bash
ls -lh viroforge/data/viral_genomes.db
```

### Collection Not Found
List available collections:
```bash
python scripts/generate_fastq_dataset.py --list-collections
```

## Phase 5 Enhancements (Completed)

✅ **Advanced VLP Simulation**: Full enrichment framework with 5 protocols
✅ **Contamination Simulation**: Type-specific contaminants at 3 levels (clean, realistic, heavy)
✅ **Size-Based Filtration**: Virion size estimation and protocol-specific retention
✅ **Contamination Reduction**: Type-specific reduction (host DNA, rRNA, bacteria, PhiX)
✅ **Literature Validation**: All parameters validated against published studies

## Future Enhancement Opportunities

1. **Amplification Bias**: Add RdAB/MDA amplification simulation to workflow
2. **Custom Abundance Distributions**: Support user-defined abundance profiles
3. **Multi-Sample Generation**: Generate related samples for longitudinal studies
4. **Quality Metrics**: Add read quality distribution controls
5. **Real Contamination Databases**: Integrate actual host genomes and SILVA rRNA
6. **Performance Optimization**: Parallel FASTQ generation for large batches

## Integration with Hecatomb

Generated datasets are designed for benchmarking Hecatomb's virome analysis pipeline:

```bash
# Generate benchmark dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmark/gut_virome \
    --coverage 10

# Run Hecatomb
hecatomb run \
    --reads benchmark/gut_virome/fastq/*_R{1,2}.fastq \
    --outdir results/gut_benchmark

# Compare results to ground truth
python scripts/evaluate_results.py \
    --results results/gut_benchmark \
    --truth benchmark/gut_virome/metadata/gut_virome_metadata.json
```

## References

- InSilicoSeq: [https://github.com/HadrienG/InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)
- ViroForge Phase 3: See `docs/BODY_SITE_COLLECTIONS.md`
- ICTV Taxonomy: See database taxonomy mappings

## Summary

Phase 4 (Enhanced in Phase 5) provides a complete infrastructure for generating realistic synthetic virome datasets with:
- ✅ Database-driven collection loading (14,423 RefSeq genomes)
- ✅ Enhanced VLP enrichment (5 protocols with size-based filtration)
- ✅ Type-specific contamination reduction (literature-validated)
- ✅ Platform-specific error profiles (NovaSeq, MiSeq, HiSeq)
- ✅ Complete ground truth tracking (viral + contaminants)
- ✅ Batch generation capabilities (6 presets)
- ✅ Multiple use case presets (VLP comparison, platform comparison, coverage series)

This enables systematic benchmarking of virome analysis tools with known composition datasets across different VLP protocols, experimental conditions, and sequencing platforms.

**Literature Validated**: All VLP parameters validated against peer-reviewed studies (Thurber et al. 2009, Shkoporov et al. 2018, Roux et al. 2016, Lim et al. 2020)
