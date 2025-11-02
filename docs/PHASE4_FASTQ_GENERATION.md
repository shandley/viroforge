# Phase 4: FASTQ Generation Workflows

**Status**: ✅ Completed
**Date**: 2025-11-01
**ViroForge Version**: 0.3.0

## Overview

Phase 4 implements a complete FASTQ generation pipeline that converts curated body site collections from the database into realistic synthetic sequencing datasets. This infrastructure enables benchmarking of virome analysis pipelines with ground truth data.

## Components

### 1. Core Generation Script

**File**: `scripts/generate_fastq_dataset.py`

Generates FASTQ files from a single body site collection with configurable parameters.

**Key Features**:
- Loads genomes directly from SQLite database
- Simulates VLP enrichment (optional)
- Integrates with InSilicoSeq for realistic error profiles
- Exports complete ground truth metadata (JSON + TSV)
- Supports NovaSeq, MiSeq, and HiSeq platforms

**Usage**:
```bash
# Generate gut virome with VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --coverage 10 \
    --platform novaseq

# List available collections
python scripts/generate_fastq_dataset.py --list-collections

# Generate without VLP (bulk metagenome)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output data/fastq/marine_bulk \
    --no-vlp \
    --coverage 10
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

3. **`vlp-comparison`**: VLP vs bulk comparison
   - Collections: Gut and Marine
   - Generates both VLP and non-VLP versions
   - Coverage: 10x

4. **`platform-comparison`**: Cross-platform comparison
   - Collection: Gut
   - Platforms: NovaSeq, MiSeq, HiSeq
   - Coverage: 10x

5. **`coverage-series`**: Coverage depth series
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

### VLP Enrichment Simulation

VLP (Virus-Like Particle) enrichment is simulated by adjusting relative abundances:

```python
if apply_vlp:
    # Increase viral recovery ~20% with variation
    abundance = abundance * (1.0 + np.random.normal(0.2, 0.05))
```

This simulates the enrichment effect of nuclease treatment that removes non-encapsidated DNA.

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

**JSON** (`metadata.json`):
```json
{
  "generation_info": {
    "timestamp": "2025-11-01T18:22:09.895652",
    "viroforge_version": "0.3.0",
    "random_seed": 42
  },
  "collection": {
    "id": 16,
    "name": "Mouse Gut Virome - Laboratory (C57BL/6)",
    "description": "Mouse gut virome from laboratory C57BL/6 mice...",
    "n_genomes": 22
  },
  "configuration": {
    "coverage": 1.0,
    "n_reads": null,
    "read_length": 150,
    "insert_size": 350,
    "platform": "miseq",
    "vlp_enrichment": true,
    "vlp_efficiency": 0.95
  },
  "genomes": [...]
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

## Future Enhancements

1. **Advanced VLP Simulation**: Integrate full ViroForge enrichment framework
2. **Amplification Bias**: Add MDA/SISPA amplification simulation
3. **Contamination Simulation**: Add non-viral contaminants at configurable levels
4. **Custom Abundance Distributions**: Support user-defined abundance profiles
5. **Multi-Sample Generation**: Generate related samples for longitudinal studies
6. **Quality Metrics**: Add read quality distribution controls

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

Phase 4 provides a complete infrastructure for generating realistic synthetic virome datasets with:
- ✅ Database-driven collection loading
- ✅ VLP enrichment simulation
- ✅ Platform-specific error profiles
- ✅ Complete ground truth tracking
- ✅ Batch generation capabilities
- ✅ Multiple use case presets

This enables systematic benchmarking of virome analysis tools with known composition datasets across different experimental conditions and sequencing platforms.
