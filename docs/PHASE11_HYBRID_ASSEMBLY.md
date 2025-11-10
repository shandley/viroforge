# Phase 11: Hybrid Assembly Support

**Date**: 2025-11-10
**Status**: In Progress
**Goal**: Enable matched short-read + long-read dataset generation for hybrid assembly benchmarking

---

## Overview

Hybrid assembly combines short-read accuracy with long-read length to produce superior assemblies. ViroForge will support generating matched datasets from the same collection with identical genome compositions.

## Architecture Design

### 1. Matched Dataset Generation

**Key Requirement**: Short and long reads must come from the **same genome composition**

**Approach**:
- Generate both platforms in a single script run
- Use the **same random seed** for both platforms
- Apply the **same VLP enrichment** (with read-type adjustment)
- Produce **identical genome lists and abundances**

### 2. Command-Line Interface

#### Option 1: Hybrid Flag (Preferred)
```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_hybrid \
    --hybrid-mode \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 \
    --depth 15 \
    --seed 42
```

**Advantages**:
- Clear intent (user wants hybrid dataset)
- Single run generates both
- Guaranteed matching composition

#### Option 2: Separate Runs with Seed Matching
```bash
# Short reads
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_short \
    --platform novaseq \
    --coverage 30 \
    --seed 42

# Long reads (SAME seed, SAME collection)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_long \
    --platform pacbio-hifi \
    --depth 15 \
    --seed 42
```

**Advantages**:
- Already works with current code
- More flexible (can generate at different times)
- Simpler implementation

**Disadvantages**:
- User must remember to use same seed
- Separate output directories
- No explicit linking metadata

### 3. Output Structure

#### Hybrid Mode Structure
```
data/gut_hybrid/
├── short_reads/
│   ├── fasta/
│   │   └── gut_virome.fasta
│   ├── fastq/
│   │   ├── gut_virome_R1.fastq
│   │   └── gut_virome_R2.fastq
│   └── metadata/
│       ├── gut_virome_metadata.json
│       ├── gut_virome_composition.tsv
│       └── gut_virome_abundances.txt
├── long_reads/
│   ├── fasta/
│   │   └── gut_virome.fasta (symlink to short_reads/fasta)
│   ├── fastq/
│   │   └── gut_virome_hifi.fastq.gz
│   └── metadata/
│       ├── gut_virome_metadata.json
│       ├── gut_virome_composition.tsv
│       └── gut_virome_ground_truth.tsv
└── hybrid_metadata.json  # Links the two datasets
```

#### Hybrid Metadata Format
```json
{
  "hybrid_dataset": true,
  "generation_timestamp": "2025-11-10T12:00:00",
  "random_seed": 42,
  "collection": {
    "id": 9,
    "name": "Human Gut Virome",
    "n_genomes": 25
  },
  "short_reads": {
    "platform": "novaseq",
    "coverage": 30,
    "read_length": 150,
    "output_dir": "short_reads/",
    "r1": "short_reads/fastq/gut_virome_R1.fastq",
    "r2": "short_reads/fastq/gut_virome_R2.fastq"
  },
  "long_reads": {
    "platform": "pacbio-hifi",
    "depth": 15,
    "read_length_mean": 15000,
    "passes": 10,
    "output_dir": "long_reads/",
    "reads": "long_reads/fastq/gut_virome_hifi.fastq.gz"
  },
  "composition_consistency": {
    "same_seed": true,
    "same_genomes": true,
    "same_abundances": true,
    "same_vlp_protocol": true
  }
}
```

### 4. Implementation Plan

#### Phase 11.1: Enhance generate_fastq_dataset.py
- Add `--hybrid-mode` flag
- Add `--short-platform` and `--long-platform` options
- Refactor to support generating both platforms in one run
- Create nested output directories
- Generate hybrid metadata JSON

#### Phase 11.2: Ensure Composition Matching
- Verify same seed produces identical genome lists
- Verify abundances match (accounting for VLP read-type adjustment)
- Add validation checks

#### Phase 11.3: Documentation
- Create hybrid assembly tutorial
- Document recommended assemblers (SPAdes hybrid, Unicycler, MaSuRCA)
- Provide benchmarking examples

#### Phase 11.4: Testing
- Add hybrid generation tests
- Verify composition consistency
- Test all platform combinations

### 5. Supported Hybrid Combinations

| Short Read | Long Read | Use Case |
|------------|-----------|----------|
| NovaSeq | PacBio HiFi | High accuracy + complete genomes |
| NovaSeq | Nanopore | Cost-effective, ultra-long scaffolds |
| MiSeq | PacBio HiFi | Moderate depth, high accuracy |
| HiSeq | Nanopore | Legacy data + long reads |

### 6. Hybrid Assembly Workflow Example

```bash
# 1. Generate hybrid dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_hybrid \
    --hybrid-mode \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 \
    --depth 15 \
    --vlp-protocol tangential_flow \
    --seed 42

# 2. Run hybrid assembler (Unicycler example)
unicycler \
    -1 data/gut_hybrid/short_reads/fastq/gut_virome_R1.fastq \
    -2 data/gut_hybrid/short_reads/fastq/gut_virome_R2.fastq \
    -l data/gut_hybrid/long_reads/fastq/gut_virome_hifi.fastq.gz \
    -o results/unicycler_hybrid

# 3. Evaluate against ground truth
# Compare with data/gut_hybrid/metadata/composition.tsv
```

### 7. Alternative: Separate Runs (Already Supported)

**Current Capability**: Users can generate matched datasets with separate runs:

```bash
# Generate short reads
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_novaseq \
    --platform novaseq \
    --coverage 30 \
    --seed 42

# Generate long reads (SAME collection, SAME seed)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut_hifi \
    --platform pacbio-hifi \
    --depth 15 \
    --seed 42

# Use both with hybrid assembler
unicycler \
    -1 data/gut_novaseq/fastq/gut_virome_R1.fastq \
    -2 data/gut_novaseq/fastq/gut_virome_R2.fastq \
    -l data/gut_hifi/fastq/gut_virome_hifi.fastq.gz \
    -o results/unicycler
```

**This already works!** The hybrid mode would just make it more convenient and explicit.

---

## Implementation Decision

### Recommended Approach: Document Existing Capability + Add Convenience Script

**Rationale**:
1. **Already works**: Users can generate matched datasets with same seed
2. **Simple**: No major code changes needed
3. **Flexible**: Users control timing and parameters

**Enhancements**:
1. Create `scripts/generate_hybrid_dataset.py` convenience script
2. Add validation that seed/collection match
3. Document in tutorial with examples
4. Optionally add `--hybrid-mode` flag in future

### Quick Implementation

#### Phase 11.1: Convenience Script (1-2 hours)
Create `scripts/generate_hybrid_dataset.py` that calls `generate_fastq_dataset.py` twice with matched parameters.

#### Phase 11.2: Tutorial (1-2 hours)
Create `docs/HYBRID_ASSEMBLY_TUTORIAL.md` with:
- How to generate matched datasets
- Supported hybrid assemblers
- Benchmarking workflows
- Examples

#### Phase 11.3: Validation Utility (1 hour)
Create `scripts/validate_hybrid_composition.py` to verify two datasets match.

---

## Next Steps

1. Create convenience script for hybrid generation
2. Create hybrid assembly tutorial
3. Add validation utility
4. Update main README

Total Estimated Time: 4-5 hours
