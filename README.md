# ViroForge

**Forging synthetic viromes for benchmarking and validation**

A comprehensive mock metavirome data generator for testing and validating virome analysis pipelines.

[![Tests](https://img.shields.io/badge/tests-178%20passing-brightgreen)](tests/)
[![Phase](https://img.shields.io/badge/Phase%204-Complete-success)](docs/IMPLEMENTATION_PLAN.md)
[![Collections](https://img.shields.io/badge/collections-8%20body%20sites-blue)](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)
[![Genomes](https://img.shields.io/badge/genomes-14%2C423%20RefSeq-blue)](docs/GENOME_DATABASE_DESIGN.md)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## Overview

ViroForge generates realistic synthetic virome sequencing datasets with complete ground truth metadata, enabling rigorous validation of QC pipelines, assembly tools, taxonomic classifiers, and analysis workflows.

**What makes ViroForge different?** It's the first simulator to model the complete virome workflow: from viral community composition through VLP enrichment, library amplification, sequencing artifacts, to final FASTQ files.

### Key Features

- **8 Curated Body Site Collections** - Literature-validated virome compositions (gut, oral, skin, respiratory, marine, soil, freshwater, mouse gut)
- **14,423 RefSeq Viral Genomes** - Complete database with ICTV taxonomy integration
- **Database-Driven FASTQ Generation** - Direct generation from curated collections with InSilicoSeq
- **Complete Ground Truth** - Taxonomic composition, abundance tables, genome-read mappings
- **VLP Enrichment Simulation** - Realistic filtration, nuclease treatment effects
- **Platform-Specific Error Models** - NovaSeq, MiSeq, HiSeq with realistic artifacts
- **Reproducible Benchmarks** - Random seeds, complete metadata, known composition
- **Production Ready** - 178 tests passing, used for real pipeline validation

---

## Why ViroForge?

### The Problem

Current virome analysis tool validation approaches are limited:
- **Physical synthetic communities**: Expensive ($10k+), time-consuming, limited complexity
- **Existing simulators**: Bacterial-focused (CAMISIM), don't model VLP enrichment
- **Real datasets**: Unknown ground truth, can't systematically test edge cases

### The Solution

ViroForge enables:
- ✅ Unlimited synthetic datasets with complete ground truth
- ✅ VLP enrichment vs bulk metagenome comparisons
- ✅ Realistic contamination profiles (host DNA, rRNA, reagent bacteria)
- ✅ Library prep and sequencing artifact modeling
- ✅ Cross-platform reproducibility testing (NovaSeq, MiSeq, etc.)
- ✅ Standardized benchmarking datasets for the community

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/shandley/viroforge.git
cd viroforge

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install biopython numpy pandas
pip install insilicoseq  # For FASTQ generation
```

### Generate FASTQ from Curated Collections (Recommended)

**Step 1: List Available Collections**

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

**Step 2: Generate FASTQ Dataset**

```bash
# Generate gut virome with VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --coverage 10 \
    --platform novaseq
```

Output:
```
output/
├── fasta/
│   └── gut_virome_adult_healthy_western_diet.fasta  # Reference genomes
├── fastq/
│   ├── gut_virome_adult_healthy_western_diet_R1.fastq  # Forward reads
│   └── gut_virome_adult_healthy_western_diet_R2.fastq  # Reverse reads
└── metadata/
    ├── gut_virome_adult_healthy_western_diet_metadata.json  # Complete ground truth
    └── gut_virome_adult_healthy_western_diet_composition.tsv  # Abundance table
```

**Step 3: Use with Hecatomb or Other Pipelines**

```bash
# Run Hecatomb on generated dataset
hecatomb run \
    --reads data/fastq/gut_virome/fastq/*_R{1,2}.fastq \
    --outdir results/gut_benchmark

# Compare results to ground truth
# See metadata/gut_virome_metadata.json for known composition
```

See [FASTQ Generation Guide](docs/PHASE4_FASTQ_GENERATION.md) for detailed documentation.

---

## Use Cases

### 1. Benchmark Virome Analysis Pipelines

```bash
# Generate gut virome benchmark dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmarks/gut_10x \
    --coverage 10 \
    --platform novaseq

# Run your pipeline
your_pipeline benchmarks/gut_10x/fastq/*_R{1,2}.fastq

# Compare results to ground truth
# See benchmarks/gut_10x/metadata/gut_virome_metadata.json
```

### 2. Compare VLP vs Bulk Metagenome

```bash
# Generate with VLP enrichment (default)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_vlp \
    --coverage 10

# Generate bulk metagenome (no VLP)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_bulk \
    --coverage 10 \
    --no-vlp

# Compare viral recovery rates between VLP and bulk
```

### 3. Cross-Platform Reproducibility Testing

```bash
# Generate datasets on different platforms
python scripts/generate_fastq_dataset.py --collection-id 9 --platform novaseq --output novaseq_gut
python scripts/generate_fastq_dataset.py --collection-id 9 --platform miseq --output miseq_gut
python scripts/generate_fastq_dataset.py --collection-id 9 --platform hiseq --output hiseq_gut

# Compare platform-specific artifacts and assembly quality
```

### 4. Batch Generation for Comprehensive Benchmarks

```bash
# Generate all 8 collections at 10x coverage
python scripts/batch_generate_fastq.py \
    --preset benchmark-standard \
    --output data/benchmark_suite

# Or run quick tests
python scripts/batch_generate_fastq.py \
    --preset quick-test \
    --output data/test_datasets
```

See [docs/PHASE4_FASTQ_GENERATION.md](docs/PHASE4_FASTQ_GENERATION.md) for complete documentation.

---

## Complete Workflow Example

ViroForge models the complete virome data generation pipeline:

```
Viral Community → Contamination → VLP Enrichment → Amplification → Sequencing → Platform Artifacts
      ↓                 ↓                ↓                ↓               ↓              ↓
   50 genomes      + Host DNA      97% viral        Length/GC      100K reads      +2.5% polyG
   Gut-specific    + Bacteria       (1.94x)          bias           Paired-end      +9% optical dups
   Log-normal      + Fungal                          (RdAB)                         +1.5% index hop
```

**Ground Truth Tracking**: Every read is traced back to its source genome throughout the entire pipeline.

---

## Pre-Built Components

### VLP Enrichment Protocols

| Protocol | Filtration | Nuclease | Use Case |
|----------|------------|----------|----------|
| `standard_vlp()` | 0.2 μm TFF | 95% | Most common protocol |
| `iron_chloride_vlp()` | FeCl3 precipitation | 98% | High purity studies |
| `ultracentrifugation_vlp()` | Density gradient | 90% | Traditional method |
| `syringe_filter_vlp()` | 0.45 μm syringe | 80% | Field-friendly |

### Amplification Methods

| Method | Length Bias | GC Bias | Use Case |
|--------|-------------|---------|----------|
| `rdab_40_cycles()` | Strong | Moderate | Standard virome protocol |
| `rdab_30_cycles()` | Moderate | Moderate | Less bias needed |
| `mda_standard()` | None | Extreme | Low-biomass samples |
| `linker_standard()` | None | Minimal | Modern protocols |
| `no_amplification()` | None | None | High-biomass control |

### Platform Profiles

| Platform | Flow Cell | PolyG Tails | Optical Dups | Index Hopping |
|----------|-----------|-------------|--------------|---------------|
| `novaseq_6000()` | Patterned | 2.5% | 9% | 1.5% |
| `nextseq_2000()` | Patterned | 2.0% | 7% | 1.0% |
| `miseq()` | Cluster | 0% | 2.5% | 0.1% |
| `hiseq_2500()` | Cluster | 0% | 4.5% | 0.2% |
| `no_artifacts()` | Ideal | 0% | 0% | 0% |

---

## Examples

ViroForge includes comprehensive examples in the `examples/` directory:

### Basic Usage
- **`create_community_example.py`** - Create viral communities
- **`create_contamination_example.py`** - Model contamination
- **`vlp_enrichment_basic.py`** - Apply VLP enrichment

### Protocol Comparisons
- **`vlp_protocol_comparison.py`** - Compare different VLP methods
- **`vlp_vs_bulk_comparison.py`** - VLP vs bulk metagenome
- **`amplification_comparison.py`** - Compare amplification methods
- **`platform_comparison.py`** - Compare sequencing platforms

### Complete Workflows
- **`complete_workflow_integrated.py`** - End-to-end pipeline (recommended starting point)
- **`cross_platform_workflow.py`** - NovaSeq vs MiSeq comparison

Run any example:
```bash
python examples/complete_workflow_integrated.py
```

---

## Output Files

### Ground Truth Metadata

Every ViroForge dataset includes complete ground truth:

**`ground_truth_composition.tsv`**
```
genome_id    taxonomy                   abundance    length    gc_content
NC_001416    Enterobacteria phage T7    0.1234      39937     0.485
NC_007458    Escherichia phage MS2      0.0567      3569      0.518
...
```

**`ground_truth_read_mapping.tsv`**
```
read_id              genome_id    genome_name              family
read_001_forward     NC_001416    Enterobacteria phage T7  Podoviridae
read_002_forward     NC_007458    Escherichia phage MS2    Leviviridae
...
```

**`pipeline_summary.txt`**
```
ViroForge Complete Workflow Summary
======================================================================

Pipeline Configuration:
  Body site:         Gut
  Viral genomes:     50
  Contamination:     Realistic
  VLP enrichment:    Standard (0.2 μm, 95% nuclease)
  Amplification:     RdAB (40 cycles)
  Platform:          NovaSeq 6000

Pipeline Stages:
1. Initial Composition
   Viral fraction:    50.0%

2. After VLP Enrichment
   Viral fraction:    97.1%
   Enrichment:        1.94x

3. After Amplification
   Viral fraction:    100.0%

4. Platform Artifacts
   PolyG tails:       2.5%
   Optical dups:      9.0%
   Index hopping:     1.5%
```

---

## Documentation

### User Documentation
- **[FASTQ Generation Guide](docs/PHASE4_FASTQ_GENERATION.md)** - Generate datasets from curated collections
- **[Collection Implementation Guide](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)** - How collections were curated
- **[Database Design](docs/GENOME_DATABASE_DESIGN.md)** - RefSeq genome database schema
- **[User Guide](docs/USER_GUIDE.md)** - Comprehensive usage guide (legacy API)
- **[Tutorial](docs/TUTORIAL.md)** - Step-by-step walkthrough (legacy API)

### Body Site Collections
- **[Collection Overview](docs/BODY_SITE_COLLECTIONS.md)** - Original curation plans
- **[Gut Virome Curation](docs/GUT_VIROME_CURATION.md)** - Detailed gut virome plan
- **[Curation Workflow](docs/BODY_SITE_CURATION_WORKFLOW.md)** - Technical workflow

### Developer Documentation
- **[Design Rationale](docs/DESIGN_RATIONALE.md)** - Design decisions and literature review
- **[Implementation Plan](docs/IMPLEMENTATION_PLAN.md)** - Phased development roadmap
- **[VLP Biology](docs/VLP_ENRICHMENT_BIOLOGY.md)** - VLP enrichment biology guide
- **[API Reference](docs/API.md)** - Python API documentation (legacy)

### Scripts Documentation
- **[Exploration Tools](scripts/README_EXPLORATION_TOOLS.md)** - Database exploration utilities
- **[Helper Utilities](scripts/README_HELPER_UTILITIES.md)** - Supporting scripts

---

## Project Status

**Current Version**: 0.3.0

**Phase 4: Complete** | **Production Ready**

### Completed Phases

**Phase 1: Core Simulator (Complete)**
- Viral community composition (5 body sites, 3 abundance distributions)
- Contamination profiles (realistic host, bacterial, fungal DNA)
- FASTQ generation with ground truth tracking
- Comprehensive validation framework
- 158 tests passing

**Phase 2: Virome-Specific Features (Complete)**
- VLP Enrichment Framework (40 tests)
  - Filtration models, nuclease treatment
  - 4 pre-defined protocols
- Amplification Bias Framework (31 tests)
  - RdAB, MDA, linker amplification
  - Length and GC bias modeling
- Platform Artifact Framework (33 tests)
  - PolyG tails, optical duplicates, index hopping
  - 5 platform profiles (NovaSeq, MiSeq, HiSeq)
- Integration & Workflows (20 tests)
  - End-to-end pipelines, cross-platform comparisons

**Phase 3: Genome Database & Collections (Complete)**
- RefSeq viral genome database (14,423 genomes)
- ICTV taxonomy integration (53.9% coverage)
- 8 curated body site collections (1,198 genomes):
  - Gut, Oral, Skin, Respiratory (human)
  - Marine, Soil, Freshwater (environmental)
  - Mouse Gut (model organism)
- Literature-validated compositions
- Automated curation workflows

**Phase 4: FASTQ Generation (Complete)**
- Database-driven FASTQ generation
- InSilicoSeq integration for realistic reads
- VLP enrichment simulation
- Platform-specific error models
- Complete ground truth metadata export
- Batch generation with presets
- Comprehensive documentation

### Test Coverage

```
Total Tests:        178 passing (100%)
Unit Tests:         158
Integration Tests:  20
Collections:        8 body sites
Genomes:            14,423 RefSeq viral genomes
Coverage:           Comprehensive
```

### Implementation Status

- **Phase 1** - Core functionality (Complete)
- **Phase 2** - Virome-specific features (Complete)
- **Phase 3** - Genome database & collections (Complete)
- **Phase 4** - FASTQ generation (Complete)
- **Phase 5** - Publication & community release (Planned)

---

## Literature Validation

All ViroForge parameters are validated against peer-reviewed literature:

**VLP Enrichment**
- Shkoporov & Hill (2019) *Nat Rev Microbiol* - VLP protocols
- Zolfo et al. (2019) *Microbiome* - ViromeQC metrics

**Amplification Bias**
- Kim et al. (2013) *Nat Methods* - Amplification bias characterization
- Marine et al. (2014) *PeerJ* - Transposase-based protocols
- Duhaime et al. (2012) *Environ Microbiol* - MDA artifacts

**Platform Artifacts**
- Costello et al. (2018) *BMC Genomics* - Index swapping
- Chen et al. (2017) Illumina Technical Note - NovaSeq chemistry
- Sinha et al. (2017) *Genome Res* - Index switching rates

See `docs/VLP_ENRICHMENT_BIOLOGY.md` for detailed literature review.

---

## Contributing

We welcome contributions! ViroForge is open-source (MIT License) and community-driven.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Write tests for your changes
4. Ensure all tests pass (`pytest tests/`)
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

### Development Setup

```bash
# Clone and install in development mode
git clone https://github.com/shandley/viroforge.git
cd viroforge
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run specific test suite
pytest tests/test_enrichment.py -v
pytest tests/test_integration_workflow.py -v
```

---

## Citation

### Software

```bibtex
@software{viroforge2025,
  title = {ViroForge: A Synthetic Virome Data Generator},
  author = {Handley, Scott and contributors},
  year = {2025},
  url = {https://github.com/shandley/viroforge},
  version = {0.2.0-dev}
}
```

### Publication

Manuscript in preparation. Check back soon for preprint!

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

ViroForge was developed to support virome analysis pipeline validation and benchmarking in the microbiome research community.

### Related Projects

- **[Hecatomb](https://github.com/shandley/hecatomb)** - Viral metagenome assembly pipeline
- **[ViromeQC](https://github.com/SegataLab/viromeqc)** - Virome enrichment assessment
- **[CAMISIM](https://github.com/CAMI-challenge/CAMISIM)** - Metagenome simulator
- **[InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)** - Illumina read simulator (integrated)

### Team

**Principal Investigator**: Scott Handley
**Institution**: Washington University in St. Louis
**Lab Website**: [Handley Lab](https://www.handleylab.org)

---

## Support

- 📖 **Documentation**: Check `docs/` directory
- 💬 **Questions**: Open a [GitHub Discussion](https://github.com/shandley/viroforge/discussions)
- 🐛 **Bug Reports**: Open a [GitHub Issue](https://github.com/shandley/viroforge/issues)
- 📧 **Email**: scott.handley@wustl.edu

---

## Status

**Production Ready** - Phase 4 Complete

ViroForge is ready for use in benchmarking studies. All core functionality is complete and thoroughly tested:

- 8 curated body site collections with literature-validated compositions
- 14,423 RefSeq viral genomes with ICTV taxonomy
- Database-driven FASTQ generation with complete ground truth
- Platform-specific error models (NovaSeq, MiSeq, HiSeq)
- VLP enrichment simulation
- Comprehensive documentation

**Last Updated**: 2025-11-01
