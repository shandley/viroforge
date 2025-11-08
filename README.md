# ViroForge

**Forging synthetic viromes for benchmarking and validation**

A comprehensive mock metavirome data generator for testing and validating virome analysis pipelines.

[![Tests](https://img.shields.io/badge/tests-28%20passing-brightgreen)](tests/)
[![Phase](https://img.shields.io/badge/Phase%205-Complete-success)](lab-notebook/sessions/2025-11/)
[![Collections](https://img.shields.io/badge/collections-8%20body%20sites-blue)](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)
[![Genomes](https://img.shields.io/badge/genomes-14%2C423%20RefSeq-blue)](docs/GENOME_DATABASE_DESIGN.md)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## Overview

ViroForge generates realistic synthetic virome sequencing datasets with complete ground truth metadata, enabling rigorous validation of QC pipelines, assembly tools, taxonomic classifiers, and analysis workflows.

**What makes ViroForge different?** It's the first simulator to model VLP enrichment with size-based filtration and type-specific contamination reduction, providing realistic virome datasets for benchmarking.

### Key Features

- **8 Curated Body Site Collections** - Literature-validated virome compositions (gut, oral, skin, respiratory, marine, soil, freshwater, mouse gut)
- **14,423 RefSeq Viral Genomes** - Complete database with ICTV taxonomy integration (53.9% coverage)
- **Database-Driven FASTQ Generation** - Direct generation from curated collections with InSilicoSeq
- **Enhanced VLP Enrichment Modeling** - Size-based filtration with 5 protocols (tangential flow, syringe, ultracentrifugation, Norgen, bulk)
- **Type-Specific Contamination Reduction** - Host DNA, rRNA, bacteria, PhiX with protocol-dependent efficiency
- **Complete Ground Truth** - Taxonomic composition, abundance tables, genome-read mappings, contaminant sequences
- **Platform-Specific Error Models** - NovaSeq, MiSeq, HiSeq with realistic artifacts
- **Reproducible Benchmarks** - Random seeds, complete metadata, known composition
- **Production Ready** - Comprehensive testing, literature-validated parameters

---

## Why ViroForge?

### The Problem

Current virome analysis tool validation approaches are limited:
- **Physical synthetic communities**: Expensive ($10k+), time-consuming, limited complexity
- **Existing simulators**: Bacterial-focused (CAMISIM), don't model VLP enrichment
- **Real datasets**: Unknown ground truth, can't systematically test edge cases

### The Solution

ViroForge enables:
- Unlimited synthetic datasets with complete ground truth
- VLP enrichment vs bulk metagenome comparisons
- Realistic contamination profiles with protocol-dependent reduction
- Size-based viral filtration modeling
- Cross-platform reproducibility testing (NovaSeq, MiSeq, HiSeq)
- Standardized benchmarking datasets for the community

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

### Generate FASTQ from Curated Collections

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

**Step 2: Generate FASTQ Dataset with VLP Enrichment**

```bash
# Generate gut virome with tangential flow VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --coverage 10 \
    --platform novaseq \
    --vlp-protocol tangential_flow
```

Output:
```
output/
├── fasta/
│   └── collection_9.fasta  # Reference genomes with abundances
├── fastq/
│   ├── collection_9_R1.fastq  # Forward reads
│   └── collection_9_R2.fastq  # Reverse reads
└── metadata/
    ├── collection_9_metadata.json  # Complete ground truth (viral + contaminants)
    ├── collection_9_composition.tsv  # Abundance table
    └── collection_9_abundances.txt  # ISS abundance file
```

**Step 3: Use with Hecatomb or Other Pipelines**

```bash
# Run Hecatomb on generated dataset
hecatomb run \
    --reads data/fastq/gut_virome/fastq/*_R{1,2}.fastq \
    --outdir results/gut_benchmark

# Compare results to ground truth
# See metadata/collection_9_metadata.json for known composition
```

See [FASTQ Generation Guide](scripts/README_FASTQ_GENERATION.md) for detailed documentation.

---

## Use Cases

### 1. Benchmark Virome Analysis Pipelines

```bash
# Generate gut virome benchmark dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output benchmarks/gut_10x \
    --coverage 10 \
    --vlp-protocol tangential_flow

# Run your pipeline
your_pipeline benchmarks/gut_10x/fastq/*_R{1,2}.fastq

# Compare results to ground truth
# See benchmarks/gut_10x/metadata/collection_9_metadata.json
```

### 2. Compare VLP Protocols

```bash
# Generate datasets with different VLP protocols
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output data/vlp_comparison

# Compares: tangential_flow, syringe, ultracentrifugation, norgen, bulk
```

### 3. VLP vs Bulk Metagenome Comparison

```bash
# Generate with VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_vlp \
    --coverage 10 \
    --vlp-protocol tangential_flow

# Generate bulk metagenome (no VLP)
python scripts/generate_fastq_dataset.py \
    --collection-id 13 \
    --output marine_bulk \
    --coverage 10 \
    --no-vlp

# Compare viral recovery rates and contamination levels
```

### 4. Cross-Platform Reproducibility Testing

```bash
# Generate datasets on different platforms
python scripts/generate_fastq_dataset.py --collection-id 9 --platform novaseq --output novaseq_gut
python scripts/generate_fastq_dataset.py --collection-id 9 --platform miseq --output miseq_gut
python scripts/generate_fastq_dataset.py --collection-id 9 --platform hiseq --output hiseq_gut

# Compare platform-specific artifacts and assembly quality
```

### 5. Batch Generation for Comprehensive Benchmarks

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

See [scripts/README_FASTQ_GENERATION.md](scripts/README_FASTQ_GENERATION.md) for complete documentation.

---

## VLP Enrichment Protocols

ViroForge models 5 VLP enrichment protocols with realistic size-based filtration and contamination reduction:

| Protocol | Method | Filtration | Contamination Reduction | Use Case |
|----------|--------|------------|-------------------------|----------|
| `tangential_flow` | 0.2 μm TFF | Size-based (sigmoid) | 93.5% | Most common protocol |
| `syringe` | 0.22 μm syringe | Size-based (step) | 90.2% | Field-friendly |
| `ultracentrifugation` | Density gradient | Mild size bias | 91.8% | Traditional method |
| `norgen` | Column-based | Modest size bias | 89.1% | Commercial kit |
| `none` (--no-vlp) | Bulk metagenome | None | 0% | Control/comparison |

**Features:**
- Virion size estimation from genome length and type (dsDNA, ssDNA, ssRNA, dsRNA)
- Protocol-specific filtration curves
- Type-specific contamination reduction:
  - Host DNA: 95-99% (DNase treatment)
  - rRNA: 90-98% (size-based removal)
  - Bacteria: 85-95% (filtration)
  - PhiX: 0-20% (treated as small virus)
- Literature-validated parameters (Lim et al. 2020, Thurber et al. 2009, Reyes et al. 2012)

See [VLP Integration Guide](docs/VLP_CONTAMINATION_INTEGRATION.md) for details.

---

## Output Files

### Ground Truth Metadata

Every ViroForge dataset includes complete ground truth metadata:

**`metadata.json`** - Complete information including:
- Collection metadata (ID, name, environment)
- Configuration (coverage, platform, VLP protocol, contamination level)
- All sequences (viral genomes + contaminants):
  - Genome ID, name, type (viral/contaminant)
  - Length, GC content
  - Relative abundance
  - Taxonomy (family, genus, species)
- Enrichment statistics (if VLP applied)
- Generation timestamp and ViroForge version

**`composition.tsv`** - Tab-separated abundance table
```
genome_id    genome_name              length    gc_content    relative_abundance    family
NC_001416    Enterobacteria phage T7  39937     48.5          0.1234                Podoviridae
NC_007458    Escherichia phage MS2    3569      51.8          0.0567                Leviviridae
...
```

**`abundances.txt`** - InSilicoSeq abundance file

See [Phase 4 Documentation](docs/PHASE4_FASTQ_GENERATION.md) for output format details.

---

## Documentation

### User Guides
- **[Quick Start](docs/QUICKSTART.md)** - Quick start for script-based workflow
- **[FASTQ Generation Guide](scripts/README_FASTQ_GENERATION.md)** - Generate datasets from collections
- **[Database Exploration Tools](scripts/README_EXPLORATION_TOOLS.md)** - Explore genome database

### Technical Documentation
- **[Collection Implementation Guide](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)** - Collection curation rationale
- **[Genome Database Design](docs/GENOME_DATABASE_DESIGN.md)** - Database schema
- **[Phase 4: FASTQ Generation](docs/PHASE4_FASTQ_GENERATION.md)** - FASTQ workflow
- **[Phase 5: VLP Integration](docs/PHASE5_TASK2_FASTQ_INTEGRATION.md)** - VLP enrichment integration
- **[Phase 5: Validation Report](docs/PHASE5_TASK3_VALIDATION_REPORT.md)** - Comprehensive validation
- **[VLP Integration Guide](docs/VLP_CONTAMINATION_INTEGRATION.md)** - Usage guide
- **[Validation Test Suite](docs/VALIDATION_TEST_SUITE.md)** - Pipeline validation

### API Documentation (Legacy)
- **[Tutorial](docs/TUTORIAL.md)** - Step-by-step programmatic API tutorial
- **[User Guide](docs/USER_GUIDE.md)** - Comprehensive programmatic API guide
- **[API Reference](docs/API.md)** - Python API documentation

**Note**: Legacy API (Phases 1-2) allows programmatic access to individual components. Current workflow uses script-based interface with curated collections.

---

## Project Status

**Current Version**: 0.4.0

**Phase 5: Complete** | **Production Ready**

### Completed Phases

**Phase 1: Core Simulator (Complete)**
- Viral community composition
- Contamination profiles (host DNA, bacteria, rRNA, PhiX)
- FASTQ generation with ground truth tracking
- Comprehensive validation framework

**Phase 2: Virome-Specific Features (Complete)**
- VLP enrichment framework (basic implementation)
- Amplification bias framework (RdAB, MDA, linker)
- Platform artifact framework (polyG, optical duplicates, index hopping)
- Integration and complete workflows

**Phase 3: Genome Database & Collections (Complete)**
- RefSeq viral genome database (14,423 genomes)
- ICTV taxonomy integration (53.9% coverage)
- 8 curated body site collections (1,198 genomes)
- Literature-validated compositions
- Automated curation workflows

**Phase 4: FASTQ Generation (Complete)**
- Database-driven FASTQ generation
- InSilicoSeq integration for realistic reads
- Platform-specific error models
- Complete ground truth metadata export
- Batch generation with presets

**Phase 5: Enhanced VLP Modeling (Complete)**
- Size-based filtration modeling (virion diameter from genome properties)
- 5 VLP protocols with literature-validated parameters
- Type-specific contamination reduction (host DNA, rRNA, bacteria, PhiX)
- Protocol-dependent efficiency modeling
- Comprehensive testing and validation
- Complete integration with FASTQ generation

### Test Coverage

```
Unit Tests:             16/16 passing (VLP/contamination)
Integration Tests:      12 tests (3 dry-run passing, 9 require ISS)
Collection Validation:  8/8 collections validated
Genomes:                14,423 RefSeq viral genomes
Literature Validation:  5/5 metrics validated
```

---

## Literature Validation

All ViroForge parameters are validated against peer-reviewed literature:

**VLP Enrichment and Contamination Reduction**
- Lim et al. (2020) - VLP protocol comparison
- Thurber et al. (2009) - Ultracentrifugation efficiency
- Reyes et al. (2012) - DNase treatment efficiency
- Kim et al. (2015) - Host DNA contamination reduction
- Solonenko et al. (2013) - VLP enrichment factors

**Virion Size Relationships**
- Cui et al. (2014) - Genome length to virion size
- Nasir et al. (2017) - Virion size distributions
- Danovaro et al. (2011) - Marine viral particle sizes

See [Phase 5 Validation Report](docs/PHASE5_TASK3_VALIDATION_REPORT.md) for detailed validation results.

---

## Contributing

We welcome contributions! ViroForge is open-source (MIT License) and community-driven.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Write tests for your changes
4. Ensure all tests pass (`pytest tests/ -v`)
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

# Run specific test suites
pytest tests/test_vlp_contamination.py -v
pytest tests/test_fastq_integration.py -v --tb=short
```

---

## Citation

### Software

```bibtex
@software{viroforge2025,
  title = {ViroForge: A Synthetic Virome Data Generator with Enhanced VLP Modeling},
  author = {Handley, Scott and contributors},
  year = {2025},
  url = {https://github.com/shandley/viroforge},
  version = {0.4.0}
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

- **Documentation**: Check `docs/` directory
- **Questions**: Open a [GitHub Discussion](https://github.com/shandley/viroforge/discussions)
- **Bug Reports**: Open a [GitHub Issue](https://github.com/shandley/viroforge/issues)
- **Email**: scott.handley@wustl.edu

---

## Status

**Production Ready** - Phase 5 Complete

ViroForge is ready for use in benchmarking studies. All core functionality is complete and thoroughly tested:

- 8 curated body site collections with literature-validated compositions
- 14,423 RefSeq viral genomes with ICTV taxonomy
- Enhanced VLP enrichment modeling with size-based filtration
- Type-specific contamination reduction
- Database-driven FASTQ generation with complete ground truth
- Platform-specific error models (NovaSeq, MiSeq, HiSeq)
- Comprehensive documentation and validation

**Last Updated**: 2025-11-08
