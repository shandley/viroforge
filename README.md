# ViroForge

**Forging synthetic viromes for benchmarking and validation**

A comprehensive mock metavirome data generator for testing and validating virome analysis pipelines, now with **RNA virome support** and critical disease/environmental collections.

[![Tests](https://img.shields.io/badge/tests-70%2B%20passing-brightgreen)](tests/)
[![Phase](https://img.shields.io/badge/Phase%208.2-Complete-success)](ROADMAP.md)
[![Collections](https://img.shields.io/badge/collections-23%20curated-blue)](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)
[![Genomes](https://img.shields.io/badge/genomes-14%2C423%20RefSeq-blue)](docs/GENOME_DATABASE_DESIGN.md)
[![Taxonomy](https://img.shields.io/badge/taxonomy-57.1%25%20ICTV-blue)](docs/TAXONOMY_BUG_FIX.md)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## 🆕 What's New in v0.6.0

### Phase 8: RNA Virome Workflow (November 2025)
- **✨ RNA virome support** - Complete workflow for RNA virus sequencing
  - Reverse transcription modeling with virus-type specific efficiency (40-90%)
  - rRNA depletion (Ribo-Zero/RiboMinus) with 10-20x viral enrichment
  - RNA degradation and fragmentation modeling
- **🧬 3 RNA virus collections** - Respiratory, arbovirus, fecal RNA viromes
- **🔬 RNA-specific contamination** - Realistic rRNA profiles (90% → 10% post-depletion)

### Phase 7: Critical Collections & Taxonomy Fix (November 2025)
- **🏥 Disease state collections** - IBD, HIV+, CF respiratory viromes
- **🌊 Wastewater virome** - Epidemiological surveillance applications
- **🐛 Major taxonomy bug fixed** - Enhanced fuzzy matching fixed 469 genomes (7.1% of unmatched)
  - **CRITICAL**: HIV+ collection now includes herpesviruses (EBV, KSHV)
  - **MAJOR**: Fecal RNA collection +81% size with rotavirus/norovirus
- **📊 23 total collections** - From 8 to 23 curated virome collections

See [Taxonomy Bug Fix Documentation](docs/TAXONOMY_BUG_FIX.md) for complete details.

---

## Overview

ViroForge generates realistic synthetic virome sequencing datasets with complete ground truth metadata, enabling rigorous validation of QC pipelines, assembly tools, taxonomic classifiers, and analysis workflows.

**What makes ViroForge unique?**
- First simulator to model VLP enrichment with size-based filtration
- **NEW**: First simulator with RNA virome workflow (RT, rRNA depletion, degradation)
- Complete ground truth for both DNA and RNA viromes
- Disease state and environmental surveillance collections
- Enhanced taxonomy with fuzzy matching for strain-specific names

### Key Features

- **23 Curated Virome Collections** - Literature-validated compositions:
  - **Host-associated (15)**: Healthy gut/oral/skin/respiratory, disease states (IBD, HIV+, CF), VLP comparisons
  - **Environmental (5)**: Marine, soil, freshwater, wastewater
  - **RNA viromes (3)**: Respiratory RNA, arbovirus, fecal RNA
- **14,423 RefSeq Viral Genomes** - Complete database with enhanced ICTV taxonomy (57.1% coverage after fix)
- **DNA & RNA Virome Workflows** - Complete support for both molecule types
  - DNA: VLP enrichment, amplification bias, sequencing artifacts
  - RNA: Reverse transcription, rRNA depletion (Ribo-Zero), RNA degradation
- **Enhanced VLP Enrichment Modeling** - Size-based filtration with 5 protocols
- **Type-Specific Contamination Reduction** - DNA and RNA-specific profiles
- **Complete Ground Truth** - Taxonomic composition, abundance tables, workflow statistics
- **Platform-Specific Error Models** - NovaSeq, MiSeq, HiSeq with realistic artifacts
- **Reproducible Benchmarks** - Random seeds, complete metadata, known composition
- **Production Ready** - 70+ comprehensive tests, literature-validated parameters

---

## Why ViroForge?

### The Problem

Current virome analysis tool validation approaches are limited:
- **Physical synthetic communities**: Expensive ($10k+), time-consuming, limited complexity
- **Existing simulators**: Bacterial-focused (CAMISIM), don't model VLP enrichment or RNA workflows
- **Real datasets**: Unknown ground truth, can't systematically test edge cases
- **RNA viromes**: No existing tools model RT efficiency, rRNA depletion, or RNA degradation

### The Solution

ViroForge enables:
- Unlimited synthetic datasets with complete ground truth (DNA and RNA)
- VLP enrichment vs bulk metagenome comparisons
- RNA virome workflows with Ribo-Zero modeling (90% rRNA → 10%)
- Realistic contamination profiles with protocol-dependent reduction
- Disease state and environmental surveillance benchmarking
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

### Generate DNA Virome Dataset

**Step 1: List Available Collections**

```bash
python scripts/generate_fastq_dataset.py --list-collections
```

**Step 2: Generate FASTQ Dataset with VLP Enrichment**

```bash
# Generate gut virome with tangential flow VLP enrichment
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/fastq/gut_virome \
    --molecule-type dna \
    --coverage 10 \
    --platform novaseq \
    --vlp-protocol tangential_flow
```

### Generate RNA Virome Dataset ✨ NEW

```bash
# Generate respiratory RNA virome with Ribo-Zero depletion
python scripts/generate_fastq_dataset.py \
    --collection-id 21 \
    --output data/fastq/respiratory_rna \
    --molecule-type rna \
    --rna-primer random_hexamer \
    --rna-depletion ribo_zero \
    --coverage 10 \
    --platform novaseq
```

**RNA workflow includes:**
- Reverse transcription with virus-type specific efficiency
- rRNA depletion (90% → 10% with Ribo-Zero)
- RNA degradation and fragmentation
- RNA-specific contamination profiles

Output:
```
output/
├── fasta/
│   └── collection_21.fasta  # Reference genomes
├── fastq/
│   ├── collection_21_R1.fastq  # Forward reads
│   └── collection_21_R2.fastq  # Reverse reads
└── metadata/
    ├── collection_21_metadata.json  # Complete ground truth + RNA workflow stats
    ├── collection_21_composition.tsv  # Abundance table
    └── collection_21_abundances.txt  # ISS abundance file
```

See [FASTQ Generation Guide](scripts/README_FASTQ_GENERATION.md) for detailed documentation.

---

## Use Cases

### 1. Benchmark RNA Virome Analysis Pipelines ✨ NEW

```bash
# Generate respiratory RNA virome with complete workflow
python scripts/generate_fastq_dataset.py \
    --collection-id 21 \
    --output benchmarks/respiratory_rna_10x \
    --molecule-type rna \
    --rna-depletion ribo_zero \
    --coverage 10

# Run your pipeline
your_pipeline benchmarks/respiratory_rna_10x/fastq/*_R{1,2}.fastq

# Compare to ground truth (includes RT efficiency, rRNA depletion stats)
```

### 2. Compare Ribo-Zero Efficiency ✨ NEW

```bash
# Generate with Ribo-Zero
python scripts/generate_fastq_dataset.py \
    --collection-id 21 \
    --output rna_with_ribozero \
    --molecule-type rna \
    --rna-depletion ribo_zero

# Generate without Ribo-Zero (failed depletion)
python scripts/generate_fastq_dataset.py \
    --collection-id 21 \
    --output rna_no_ribozero \
    --molecule-type rna \
    --rna-depletion none

# Compare: 90% rRNA vs 10% rRNA, 10-20x viral enrichment difference
```

### 3. Benchmark Disease State Detection

```bash
# Generate healthy vs disease state collections
python scripts/generate_fastq_dataset.py --collection-id 9  --output healthy_gut  # Healthy
python scripts/generate_fastq_dataset.py --collection-id 18 --output ibd_gut      # IBD
python scripts/generate_fastq_dataset.py --collection-id 19 --output hiv_gut      # HIV+

# Test pipeline's ability to detect virome dysbiosis
```

### 4. Wastewater Surveillance Benchmarking ✨ NEW

```bash
# Generate wastewater virome for epidemiological surveillance
python scripts/generate_fastq_dataset.py \
    --collection-id 17 \
    --output wastewater_surveillance \
    --coverage 10

# Test pathogen detection (SARS-CoV-2, rotavirus, norovirus, etc.)
```

### 5. VLP vs Bulk Metagenome Comparison

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

---

## Virome Collections

ViroForge includes **23 curated virome collections** representing diverse viral ecosystems:

### Original Collections (1-8)
1. **Healthy Human Gut** (134 genomes) - Western diet, adult
2. **Healthy Human Skin** (42 genomes) - Sebaceous sites
3. **Healthy Human Oral** (67 genomes) - Saliva
4. **Healthy Human Urogenital** (31 genomes) - Vaginal virome
5. **Healthy Human Respiratory** (58 genomes) - Nasopharynx
6. **Marine Virome** (78 genomes) - Coastal surface water
7. **Soil Virome** (82 genomes) - Agricultural
8. **Freshwater Virome** (71 genomes) - Lake surface water

### VLP Comparison Collections (9-15)
9-15. **VLP Protocol Comparisons** - Test different enrichment methods
- Baseline, high/low bacterial lysis, prophage induction, eukaryotic virus shedding

### Amplification Comparison (16)
16. **Pre-Amplification Control** (100 genomes) - Method comparison baseline

### Critical Collections - Phase 7 ✨ NEW (17-20)
17. **Wastewater Virome** (352 genomes) - Epidemiological surveillance
    - Enteric viruses, bacteriophages, emerging pathogens (SARS-CoV-2, mpox)
18. **IBD Gut Virome** (90 genomes) - Inflammatory bowel disease dysbiosis
    - Reduced diversity, altered Caudovirales, increased temperate phages
19. **HIV+ Gut Virome** (55 genomes) - HIV-associated gut dysbiosis
    - **CRITICAL FIX**: Now includes human herpesviruses (EBV, KSHV, HSV-1, VZV)
20. **CF Respiratory Virome** (81 genomes) - Cystic fibrosis lung
    - Pseudomonas/Staphylococcus phages, respiratory viruses (influenza, RSV)

### RNA Virome Collections - Phase 8 ✨ NEW (21-23)
21. **Human Respiratory RNA Virome** (56 genomes) - RNA respiratory viruses
    - Influenza, RSV, coronaviruses, rhinoviruses, enteroviruses
22. **Arbovirus Environmental** (39 genomes) - Mosquito-associated RNA viruses
    - Flaviviruses, alphaviruses, bunyaviruses
23. **Fecal RNA Virome** (58 genomes) - Enteric RNA viruses
    - **MAJOR FIX**: Now includes rotavirus (12) and norovirus (15) - +81% size

**Total**: 1,444 genomes across 23 diverse environments, all with literature-validated compositions.

See [Collection Implementation Guide](docs/COLLECTION_IMPLEMENTATION_GUIDE.md) for detailed curation rationale.

---

## RNA Virome Workflow ✨ NEW

ViroForge is the first simulator to model complete RNA virome workflows:

### Reverse Transcription
- **Virus-type specific efficiency**:
  - ssRNA+ (positive sense): 70-90% (norovirus, poliovirus, coronaviruses)
  - ssRNA- (negative sense): 50-70% (influenza, RSV, paramyxoviruses)
  - dsRNA: 40-80% (rotavirus, reovirus)
- **RT artifacts**: Template switching (~2%), 5'/3' truncation (~15%)
- **Primer types**: Random hexamer (default), random octamer, oligo-dT, specific

### rRNA Depletion (Ribo-Zero/RiboMinus)
- **Critical for RNA viromes**: 80-95% of RNA is rRNA (vs ~5% for DNA)
- **Depletion efficiency**: 90-95% removal (Ribo-Zero), 85-90% (RiboMinus)
- **Viral enrichment**: 10-20x increase in viral reads
- **Before depletion**: 90% rRNA, 1% viral (unusable)
- **After depletion**: 10% rRNA, 20% viral (excellent)

### RNA Degradation
- **10-100x faster than DNA**: RNase contamination, chemical instability
- **Fragmentation**: 2-4 fragments per degraded sequence
- **5'/3' bias**: 5' end more degraded, uneven coverage

### RNA-Specific Contamination
- **Host RNA**: 90% rRNA before depletion → 10% after
- **Bacterial RNA**: 16S/23S rRNA + mRNA from microbiome
- **Dramatically different from DNA**: >10x more contamination without Ribo-Zero

**Command-line flags**:
```bash
--molecule-type {dna,rna}                    # Select workflow
--rna-primer {random_hexamer,random_octamer,oligo_dt,specific}
--rna-depletion {ribo_zero,ribominus,none}   # rRNA depletion method
```

See [RNA Virome Workflow Documentation](viroforge/workflows/rna_virome.py) for technical details.

---

## VLP Enrichment Protocols

ViroForge models 5 VLP enrichment protocols with realistic size-based filtration and contamination reduction:

| Protocol | Method | Filtration | Contamination Reduction | Viral Recovery | Use Case |
|----------|--------|------------|-------------------------|----------------|----------|
| `tangential_flow` | 0.2 μm TFF | Size-based (sigmoid) | 91.2% | 85% | Highest purity |
| `ultracentrifugation` | Density gradient | Minimal size bias | 88.4% | 90% | Highest recovery |
| `norgen` | Column-based | Modest size bias | 87.1% | 70% | Convenient |
| `syringe` | 0.22 μm syringe | Size-based (step) | 85.7% | 60% | Field-friendly |
| `none` (--no-vlp) | Bulk metagenome | None | 0% | 100% | Control/comparison |

**Features:**
- Virion size estimation from genome length and type (dsDNA, ssDNA, ssRNA, dsRNA)
- Protocol-specific filtration curves
- Type-specific contamination reduction:
  - Host DNA/RNA: 95-99% (DNase treatment)
  - rRNA: 90-98% (size-based removal + Ribo-Zero for RNA)
  - Bacteria: 85-95% (filtration)
  - PhiX: 0-20% (treated as small virus)

### Contamination Levels

| Level | DNA Initial | RNA Initial (no Ribo-Zero) | After VLP (TFF) | Use Case |
|-------|------------|---------------------------|-----------------|----------|
| **Clean** | ~0.7% | ~92% | ~0.1% / ~5% | Optimal prep |
| **Realistic** | ~7.4% | ~95% | ~0.7% / ~10% | Typical (default) |
| **Heavy** | ~27.1% | ~97% | ~2.0% / ~20% | Failed prep |

---

## Taxonomy Bug Fix & Enhancement 🐛

### Problem Discovered (November 2025)

**46% of database (6,651/14,423 genomes) had `family='Unknown'`** due to mismatches between:
- **RefSeq**: Strain-specific names (e.g., "Influenza A virus (A/California/07/2009(H1N1))")
- **ICTV**: General species names (e.g., "influenza A virus")

### Impact on Collections

**CRITICAL - Collection 19 (HIV+ Gut)**:
- Had **ZERO herpesviruses** (scientifically invalid - HIV+ patients show herpesvirus reactivation)
- Fixed: 0 → 6 human herpesviruses (EBV, KSHV, HSV-1, VZV, HHV-6B, HHV-7)

**MAJOR - Collection 23 (Fecal RNA)**:
- Only 32 genomes, missing rotavirus/norovirus (primary enteric pathogens)
- Fixed: 32 → 58 genomes (**+81%**), now has 15 norovirus + 12 rotavirus

**Collections 17 (Wastewater) & 20 (CF Respiratory)**: Also fixed

### Solution Applied

Enhanced `scripts/fix_taxonomy_unmatched.py` with:
1. **Pattern-based family matching** (20+ virus families)
2. **Improved normalization** (remove "type", "strain", trailing numbers)
3. **Fuzzy matching** for strain-specific nomenclature

**Results**: Fixed **469 genomes (7.1% of unmatched)**

**Current taxonomy coverage**: 8,241/14,423 genomes (57.1%) assigned

See [Taxonomy Bug Fix Documentation](docs/TAXONOMY_BUG_FIX.md) for complete details, lessons learned, and prevention strategies.

---

## Output Files

### Ground Truth Metadata

Every ViroForge dataset includes complete ground truth metadata:

**`metadata.json`** - Complete information including:
- Collection metadata (ID, name, environment)
- Configuration (coverage, platform, VLP protocol, contamination level, **molecule type**)
- All sequences (viral genomes + contaminants):
  - Genome ID, name, type (viral/contaminant)
  - Length, GC content, relative abundance
  - Taxonomy (family, genus, species)
- Enrichment statistics (VLP and/or **RNA workflow stats**)
- **RNA-specific stats** (if applicable):
  - RT efficiency by virus type
  - rRNA depletion efficiency and viral enrichment
  - RNA degradation statistics
  - Overall recovery rate

**`composition.tsv`** - Tab-separated abundance table

**`abundances.txt`** - InSilicoSeq abundance file

---

## Documentation

### Getting Started
- **[Quick Start](docs/QUICKSTART.md)** - Generate your first dataset in 5 minutes
- **[FASTQ Generation Guide](scripts/README_FASTQ_GENERATION.md)** - Complete command-line reference
- **[VLP Protocol Comparison Tutorial](docs/VLP_PROTOCOL_COMPARISON_TUTORIAL.md)** - Compare all 5 VLP protocols
- **[RNA Virome Workflow](viroforge/workflows/rna_virome.py)** - RNA virome technical details ✨ NEW

### Critical Documentation ⚠️
- **[Taxonomy Bug Fix](docs/TAXONOMY_BUG_FIX.md)** - **READ THIS** if working with taxonomy ✨ NEW
- **[Collection Implementation Guide](docs/COLLECTION_IMPLEMENTATION_GUIDE.md)** - All 23 collections documented

### Technical Documentation
- **[Genome Database Design](docs/GENOME_DATABASE_DESIGN.md)** - Database schema
- **[Phase 4: FASTQ Generation](docs/PHASE4_FASTQ_GENERATION.md)** - Complete FASTQ workflow
- **[Phase 5: VLP Integration](docs/PHASE5_TASK2_FASTQ_INTEGRATION.md)** - VLP enrichment integration
- **[VLP Integration Guide](docs/VLP_CONTAMINATION_INTEGRATION.md)** - Technical implementation
- **[Validation Test Suite](docs/VALIDATION_TEST_SUITE.md)** - Pipeline validation

---

## Project Status

**Current Version**: 0.6.0

**Phase 8.2: Complete** | **Production Ready**

### Completed Phases

**Phase 1-2: Core Simulator (Complete)**
- Viral community composition
- Contamination profiles
- FASTQ generation with ground truth
- Amplification bias and platform artifacts

**Phase 3-4: Database & FASTQ Generation (Complete)**
- RefSeq viral genome database (14,423 genomes)
- ICTV taxonomy integration
- Database-driven FASTQ generation
- InSilicoSeq integration

**Phase 5: Enhanced VLP Modeling (Complete)**
- Size-based filtration modeling
- 5 VLP protocols with literature validation
- Type-specific contamination reduction

**Phase 6: Amplification Bias Integration (Complete)**
- RdAB, MDA, Linker amplification methods
- Amplification comparison tools

**Phase 7: Critical Collections (Complete)** ✨
- Wastewater virome (epidemiological surveillance)
- Disease state collections (IBD, HIV+, CF)
- Progressive dysbiosis modeling
- **Taxonomy bug discovery and fix (469 genomes fixed)**

**Phase 8: RNA Virome Workflow (Complete)** ✨
- **Phase 8.1**: 3 RNA virome collections (respiratory, arbovirus, fecal)
- **Phase 8.2**: Complete RNA workflow implementation
  - Reverse transcription with virus-type specific efficiency
  - rRNA depletion (Ribo-Zero/RiboMinus) modeling
  - RNA degradation and fragmentation
  - RNA-specific contamination profiles
  - Full integration with FASTQ generation
  - Comprehensive test suite (70+ tests)

### Test Coverage

```
Unit Tests:               40+ RNA workflow tests
Contamination Tests:      30+ RNA contamination tests
Integration Tests:        12 FASTQ generation tests
VLP Tests:                16 VLP/contamination tests
Collection Validation:    23/23 collections validated
Genomes:                  14,423 RefSeq viral genomes
Taxonomy Coverage:        57.1% (8,241 genomes assigned)
Literature Validation:    5/5 metrics validated
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

**RNA Virome Workflows** ✨ NEW
- Reverse transcription efficiency by virus type
- rRNA contamination levels (Qin et al. 2010, Greninger et al. 2015)
- Ribo-Zero depletion efficiency (Illumina technical documentation)
- RNA degradation rates (Fleige & Pfaffl 2006)

**Virion Size Relationships**
- Cui et al. (2014) - Genome length to virion size
- Nasir et al. (2017) - Virion size distributions
- Danovaro et al. (2011) - Marine viral particle sizes

See [Phase 5 Validation Report](docs/PHASE5_TASK3_VALIDATION_REPORT.md) for detailed validation results.

---

## Roadmap

### Completed (v0.1.0 - v0.6.0)
- ✅ Core simulator with contamination modeling
- ✅ 14,423 RefSeq viral genomes with ICTV taxonomy
- ✅ 23 curated virome collections
- ✅ Enhanced VLP enrichment (5 protocols)
- ✅ Amplification bias integration
- ✅ Critical disease/environmental collections
- ✅ **RNA virome workflow (RT, rRNA depletion, degradation)**
- ✅ **Taxonomy bug fix (469 genomes)**

### Planned (v0.7.0+)
- **Phase 9**: Additional host-associated collections (5+)
  - Blood/plasma, ocular, lung, urinary viromes
- **Phase 10**: Long-read sequencing support (PacBio HiFi, Nanopore)
- **Phase 11**: Temporal dynamics modeling
- **Phase 12**: Animal model collections (zebrafish, pig, chicken, primate)
- **Phase 13**: Environmental diversity (hot spring, hypersaline, hospital, plant)

See [ROADMAP.md](ROADMAP.md) for detailed development plans.

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

# Run RNA workflow tests
pytest tests/test_rna_workflow.py -v
pytest tests/test_rna_contamination.py -v
```

---

## Citation

If you use ViroForge in your research, please cite:

**Handley, Scott, and ViroForge Contributors. (2025). ViroForge: A Synthetic Virome Data Generator with RNA Workflow Support and Enhanced VLP Modeling (Version 0.6.0) [Computer software]. https://github.com/shandley/viroforge**

### BibTeX

```bibtex
@software{viroforge2025,
  title = {ViroForge: A Synthetic Virome Data Generator with RNA Workflow Support and Enhanced VLP Modeling},
  author = {Handley, Scott and contributors},
  year = {2025},
  url = {https://github.com/shandley/viroforge},
  version = {0.6.0},
  note = {Includes DNA and RNA virome workflows, 23 curated collections, and enhanced taxonomy}
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
**Development**: Claude Code (Anthropic) - Phase 7-8 RNA workflow & taxonomy fixes

---

## Support

- **Documentation**: Check `docs/` directory and `.claude/claude.md`
- **Questions**: Open a [GitHub Discussion](https://github.com/shandley/viroforge/discussions)
- **Bug Reports**: Open a [GitHub Issue](https://github.com/shandley/viroforge/issues)
- **Email**: scott.handley@wustl.edu

---

## Status

**Production Ready** - Phase 8.2 Complete

ViroForge is ready for use in benchmarking studies. All core functionality is complete and thoroughly tested:

- **23 curated virome collections** spanning healthy, disease, and environmental samples
- **14,423 RefSeq viral genomes** with enhanced ICTV taxonomy (57.1% coverage)
- **DNA virome workflow**: VLP enrichment (5 protocols), amplification bias, sequencing artifacts
- **RNA virome workflow**: Reverse transcription, rRNA depletion (Ribo-Zero), RNA degradation ✨ NEW
- **Enhanced taxonomy**: Fuzzy matching fixed 469 genomes, critical for HIV+ and RNA collections
- **Complete ground truth**: Taxonomic composition, abundance, workflow statistics
- **Platform support**: NovaSeq, MiSeq, HiSeq with realistic error models
- **Comprehensive testing**: 70+ tests covering all workflow components

**Last Updated**: 2025-11-09
