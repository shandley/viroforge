# ViroForge

**Forging synthetic viromes for benchmarking and validation**

A comprehensive mock metavirome data generator for testing and validating virome analysis pipelines.

---

## Overview

ViroForge generates realistic synthetic virome sequencing datasets with complete ground truth metadata, enabling rigorous validation of QC pipelines, assembly tools, taxonomic classifiers, and analysis workflows.

### Key Features

- 🧬 **Virome-specific simulation** - Models VLP enrichment, viral contamination patterns
- 🔬 **Library prep biases** - RdAB amplification, mechanical shearing, fragmentation
- 📊 **Sequencing artifacts** - NovaSeq polyG tails, optical duplicates, quality profiles
- 🌍 **Body-site profiles** - Gut, oral, skin, respiratory-specific viral compositions
- 📈 **Complete ground truth** - Taxonomic composition, abundance tables, read mappings
- 🎯 **Pre-built scenarios** - Ready-to-use profiles for common use cases
- ⚙️ **Highly configurable** - Full control over all simulation parameters

---

## Why ViroForge?

### The Problem

Current virome analysis tool validation approaches are limited:
- **Physical synthetic communities**: Expensive ($10k+), time-consuming, limited complexity
- **Existing simulators**: Bacterial-focused (CAMISIM), don't model VLP enrichment
- **Real datasets**: Unknown ground truth, can't control parameters

### The Solution

ViroForge fills the gap by providing:
- ✅ Unlimited synthetic datasets with complete ground truth
- ✅ VLP enrichment vs bulk metagenome comparisons
- ✅ Realistic contamination profiles (host DNA, rRNA, PhiX, reagent bacteria)
- ✅ Library prep and sequencing artifact modeling
- ✅ Standardized benchmarking datasets for the community

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/shandley/viroforge.git
cd viroforge

# Install dependencies
pip install -e .
```

### Generate Your First Mock Dataset

```bash
# Use a pre-built scenario
viroforge create \
  --profile gut_virome_vlp \
  --reads 10M \
  --output my_first_mock/

# Outputs:
# - my_first_mock/reads_R1.fastq.gz
# - my_first_mock/reads_R2.fastq.gz
# - my_first_mock/ground_truth_abundance.tsv
# - my_first_mock/ground_truth_taxonomy.tsv
# - my_first_mock/ground_truth_reads.tsv
# - my_first_mock/qc_metrics.json
```

---

## Use Cases

### 1. Validate QC Pipelines

```bash
# Generate dataset with known contamination
viroforge create --profile failed_vlp --host-dna 15.0 --output test_qc/

# Run through your QC pipeline
snakemake --use-conda test_qc/

# Verify: Should flag FAIL for high host contamination
```

### 2. Benchmark Analysis Tools

```bash
# Generate standardized test dataset
viroforge create --profile gut_virome_diverse --output benchmark/

# Run multiple classifiers
kraken2 --db viral benchmark/reads_R1.fastq.gz
kaiju -t nodes.dmp -f viral.fmi benchmark/reads_R1.fastq.gz

# Compare to ground truth
viroforge evaluate benchmark/ classifier_output.txt
```

### 3. Protocol Comparison

```bash
# Same community, different methods
viroforge create --profile gut_virome --vlp-enrichment true --output vlp/
viroforge create --profile gut_virome --vlp-enrichment false --output bulk/

# Compare viral recovery, contamination, diversity
```

### 4. Power Analysis

```bash
# Test different sequencing depths
for depth in 1M 5M 10M 50M; do
  viroforge create --profile gut_virome --reads $depth --output depth_${depth}/
done

# Determine minimum reads needed for your analysis
```

---

## Pre-Built Scenarios

ViroForge includes ready-to-use scenarios:

| Scenario | Description | Use Case |
|----------|-------------|----------|
| `gut_virome_clean` | Clean VLP-enriched gut virome | Algorithm development |
| `gut_virome_realistic` | Moderate contamination (2% host, 5% rRNA) | Realistic testing |
| `failed_vlp` | High contamination (15% host, 20% rRNA) | QC validation |
| `oral_virome` | Oral-specific viral composition | Body-site analysis |
| `low_biomass_skin` | Low viral load, high reagent contamination | Low-biomass samples |
| `novaseq_polyg_heavy` | Heavy polyG artifacts | Platform artifact testing |
| `optical_dup_high` | High optical duplicate rate | Duplicate removal testing |

---

## Configuration

### Using Pre-Built Profiles (Easy)

```bash
viroforge create --profile gut_virome_vlp --output my_data/
```

### Custom Configuration (Advanced)

Create `my_config.yaml`:

```yaml
community:
  source: gut_virome_database
  n_species: 50
  diversity: high
  abundance_distribution: lognormal

contamination:
  host_dna: 2.0  # percent
  rrna: 5.0
  phix: 0.1
  reagent_bacteria: 0.5

library_prep:
  method: vlp_rdab
  vlp_enrichment: true
  pcr_cycles: 40
  fragmentation: mechanical_shearing

sequencing:
  platform: novaseq
  read_length: 150
  insert_size: 300
  coverage: 10x
```

Run with custom config:

```bash
viroforge create --config my_config.yaml --output my_data/
```

---

## Output Files

### Sequencing Reads
- `reads_R1.fastq.gz` - Forward reads
- `reads_R2.fastq.gz` - Reverse reads

### Ground Truth Metadata
- `ground_truth_abundance.tsv` - True viral abundance matrix
- `ground_truth_taxonomy.tsv` - Complete taxonomic assignments
- `ground_truth_reads.tsv` - Read origin mapping (read → genome → taxonomy)
- `ground_truth_assembly.tsv` - Expected assembly outcomes

### QC Metrics
- `qc_metrics.json` - Expected QC values (ViromeQC score, contamination %)

---

## Documentation

- **[Design Rationale](docs/DESIGN_RATIONALE.md)** - Complete design thinking and literature review
- **[Installation Guide](docs/INSTALLATION.md)** - Detailed installation instructions
- **[User Guide](docs/USER_GUIDE.md)** - Comprehensive usage documentation
- **[Tutorial](docs/TUTORIAL.md)** - Step-by-step examples
- **[API Reference](docs/API.md)** - Programmatic usage

---

## Project Status

**Current Version**: 0.1.0-dev (Phase 1: ~80% Complete)

### ✅ Phase 1 Complete
- Community composition (5 body sites, 3 abundance distributions)
- Contamination profiles (4 levels: clean, realistic, heavy, failed)
- FASTQ generation (InSilicoSeq integration)
- Validation framework (prevents FASTQ quality issues)
- Ground truth tracking (complete metadata)
- Comprehensive testing (100% pass rate)

### 🚧 Phase 2 Starting (Virome-Specific Features)
**Next 12 weeks**: Implementing core virome biology
- VLP enrichment framework (Weeks 1-3)
- Amplification bias (RdAB, MDA) (Weeks 4-6)
- Platform artifacts (NovaSeq polyG, optical dups) (Weeks 7-8)
- Integration & validation (Weeks 9-12)

See `docs/IMPLEMENTATION_PLAN.md` for detailed Phase 2 plan.

### Roadmap
- **Phase 1** ✅ Core functionality (80% complete)
- **Phase 2** 🚧 Virome-specific features (starting)
- **Phase 3** ⏳ Validation and refinement
- **Phase 4** ⏳ Publication and release

---

## Contributing

We welcome contributions! This is an open-source project under the MIT License.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## Citation

**Software**:
```
ViroForge: A synthetic virome data generator
Scott Handley Lab, Washington University in St. Louis
https://github.com/shandley/viroforge
```

**Publication**: Coming soon (manuscript in preparation)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

Developed to support virome analysis pipeline validation and benchmarking.

### Related Projects
- [lab-virome-QC](https://github.com/shandley/lab-virome-QC) - VLP-enriched virome QC pipeline
- [ViromeQC](https://github.com/SegataLab/viromeqc) - Virome enrichment assessment
- [CAMISIM](https://github.com/CAMI-challenge/CAMISIM) - Metagenome simulator (inspiration)
- [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) - Illumina read simulator (integrated)

### Funding
- [Add funding sources]

---

## Contact

**Principal Investigator**: Scott Handley
**Email**: scott.handley@wustl.edu
**Lab Website**: [Add lab website]
**GitHub Issues**: https://github.com/shandley/viroforge/issues

---

## Status

⚠️ **Under Active Development** - This project is in early development (Phase 1). APIs and features are subject to change.

**Last Updated**: 2025-01-30
