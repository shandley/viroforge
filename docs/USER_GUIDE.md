# ViroForge User Guide

**Comprehensive Guide to Synthetic Virome Data Generation**

This guide provides in-depth documentation of all ViroForge features, parameters, and best practices.

---

## Table of Contents

1. [Installation](#installation)
2. [Core Concepts](#core-concepts)
3. [Viral Communities](#viral-communities)
4. [Contamination Profiles](#contamination-profiles)
5. [VLP Enrichment](#vlp-enrichment)
6. [Amplification Bias](#amplification-bias)
7. [Platform Artifacts](#platform-artifacts)
8. [Complete Workflows](#complete-workflows)
9. [Ground Truth Files](#ground-truth-files)
10. [Best Practices](#best-practices)
11. [Troubleshooting](#troubleshooting)
12. [Advanced Topics](#advanced-topics)

---

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Virtual environment (recommended)

### Basic Installation

```bash
# Clone repository
git clone https://github.com/shandley/viroforge.git
cd viroforge

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install ViroForge
pip install -e .
```

### Development Installation

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Run tests to verify installation
pytest tests/ -v
```

### Optional Dependencies

```bash
# For FASTQ generation (InSilicoSeq)
conda install -c bioconda insilicoseq

# For visualization
pip install matplotlib seaborn
```

---

## Core Concepts

### The ViroForge Workflow

ViroForge models the complete virome data generation pipeline:

```
1. Viral Community      → Body-site specific viral composition
2. Contamination        → Host DNA, bacterial DNA, reagent contamination
3. Mock Composition     → Combines viral + contamination (e.g., 50/50)
4. VLP Enrichment       → Filtration + nuclease treatment (→ 97% viral)
5. Amplification Bias   → RdAB/MDA/Linker amplification
6. Read Generation      → Paired-end Illumina reads
7. Platform Artifacts   → PolyG tails, optical duplicates, index hopping
8. FASTQ Output         → Final sequencing files + ground truth
```

### Key Principles

**1. Ground Truth Tracking**
- Every read traced back to source genome
- Complete metadata preserved throughout pipeline
- Enables precise validation of analysis tools

**2. Literature Validation**
- All parameters from peer-reviewed literature
- Realistic artifact rates
- Biologically plausible viral compositions

**3. Reproducibility**
- Random seeds for deterministic output
- Complete parameter documentation
- Version-controlled parameters

**4. Modularity**
- Independent components
- Mix and match protocols
- Custom parameter configuration

---

## Viral Communities

### Creating Viral Communities

#### Pre-Built Body Site Profiles

```python
from viroforge.core.community import create_body_site_profile

# Gut virome
gut = create_body_site_profile(
    body_site='gut',
    n_genomes=50,
    random_seed=42
)

# Available body sites:
# - 'gut': Gut-associated viruses (Caudovirales-dominated)
# - 'oral': Oral cavity viruses (Siphoviridae-rich)
# - 'skin': Skin-associated viruses (Podoviridae)
# - 'respiratory': Respiratory tract viruses
# - 'diverse': Maximum diversity (all families)
```

#### Custom Viral Community

```python
from viroforge.core.community import ViralCommunity, ViralGenome

# Create custom community
community = ViralCommunity(name='my_custom_virome')

# Add individual genomes
genome1 = ViralGenome(
    genome_id='NC_001416',
    organism='Enterobacteria phage T7',
    family='Podoviridae',
    sequence='ATCG...',  # Full genome sequence
    abundance=0.15       # Relative abundance (0-1)
)
community.add_genome(genome1)

# Add more genomes...
# Total abundances should sum to 1.0
```

### Abundance Distributions

ViroForge supports three abundance distribution models:

**1. Log-Normal Distribution** (Default - Most Realistic)
```python
community = create_body_site_profile(
    body_site='gut',
    n_genomes=50,
    abundance_distribution='lognormal',
    lognormal_mean=0.0,    # Mean of log-transformed abundances
    lognormal_sigma=1.0,   # Spread of distribution
    random_seed=42
)
```
- Mimics natural viral communities
- Few dominant viruses, many rare ones
- Literature: [Power et al. 2015]

**2. Power Law Distribution**
```python
community = create_body_site_profile(
    body_site='gut',
    n_genomes=50,
    abundance_distribution='power_law',
    power_law_alpha=1.5,   # Exponent (higher = more uneven)
    random_seed=42
)
```
- Even more uneven than log-normal
- Heavy-tailed distribution
- Good for extreme diversity

**3. Even Distribution**
```python
community = create_body_site_profile(
    body_site='gut',
    n_genomes=50,
    abundance_distribution='even',
    random_seed=42
)
```
- All genomes equal abundance (1/n each)
- Useful for controlled testing
- Not biologically realistic

### Community Statistics

```python
# Get summary statistics
stats = community.get_summary_stats()

print(f"Total genomes: {stats['total_genomes']}")
print(f"Families: {stats['families']}")
print(f"Mean length: {stats['mean_length']:,.0f} bp")
print(f"Mean GC: {stats['mean_gc']*100:.1f}%")
print(f"Total abundance: {stats['total_abundance']:.4f}")  # Should be 1.0

# Access individual genomes
for genome in community.genomes:
    print(f"{genome.organism}: {genome.abundance:.4f}")
```

---

## Contamination Profiles

### Pre-Defined Contamination Levels

```python
from viroforge.core.contamination import create_contamination_profile

# Clean VLP enrichment (minimal contamination)
clean = create_contamination_profile('clean', random_seed=42)
# ~2% total: 1% host, 0.5% bacterial, 0.5% other

# Realistic VLP enrichment (typical)
realistic = create_contamination_profile('realistic', random_seed=42)
# ~5% total: 2% host, 2% bacterial, 1% other

# Heavy contamination (sub-optimal)
heavy = create_contamination_profile('heavy', random_seed=42)
# ~10% total: 5% host, 3% bacterial, 2% other

# Failed VLP enrichment (should flag QC)
failed = create_contamination_profile('failed', random_seed=42)
# ~35% total: 15% host, 15% bacterial, 5% other
```

### Contamination Sources

**1. Host DNA**
- Human genomic DNA (default)
- Can specify organism: `host_organism='mouse'`, `'rat'`, `'human'`
- Removed by nuclease treatment + filtration

**2. Bacterial DNA**
- Gut microbiome bacteria (if body_site='gut')
- Reagent contamination (Pseudomonas, etc.)
- Partially removed by VLP enrichment

**3. Fungal DNA**
- Environmental fungi
- Reagent contamination
- Partially removed by VLP enrichment

**4. PhiX Control Spike-In**
- Illumina sequencing control
- Typical: 0.1-1%
- Should be detected/removed in QC

### Custom Contamination

```python
from viroforge.core.contamination import ContaminationProfile, ContaminantGenome

# Create custom profile
contamination = ContaminationProfile(name='my_custom_contam')

# Add host DNA
host = ContaminantGenome(
    genome_id='host_chr1',
    organism='Homo sapiens',
    sequence='ATCG...',
    abundance=0.05,  # 5%
    contaminant_type='host'
)
contamination.add_contaminant(host)

# Add more contaminants...
```

---

## VLP Enrichment

VLP (Virus-Like Particle) enrichment is THE defining feature of viromics.

### Pre-Built VLP Protocols

```python
from viroforge.enrichment import (
    standard_vlp,
    iron_chloride_vlp,
    ultracentrifugation_vlp,
    syringe_filter_vlp
)

# Standard tangential flow filtration
standard = standard_vlp()
# 0.2 μm TFF, 95% nuclease efficiency, 80% recovery

# Iron chloride precipitation (Conceição-Neto et al.)
iron = iron_chloride_vlp()
# FeCl3 precipitation, 98% nuclease efficiency, 60% recovery

# Ultracentrifugation (traditional method)
ultra = ultracentrifugation_vlp()
# Density gradient, 90% nuclease efficiency, 50% recovery

# Syringe filter (field-friendly)
syringe = syringe_filter_vlp()
# 0.45 μm syringe, 80% nuclease efficiency, 85% recovery
```

### Custom VLP Protocol

```python
from viroforge.enrichment import VLPEnrichment

# Create custom protocol
custom_vlp = VLPEnrichment(
    # Filtration parameters
    filtration_method='tangential_flow',  # or 'syringe', 'ultracentrifugation', 'none'
    filtration_cutoff_um=0.2,            # Filter pore size (μm)
    prefiltration_cutoff_um=0.45,        # Pre-filter (optional, None to disable)

    # Nuclease treatment
    nuclease_treatment=True,              # Enable nuclease treatment
    nuclease_efficiency=0.95,             # Efficiency (0-1)

    # Size-based retention
    size_retention_curve='sigmoid',       # or 'step' function
    stochastic_variation=0.2,            # Random variation (0-1)

    # Reproducibility
    random_seed=42
)

# Apply to composition
custom_vlp.apply(composition)
```

### VLP Parameters Explained

**Filtration Method**
- `'tangential_flow'`: Most common, continuous flow (TFF)
- `'syringe'`: Simple syringe filter, field-friendly
- `'ultracentrifugation'`: Density gradient, traditional
- `'none'`: No filtration (bulk metagenome)

**Filtration Cutoff**
- `0.1 μm`: Removes bacteria, keeps large viruses only
- `0.2 μm`: Standard virome cutoff, keeps most viruses
- `0.45 μm`: Loose cutoff, keeps bacteria-sized particles
- Smaller = higher purity, lower yield

**Nuclease Efficiency**
- `0.95`: Excellent nuclease treatment (standard)
- `0.80`: Good nuclease treatment
- `0.60`: Sub-optimal nuclease treatment
- `0.30`: Failed nuclease treatment (should flag QC)
- Removes free DNA/RNA, not virion-packaged genomes

**Size Retention Curve**
- `'sigmoid'`: Gradual size-based retention (realistic)
- `'step'`: Sharp cutoff at filter size (idealized)

**Stochastic Variation**
- `0.0`: Perfect reproducibility (idealized)
- `0.2`: Normal variation (realistic)
- `0.5`: High variation (poor technique)
- Models prep-to-prep variability

### Modeling Failed VLP Enrichment

```python
# Failed enrichment (for QC validation)
failed_vlp = VLPEnrichment(
    filtration_cutoff_um=0.45,      # Too large
    nuclease_efficiency=0.3,        # Poor nuclease
    stochastic_variation=0.5,       # High variability
    random_seed=42
)

failed_vlp.apply(composition)
# Result: <80% viral (should flag ViromeQC)
```

---

## Amplification Bias

Most virome protocols require PCR amplification, which introduces systematic bias.

### Pre-Built Amplification Methods

```python
from viroforge.amplification import (
    no_amplification,
    rdab_30_cycles,
    rdab_40_cycles,
    mda_standard,
    mda_overnight,
    linker_standard
)

# No amplification (high biomass control)
none = no_amplification()

# RdAB amplification (most common virome method)
rdab30 = rdab_30_cycles()  # Less bias
rdab40 = rdab_40_cycles()  # Standard protocol

# MDA amplification (low biomass samples)
mda = mda_standard()        # 4-16 hours
mda_long = mda_overnight()  # Overnight (more bias)

# Linker amplification (modern, minimal bias)
linker = linker_standard()
```

### RdAB Amplification Details

```python
from viroforge.amplification import RdABAmplification

# Custom RdAB protocol
rdab = RdABAmplification(
    cycles=40,              # PCR cycles (20-45)
    efficiency=0.90,        # Amplification efficiency per cycle
    length_bias_strength=1.0,    # Length bias intensity
    gc_bias_strength=1.0,        # GC bias intensity
    gc_optimal=0.5,              # Optimal GC content
    gc_tolerance=0.15,           # GC tolerance window
    random_seed=42
)
```

**Bias Models:**
- **Length Bias**: Short genomes amplify better than long
  - `exp(-0.015 * length_kb * strength)`
  - Stronger with more cycles
- **GC Bias**: Extreme GC amplifies poorly
  - `exp(-((gc - optimal) / tolerance)^2 * strength)`
  - Quadratic penalty for extreme GC

### MDA Amplification Details

```python
from viroforge.amplification import MDAAmplification

# Custom MDA protocol
mda = MDAAmplification(
    duration_hours=8,            # Reaction time (4-24 hours)
    gc_bias_strength=2.0,        # GC bias (stronger than RdAB)
    gc_optimal=0.5,              # Optimal GC
    gc_tolerance=0.10,           # Tighter tolerance
    stochasticity=0.3,           # Random variation (high)
    chimera_rate=0.05,           # Chimera formation (5%)
    random_seed=42
)
```

**MDA Characteristics:**
- **Extreme GC Bias**: 10-1000x range in final abundances
- **High Stochasticity**: Large random variations
- **No Length Bias**: All fragments amplified equally
- **Chimera Formation**: Template switching artifacts

### Linker Amplification Details

```python
from viroforge.amplification import LinkerAmplification

# Custom linker protocol
linker = LinkerAmplification(
    cycles=20,               # Fewer cycles needed
    efficiency=0.92,         # High efficiency
    gc_bias_strength=0.2,    # Minimal GC bias
    gc_optimal=0.5,
    gc_tolerance=0.20,       # Wide tolerance
    random_seed=42
)
```

**Linker Advantages:**
- **No Length Bias**: Adapters on all fragments
- **Minimal GC Bias**: Modern polymerases
- **Fewer Cycles**: Less bias accumulation
- **Modern Method**: Increasingly common

---

## Platform Artifacts

Different Illumina platforms introduce different artifacts.

### Pre-Built Platform Profiles

```python
from viroforge.artifacts import (
    novaseq_6000,
    nextseq_2000,
    miseq,
    hiseq_2500,
    no_artifacts
)

# NovaSeq 6000 (patterned flow cell)
novaseq = novaseq_6000()
# PolyG: 2.5%, OpticalDup: 9%, IndexHop: 1.5%

# NextSeq 2000 (patterned flow cell)
nextseq = nextseq_2000()
# PolyG: 2.0%, OpticalDup: 7%, IndexHop: 1.0%

# MiSeq (cluster-based flow cell)
miseq_platform = miseq()
# PolyG: 0%, OpticalDup: 2.5%, IndexHop: 0.1%

# HiSeq 2500 (cluster-based, legacy)
hiseq = hiseq_2500()
# PolyG: 0%, OpticalDup: 4.5%, IndexHop: 0.2%

# Ideal (no artifacts - control)
ideal = no_artifacts()
```

### Platform Comparison

| Platform | Flow Cell | PolyG | Opt Dup | Index Hop | Throughput |
|----------|-----------|-------|---------|-----------|------------|
| NovaSeq 6000 | Patterned | ✓ | High | High | Very High |
| NextSeq 2000 | Patterned | ✓ | Medium | Medium | High |
| MiSeq | Cluster | ✗ | Low | Low | Low |
| HiSeq 2500 | Cluster | ✗ | Medium | Low | Medium |

**Key Insight**: Patterned flow cells (NovaSeq, NextSeq) have polyG tails. Cluster-based flow cells (MiSeq, HiSeq) do NOT.

### Custom Platform Configuration

```python
from viroforge.artifacts import (
    PlatformProfile,
    PolyGTailArtifact,
    OpticalDuplicateArtifact,
    IndexHoppingArtifact
)

# Create custom platform
custom_platform = PlatformProfile(
    name='My Custom Platform',
    flow_cell_type='patterned',
    artifacts=[
        PolyGTailArtifact(
            frequency=0.025,     # 2.5% of reads
            min_length=10,       # Minimum polyG length
            max_length=50,       # Maximum polyG length
            r1_rate=0.3,         # R1 affected rate (relative)
            r2_rate=1.0,         # R2 affected rate (relative, R2 > R1)
            random_seed=42
        ),
        OpticalDuplicateArtifact(
            rate=0.09,              # 9% duplication rate
            proximity_threshold=100, # Distance for duplication
            tile_size=10000,        # Flow cell tile size
            random_seed=42
        ),
        IndexHoppingArtifact(
            rate=0.015,          # 1.5% hopping rate
            n_indexes=96,        # Number of sample indexes
            random_seed=42
        )
    ]
)
```

### Artifact Descriptions

**1. PolyG Tails**
- **Cause**: Incomplete fluorophore quenching in patterned flow cells
- **Effect**: Homopolymer G runs added to read ends
- **Typical**: 1-3% of reads, 10-50bp tails
- **R2 > R1**: R2 more affected than R1
- **Only in**: NovaSeq, NextSeq (patterned flow cells)
- **Remove with**: Fastp, Cutadapt with polyG trimming

**2. Optical Duplicates**
- **Cause**: Adjacent cluster signal bleed
- **Effect**: Duplicate read pairs with nearby coordinates
- **Typical**: 2.5-9% depending on platform
- **All platforms**: But rate varies
- **Remove with**: Picard MarkDuplicates, samtools markdup

**3. Index Hopping**
- **Cause**: Free adapter barcode mixing in multiplexed libraries
- **Effect**: Reads misassigned to wrong sample
- **Typical**: 0.1-2% depending on platform
- **Higher in**: Patterned flow cells
- **Mitigate**: Unique dual indexes, higher cluster density

---

## Complete Workflows

### Basic Workflow

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles

# 1. Create composition
community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
contamination = create_contamination_profile('realistic', random_seed=42)
composition = MockViromeComposition(
    name='basic_workflow',
    viral_community=community,
    contamination_profile=contamination,
    viral_fraction=0.5
)

# 2. Apply VLP enrichment
vlp = standard_vlp()
vlp.apply(composition)

# 3. Apply amplification
amplification = rdab_40_cycles()
amplification.apply(composition)

# 4. Ready for read generation
print(f"Final viral fraction: {composition.viral_fraction*100:.1f}%")
```

### Complete End-to-End Pipeline

See `examples/complete_workflow_integrated.py` for full FASTQ generation workflow.

---

## Ground Truth Files

### Composition Table

**`ground_truth_composition.tsv`**

| genome_id | organism | family | abundance | length | gc_content |
|-----------|----------|--------|-----------|--------|------------|
| NC_001416 | Enterobacteria phage T7 | Podoviridae | 0.1234 | 39937 | 0.485 |
| NC_007458 | Escherichia phage MS2 | Leviviridae | 0.0567 | 3569 | 0.518 |

**Usage:**
```python
import pandas as pd
truth = pd.read_csv('ground_truth_composition.tsv', sep='\t')
print(truth.head())
```

### Read Mapping Table

**`ground_truth_read_mapping.tsv`**

| read_id | genome_id | genome_name | family |
|---------|-----------|-------------|--------|
| read_001_forward | NC_001416 | Enterobacteria phage T7 | Podoviridae |
| read_002_forward | NC_007458 | Escherichia phage MS2 | Leviviridae |

**Usage:**
```python
mapping = pd.read_csv('ground_truth_read_mapping.tsv', sep='\t')
# Validate classifier output
predicted = pd.read_csv('classifier_output.txt', sep='\t')
merged = mapping.merge(predicted, on='read_id')
accuracy = (merged['true_family'] == merged['predicted_family']).mean()
print(f"Classification accuracy: {accuracy*100:.1f}%")
```

---

## Best Practices

### 1. Always Use Random Seeds

```python
# Good: Reproducible
community = create_body_site_profile('gut', n_genomes=50, random_seed=42)

# Bad: Not reproducible
community = create_body_site_profile('gut', n_genomes=50)
```

### 2. Match Body Site to Research Question

```python
# Studying gut viromes? Use gut profile
gut_community = create_body_site_profile('gut', ...)

# Studying oral viromes? Use oral profile
oral_community = create_body_site_profile('oral', ...)
```

### 3. Use Realistic Contamination

```python
# For testing VLP success
realistic_contam = create_contamination_profile('realistic')

# For testing QC flags
failed_contam = create_contamination_profile('failed')
```

### 4. Document Your Parameters

```python
# Save parameters for reproducibility
import json

params = {
    'body_site': 'gut',
    'n_genomes': 50,
    'contamination_level': 'realistic',
    'vlp_protocol': 'standard',
    'amplification': 'rdab_40',
    'platform': 'novaseq',
    'random_seed': 42
}

with open('parameters.json', 'w') as f:
    json.dump(params, f, indent=2)
```

### 5. Validate Abundance Sums

```python
# Total abundance should be 1.0
total = composition.get_total_abundance()
assert abs(total - 1.0) < 0.01, f"Abundance sum error: {total}"
```

### 6. Export Ground Truth

```python
# Always export ground truth for validation
composition.export_composition('ground_truth_composition.tsv')
# Export read mappings after read generation
```

---

## Troubleshooting

### Common Issues

**1. ImportError: No module named 'viroforge'**
```bash
# Solution: Install in development mode
pip install -e .
```

**2. Total abundance != 1.0**
```python
# Check abundance distribution
total = composition.get_total_abundance()
print(f"Total: {total}")

# Renormalize if needed
composition.renormalize_abundances()
```

**3. VLP enrichment not improving viral fraction**
```python
# Check initial viral fraction
print(f"Before VLP: {composition.viral_fraction}")

# Verify VLP was applied
vlp = standard_vlp()
vlp.apply(composition)
print(f"After VLP: {composition.viral_fraction}")

# Should increase from 0.5 to ~0.97
```

**4. No polyG tails in reads**
```python
# Check platform type
platform = novaseq_6000()
print(f"Flow cell: {platform.flow_cell_type}")

# MiSeq has NO polyG (cluster-based)
# Only patterned flow cells have polyG
```

**5. Tests failing**
```bash
# Update dependencies
pip install -e ".[dev]" --upgrade

# Run specific test
pytest tests/test_enrichment.py -v

# Check Python version
python --version  # Should be 3.8+
```

---

## Advanced Topics

### Custom Genome Databases

```python
from viroforge.core.community import ViralCommunity, ViralGenome
from Bio import SeqIO

# Load genomes from FASTA
community = ViralCommunity(name='custom')

for record in SeqIO.parse('my_genomes.fasta', 'fasta'):
    genome = ViralGenome(
        genome_id=record.id,
        organism=record.description,
        family='Unknown',
        sequence=str(record.seq),
        abundance=0.0  # Will set later
    )
    community.add_genome(genome)

# Set abundances (must sum to 1.0)
n = len(community.genomes)
for genome in community.genomes:
    genome.abundance = 1.0 / n
```

### Batch Processing

```python
# Generate multiple datasets
for seed in range(10):
    community = create_body_site_profile('gut', n_genomes=50, random_seed=seed)
    contamination = create_contamination_profile('realistic', random_seed=seed)

    composition = MockViromeComposition(
        name=f'batch_{seed}',
        viral_community=community,
        contamination_profile=contamination,
        viral_fraction=0.5
    )

    vlp = standard_vlp()
    vlp.apply(composition)

    # Generate reads...
    print(f"Dataset {seed} complete: {composition.viral_fraction*100:.1f}% viral")
```

### Parameter Sweeps

```python
# Test different VLP enrichment efficiencies
efficiencies = [0.3, 0.5, 0.7, 0.9, 0.95]

results = []
for eff in efficiencies:
    composition = MockViromeComposition(...)

    vlp = VLPEnrichment(nuclease_efficiency=eff, random_seed=42)
    vlp.apply(composition)

    results.append({
        'efficiency': eff,
        'viral_fraction': composition.viral_fraction
    })

import pandas as pd
df = pd.DataFrame(results)
print(df)
```

### Custom Artifact Rates

```python
from viroforge.artifacts import PlatformProfile, PolyGTailArtifact

# Test different artifact rates
for polyg_rate in [0.01, 0.025, 0.05, 0.10]:
    platform = PlatformProfile(
        name=f'Test_{polyg_rate*100:.0f}pct_polyG',
        flow_cell_type='patterned',
        artifacts=[
            PolyGTailArtifact(frequency=polyg_rate, random_seed=42)
        ]
    )

    # Apply to reads and count artifacts...
```

---

## Additional Resources

- **Tutorial**: `docs/TUTORIAL.md` - Step-by-step guide
- **API Reference**: `docs/API.md` - Complete API documentation
- **Examples**: `examples/` - Working code examples
- **Design Rationale**: `docs/DESIGN_RATIONALE.md` - Why we made certain choices
- **VLP Biology**: `docs/VLP_ENRICHMENT_BIOLOGY.md` - Detailed VLP background
- **Lab Notebook**: `lab-notebook/INDEX.md` - Development history

---

## Getting Help

- 📖 **Documentation**: Check `docs/` directory
- 💬 **Discussions**: https://github.com/shandley/viroforge/discussions
- 🐛 **Bug Reports**: https://github.com/shandley/viroforge/issues
- 📧 **Email**: scott.handley@wustl.edu

---

**Happy Forging!** 🔨🦠
