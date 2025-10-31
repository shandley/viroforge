# ViroForge Mock Virome Creation Workflow

**Purpose**: This document describes the complete workflow for creating mock virome datasets, including when and how different components are added, modified, and processed.

**Last Updated**: 2025-01-31
**Status**: Phase 1 Complete, Phase 2 in Development

---

## Table of Contents

1. [Overview](#overview)
2. [Phase 1 Workflow (Current)](#phase-1-workflow-current)
3. [Phase 2 Workflow (In Development)](#phase-2-workflow-in-development)
4. [Component Details](#component-details)
5. [Order of Operations](#order-of-operations)
6. [Examples](#examples)

---

## Overview

ViroForge creates realistic mock virome datasets through a multi-stage pipeline that mimics the biological and technical processes of real virome sequencing. Understanding the order of operations is critical because **each stage modifies either genome abundances or read sequences**.

### Key Principle: Two Types of Modifications

| Modification Type | What Changes | When Applied | Examples |
|-------------------|--------------|--------------|----------|
| **Abundance-Modifying** | Relative abundance of genomes | BEFORE read generation | VLP enrichment, amplification bias, library prep bias |
| **Read-Modifying** | The reads themselves | AFTER read generation | PolyG tails, optical duplicates, quality scores |

---

## Phase 1 Workflow (Current)

**Status**: ✅ Fully Implemented and Tested

### Pipeline Stages

```
┌─────────────────────────────────────────────────────────────┐
│ Stage 1: VIRAL COMMUNITY CREATION                           │
├─────────────────────────────────────────────────────────────┤
│ • Sample or create viral genomes                            │
│ • Assign taxonomic information                              │
│ • Apply abundance distribution (lognormal, powerlaw, even)  │
│                                                              │
│ Output: ViralCommunity object                               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 2: CONTAMINATION PROFILE CREATION                     │
├─────────────────────────────────────────────────────────────┤
│ • Add host DNA (human, mouse, etc.)                         │
│ • Add rRNA contamination (bacterial, archaeal, eukaryotic)  │
│ • Add reagent bacteria (kit contamination)                  │
│ • Add PhiX174 control spike-in                              │
│                                                              │
│ Output: ContaminationProfile object                         │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 3: COMPOSITION ASSEMBLY                               │
├─────────────────────────────────────────────────────────────┤
│ • Combine viral community + contamination profile           │
│ • Normalize abundances to target viral_fraction             │
│ • Example: 85% viral, 15% contamination                     │
│                                                              │
│ Output: MockViromeComposition object                        │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 4: READ GENERATION (InSilicoSeq)                      │
├─────────────────────────────────────────────────────────────┤
│ • Write genomes + abundances to temp files                  │
│ • Call InSilicoSeq with platform error model                │
│ • Generate paired-end FASTQ reads                           │
│ • Validate every FASTQ record (seq/qual length match)       │
│                                                              │
│ Output: R1/R2 FASTQ files + ground truth metadata           │
└─────────────────────────────────────────────────────────────┘
```

### Code Example (Phase 1)

```python
from viroforge.core import create_body_site_profile, create_contamination_profile
from viroforge.utils import MockViromeComposition
from viroforge.simulators import generate_reads

# Stage 1: Create viral community
viral_community = create_body_site_profile(
    body_site='gut',
    n_genomes=50,
    random_seed=42
)

# Stage 2: Create contamination profile
contamination = create_contamination_profile(
    profile='realistic',  # 2% host, 5% rRNA, 0.5% bacteria, 0.1% PhiX
    host_organism='human',
    random_seed=42
)

# Stage 3: Combine into composition
composition = MockViromeComposition(
    name='gut_virome_realistic',
    viral_community=viral_community,
    contamination_profile=contamination,
    viral_fraction=0.85  # 85% viral, 15% contamination
)

# Stage 4: Generate reads
output = generate_reads(
    composition=composition,
    output_prefix='my_dataset',
    n_reads=1_000_000,
    model='NovaSeq',
    random_seed=42
)

# Output: my_dataset_R1.fastq.gz, my_dataset_R2.fastq.gz, ground_truth.tsv
```

---

## Phase 2 Workflow (In Development)

**Status**: 🚧 Planned for Implementation (Weeks 1-12)

### Extended Pipeline with Virome-Specific Features

```
┌─────────────────────────────────────────────────────────────┐
│ Stage 1: VIRAL COMMUNITY CREATION                           │
│ (Same as Phase 1)                                           │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 2: CONTAMINATION PROFILE CREATION                     │
│ (Same as Phase 1)                                           │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 3: INITIAL COMPOSITION (Pre-Enrichment)               │
├─────────────────────────────────────────────────────────────┤
│ • Combine viral + contamination (BULK METAGENOME state)     │
│ • Represents sample BEFORE any lab processing               │
│ • Example: 50% viral, 30% host, 15% rRNA, 5% bacteria       │
│                                                              │
│ Output: MockViromeComposition (pre-enrichment)              │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 4: VLP ENRICHMENT (ABUNDANCE-MODIFYING) 🆕             │
├─────────────────────────────────────────────────────────────┤
│ • Apply size-based filtration (0.2-0.45 μm)                 │
│ • Apply nuclease treatment (removes free DNA)               │
│ • Apply family-specific enrichment factors                  │
│ • Dramatically increases viral fraction                     │
│                                                              │
│ MODIFIES: Genome abundances                                 │
│ Example: 50% viral → 95% viral, host 30% → 0.5%            │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 5: AMPLIFICATION BIAS (ABUNDANCE-MODIFYING) 🆕         │
├─────────────────────────────────────────────────────────────┤
│ • Apply method-specific bias (RdAB, MDA, Linker, or none)  │
│ • Length-dependent amplification (short favored)            │
│ • GC-dependent amplification                                │
│ • Cycle-dependent scaling                                   │
│                                                              │
│ MODIFIES: Genome abundances                                 │
│ Example: Short genomes 10x enriched after 40 PCR cycles    │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 6: LIBRARY PREP BIAS (ABUNDANCE-MODIFYING) 🆕          │
├─────────────────────────────────────────────────────────────┤
│ • Apply fragmentation method bias                           │
│   - Mechanical shearing: minimal bias                       │
│   - Tagmentation (Nextera): GC-dependent bias               │
│   - Enzymatic: moderate bias                                │
│                                                              │
│ MODIFIES: Genome abundances (within-genome coverage)        │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 7: READ GENERATION (InSilicoSeq)                      │
│ (Same as Phase 1)                                           │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage 8: PLATFORM ARTIFACTS (READ-MODIFYING) 🆕              │
├─────────────────────────────────────────────────────────────┤
│ • Add polyG tails (NovaSeq 2-channel chemistry)             │
│ • Add optical duplicates (patterned flow cells)             │
│ • Add index hopping (multiplexed samples)                   │
│                                                              │
│ MODIFIES: Read sequences themselves                         │
│ Example: 20% of reads get polyG tails added                 │
└─────────────────────────────────────────────────────────────┘
                            ↓
                   FINAL FASTQ OUTPUT
```

### Code Example (Phase 2)

```python
from viroforge.core import create_body_site_profile, create_contamination_profile
from viroforge.utils import MockViromeComposition
from viroforge.enrichment import VLPEnrichment  # 🆕 Phase 2
from viroforge.amplification import RdABAmplification  # 🆕 Phase 2
from viroforge.library_prep import MechanicalShearing  # 🆕 Phase 2
from viroforge.simulators import generate_reads

# Stages 1-2: Create base components (same as Phase 1)
viral_community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
contamination = create_contamination_profile('realistic', random_seed=42)

# Stage 3: Create BULK METAGENOME composition (pre-enrichment)
composition = MockViromeComposition(
    name='gut_virome',
    viral_community=viral_community,
    contamination_profile=contamination,
    viral_fraction=0.50  # 50% viral BEFORE enrichment (bulk state)
)

# Stage 4: Apply VLP enrichment 🆕
vlp = VLPEnrichment(
    filtration_cutoff_um=0.2,
    nuclease_efficiency=0.95
)
vlp.apply(composition)  # Now ~95% viral after enrichment!

# Stage 5: Apply amplification bias 🆕
rdab = RdABAmplification(cycles=40)
rdab.apply(composition)  # Short genomes now over-represented

# Stage 6: Apply library prep bias 🆕
shearing = MechanicalShearing()
shearing.apply(composition)  # Minimal bias for mechanical shearing

# Stage 7: Generate reads
output = generate_reads(
    composition=composition,
    output_prefix='my_vlp_dataset',
    n_reads=1_000_000,
    model='NovaSeq',
    random_seed=42
)

# Stage 8: Platform artifacts 🆕
# (This happens automatically inside generate_reads or as post-processing)
# - PolyG tails added to 20% of reads
# - Optical duplicates created at 5% rate
```

---

## Component Details

### 1. Viral Community

**What it contains:**
- Collection of viral genomes (synthetic or real)
- Taxonomic assignments (family, genus, species)
- Relative abundances following specified distribution
- Genome sequences (or placeholders)

**Sources:**
- Synthetic genomes (current default): Random sequences of realistic lengths
- Real genomes (optional): Sampled from FASTA databases (RefSeq, IMG/VR, custom)

**Abundance distributions:**
- **Log-normal** (realistic): Most viromes, some dominant species, long tail
- **Power-law** (highly uneven): Few hyper-dominant species, steep decline
- **Even** (synthetic control): All species equal abundance

### 2. Contamination Profile

**What it contains:**
- Host DNA fragments (human, mouse, or custom organism)
- rRNA sequences (bacterial 16S/23S, archaeal, eukaryotic 18S/28S)
- Reagent bacteria genomes (kit contamination: *Delftia*, *Ralstonia*, etc.)
- PhiX174 spike-in control

**Contamination levels:**

| Level | Host DNA | rRNA | Reagent Bacteria | PhiX | Total Contam | Use Case |
|-------|----------|------|------------------|------|--------------|----------|
| **clean** | 0.1% | 0.5% | 0.1% | 0.1% | 0.8% | Excellent VLP enrichment |
| **realistic** | 2.0% | 5.0% | 0.5% | 0.1% | 7.6% | Typical VLP-enriched virome |
| **heavy** | 10.0% | 15.0% | 2.0% | 0.1% | 27.1% | Poor VLP enrichment |
| **failed** | 15.0% | 20.0% | 4.0% | 0.1% | 39.1% | Complete VLP failure |

### 3. VLP Enrichment (Phase 2)

**What it does:**
- Models physical separation of virus-like particles from cells/debris
- Applies size-based retention (0.2-0.45 μm filtration)
- Applies nuclease treatment efficiency (removes free DNA)
- Applies family-specific enrichment factors (from ViromeQC literature)

**Effect on composition:**
- Viral fraction increases dramatically (50% → 95%)
- Host DNA decreases 10-100x (30% → 0.3-3%)
- rRNA decreases 5-20x (15% → 0.75-3%)
- Reagent bacteria partially removed
- Small viruses slightly enriched, large viruses slightly depleted

**See**: [VLP_ENRICHMENT_BIOLOGY.md](VLP_ENRICHMENT_BIOLOGY.md) for detailed explanation

### 4. Amplification Bias (Phase 2)

**What it does:**
- Models PCR or isothermal amplification effects on abundances
- Length-dependent bias (short genomes amplify better)
- GC-dependent bias (moderate GC optimal)
- Method-specific characteristics (RdAB, MDA, Linker, None)

**Methods:**

| Method | Length Bias | GC Bias | Stochasticity | Use Case |
|--------|-------------|---------|---------------|----------|
| **RdAB** | Strong | Moderate | Low | Standard virome protocol (40 cycles) |
| **MDA** | Minimal | Extreme | High | Low-biomass samples, marine viromes |
| **Linker** | Minimal | Weak | Low | Alternative amplification method |
| **None** | None | None | None | High-biomass samples, no amplification |

### 5. Library Prep Bias (Phase 2)

**What it does:**
- Models fragmentation method effects on sequence representation
- GC-dependent coverage for tagmentation
- Uniform coverage for mechanical shearing

**Methods:**
- **Mechanical shearing** (Covaris, sonication): Minimal bias, uniform coverage
- **Tagmentation** (Nextera, Tn5): GC bias, AT-rich underrepresented
- **Enzymatic** (NEBNext, KAPA): Intermediate bias

### 6. Platform Artifacts (Phase 2)

**What it does:**
- Adds platform-specific sequencing artifacts to reads
- Does NOT change genome abundances, modifies read sequences

**Artifacts:**

| Platform | PolyG Tails | Optical Duplicates | Index Hopping |
|----------|-------------|-------------------|---------------|
| **NovaSeq** | Yes (20%) | Yes (5%) | Yes (0.5%) |
| **MiSeq** | No | Low (1%) | Minimal |
| **HiSeq** | No | Low (1%) | Minimal |

---

## Order of Operations

### Critical Sequence Understanding

The order matters because **earlier stages affect what later stages act upon**:

```
BIOLOGICAL SAMPLE
        ↓
[Contamination is PRESENT]
        ↓
VLP ENRICHMENT ← Reduces contamination, enriches viruses
        ↓
ENRICHED SAMPLE (mostly viral)
        ↓
AMPLIFICATION ← Biases viral abundances
        ↓
AMPLIFIED DNA (biased abundances)
        ↓
LIBRARY PREP ← Fragments DNA, may introduce bias
        ↓
LIBRARY (ready for sequencing)
        ↓
SEQUENCING ← Adds platform artifacts
        ↓
RAW FASTQ READS
```

### Example: Following Host DNA Through Pipeline

```
Initial Sample (bulk metagenome):
├── Host DNA: 30% of total DNA

After VLP Enrichment (0.2 μm filtration + 95% nuclease):
├── Most host cells removed (too large)
├── Free host DNA degraded (95% efficiency)
├── Small host DNA fragments remain
└── Host DNA: 0.5% of total DNA (60x reduction)

After Amplification (RdAB 40 cycles):
├── Host DNA fragments amplify
├── But viral genomes may amplify better (shorter on average)
└── Host DNA: 0.3% of total DNA (viral preferentially amplified)

After Library Prep (mechanical shearing):
├── Host DNA fragmented uniformly
└── Host DNA: 0.3% of total DNA (no change)

After Sequencing (NovaSeq):
├── Some host DNA reads get polyG tails (20% of all reads)
└── Host DNA: 0.3% of reads in FASTQ
```

### Example: Comparing VLP vs Bulk from Same Starting Community

```python
# Create one starting community
base_viral_community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
base_contamination = create_contamination_profile('realistic', random_seed=42)

# VLP-ENRICHED WORKFLOW
vlp_composition = MockViromeComposition(
    name='vlp_enriched',
    viral_community=base_viral_community.copy(),  # Same community
    contamination_profile=base_contamination.copy(),  # Same contamination
    viral_fraction=0.50  # Start as bulk (50% viral)
)

vlp_enrichment = VLPEnrichment(filtration_cutoff_um=0.2, nuclease_efficiency=0.95)
vlp_enrichment.apply(vlp_composition)  # Now ~95% viral

vlp_output = generate_reads(vlp_composition, 'vlp_dataset', n_reads=1_000_000)

# BULK METAGENOME WORKFLOW (NO ENRICHMENT)
bulk_composition = MockViromeComposition(
    name='bulk_metagenome',
    viral_community=base_viral_community.copy(),  # Same community
    contamination_profile=base_contamination.copy(),  # Same contamination
    viral_fraction=0.50  # Stays as bulk (50% viral)
)

# NO VLP enrichment applied!
bulk_output = generate_reads(bulk_composition, 'bulk_dataset', n_reads=1_000_000)

# RESULT:
# - Both datasets have SAME viral species
# - But VERY DIFFERENT abundances
# - VLP: 95% viral reads, Bulk: 50% viral reads
# - VLP: 0.5% host DNA, Bulk: 30% host DNA
```

---

## Examples

### Example 1: Phase 1 Simple Dataset

```python
from viroforge.simulators import quick_generate

# One-liner - creates gut virome with realistic contamination
output = quick_generate(
    body_site='gut',
    contamination_level='realistic',
    n_reads=1_000_000,
    random_seed=42
)

# Output:
# - reads_R1.fastq.gz: Forward reads
# - reads_R2.fastq.gz: Reverse reads
# - ground_truth_abundance.tsv: True abundances
# - ground_truth_taxonomy.tsv: Taxonomic assignments
```

### Example 2: Phase 1 Custom Contamination

```python
from viroforge.core import create_body_site_profile, create_contamination_profile
from viroforge.utils import MockViromeComposition
from viroforge.simulators import generate_reads

# Create custom contamination
contamination = create_contamination_profile(
    profile='realistic',
    host_dna_pct=5.0,      # Override: 5% instead of 2%
    rrna_pct=10.0,         # Override: 10% instead of 5%
    phix_pct=0.5           # Override: 0.5% instead of 0.1%
)

# Create composition
composition = MockViromeComposition(
    name='custom_contamination',
    viral_community=create_body_site_profile('gut', n_genomes=50),
    contamination_profile=contamination,
    viral_fraction=0.85
)

# Generate reads
output = generate_reads(composition, 'custom_dataset', n_reads=1_000_000)
```

### Example 3: Phase 2 VLP vs Bulk Comparison (Planned)

```python
from viroforge.enrichment import VLPEnrichment, no_enrichment

# Create base composition (same for both)
base_community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
base_contamination = create_contamination_profile('realistic', random_seed=42)

# VLP-ENRICHED
vlp_comp = MockViromeComposition(
    'vlp', base_community.copy(), base_contamination.copy(), viral_fraction=0.50
)
VLPEnrichment(filtration_cutoff_um=0.2).apply(vlp_comp)
vlp_output = generate_reads(vlp_comp, 'vlp_dataset', n_reads=1_000_000)

# BULK METAGENOME
bulk_comp = MockViromeComposition(
    'bulk', base_community.copy(), base_contamination.copy(), viral_fraction=0.50
)
# No enrichment applied
bulk_output = generate_reads(bulk_comp, 'bulk_dataset', n_reads=1_000_000)

# Compare results to test if your pipeline handles them differently
```

### Example 4: Phase 2 Complete Workflow (Planned)

```python
from viroforge.enrichment import VLPEnrichment
from viroforge.amplification import RdABAmplification
from viroforge.library_prep import MechanicalShearing
from viroforge.platform import NovaSeqPlatform

# Create base composition
composition = MockViromeComposition(
    name='complete_workflow',
    viral_community=create_body_site_profile('gut', 50),
    contamination_profile=create_contamination_profile('realistic'),
    viral_fraction=0.50  # Bulk starting state
)

# Apply full processing pipeline
VLPEnrichment(filtration_cutoff_um=0.2, nuclease_efficiency=0.95).apply(composition)
RdABAmplification(cycles=40).apply(composition)
MechanicalShearing().apply(composition)

# Generate reads with platform artifacts
output = generate_reads(
    composition,
    'complete_dataset',
    n_reads=1_000_000,
    model='NovaSeq',
    add_artifacts=True  # PolyG tails, optical duplicates
)
```

---

## Summary

### Key Takeaways

1. **Order matters**: Contamination → Enrichment → Amplification → Library Prep → Sequencing → Artifacts

2. **Two types of modifications**:
   - Abundance-modifying (changes genome ratios)
   - Read-modifying (changes read sequences)

3. **Phase 1 is complete and working**: Viral community + contamination + read generation

4. **Phase 2 adds virome-specific features**: VLP enrichment, amplification bias, platform artifacts

5. **VLP enrichment is the key differentiator**: Transforms bulk metagenome into virome

6. **Flexibility at every stage**: Users can customize all parameters or use pre-defined templates

---

## See Also

- [VLP_ENRICHMENT_BIOLOGY.md](VLP_ENRICHMENT_BIOLOGY.md) - Detailed explanation of VLP enrichment biology
- [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) - Phase 2 implementation timeline
- [DESIGN_RATIONALE.md](DESIGN_RATIONALE.md) - Original design document
- [examples/](../examples/) - Code examples for all use cases

---

**Questions or feedback?** Open an issue at https://github.com/shandley/viroforge/issues
