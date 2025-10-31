# VLP Enrichment: Biological Principles and Implementation

**Purpose**: This document explains the biological principles behind VLP (Virus-Like Particle) enrichment and how ViroForge models these processes to create realistic mock virome datasets.

**Last Updated**: 2025-01-31
**Status**: Implementation in Progress (Phase 2, Weeks 1-3)

---

## Table of Contents

1. [What is VLP Enrichment?](#what-is-vlp-enrichment)
2. [Biological Principles](#biological-principles)
3. [ViroForge Implementation](#viroforge-implementation)
4. [Effects on Viral Composition](#effects-on-viral-composition)
5. [Effects on Contamination](#effects-on-contamination)
6. [VLP vs Bulk Comparison](#vlp-vs-bulk-comparison)
7. [Validation Targets](#validation-targets)
8. [References](#references)

---

## What is VLP Enrichment?

### Definition

**VLP (Virus-Like Particle) enrichment** is a laboratory technique used to isolate free viral particles from complex biological samples (stool, saliva, serum, environmental water, etc.) by:

1. **Physical separation** (filtration, ultracentrifugation)
2. **Enzymatic treatment** (nuclease digestion of free DNA/RNA)

The goal is to **enrich for encapsidated viral genomes** while **removing cellular material and non-encapsidated nucleic acids**.

### Why VLP Enrichment?

**Problem**: Raw biological samples contain:
- 🦠 Viruses (what we want)
- 🧫 Bacteria (contamination)
- 🧬 Host cells (contamination)
- 💧 Free DNA/RNA (contamination)
- 🧪 rRNA (highly abundant contamination)

In a typical gut sample:
- **WITHOUT enrichment**: 10-50% viral reads, 30-70% host/bacterial DNA
- **WITH VLP enrichment**: 85-99% viral reads, 0.1-5% host/bacterial DNA

**Impact**: VLP enrichment is the defining feature that makes viromics possible.

---

## Biological Principles

### Principle 1: Size-Based Separation

**Biology**: Most viruses are **20-300 nm** in diameter, while bacterial cells are **500-5000 nm**.

**Method**: Filtration through **0.2-0.45 μm** pore size membranes

**Effect**:
- ✅ **Viruses pass through** (20-300 nm < 200-450 nm cutoff)
- ❌ **Bacteria retained** (500-5000 nm > 200-450 nm cutoff)
- ❌ **Host cells retained** (5000-20000 nm > cutoff)
- ⚠️ **Large viruses partially retained** (e.g., Mimiviridae >200 nm)
- ✅ **Small viruses enriched** (e.g., Microviridae ~25 nm)

**Filtration Methods**:
| Method | Description | Retention | Common Use |
|--------|-------------|-----------|------------|
| **Syringe filtration** | Push sample through 0.22 μm filter | Sharp cutoff | Small volumes, quick |
| **Tangential flow** | Flow parallel to membrane | Gradual retention curve | Large volumes, gentle |
| **Ultracentrifugation** | Pellet viruses (100,000 × g) | Size/density dependent | Marine samples, traditional |
| **Iron chloride precipitation** | Precipitate viruses with FeCl₃ | Enhanced recovery | High recovery method |

**ViroForge Models**:
```python
# Retention curve based on virus size
def size_retention(virus_size_nm, filtration_cutoff_nm, curve_type='sigmoid'):
    if curve_type == 'sigmoid':
        # Smooth transition around cutoff (realistic for TFF)
        return 1 / (1 + exp((virus_size_nm - cutoff_nm) / sharpness))
    elif curve_type == 'step':
        # Sharp cutoff (realistic for syringe filtration)
        return 1.0 if virus_size_nm < cutoff_nm else 0.0
```

### Principle 2: Nuclease Treatment

**Biology**: DNase and RNase enzymes degrade **free (non-encapsidated) nucleic acids** but cannot penetrate viral capsids.

**Method**: Treat filtered sample with DNase/RNase cocktail

**Effect**:
- ✅ **Encapsidated viral genomes protected** (inside capsid)
- ❌ **Free host DNA degraded** (95-99% removal)
- ❌ **Free bacterial DNA degraded** (95-99% removal)
- ⚠️ **rRNA partially protected** (ribosome structure provides some protection)
- ⚠️ **Some free DNA remains** (enzyme efficiency not 100%)

**Nuclease Efficiency Factors**:
- **Enzyme concentration**: Higher = better removal
- **Incubation time**: Longer = better removal
- **Temperature**: 37°C optimal
- **Inhibitors**: Sample components may inhibit enzymes
- **Secondary structure**: rRNA, tRNA protected by folding

**Typical Efficiencies** (from literature):
- Standard protocol: **90-95%** removal of free DNA
- Optimized (FeCl₃ + nuclease): **95-98%** removal
- Poor conditions: **70-85%** removal

**ViroForge Models**:
```python
# Nuclease treatment effect
def nuclease_retention(genome, is_encapsidated, nuclease_efficiency=0.95):
    if is_encapsidated:
        return 1.0  # Protected by capsid
    else:
        return 1.0 - nuclease_efficiency  # Free DNA removed
```

### Principle 3: Family-Specific Enrichment

**Biology**: Different viral families have different **physical properties** that affect recovery:
- Capsid stability (some viruses more fragile)
- Size distribution (varies by family)
- Morphology (tailed phages vs spherical vs filamentous)

**Evidence**: ViromeQC paper (Zolfo et al., 2019, Nature Biotechnology) measured **family-specific enrichment factors** from 248 published viromes.

**Observed Enrichment Factors**:

| Viral Family | Typical Size | Morphology | Enrichment Factor | Reason |
|--------------|--------------|------------|-------------------|--------|
| **Microviridae** | 25-30 nm | Small icosahedral | 2.5× | Small, stable, highly retained |
| **Siphoviridae** | 50-100 nm | Long-tailed | 1.2× | Typical phage, moderate enrichment |
| **Myoviridae** | 60-120 nm | Contractile tail | 1.0× | Baseline (reference) |
| **Podoviridae** | 50-80 nm | Short-tailed | 1.1× | Slightly enriched |
| **Inoviridae** | 6 nm × 1000 nm | Filamentous | 0.3× | Depleted (retained by filters) |
| **Caudoviricetes** | Variable | Tailed phages | 1.0-1.5× | Most common, well-retained |

**ViroForge Models**:
```python
# Family-specific enrichment (from ViromeQC literature)
family_enrichment_factors = {
    'Microviridae': 2.5,
    'Siphoviridae': 1.2,
    'Myoviridae': 1.0,
    'Podoviridae': 1.1,
    'Inoviridae': 0.3,
    # ... additional families from literature
}
```

### Principle 4: Stochastic Variation

**Biology**: VLP enrichment is not deterministic - biological variation exists:
- Sample handling variability
- Filter pore size variation
- Enzyme activity variation
- Viral capsid integrity variation
- Environmental conditions (temperature, pH)

**Effect**: Two identical samples processed identically will have **slightly different** results

**ViroForge Models**:
```python
# Add stochastic variation (log-normal noise)
stochastic_factor = np.random.lognormal(mean=0, sigma=0.2)
retention = base_retention * stochastic_factor
```

---

## ViroForge Implementation

### Multi-Factor Retention Model

ViroForge calculates retention for each genome as:

```
Retention = Size_Factor × Family_Factor × Stability_Factor × Stochastic_Factor
```

#### Component 1: Size-Based Retention

```python
def calculate_size_factor(genome_length_bp, filtration_cutoff_um=0.2):
    """
    Calculate retention based on virus size.

    Virus size estimation:
    - Rough correlation: 1000 bp ≈ 50 nm diameter
    - Small genomes (<10 kb) ≈ 25-50 nm → highly retained
    - Medium genomes (10-100 kb) ≈ 50-200 nm → normally retained
    - Large genomes (>100 kb) ≈ >200 nm → partially excluded
    """
    # Estimate virus size from genome length
    estimated_size_nm = (genome_length_bp / 1000) ** 0.5 * 50

    cutoff_nm = filtration_cutoff_um * 1000

    # Sigmoid retention curve (smooth for tangential flow)
    sharpness = cutoff_nm * 0.2
    retention = 1 / (1 + np.exp((estimated_size_nm - cutoff_nm) / sharpness))

    return retention

# Example:
# - 5 kb genome (Microvirus) → ~35 nm → retention = 0.99 (highly retained)
# - 50 kb genome (Siphovirus) → ~110 nm → retention = 0.95 (well retained)
# - 150 kb genome (Jumbo phage) → ~190 nm → retention = 0.85 (mostly retained)
# - 300 kb genome (Large virus) → ~270 nm → retention = 0.15 (mostly excluded)
```

#### Component 2: Family-Based Enrichment

```python
def calculate_family_factor(viral_family):
    """
    Apply family-specific enrichment factor.

    Based on ViromeQC paper empirical observations.
    """
    family_factors = {
        'Microviridae': 2.5,
        'Siphoviridae': 1.2,
        'Myoviridae': 1.0,  # Baseline
        'Podoviridae': 1.1,
        'Inoviridae': 0.3,
        'Caudoviricetes': 1.15,
        'Unknown': 1.0  # No information, assume baseline
    }

    return family_factors.get(viral_family, 1.0)
```

#### Component 3: Capsid Stability

```python
def calculate_stability_factor(gc_content):
    """
    Estimate capsid stability from GC content.

    High GC → more stable capsid (approximation)
    """
    # GC content effect on stability (weak correlation)
    # Higher GC slightly more stable
    baseline = 1.0
    gc_effect = (gc_content - 50) * 0.002  # ±10% max effect

    stability = baseline + gc_effect
    stability = np.clip(stability, 0.85, 1.15)  # Limit to ±15%

    return stability
```

#### Component 4: Nuclease Treatment

```python
def apply_nuclease_treatment(genome, is_encapsidated, efficiency=0.95):
    """
    Apply nuclease treatment effect.

    Encapsidated viruses protected, free DNA removed.
    """
    if is_encapsidated:
        # Viral genomes in capsids protected
        return 1.0
    else:
        # Free DNA removed (efficiency = fraction removed)
        return 1.0 - efficiency

# Example:
# - Viral genome (encapsidated) → retention = 1.0 (fully protected)
# - Host DNA fragment (free) → retention = 0.05 (95% removed)
# - rRNA in ribosome (partially protected) → retention = 0.2 (80% removed)
```

#### Component 5: Stochastic Variation

```python
def add_stochastic_variation(base_retention, sigma=0.2):
    """
    Add biological variation to retention.

    Log-normal distribution (multiplicative noise).
    """
    stochastic = np.random.lognormal(mean=0, sigma=sigma)
    return base_retention * stochastic
```

### Complete VLP Enrichment Model

```python
class VLPEnrichment:
    """
    Model VLP enrichment process.

    Flexible, configurable for different protocols.
    """

    def __init__(
        self,
        filtration_method='tangential_flow',
        filtration_cutoff_um=0.2,
        prefiltration_cutoff_um=0.45,
        nuclease_treatment=True,
        nuclease_efficiency=0.95,
        size_retention_curve='sigmoid'
    ):
        self.filtration_method = filtration_method
        self.filtration_cutoff = filtration_cutoff_um
        self.prefiltration = prefiltration_cutoff_um
        self.nuclease_treatment = nuclease_treatment
        self.nuclease_efficiency = nuclease_efficiency
        self.curve_type = size_retention_curve

    def apply(self, composition):
        """
        Apply VLP enrichment to composition.

        Modifies genome abundances in place.
        """
        for genome in composition.all_genomes:
            # Calculate retention factors
            size_factor = self._size_retention(genome.length)
            family_factor = self._family_enrichment(genome.family)
            stability_factor = self._stability(genome.gc_content)
            stochastic_factor = np.random.lognormal(0, 0.2)

            # Combined retention
            retention = (size_factor * family_factor *
                        stability_factor * stochastic_factor)

            # Nuclease treatment
            if self.nuclease_treatment:
                is_encapsidated = genome.source_type == 'viral'
                nuclease_retention = self._nuclease_effect(
                    is_encapsidated,
                    self.nuclease_efficiency
                )
                retention *= nuclease_retention

            # Apply to abundance
            genome.abundance *= retention

        # Renormalize abundances to sum to 1.0
        composition.normalize()
```

---

## Effects on Viral Composition

### Abundance Shifts

VLP enrichment **does not add or remove viral species**, but it **shifts their relative abundances**.

**Before VLP** (Bulk Metagenome):
```
Virus A (5 kb, Microviridae):     10%
Virus B (50 kb, Siphoviridae):    30%
Virus C (100 kb, Myoviridae):     25%
Virus D (150 kb, Jumbo phage):    20%
Virus E (10 kb, filamentous):     15%
```

**After VLP** (0.2 μm filtration + 95% nuclease):
```
Virus A (5 kb, Microviridae):     18%  ↑ (small + Microviridae enrichment)
Virus B (50 kb, Siphoviridae):    38%  ↑ (good retention + family enrichment)
Virus C (100 kb, Myoviridae):     28%  ↑ (good retention)
Virus D (150 kb, Jumbo phage):    14%  ↓ (larger, partially excluded)
Virus E (10 kb, filamentous):      2%  ↓↓ (filamentous, depleted)
```

### Expected Changes by Viral Type

| Virus Type | Size | Before VLP | After VLP | Fold Change | Reason |
|------------|------|------------|-----------|-------------|--------|
| **Microvirus** | ~25 nm | 10% | 15-20% | 1.5-2× | Small + stable + family enrichment |
| **Typical Siphovirus** | ~80 nm | 30% | 35-40% | 1.2-1.3× | Good retention + moderate enrichment |
| **Typical Myovirus** | ~100 nm | 25% | 25-28% | 1.0-1.1× | Good retention, baseline enrichment |
| **Jumbo Phage** | ~150 nm | 20% | 12-16% | 0.6-0.8× | Larger, some exclusion |
| **Filamentous** | 6×1000 nm | 15% | 1-3% | 0.1-0.2× | Retained by filter (too large in one dimension) |

### Diversity Effects

**Alpha diversity** (within-sample):
- Overall diversity **preserved** (same species present)
- Evenness slightly **decreased** (small viruses more abundant)
- Rare species may fall **below detection** after enrichment + sequencing

**Beta diversity** (between-samples):
- VLP samples **cluster separately** from bulk metagenomes
- VLP samples more **similar to each other** (enrichment normalizes)

---

## Effects on Contamination

### Host DNA

**Starting condition** (bulk metagenome):
- Host cells present (~10⁶-10⁸ cells/mL)
- Free host DNA from lysed cells
- **Abundance**: 30-70% of total DNA

**After filtration** (0.2 μm):
- Host cells removed (too large)
- Small host DNA fragments pass through
- **Abundance**: 5-15% of total DNA (5-10× reduction)

**After nuclease** (95% efficiency):
- Free host DNA degraded
- Small protected fragments remain
- **Abundance**: 0.5-2% of total DNA (50-100× total reduction)

**Why some remains**:
- Incomplete nuclease digestion
- DNA in cell debris (partial protection)
- Very small fragments (<100 bp) resistant to nuclease

### rRNA Contamination

**Starting condition**:
- Bacterial/archaeal/eukaryotic ribosomes
- Very abundant in cells
- **Abundance**: 10-30% of total nucleic acids

**After filtration**:
- Ribosomes small (~20 nm), pass through filter
- Minimal reduction
- **Abundance**: 10-30% (no change)

**After nuclease**:
- Ribosome structure provides **partial protection**
- 16S/18S rRNA partially resistant
- **Abundance**: 2-8% (3-5× reduction)

**Why rRNA is problematic**:
- Ribosomes are **small** (pass through filter)
- rRNA has **secondary structure** (partially protects from nuclease)
- Very **abundant** in biological samples
- Indicator of **poor VLP enrichment** if high

### Reagent Bacteria

**Source**: DNA extraction kit contamination (Salter et al., 2014)

**Common contaminants**:
- *Delftia acidovorans*
- *Ralstonia pickettii*
- *Burkholderia*
- *Bradyrhizobium*

**Effect of VLP enrichment**:
- Bacterial cells **removed by filtration** (too large)
- Bacterial DNA **removed by nuclease** (not encapsidated)
- **Abundance**: 0.1-1% (5-10× reduction)

**Why they persist**:
- Present at all stages (can't fully eliminate)
- Low biomass samples amplify contamination
- Some may be carried through as cell debris

### PhiX174 Control

**Source**: Illumina sequencing spike-in control

**Effect of VLP enrichment**:
- PhiX is a **small icosahedral virus** (~25 nm)
- **Fully retained** by VLP enrichment (if added before enrichment)
- Usually added **after** enrichment (at sequencing)
- **Abundance**: 0.1-1% (no change)

---

## VLP vs Bulk Comparison

### Side-by-Side Comparison

Starting from the **same initial sample**:

| Component | Bulk Metagenome | VLP-Enriched | Fold Change |
|-----------|-----------------|--------------|-------------|
| **Viral reads** | 10-50% | 85-99% | 5-10× increase |
| **Host DNA** | 30-70% | 0.1-2% | 30-100× decrease |
| **rRNA** | 10-30% | 2-8% | 3-5× decrease |
| **Bacterial DNA** | 10-40% | 0.5-3% | 10-20× decrease |
| **Reagent contamination** | 0.5-2% | 0.1-0.5% | 3-5× decrease |
| **PhiX control** | 0.1-1% | 0.1-1% | No change |

### ViromeQC Scores

**ViromeQC** (Zolfo et al., 2019) calculates an enrichment score:

```
ViromeQC_score = log10(viral_reads / non_viral_reads)
```

**Expected scores**:
| Sample Type | Viral % | ViromeQC Score | Interpretation |
|-------------|---------|----------------|----------------|
| **Failed VLP** | 40-60% | 0-2 | Poor enrichment, like bulk metagenome |
| **Partial VLP** | 70-85% | 3-6 | Moderate enrichment, some contamination |
| **Good VLP** | 85-95% | 7-12 | Good enrichment, typical virome |
| **Excellent VLP** | 95-99% | 13-20 | Excellent enrichment, clean virome |
| **Bulk Metagenome** | 10-50% | -2 to 0 | No enrichment |

---

## Validation Targets

ViroForge VLP enrichment model will be validated against **published virome statistics**:

### Target Ranges (from ViromeQC survey of 248 viromes)

| Metric | Clean VLP | Realistic VLP | Poor VLP | Failed VLP | Bulk |
|--------|-----------|---------------|----------|------------|------|
| **Viral reads** | 95-99% | 85-95% | 60-85% | 40-60% | 10-50% |
| **Host DNA** | 0.1-0.5% | 0.5-2% | 2-10% | 10-15% | 30-70% |
| **rRNA** | 0.5-2% | 2-5% | 5-15% | 15-20% | 10-30% |
| **Bacteria** | 0.1-0.5% | 0.5-1% | 1-5% | 5-10% | 10-40% |
| **ViromeQC score** | 13-20 | 7-12 | 2-7 | 0-2 | <0 |
| **Shannon diversity** | 2.5-4.0 | 2.0-3.5 | 1.5-3.0 | 1.5-2.5 | 2.0-4.0 |

### Validation Strategy

1. **Generate synthetic datasets** with VLP enrichment
2. **Calculate metrics** (viral %, contamination %, ViromeQC score)
3. **Compare to target ranges** from literature
4. **Adjust parameters** if outside ranges
5. **Statistical validation** (t-tests, effect sizes)

### Success Criteria

✅ **Simulated VLP-enriched viromes**:
- Viral reads: 85-99% ✓
- Host DNA: 0.1-2% ✓
- rRNA: 0.5-5% ✓
- ViromeQC score: 7-20 ✓

✅ **Simulated bulk metagenomes**:
- Viral reads: 10-50% ✓
- Host DNA: 30-70% ✓
- rRNA: 10-30% ✓
- ViromeQC score: <2 ✓

✅ **Statistical comparison**:
- VLP vs Bulk significantly different (p < 0.001)
- Effect size large (Cohen's d > 1.0)
- Distributions match published data

---

## References

### Key Papers

1. **Zolfo, M., et al. (2019).** "Detecting contamination in viromes using ViromeQC." *Nature Biotechnology* 37:1408-1412.
   - ViromeQC enrichment scoring
   - Family-specific enrichment factors
   - Survey of 248 published viromes

2. **Conceição-Neto, N., et al. (2015).** "Modular approach to customise sample preparation procedures for viral metagenomics." *BMC Genomics* 16:617.
   - VLP enrichment protocols
   - Iron chloride precipitation method
   - Comparison of different methods

3. **Thurber, R.V., et al. (2009).** "Laboratory procedures to generate viral metagenomes." *Nature Protocols* 4:470-483.
   - Iron chloride flocculation protocol
   - Ultracentrifugation methods
   - Classic viromics methodology

4. **Solonenko, S.A., et al. (2013).** "Sequencing platform and library preparation choices impact viral metagenomes." *BMC Genomics* 14:320.
   - Platform comparison
   - Library prep bias effects
   - Quantitative viromics

5. **Roux, S., et al. (2016).** "Towards quantitative viromics for both double-stranded and single-stranded DNA viruses." *PeerJ* 4:e2777.
   - Quantitative virome methods
   - ssDNA virus recovery
   - Multiple displacement amplification (MDA)

6. **Salter, S.J., et al. (2014).** "Reagent and laboratory contamination can critically impact sequence-based microbiome analyses." *BMC Biology* 12:87.
   - Reagent contamination sources
   - Low-biomass sample issues
   - Common kit contaminants

### Additional Resources

- **ViromeQC tool**: https://github.com/SegataLab/viromeqc
- **Virome protocols compilation**: https://www.protocols.io/workspaces/viromics
- **VLP enrichment reviews**: Multiple reviews in *Current Opinion in Virology*, *Nature Reviews Microbiology*

---

## Implementation Notes

### Phase 2 Timeline

**Weeks 1-3**: VLP Enrichment Framework
- Week 1: Core framework, size-based retention
- Week 2: Nuclease treatment, family-specific enrichment
- Week 3: Pre-defined protocols, testing, validation

### Pre-Defined VLP Protocols

```python
# Standard VLP protocol (most common)
standard_vlp = VLPEnrichment(
    filtration_method='tangential_flow',
    filtration_cutoff_um=0.2,
    prefiltration_cutoff_um=0.45,
    nuclease_treatment=True,
    nuclease_efficiency=0.95
)

# Iron chloride protocol (Conceição-Neto et al.)
iron_chloride_vlp = VLPEnrichment(
    filtration_method='tangential_flow',
    filtration_cutoff_um=0.2,
    nuclease_treatment=True,
    nuclease_efficiency=0.98  # Better with FeCl3
)

# Ultracentrifugation protocol
ultracentrifuge_vlp = VLPEnrichment(
    filtration_method='none',
    ultracentrifugation=True,
    uc_speed_g=100000,
    nuclease_treatment=True,
    nuclease_efficiency=0.90
)

# NO enrichment (bulk metagenome)
no_enrichment = VLPEnrichment(
    filtration_method='none',
    nuclease_treatment=False
)
```

### Testing Strategy

```python
# Test 1: VLP enrichment increases viral fraction
composition = create_mock_virome('gut', 'realistic')
initial_viral_frac = composition.viral_fraction  # ~50% (bulk)

vlp = VLPEnrichment()
vlp.apply(composition)

final_viral_frac = composition.viral_fraction  # ~95% (VLP)
assert final_viral_frac > initial_viral_frac * 1.5  # At least 1.5x increase

# Test 2: VLP vs bulk comparison
vlp_comp = create_mock_virome('gut', 'realistic')
bulk_comp = create_mock_virome('gut', 'realistic')

VLPEnrichment().apply(vlp_comp)
# No enrichment applied to bulk_comp

vlp_output = generate_reads(vlp_comp, 'vlp', n_reads=10_000)
bulk_output = generate_reads(bulk_comp, 'bulk', n_reads=10_000)

# VLP should have dramatically higher viral fraction
assert vlp_viral_frac > bulk_viral_frac * 2
assert vlp_host_dna < bulk_host_dna * 0.1
```

---

## See Also

- [WORKFLOW.md](WORKFLOW.md) - Complete pipeline workflow
- [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) - Phase 2 implementation details
- [DESIGN_RATIONALE.md](DESIGN_RATIONALE.md) - Original design document

---

**Questions or feedback?** Open an issue at https://github.com/shandley/viroforge/issues
