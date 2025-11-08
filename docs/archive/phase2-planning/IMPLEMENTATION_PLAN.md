# ViroForge Implementation Plan: Community-Focused Phase 2

**Date**: 2025-01-30
**Approach**: Focused depth with broad applicability
**Timeline**: 2-3 months
**Goal**: Fully developed, tested, robust, and lab-agnostic virome simulator

---

## Strategic Direction

### User Requirements
✅ Focused implementation (not trying to do everything)
✅ Fully developed and tested
✅ Lab-agnostic (not specific to one workflow)
✅ Useful for broader virome community
✅ Robust and production-ready

### Modified Approach: "Focused Breadth"

**Philosophy**: Implement core virome-specific features **deeply and flexibly** rather than shallowly and broadly.

**Instead of**:
- ❌ Your lab's specific VLP protocol → One fixed implementation
- ❌ Your lab's specific RdAB protocol → One fixed implementation
- ❌ Your lab's specific NovaSeq setup → One fixed implementation

**We implement**:
- ✅ **VLP enrichment framework** → Configurable for any VLP protocol
- ✅ **Amplification bias framework** → Supports RdAB, MDA, linker-amplified, none
- ✅ **Platform artifact framework** → Supports NovaSeq, MiSeq, HiSeq, NextSeq
- ✅ **Library prep framework** → Supports shearing, tagmentation, enzymatic

**Key Difference**: Build flexible, parameterized implementations that work for many labs, not just one.

---

## Core Features to Implement (Tier 1)

### 1. VLP Enrichment Framework (Weeks 1-3)

**Goal**: Model VLP enrichment in a way that works for ANY lab's protocol

#### Flexible Design

```python
class VLPEnrichment:
    """
    Configurable VLP enrichment modeling.

    Works for different protocols:
    - Standard VLP (0.2 μm filtration + nuclease)
    - Ultracentrifugation-based
    - Tangential flow filtration
    - Commercial kits (Norgen, Zymo, etc.)
    """

    def __init__(
        self,
        filtration_method='tangential_flow',  # or 'syringe', 'centrifugal', 'none'
        filtration_cutoff_um=0.2,  # Configurable
        prefiltration_cutoff_um=0.45,  # Optional pre-filter
        nuclease_treatment=True,
        nuclease_efficiency=0.95,  # Configurable
        size_retention_curve='sigmoid',  # or 'step', 'linear'
        ultracentrifugation=False,
        uc_speed_g=100000,  # If using UC
    ):
        self.filtration_method = filtration_method
        self.filtration_cutoff = filtration_cutoff_um
        # ... store all parameters

    def apply(self, composition):
        """
        Apply VLP enrichment to composition.

        Adjusts genome abundances based on:
        - Genome size (retention curve)
        - Encapsidation status
        - Nuclease treatment efficiency
        - Method-specific biases
        """
        for genome in composition.all_genomes:
            # Size-based retention
            retention = self._calculate_retention(genome)

            # Nuclease treatment (if enabled)
            if self.nuclease_treatment and not genome.is_encapsidated:
                retention *= (1 - self.nuclease_efficiency)

            # Method-specific adjustments
            retention *= self._method_specific_bias(genome)

            genome.abundance *= retention

        composition.normalize()

    def _calculate_retention(self, genome):
        """Calculate size-based retention."""
        if self.size_retention_curve == 'sigmoid':
            # Smooth transition around cutoff
            return self._sigmoid_retention(genome.size_nm, self.filtration_cutoff)
        elif self.size_retention_curve == 'step':
            # Sharp cutoff
            return 1.0 if genome.size_nm < self.filtration_cutoff else 0.0
        # ... other curve types

    def _method_specific_bias(self, genome):
        """Method-specific recovery biases."""
        if self.filtration_method == 'tangential_flow':
            # TFF tends to have better large particle retention
            return self._tff_bias(genome)
        elif self.filtration_method == 'syringe':
            # Syringe filtration more binary
            return self._syringe_bias(genome)
        # ... other methods
```

#### Pre-defined Protocols (for convenience)

```python
# Standard VLP protocol (most common)
standard_vlp = VLPEnrichment(
    filtration_method='tangential_flow',
    filtration_cutoff_um=0.2,
    prefiltration_cutoff_um=0.45,
    nuclease_treatment=True,
    nuclease_efficiency=0.95,
)

# Iron chloride precipitation + filtration (Conceição-Neto et al.)
iron_chloride_vlp = VLPEnrichment(
    filtration_method='tangential_flow',
    filtration_cutoff_um=0.2,
    prefiltration_cutoff_um=0.8,
    nuclease_treatment=True,
    nuclease_efficiency=0.98,  # Higher with FeCl3
)

# Ultracentrifugation-based
ultracentrifuge_vlp = VLPEnrichment(
    filtration_method='none',
    ultracentrifugation=True,
    uc_speed_g=100000,
    nuclease_treatment=True,
    nuclease_efficiency=0.90,
)

# Commercial kit (Norgen)
norgen_vlp = VLPEnrichment(
    filtration_method='kit',
    kit_type='norgen',
    nuclease_treatment=False,  # Some kits skip this
)

# NO enrichment (bulk metagenome)
no_enrichment = VLPEnrichment(
    filtration_method='none',
    nuclease_treatment=False,
)
```

#### Testing Strategy

**Test with multiple protocols**:
1. Standard VLP (0.2 μm TFF + nuclease)
2. Bulk metagenome (no enrichment)
3. Ultracentrifugation protocol
4. Commercial kit protocol

**Validation**:
- VLP should have 90-99% viral reads
- Bulk should have 10-50% viral reads (depending on sample)
- Results should match literature ranges

---

### 2. Amplification Bias Framework (Weeks 4-6)

**Goal**: Support multiple amplification methods used in the field

#### Flexible Design

```python
class AmplificationMethod:
    """Base class for amplification methods."""

    def apply(self, composition):
        """Apply method-specific bias to composition."""
        raise NotImplementedError


class RdABAmplification(AmplificationMethod):
    """
    RdAB (Random RT + Sequenase + PCR) amplification.

    Used by: Many virome labs
    Characteristics:
    - Length-dependent bias (exponential, favors short)
    - GC-dependent bias (moderate GC optimal)
    - Cycle-dependent scaling
    """

    def __init__(
        self,
        cycles=40,
        length_bias_strength=1.0,  # Configurable intensity
        gc_bias_strength=1.0,
        optimal_gc=0.50,
        gc_tolerance=0.15,
    ):
        self.cycles = cycles
        self.length_bias_strength = length_bias_strength
        self.gc_bias_strength = gc_bias_strength
        self.optimal_gc = optimal_gc
        self.gc_tolerance = gc_tolerance

    def apply(self, composition):
        for genome in composition.all_genomes:
            # Length bias (exponential advantage for short genomes)
            length_factor = self._length_efficiency(genome.length)

            # GC bias (quadratic, peak at optimal_gc)
            gc_factor = self._gc_efficiency(genome.gc_content)

            # Combined bias
            genome.abundance *= (length_factor ** self.cycles) * gc_factor

        composition.normalize()


class MDAAmplification(AmplificationMethod):
    """
    Multiple Displacement Amplification (φ29 polymerase).

    Used by: Low-biomass virome studies
    Characteristics:
    - EXTREME GC bias (10-1000x)
    - Random, highly stochastic
    - Chimera formation (5-30%)
    """

    def __init__(
        self,
        amplification_time_hours=4,
        gc_bias_strength=2.0,  # Stronger than PCR
        stochasticity=0.3,  # Random variation
        chimera_rate=0.15,
    ):
        self.amplification_time = amplification_time_hours
        self.gc_bias_strength = gc_bias_strength
        self.stochasticity = stochasticity
        self.chimera_rate = chimera_rate

    def apply(self, composition):
        for genome in composition.all_genomes:
            # Extreme GC bias
            gc_factor = self._extreme_gc_bias(genome.gc_content)

            # Random stochastic variation
            stochastic_factor = np.random.lognormal(0, self.stochasticity)

            # Combined
            genome.abundance *= gc_factor * stochastic_factor

        composition.normalize()

        # Generate chimeras (if generating reads)
        # ... chimera generation logic


class LinkerAmplification(AmplificationMethod):
    """
    Linker-based amplification.

    Used by: Some viral metagenomics protocols
    Characteristics:
    - Less bias than RdAB or MDA
    - No length bias (linkers on all fragments)
    - Moderate GC bias
    """

    def __init__(
        self,
        cycles=25,  # Usually fewer cycles than RdAB
        gc_bias_strength=0.5,  # Weaker than RdAB
    ):
        self.cycles = cycles
        self.gc_bias_strength = gc_bias_strength


class NoAmplification(AmplificationMethod):
    """No amplification (high-biomass samples)."""

    def apply(self, composition):
        # No bias
        pass
```

#### Pre-defined Amplification Protocols

```python
# Common amplification methods
amplification_methods = {
    'rdab_40': RdABAmplification(cycles=40),
    'rdab_30': RdABAmplification(cycles=30),
    'mda_standard': MDAAmplification(amplification_time_hours=4),
    'mda_overnight': MDAAmplification(amplification_time_hours=16),
    'linker': LinkerAmplification(cycles=25),
    'none': NoAmplification(),
}
```

#### Testing Strategy

**Test scenarios**:
1. RdAB 40-cycle (strong bias)
2. RdAB 30-cycle (moderate bias)
3. MDA (extreme bias)
4. Linker (minimal bias)
5. No amplification (no bias)

**Validation**:
- Length bias should scale with cycle number
- MDA should show extreme GC bias
- No amplification should preserve input abundances

---

### 3. Platform Artifact Framework (Weeks 7-8)

**Goal**: Support different sequencing platforms and their artifacts

#### Flexible Design

```python
class SequencingPlatform:
    """Base class for sequencing platforms."""

    def add_artifacts(self, reads):
        """Add platform-specific artifacts to reads."""
        raise NotImplementedError


class NovaSeqPlatform(SequencingPlatform):
    """
    NovaSeq (2-channel chemistry).

    Artifacts:
    - PolyG tails (20% of reads)
    - Optical duplicates (5% rate)
    - Index hopping (0.1-1%)
    """

    def __init__(
        self,
        polyg_rate=0.2,
        polyg_min_length=10,
        polyg_max_length=100,
        optical_duplicate_rate=0.05,
        index_hopping_rate=0.005,
    ):
        self.polyg_rate = polyg_rate
        self.polyg_min_length = polyg_min_length
        self.polyg_max_length = polyg_max_length
        self.optical_duplicate_rate = optical_duplicate_rate
        self.index_hopping_rate = index_hopping_rate

    def add_artifacts(self, reads):
        # Add polyG tails
        reads = self._add_polyg_tails(reads)

        # Add optical duplicates
        reads = self._add_optical_duplicates(reads)

        # Add index hopping (if multiplexed)
        if self.multiplexed:
            reads = self._add_index_hopping(reads)

        return reads


class MiSeqPlatform(SequencingPlatform):
    """
    MiSeq (4-channel chemistry).

    Artifacts:
    - No polyG tails (different chemistry)
    - Lower optical duplicates (random flow cell)
    - Different quality profile
    """

    def __init__(
        self,
        optical_duplicate_rate=0.01,  # Lower than NovaSeq
    ):
        self.optical_duplicate_rate = optical_duplicate_rate

    def add_artifacts(self, reads):
        # Only optical duplicates (no polyG)
        reads = self._add_optical_duplicates(reads)
        return reads


class HiSeqPlatform(SequencingPlatform):
    """
    HiSeq (4-channel chemistry, random flow cell).

    Similar to MiSeq but different read lengths/quality.
    """
    pass


class NextSeqPlatform(SequencingPlatform):
    """
    NextSeq (2-channel chemistry, like NovaSeq).

    Similar polyG issues to NovaSeq.
    """
    pass
```

#### Testing Strategy

**Test platforms**:
1. NovaSeq (with polyG + optical dups)
2. MiSeq (no polyG, low optical dups)
3. HiSeq (no polyG, low optical dups)

**Validation**:
- NovaSeq should have ~20% reads with polyG
- MiSeq should have 0% reads with polyG
- Optical duplicate rates should match expected values

---

### 4. Library Preparation Framework (Weeks 7-8)

**Goal**: Support different fragmentation methods

```python
class LibraryPrepMethod:
    """Base class for library prep methods."""

    def fragment_genome(self, genome, target_insert_size):
        """Fragment genome with method-specific bias."""
        raise NotImplementedError


class MechanicalShearing(LibraryPrepMethod):
    """
    Mechanical shearing (Covaris, sonication).

    Characteristics:
    - Relatively uniform
    - Minimal sequence bias
    - Insert size normally distributed
    """

    def fragment_genome(self, genome, target_insert_size):
        # Uniform fragmentation
        fragments = self._uniform_fragments(genome, target_insert_size)
        return fragments


class Tagmentation(LibraryPrepMethod):
    """
    Tagmentation (Nextera, Tn5 transposase).

    Characteristics:
    - GC-dependent fragmentation
    - AT-rich regions underrepresented
    - Sequence motif preferences
    """

    def __init__(
        self,
        gc_bias_strength=0.8,
    ):
        self.gc_bias_strength = gc_bias_strength

    def fragment_genome(self, genome, target_insert_size):
        # GC-biased fragmentation
        fragments = self._gc_biased_fragments(genome, target_insert_size)
        return fragments


class EnzymaticFragmentation(LibraryPrepMethod):
    """
    Enzymatic fragmentation (NEBNext, KAPA).

    Characteristics:
    - More uniform than tagmentation
    - Less bias than mechanical shearing
    - Some sequence preferences
    """
    pass
```

---

## Implementation Architecture

### Modular, Composable Design

```python
from viroforge.enrichment import VLPEnrichment, no_enrichment
from viroforge.amplification import RdABAmplification, MDAAmplification, NoAmplification
from viroforge.library_prep import MechanicalShearing, Tagmentation
from viroforge.platform import NovaSeqPlatform, MiSeqPlatform

# Example 1: Your lab's workflow
workflow = ViromePipeline(
    enrichment=VLPEnrichment(
        filtration_cutoff_um=0.2,
        nuclease_efficiency=0.95,
    ),
    amplification=RdABAmplification(cycles=40),
    library_prep=MechanicalShearing(),
    platform=NovaSeqPlatform(polyg_rate=0.2),
)

# Example 2: MDA-based workflow (common in marine viromics)
marine_workflow = ViromePipeline(
    enrichment=VLPEnrichment(
        filtration_cutoff_um=0.2,
        nuclease_efficiency=0.98,
    ),
    amplification=MDAAmplification(amplification_time_hours=16),
    library_prep=Tagmentation(),
    platform=MiSeqPlatform(),
)

# Example 3: Bulk metagenome (no enrichment, no amplification)
bulk_metagenome = ViromePipeline(
    enrichment=no_enrichment,
    amplification=NoAmplification(),
    library_prep=MechanicalShearing(),
    platform=NovaSeqPlatform(),
)

# Example 4: High-biomass sample (no amplification needed)
high_biomass_virome = ViromePipeline(
    enrichment=VLPEnrichment(),
    amplification=NoAmplification(),  # Enough DNA, skip amplification
    library_prep=Tagmentation(),
    platform=NovaSeqPlatform(),
)

# Generate reads
composition = create_mock_virome(...)
workflow.process(composition)
output = generate_reads(composition, ...)
```

### Pre-defined Workflow Templates

```python
# Pre-defined workflows for convenience
WORKFLOWS = {
    # VLP workflows
    'vlp_rdab_novaseq': ViromePipeline(
        VLPEnrichment(), RdABAmplification(40),
        MechanicalShearing(), NovaSeqPlatform()
    ),
    'vlp_mda_miseq': ViromePipeline(
        VLPEnrichment(), MDAAmplification(),
        Tagmentation(), MiSeqPlatform()
    ),

    # Bulk metagenome workflows
    'bulk_no_amp': ViromePipeline(
        no_enrichment, NoAmplification(),
        MechanicalShearing(), NovaSeqPlatform()
    ),

    # Marine virome workflows
    'marine_vlp_mda': ViromePipeline(
        VLPEnrichment(nuclease_efficiency=0.98),
        MDAAmplification(amplification_time_hours=16),
        Tagmentation(), MiSeqPlatform()
    ),
}

# Easy to use
output = quick_generate(
    body_site='gut',
    workflow='vlp_rdab_novaseq',  # Pre-defined template
    n_reads=1_000_000,
)
```

---

## Documentation Strategy

### 1. Method Comparison Guide

Create documentation showing how different methods affect results:

```markdown
# Method Comparison Guide

## VLP vs Bulk

Generate two datasets from same composition:

```python
composition = create_gut_virome(n_genomes=50)

# VLP-enriched
vlp_output = generate_reads(
    composition, workflow='vlp_rdab_novaseq', n_reads=1M
)

# Bulk metagenome
bulk_output = generate_reads(
    composition, workflow='bulk_no_amp', n_reads=1M
)

# Compare:
# - VLP should have 90-99% viral reads
# - Bulk should have 10-50% viral reads
# - Host contamination drastically different
```

## RdAB vs MDA

```python
# RdAB (your lab)
rdab_output = generate_reads(
    composition,
    amplification=RdABAmplification(cycles=40),
    ...
)

# MDA (alternative)
mda_output = generate_reads(
    composition,
    amplification=MDAAmplification(),
    ...
)

# Compare:
# - RdAB: moderate length/GC bias
# - MDA: extreme GC bias, high stochasticity
```
```

### 2. Protocol Gallery

Document real-world protocols from literature:

```markdown
# Published Virome Protocols

## Protocol 1: Standard Gut Virome (Conceição-Neto et al. 2015)
- VLP enrichment: 0.2 μm TFF + nuclease
- Amplification: Linker-based
- Library: Nextera tagmentation
- Platform: MiSeq 2x250

## Protocol 2: Marine Virome (Brum et al. 2015)
- VLP enrichment: 0.22 μm filtration + FeCl3 precipitation
- Amplification: MDA (16h)
- Library: Nextera
- Platform: HiSeq 2x150

## Protocol 3: Low-Biomass Clinical (YOUR LAB)
- VLP enrichment: 0.2 μm TFF + nuclease
- Amplification: RdAB (40 cycles)
- Library: Mechanical shearing
- Platform: NovaSeq 2x150

... etc for 10+ published protocols
```

### 3. Tutorial Notebooks

Jupyter notebooks demonstrating different use cases:

1. **Tutorial 1**: VLP vs Bulk Comparison
2. **Tutorial 2**: Amplification Method Comparison
3. **Tutorial 3**: Platform Artifact Analysis
4. **Tutorial 4**: Custom Workflow Design
5. **Tutorial 5**: Validating Your QC Pipeline

---

## Testing & Validation Strategy

### 1. Unit Tests (Per Component)

```python
def test_vlp_enrichment():
    """Test VLP enrichment reduces bacterial DNA."""
    composition = create_composition(
        viral_genomes=50,
        bacterial_genomes=100,
    )

    vlp = VLPEnrichment()
    vlp.apply(composition)

    viral_frac = composition.viral_fraction
    assert 0.90 < viral_frac < 0.99  # VLP should be 90-99% viral


def test_rdab_length_bias():
    """Test RdAB favors short genomes."""
    composition = create_composition()

    short_genome = [g for g in composition if g.length < 10000][0]
    long_genome = [g for g in composition if g.length > 50000][0]

    initial_ratio = short_genome.abundance / long_genome.abundance

    rdab = RdABAmplification(cycles=40)
    rdab.apply(composition)

    final_ratio = short_genome.abundance / long_genome.abundance

    assert final_ratio > initial_ratio * 5  # Short should be >5x enriched
```

### 2. Integration Tests (Full Workflows)

```python
def test_vlp_rdab_novaseq_workflow():
    """Test complete workflow generates expected output."""
    composition = create_gut_virome(n_genomes=50)

    output = generate_reads(
        composition,
        workflow='vlp_rdab_novaseq',
        n_reads=10_000,
    )

    # Validate output
    assert output['r1'].exists()
    assert output['r2'].exists()

    # Check viral fraction (should be high for VLP)
    gt = pd.read_csv(output['ground_truth'])
    viral_abundance = gt[gt['source'] == 'viral']['abundance'].sum()
    assert 0.90 < viral_abundance < 0.99

    # Check for artifacts
    reads = parse_fastq(output['r1'])
    polyg_reads = count_polyg_tails(reads)
    assert 0.15 < polyg_reads / len(reads) < 0.25  # ~20% with polyG
```

### 3. Comparison Tests (Method Validation)

```python
def test_vlp_vs_bulk_comparison():
    """Test VLP enrichment dramatically increases viral fraction."""
    composition = create_gut_virome(n_genomes=50)

    # VLP
    vlp_output = generate_reads(
        composition, workflow='vlp_rdab_novaseq', n_reads=10_000
    )

    # Bulk
    bulk_output = generate_reads(
        composition, workflow='bulk_no_amp', n_reads=10_000
    )

    # Compare
    vlp_gt = pd.read_csv(vlp_output['ground_truth'])
    bulk_gt = pd.read_csv(bulk_output['ground_truth'])

    vlp_viral = vlp_gt[vlp_gt['source'] == 'viral']['abundance'].sum()
    bulk_viral = bulk_gt[bulk_gt['source'] == 'viral']['abundance'].sum()

    # VLP should be 2-10x higher viral fraction
    assert vlp_viral > bulk_viral * 2
    assert vlp_viral > 0.90
    assert bulk_viral < 0.50
```

### 4. Literature Validation

Compare simulated data to published virome statistics:

```python
def test_gut_virome_matches_literature():
    """
    Test simulated gut virome matches published statistics.

    Literature ranges (from ViromeQC survey):
    - VLP-enriched gut viromes: 85-99% viral
    - Bacterial contamination: 0.1-10%
    - Host DNA: 0.01-5%
    """
    composition = create_gut_virome(
        n_genomes=50,
        contamination='realistic',
    )

    output = generate_reads(
        composition, workflow='vlp_rdab_novaseq', n_reads=100_000
    )

    gt = pd.read_csv(output['ground_truth'])

    viral = gt[gt['source'] == 'viral']['abundance'].sum()
    bacteria = gt[gt['genome_type'] == 'reagent_bacteria']['abundance'].sum()
    host = gt[gt['genome_type'] == 'host_dna']['abundance'].sum()

    # Should match literature ranges
    assert 0.85 < viral < 0.99
    assert 0.001 < bacteria < 0.10
    assert 0.0001 < host < 0.05
```

---

## Timeline & Milestones

### Week 1-3: VLP Enrichment
- **Week 1**: Core framework, size-based retention
- **Week 2**: Nuclease treatment, method-specific bias
- **Week 3**: Pre-defined protocols, testing, documentation

**Milestone 1**: Can generate VLP vs bulk comparison datasets

### Week 4-6: Amplification Bias
- **Week 4**: RdAB implementation (length + GC bias)
- **Week 5**: MDA implementation (extreme bias + chimeras)
- **Week 6**: Linker amplification, no-amp, testing

**Milestone 2**: Can compare amplification methods

### Week 7-8: Platform Artifacts & Library Prep
- **Week 7**: NovaSeq artifacts (polyG + optical dups)
- **Week 7.5**: MiSeq/HiSeq (no polyG)
- **Week 8**: Tagmentation bias, testing

**Milestone 3**: Can model complete workflows end-to-end

### Week 9-10: Integration & Validation
- **Week 9**: Integration testing, workflow templates
- **Week 10**: Literature validation, comparison tests

**Milestone 4**: All tests passing, validated against literature

### Week 11-12: Documentation & Examples
- **Week 11**: Tutorial notebooks, protocol gallery
- **Week 12**: Method comparison guide, manuscript outline

**Milestone 5**: Complete documentation, ready for publication

---

## Success Criteria

### Technical Criteria

✅ **Modularity**: Each component independently testable
✅ **Flexibility**: Support multiple methods per component
✅ **Validation**: Match literature ranges for all methods
✅ **Documentation**: Clear examples for each workflow
✅ **Testing**: >90% code coverage, all integration tests passing

### Community Criteria

✅ **Lab-agnostic**: Works for many labs' protocols
✅ **Method coverage**: Covers most common virome methods
✅ **Reproducibility**: Same parameters → same output
✅ **Extensibility**: Easy to add new methods
✅ **Usability**: Simple API, pre-defined templates

### Scientific Criteria

✅ **Realistic**: Matches published virome statistics
✅ **Validated**: Comparison to real data shows similarity
✅ **Comprehensive**: Covers major sources of bias
✅ **Well-documented**: Methods clearly described
✅ **Publishable**: Strong enough for Nature Biotechnology

---

## Deliverables

### Code
1. `viroforge/enrichment.py` - VLP enrichment framework
2. `viroforge/amplification.py` - Amplification bias framework
3. `viroforge/library_prep.py` - Library prep methods
4. `viroforge/platform.py` - Platform artifact framework
5. `viroforge/workflows.py` - Pre-defined workflow templates
6. Updated `simulators/illumina.py` - Integration with new modules

### Tests
1. `tests/test_enrichment.py` - VLP enrichment tests
2. `tests/test_amplification.py` - Amplification bias tests
3. `tests/test_platform.py` - Platform artifact tests
4. `tests/test_workflows.py` - Integration tests
5. `tests/test_literature_validation.py` - Literature comparison

### Documentation
1. `docs/VLP_ENRICHMENT.md` - VLP enrichment guide
2. `docs/AMPLIFICATION_METHODS.md` - Amplification comparison
3. `docs/PLATFORM_ARTIFACTS.md` - Platform artifact guide
4. `docs/WORKFLOW_GALLERY.md` - Published protocol gallery
5. `docs/METHOD_COMPARISON.md` - Method comparison guide

### Tutorials
1. `tutorials/01_vlp_vs_bulk.ipynb` - VLP comparison
2. `tutorials/02_amplification_methods.ipynb` - Amplification comparison
3. `tutorials/03_platform_artifacts.ipynb` - Artifact analysis
4. `tutorials/04_custom_workflows.ipynb` - Custom workflow design
5. `tutorials/05_qc_pipeline_validation.ipynb` - QC validation

---

## Comparison to Original Options

### This Plan vs Original Options

| Aspect | Option A | Option B | Option C | **THIS PLAN** |
|--------|----------|----------|----------|---------------|
| **Timeline** | 1 month | 3 months | 2 months | **2-3 months** |
| **VLP Enrichment** | Basic | Comprehensive | Your protocol only | **Flexible framework** |
| **Amplification** | None | RdAB + MDA | RdAB only | **RdAB + MDA + Linker** |
| **Platforms** | NovaSeq | All platforms | NovaSeq | **All platforms** |
| **Library Prep** | None | All methods | Shearing | **Shearing + Tagmentation** |
| **Lab-specific** | Yes | No | **YES** | **NO (agnostic)** |
| **Community-useful** | Limited | **YES** | Limited | **YES** |
| **Publication** | Good | Excellent | Very Good | **Excellent** |
| **Validation** | Minimal | Comprehensive | Deep (one workflow) | **Comprehensive (all methods)** |

**Key Difference**: This plan implements the same core features as Option B (comprehensive) but with the focused, deep approach of Option C (tested and robust).

---

## Why This Approach is Best

### 1. Community Impact
- ✅ Works for VLP-enriched viromes (most common)
- ✅ Works for bulk metagenomes (important comparison)
- ✅ Works for MDA-based workflows (marine, low-biomass)
- ✅ Works for different platforms (NovaSeq, MiSeq, HiSeq)
- ✅ Not tied to any one lab's specific protocols

### 2. Scientific Rigor
- ✅ Multiple methods per component (shows breadth)
- ✅ Validated against literature (shows accuracy)
- ✅ Comprehensive testing (shows robustness)
- ✅ Well-documented (shows transparency)

### 3. Publication Potential
- ✅ Comprehensive enough for high-impact journal
- ✅ Novel enough (no other virome simulator like this)
- ✅ Useful enough (community will actually use it)
- ✅ Validated enough (can demonstrate accuracy)

### 4. Future-Proof
- ✅ Modular design (easy to add new methods)
- ✅ Flexible parameters (customizable for any protocol)
- ✅ Well-tested (changes won't break things)
- ✅ Community-ready (others can contribute)

---

## Next Steps

### Immediate (This Week)
1. Review this implementation plan
2. Approve approach and timeline
3. Decide on any modifications

### Next Week
1. Begin VLP enrichment framework implementation
2. Set up testing infrastructure
3. Create initial documentation structure

### Ongoing
1. Weekly progress updates
2. Continuous testing and validation
3. Literature review for parameter validation

---

## Questions for You

1. **Timeline acceptable?** 2-3 months for complete implementation?
2. **Method coverage sufficient?** VLP enrichment, RdAB, MDA, NovaSeq artifacts?
3. **Testing approach okay?** Unit tests + integration tests + literature validation?
4. **Documentation plan good?** Tutorials, protocol gallery, method comparison?
5. **Any specific methods to prioritize?** Or any to deprioritize?

---

**Ready to start implementation when you approve this plan!**
