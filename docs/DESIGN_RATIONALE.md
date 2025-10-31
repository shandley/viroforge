# ViroForge Design Rationale and Development Plan

**Project**: ViroForge - Mock Metavirome Data Generator
**Created**: 2025-01-30
**Purpose**: Comprehensive design document capturing the thinking, research, and rationale behind ViroForge

---

## Table of Contents

1. [Project Origin and Motivation](#project-origin-and-motivation)
2. [Literature Review and Gap Analysis](#literature-review-and-gap-analysis)
3. [System Architecture](#system-architecture)
4. [Technical Feasibility Assessment](#technical-feasibility-assessment)
5. [Use Cases and Applications](#use-cases-and-applications)
6. [Development Roadmap](#development-roadmap)
7. [Community Impact and Publication Strategy](#community-impact-and-publication-strategy)
8. [Integration with Existing Ecosystem](#integration-with-existing-ecosystem)

---

## Project Origin and Motivation

### Background Context

This project emerged from the need to validate a virome QC pipeline (lab-virome-QC) developed for:
- VLP (Virus-Like Particle) enriched samples
- RdAB amplification protocol (Random reverse transcription + Sequenase + 40-cycle PCR)
- Illumina NovaSeq sequencing (2-channel chemistry)
- Mechanical shearing library preparation

### The Problem

The lab needs to:
1. **Validate QC pipeline performance** - How do we know our QC flags are correct?
2. **Test downstream analysis tools** - Read-based and assembly-based classification
3. **Benchmark methods** - Compare VLP vs bulk, different protocols, different platforms
4. **Train the team** - Learn collaborative software development through a meaningful project

### The Solution

Create a comprehensive mock metavirome data generator that produces:
- Realistic sequencing reads with known ground truth
- Virome-specific contamination profiles
- Library prep and sequencing artifacts
- Multiple body-site scenarios
- Complete metadata for validation

---

## Literature Review and Gap Analysis

### Existing Metagenomic Simulators

**CAMISIM** (Fritz et al., 2019, Microbiome)
- **Strengths**: Comprehensive bacterial metagenome simulator, models abundance profiles, strain diversity, generates 2nd and 3rd generation reads
- **Limitations**: Bacterial-focused, no VLP enrichment modeling, no virus-specific contamination
- **Status**: Current gold standard for bacterial metagenome simulation
- **Citation**: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0633-6

**InSilicoSeq** (Gourlé et al., 2019, Bioinformatics)
- **Strengths**: Fast Illumina read simulation, realistic error models (HiSeq, NovaSeq), customizable from BAM files
- **Limitations**: General-purpose, no library prep biases, no contamination modeling
- **Status**: Popular for basic read simulation
- **Citation**: https://academic.oup.com/bioinformatics/article/35/3/521/5055123

**MeSS** (2024)
- **Strengths**: Recent tool (Dec 2024), metagenomic sample generation
- **Limitations**: Compared favorably to CAMISIM but still bacterial-focused
- **Status**: Emerging alternative to CAMISIM

**NeSSM, RandomReadsMG, ART, NEAT, wgsim, Mason, FASTQSim, MetaSim**
- Various general metagenomic/genomic read simulators
- None are virome-specific

### Current Virome Benchmarking Approaches

**Physical Synthetic Communities**
- Recent benchmarking studies (2023-2024) used 60-72 viral isolates
- **Advantages**: True biological ground truth
- **Disadvantages**:
  - Expensive ($10,000+)
  - Time-consuming (weeks to prepare)
  - Limited complexity (usually <100 species)
  - Not scalable
  - Lab-dependent variability
- **Citations**:
  - "Benchmarking of virome metagenomic analysis approaches" (2023)
  - Multiple 2024 benchmarking studies

**Simple Simulated Contigs**
- Many studies just fragment reference genomes
- **Advantages**: Quick, easy
- **Disadvantages**:
  - No sequencing artifacts
  - No contamination
  - No library prep biases
  - Unrealistic

**Real Published Datasets**
- Some studies use published virome data
- **Advantages**: Real biology
- **Disadvantages**:
  - Unknown ground truth
  - Can't control parameters
  - Limited availability
  - Can't test edge cases

### Identified Gaps

**No existing tool models**:
1. ✗ VLP enrichment biases (size selection, capsid stability)
2. ✗ Virome-specific contamination (PhiX, reagent bacteria, host DNA patterns)
3. ✗ RdAB or other viral amplification biases
4. ✗ NovaSeq-specific artifacts (polyG tails, optical duplicates)
5. ✗ Body-site specific viral compositions
6. ✗ VLP vs bulk comparative datasets
7. ✗ Realistic virome complexity and diversity

**Recent papers explicitly call for better simulation tools**:
- 2024 benchmarking studies note lack of standardized test datasets
- ViromeQC paper (2019) mentions need for controlled validation
- Multiple papers use inadequate simulation approaches due to lack of alternatives

---

## System Architecture

### Layered Modular Design

```
ViroForge Architecture
│
├─ Layer 1: VIRAL COMMUNITY COMPOSITION
│  ├─ Genome Sampling
│  │  ├─ RefSeq Viral database
│  │  ├─ IMG/VR (environmental viruses)
│  │  ├─ Gut Virome Database (body-site specific)
│  │  └─ ICTV taxonomy integration
│  │
│  ├─ Abundance Modeling
│  │  ├─ Log-normal distribution (realistic)
│  │  ├─ Power-law distribution (highly uneven)
│  │  ├─ Even distribution (synthetic control)
│  │  └─ Custom distributions
│  │
│  └─ Body-Site Profiles
│     ├─ Gut: crAssphage, Microviridae, Siphoviridae
│     ├─ Oral: Streptococcus/Actinomyces phages
│     ├─ Skin: Propionibacterium/Staphylococcus phages
│     ├─ Respiratory: environment-specific
│     └─ Environmental: diverse phage communities
│
├─ Layer 2: CONTAMINATION PROFILES
│  ├─ Host DNA
│  │  ├─ Human (GRCh38)
│  │  ├─ Mouse (GRCm39)
│  │  ├─ Other model organisms
│  │  └─ Abundance: 0.1% (clean) to 15% (VLP failure)
│  │
│  ├─ rRNA Contamination
│  │  ├─ Bacterial (SILVA database)
│  │  ├─ Archaeal (SILVA)
│  │  ├─ Eukaryotic (SILVA)
│  │  └─ Rfam non-coding RNAs
│  │
│  ├─ Reagent Bacteria (from literature)
│  │  ├─ Delftia acidovorans
│  │  ├─ Ralstonia spp.
│  │  ├─ Burkholderia spp.
│  │  ├─ Bradyrhizobium spp.
│  │  └─ Based on Salter et al. 2014 contamination survey
│  │
│  └─ PhiX174 Control
│     ├─ NC_001422.1 reference
│     └─ Typical 0.1-5% spike-in
│
├─ Layer 3: LIBRARY PREP SIMULATION
│  ├─ VLP Enrichment Modeling
│  │  ├─ Size selection (0.22-0.45μm filters)
│  │  ├─ Viral family enrichment factors (from ViromeQC)
│  │  ├─ Capsid stability effects
│  │  └─ VLP vs bulk comparative mode
│  │
│  ├─ RdAB Amplification Bias
│  │  ├─ PCR efficiency model: eff(L) = base_eff × e^(-L/λ)
│  │  │  where λ ≈ 1000 bp (length bias parameter)
│  │  ├─ 40-cycle exponential amplification
│  │  ├─ Random priming biases (slight GC bias)
│  │  └─ Favors fragments <1kb, under-represents >10kb
│  │
│  └─ Fragmentation
│     ├─ Mechanical shearing (random, minimal bias)
│     ├─ Fragment size distribution: mean 300bp, SD 50bp
│     └─ No sequence-specific bias (unlike tagmentation)
│
└─ Layer 4: SEQUENCING ARTIFACTS
   ├─ NovaSeq 2-Channel Chemistry
   │  ├─ PolyG tail insertion (quality-dependent)
   │  ├─ High-quality G calls when synthesis fails
   │  └─ Length: 10-100bp polyG tails at read ends
   │
   ├─ Optical Duplicates
   │  ├─ Patterned flow cell artifacts
   │  ├─ Tile-position specific patterns
   │  └─ Typical 5-15% duplicate rate
   │
   ├─ Quality Score Profiles
   │  ├─ NovaSeq error model (InSilicoSeq)
   │  ├─ Read-end quality degradation
   │  └─ Position-specific error rates
   │
   └─ Index Hopping
      ├─ Cross-sample contamination
      └─ Low level (0.1-1%)
```

### Software Components

```
viroforge/
├── core/
│   ├── community.py          # Viral community composition and sampling
│   ├── contamination.py      # Contamination profile generation
│   ├── enrichment.py         # VLP enrichment modeling
│   └── artifacts.py          # Sequencing/library prep artifacts
│
├── profiles/
│   ├── body_sites/           # Pre-defined body site compositions
│   │   ├── gut_healthy.yaml
│   │   ├── gut_dysbiotic.yaml
│   │   ├── oral_saliva.yaml
│   │   ├── skin.yaml
│   │   └── respiratory.yaml
│   │
│   ├── contaminants/         # Known contamination databases
│   │   ├── reagent_bacteria.fasta
│   │   ├── common_host_seqs.fasta
│   │   └── phix174.fasta
│   │
│   └── sequencing/           # Platform-specific error models
│       ├── novaseq.yaml
│       ├── miseq.yaml
│       └── hiseq.yaml
│
├── simulators/
│   ├── base.py              # Base simulator class
│   ├── illumina.py          # Illumina-specific (NovaSeq, MiSeq, HiSeq)
│   └── pacbio.py            # Future: long-read support
│
├── utils/
│   ├── genome_sampler.py    # Sample from viral databases
│   ├── abundance.py         # Abundance distribution modeling
│   ├── metrics.py           # Calculate ground truth metrics
│   └── validation.py        # Compare to real datasets
│
└── cli.py                    # Command-line interface
```

---

## Technical Feasibility Assessment

### Straightforward Components (✅ Low Risk)

**Viral Genome Sampling**
- Sample genomes from RefSeq viral, IMG/VR based on taxonomic profiles
- Assign abundances following specified distributions (log-normal, power-law)
- **Implementation**: Standard BioPython + NumPy sampling
- **Complexity**: Low

**Contamination Mixing**
- Spike in contaminant sequences at specified percentages
- Mix host DNA, rRNA, bacterial genomes, PhiX at defined ratios
- **Implementation**: Simple read pooling with abundance weights
- **Complexity**: Low

**Read Generation**
- Leverage InSilicoSeq or ART for base read generation with error profiles
- Both tools well-established, documented, maintained
- **Implementation**: Wrapper around existing tools
- **Complexity**: Low

**Fragment Size Modeling**
- Mechanical shearing creates random fragments following normal distribution
- Mean ~300bp, SD ~50bp for typical library prep
- **Implementation**: NumPy random sampling
- **Complexity**: Low

### Moderately Challenging Components (⚠️ Medium Risk)

**VLP Enrichment Bias**
- Model which viruses get enriched by filtration/purification
- **Data sources**:
  - ViromeQC paper enrichment factors
  - Size-based selection (0.22μm filters exclude >200nm)
  - Viral family-specific stability factors
- **Implementation**:
  - Weight viral genomes by: size × family_enrichment × stability
  - Add stochastic variation (not all viruses behave identically)
- **Complexity**: Medium (requires empirical parameterization)

**PCR Amplification Length Bias (RdAB)**
- 40-cycle PCR favors short fragments over long
- **Model**:
  ```python
  def amplification_efficiency(length_bp):
      base_eff = 0.95  # PCR efficiency per cycle
      lambda_param = 1000  # length scale (bp)
      return base_eff * np.exp(-length_bp / lambda_param)

  def amplification_copies(length_bp, cycles=40):
      eff = amplification_efficiency(length_bp)
      return (1 + eff) ** cycles
  ```
- **Implementation**:
  1. Fragment genomes first (300bp mean)
  2. Calculate amplification factor for each fragment
  3. Weight abundance by amplification factor
- **Complexity**: Medium (need to validate model)

**PolyG Tail Insertion (NovaSeq Artifact)**
- 2-channel chemistry: dark cycles called as high-quality G
- Occurs when sequencing quality drops (read ends, low complexity)
- **Implementation**:
  1. Identify read ends or low-quality regions
  2. Insert GGGGGG... sequences (10-100bp)
  3. Assign high quality scores (Q30+) to polyG
- **Complexity**: Medium (need to model when it occurs)

**Optical Duplicates (Patterned Flow Cells)**
- Duplicates share tile/position but different sequences
- Pattern: Same tile, close X/Y coordinates
- **Implementation**:
  1. Select reads for duplication (5-15% of total)
  2. Duplicate with same tile/position metadata
  3. Add slight variation in sequence (1-2 errors)
- **Complexity**: Medium (requires read metadata tracking)

### Novel Contributions (🆕 Requires Development)

**Virome-Specific Contamination Profiles**
- Curate from literature: which contaminants appear in viromes?
- Create tiered profiles: clean, moderate, heavy contamination
- **Data sources**:
  - Salter et al. 2014 (reagent contamination)
  - ViromeQC paper (contamination patterns)
  - Recent virome studies (typical contamination levels)
- **Complexity**: Medium (mostly curation, straightforward implementation)

**Body-Site Specific Viral Profiles**
- Mine published virome studies for realistic compositions
- Create representative profiles for each body site
- **Data sources**:
  - Gut Virome Database
  - Published gut virome studies (crAssphage abundance, etc.)
  - Oral, skin, respiratory virome papers
- **Complexity**: Medium (data mining + curation)

**Integrated Ground Truth Metadata**
- Track: which reads → which genomes → which taxa → which biological origin
- Enable validation of entire analysis pipelines
- **Output formats**:
  - Abundance matrices (taxa × samples)
  - Read origin mappings (read_id → genome_id → taxonomy)
  - Expected assembly graphs (which reads assemble together)
  - Expected QC metrics (ViromeQC score, host %, rRNA %)
- **Complexity**: Medium (database design + tracking)

---

## Use Cases and Applications

### 1. QC Pipeline Validation (Primary Use Case)

**Validate lab-virome-QC pipeline performance**:

```bash
# Generate test dataset with known contamination
viroforge create \
  --profile gut_virome_vlp \
  --host-dna 12.0 \
  --rrna 15.0 \
  --output test_high_contamination/

# Run through QC pipeline
snakemake --use-conda --cores 8

# Check: Pipeline should flag FAIL for high host contamination
# Verify metrics match ground truth
```

**Specific test scenarios**:

| Scenario | Purpose | Expected Outcome |
|----------|---------|------------------|
| **PolyG removal test** | Generate 20% polyG-contaminated reads | fastp removes polyG, FastQC shows clean |
| **ViromeQC validation** | Generate VLP (score 15) vs failed (score 2) | Correct pass/fail flagging |
| **Host depletion** | Generate 5% vs 15% host DNA | Different QC outcomes |
| **rRNA removal** | Generate high rRNA contamination | BBDuk removes it effectively |
| **Optical duplicates** | Generate 10% optical dupes | Clumpify removes them |

### 2. Downstream Analysis Tool Testing

**Read-based classification**:
- Know which reads should classify to which taxa
- Calculate precision, recall, F1 for Kraken2, Kaiju, VirMAP
- Test sensitivity to contamination levels

**Assembly-based analysis**:
- Know which reads should assemble together
- Verify expected contigs form
- Calculate assembly completeness (N50, L50 vs expected)
- Test binning accuracy (which contigs belong together)

**Functional annotation**:
- Know true gene content
- Validate DRAM-v, VirSorter2 annotations
- Calculate annotation accuracy

### 3. Protocol Comparison Studies

**VLP vs Bulk Metagenome**:
```bash
# Same viral community, different enrichment
viroforge create --profile gut_virome --vlp-enrichment true --output vlp_sample/
viroforge create --profile gut_virome --vlp-enrichment false --output bulk_sample/

# Compare: viral recovery, contamination levels, diversity metrics
```

**Amplification Method Comparison**:
- RdAB (40 cycles) vs reduced cycles (25)
- With vs without PCR amplification
- Compare length bias effects

**Sequencing Platform Comparison**:
- NovaSeq vs MiSeq artifacts
- Short-read vs long-read (future)

### 4. Power Analysis and Study Design

**How many reads needed?**
- Generate datasets with varying depth (1M, 5M, 10M, 50M reads)
- Assess: viral detection sensitivity, diversity recovery, assembly completeness
- Determine minimum sequencing depth for research question

**Diversity requirements**:
- Generate varying complexity (10, 50, 100, 500 species)
- Test: can tools detect rare viruses? How does complexity affect analysis?

### 5. Benchmarking and Methods Development

**Standardized test datasets**:
- Community can use same datasets for comparisons
- Enable "CAVI Challenge" (Critical Assessment of Virome Interpretation)
- Similar to CAMI challenges for metagenomics

**Training datasets for machine learning**:
- Labeled data for viral classifier training
- Known taxonomic assignments
- Known functional annotations

### 6. Educational Tool

**Teaching virome analysis concepts**:
- Students can see how parameters affect outcomes
- Visualize effects of contamination, biases
- Learn analysis pipelines with known ground truth

---

## Development Roadmap

### Phase 1: Core Functionality (Months 1-3)

**Objectives**:
- Establish basic genome sampling and contamination mixing
- Integration with InSilicoSeq for read generation
- Ground truth metadata output
- 1-2 basic scenarios working

**Deliverables**:
- [x] Project structure and repository
- [ ] `core/community.py` - Sample viral genomes from RefSeq
- [ ] `core/contamination.py` - Add host DNA, rRNA, PhiX
- [ ] `simulators/illumina.py` - Wrapper around InSilicoSeq
- [ ] `utils/genome_sampler.py` - Database sampling functions
- [ ] `utils/abundance.py` - Abundance distribution models
- [ ] CLI interface - Basic `viroforge create` command
- [ ] Ground truth outputs - Abundance tables, read mappings
- [ ] Example scenario: "gut_virome_clean"

**Success Metric**: Can generate basic mock virome dataset with contamination

### Phase 2: Virome-Specific Features (Months 4-6)

**Objectives**:
- Implement VLP enrichment modeling
- Add RdAB PCR length bias
- Create body-site specific profiles
- Implement NovaSeq artifacts

**Deliverables**:
- [ ] `core/enrichment.py` - VLP enrichment factors
- [ ] PCR amplification bias model in `core/artifacts.py`
- [ ] Body-site profiles: gut, oral, skin, respiratory
- [ ] PolyG tail insertion
- [ ] Optical duplicate generation
- [ ] 10+ pre-defined scenarios
- [ ] Comparison mode: VLP vs bulk

**Success Metric**: Can generate realistic VLP vs bulk comparative datasets

### Phase 3: Validation and Refinement (Months 7-9)

**Objectives**:
- Validate simulated data matches real virome characteristics
- Statistical comparison to published datasets
- User testing with lab analysis pipelines
- Documentation and tutorials

**Deliverables**:
- [ ] Validation scripts - Compare to real datasets
- [ ] Statistical tests - Alpha diversity, composition, quality metrics
- [ ] Integration tests - Run through lab-virome-QC
- [ ] Comprehensive documentation
- [ ] Jupyter notebook tutorials
- [ ] Example datasets on Zenodo
- [ ] Benchmarking results

**Success Metric**: Simulated data statistically indistinguishable from real viromes

### Phase 4: Publication and Release (Months 10-11)

**Objectives**:
- Prepare manuscript
- Version 1.0 release
- Community announcement

**Deliverables**:
- [ ] Methods manuscript drafted
- [ ] Pre-print on bioRxiv
- [ ] Code review and cleanup
- [ ] Version 1.0 release with DOI
- [ ] Documentation website
- [ ] Conference presentation preparation
- [ ] Social media/community announcement

**Success Metric**: Published, citable tool ready for community use

---

## Community Impact and Publication Strategy

### Why the Field Needs This

**Current limitations identified in 2024 literature**:

1. **Benchmarking studies lack standardized datasets**
   - Each study creates one-off test data
   - Not comparable across studies
   - Not shared or reusable

2. **Physical synthetic communities are not scalable**
   - Recent studies limited to 60-72 viral isolates
   - Expensive and time-consuming
   - Cannot test edge cases or extreme scenarios

3. **Existing simulators don't model virome-specific features**
   - CAMISIM is bacterial-focused
   - No VLP enrichment modeling
   - No viral contamination profiles

4. **Lack of ground truth for validation**
   - Real datasets have unknown composition
   - Cannot calculate accuracy metrics
   - Difficult to validate new methods

### ViroForge Addresses These Needs

✅ **Standardized benchmarking**: Community-wide test datasets
✅ **Scalable**: Generate unlimited datasets instantly
✅ **Virome-specific**: Models VLP enrichment, viral artifacts
✅ **Complete ground truth**: Validate entire analysis pipelines
✅ **Reproducible**: Same parameters → same dataset
✅ **Flexible**: Test edge cases, extreme scenarios

### Publication Strategy

**Target Journals** (in priority order):
1. **Nature Biotechnology** - High-impact methods (ViromeQC published here)
2. **Bioinformatics** - Standard venue for bioinformatics tools
3. **Microbiome** - Microbiome methods and applications
4. **Genome Biology** - Genomic methods and software
5. **GigaScience** - Large-scale biological data tools

**Manuscript Structure**:
1. **Introduction**: Need for virome simulation, limitations of current approaches
2. **Methods**: ViroForge architecture, algorithms, validation approach
3. **Results**:
   - Validation against real datasets
   - Benchmarking example studies
   - Comparison to physical synthetic communities
4. **Discussion**: Applications, future directions, community impact
5. **Availability**: GitHub repository, documentation, example datasets

**Pre-print Strategy**:
- bioRxiv pre-print before journal submission
- Gather community feedback
- Early adopters can start using tool
- Increases visibility and citations

**Conference Presentations**:
- ISME (International Society for Microbial Ecology)
- ASM Microbe (American Society for Microbiology)
- BOSC (Bioinformatics Open Source Conference)
- ISMB (Intelligent Systems for Molecular Biology)

**Community Engagement**:
- Twitter/X announcement thread
- Blog post on lab website
- Add to bio.tools, OMICtools registries
- Mention in virome analysis communities (forums, Slack channels)

### Expected Impact

**Citation potential**: HIGH
- Becomes standard for virome tool validation
- Similar to CAMISIM for bacterial metagenomics
- Cited by every benchmarking study
- Used for training ML models

**Adoption potential**: HIGH
- Solves real problem (no alternatives)
- Easy to use (pre-built scenarios)
- Well-documented
- Open source (MIT license)

**Long-term vision**:
- "CAVI Challenge" - Community-wide benchmarking competition
- Standard reference datasets (like HMP for microbiome)
- Integration into virome analysis pipelines
- Teaching tool for virome analysis courses

---

## Integration with Existing Ecosystem

### Upstream Integration (Data Sources)

**Viral Genome Databases**:
- **RefSeq Viral**: NCBI comprehensive viral genomes
- **IMG/VR**: Environmental viral genomes
- **Gut Virome Database**: Body-site specific
- **ICTV Taxonomy**: Proper viral classification
- **User-provided**: Custom FASTA files

**Contamination Databases**:
- **Host genomes**: Ensembl, NCBI (human, mouse, etc.)
- **SILVA**: rRNA database
- **Reagent contaminants**: Curated from literature
- **PhiX174**: Standard Illumina control

### Tool Integration (Read Generation)

**Primary**: InSilicoSeq
- Well-maintained, documented
- Realistic Illumina error models
- Fast, parallelizable
- Python integration

**Alternative**: ART (Art_Illumina)
- Widely used
- Good error models
- C++ implementation (fast)

**Future**: CAMISIM integration
- For bacterial metagenome component
- Hybrid viral + bacterial simulations

### Downstream Integration (Analysis Tools)

**QC Tools** (immediate validation):
- lab-virome-QC (user's pipeline)
- ViromeQC (enrichment scoring)
- FastQC, MultiQC (quality assessment)

**Taxonomic Classification**:
- Kraken2, Bracken (k-mer based)
- Kaiju (protein-level)
- VirMAP (viral mapping)
- Centrifuge (sequence similarity)

**Assembly Tools**:
- MEGAHIT (short-read)
- metaSPAdes (metagenome-specific)
- metaFlye (long-read)
- hybridSPAdes (hybrid)

**Viral Detection**:
- VirSorter2
- VIBRANT
- VirFinder
- DeepVirFinder
- PhaMer (transformer-based)

**Binning**:
- vRhyme (viral-specific)
- VAMB
- MetaBAT2

**Functional Annotation**:
- DRAM-v (viral metabolism)
- VirSorter2 (viral genes)
- Prokka, Prodigal (gene calling)

### Output Format Compatibility

**FASTQ** - Standard paired-end reads
```
@read_id_1 viral_genome=NC_001422|taxonomy=Phix|position=1234
ACTG...
+
IIII...
```

**Ground Truth Metadata**:

1. **Abundance table** (TSV - phyloseq compatible)
```
taxonomy,sample1,sample2,sample3
Microviridae;Phix,1500,2300,1800
Siphoviridae;crAssphage,45000,52000,48000
```

2. **Read origin mapping** (TSV)
```
read_id,genome_id,taxonomy,start_pos,end_pos,source_type
read_001,NC_001422,Phix,1234,1384,viral
read_002,GRCh38_chr1,Human,5000000,5000150,host_contaminant
```

3. **Expected QC metrics** (JSON)
```json
{
  "viromeqc_enrichment": 15.2,
  "host_dna_percent": 2.1,
  "rrna_percent": 4.8,
  "phix_percent": 0.1,
  "total_reads": 10000000,
  "viral_reads": 9300000
}
```

4. **Assembly graph expectations** (GFA)
- Which reads should assemble together
- Expected contig sequences
- Expected binning results

---

## Technical Implementation Details

### VLP Enrichment Model

**Basis**: ViromeQC enrichment factors + physical principles

```python
def calculate_vlp_enrichment(virus_genome):
    """
    Calculate VLP enrichment factor for a viral genome.

    Factors:
    1. Size: 0.22μm filter excludes particles >220nm
    2. Family: Different families have different enrichment (ViromeQC)
    3. Stability: Capsid stability affects recovery
    4. Stochastic: Add biological variation
    """
    # Size-based enrichment
    size_nm = estimate_virus_size(virus_genome.length)
    if size_nm > 220:
        size_factor = 0.01  # Mostly excluded
    elif size_nm < 50:
        size_factor = 1.5   # Highly enriched
    else:
        size_factor = 1.0   # Normal

    # Family-based enrichment (from ViromeQC paper)
    family_factors = {
        'Microviridae': 2.5,      # Small, stable
        'Siphoviridae': 1.2,      # Typical tailed phage
        'Myoviridae': 1.0,        # Typical
        'Podoviridae': 1.1,       # Slightly enriched
        'Inoviridae': 0.3,        # Filamentous, depleted
        # ... more families
    }
    family_factor = family_factors.get(virus_genome.family, 1.0)

    # Stability factor (approximate)
    gc_content = calculate_gc(virus_genome.sequence)
    stability_factor = 1.0 + (gc_content - 0.5) * 0.2  # Higher GC = more stable

    # Stochastic variation
    stochastic = np.random.lognormal(0, 0.2)

    # Combined enrichment
    enrichment = size_factor * family_factor * stability_factor * stochastic

    return enrichment
```

### PCR Amplification Bias Model

**Basis**: Exponential amplification with length-dependent efficiency

```python
def apply_pcr_bias(fragments, cycles=40, base_efficiency=0.95, lambda_bp=1000):
    """
    Model PCR amplification bias favoring short fragments.

    Args:
        fragments: List of genome fragments with abundances
        cycles: Number of PCR cycles (40 for RdAB)
        base_efficiency: Maximum PCR efficiency (0.95)
        lambda_bp: Length scale for efficiency decay (1000 bp)

    Returns:
        Fragments with adjusted abundances
    """
    for fragment in fragments:
        # Length-dependent efficiency
        length = len(fragment.sequence)
        efficiency = base_efficiency * np.exp(-length / lambda_bp)

        # Exponential amplification
        amplification_factor = (1 + efficiency) ** cycles

        # Adjust abundance
        fragment.abundance *= amplification_factor

    # Renormalize to sum to 1
    total = sum(f.abundance for f in fragments)
    for f in fragments:
        f.abundance /= total

    return fragments
```

### PolyG Tail Insertion Model

**Basis**: NovaSeq 2-channel chemistry dark cycles

```python
def insert_polyg_tails(reads, frequency=0.2, min_length=10, max_length=100):
    """
    Insert polyG tails at read ends (NovaSeq artifact).

    Args:
        reads: List of sequencing reads
        frequency: Fraction of reads affected (0.2 = 20%)
        min_length: Minimum polyG length
        max_length: Maximum polyG length

    Returns:
        Reads with polyG tails inserted
    """
    for read in reads:
        if np.random.random() < frequency:
            # PolyG length (varies)
            polyg_length = np.random.randint(min_length, max_length)

            # Insert at read end (where quality drops)
            polyg_seq = 'G' * polyg_length

            # High quality scores for polyG (that's the artifact!)
            polyg_qual = 'I' * polyg_length  # Q40

            # Append or insert
            if np.random.random() < 0.5:
                # Append to end
                read.sequence += polyg_seq
                read.quality += polyg_qual
            else:
                # Replace low-quality end
                read.sequence = read.sequence[:-polyg_length] + polyg_seq
                read.quality = read.quality[:-polyg_length] + polyg_qual

    return reads
```

### Ground Truth Tracking System

**Design**: Unique identifiers with metadata

```python
class ReadMetadata:
    """Track complete origin information for each read."""

    def __init__(self, read_id, source_genome, taxonomy, position, read_type):
        self.read_id = read_id
        self.source_genome = source_genome  # Genome ID
        self.taxonomy = taxonomy             # Full taxonomic path
        self.position = position             # Start position in genome
        self.read_type = read_type          # 'viral', 'host', 'rrna', 'phix', etc.
        self.fragment_id = None             # For assembly tracking
        self.vlp_enriched = False           # Was this enriched by VLP?
        self.pcr_amplified = False          # Was this amplified by PCR?
        self.has_polyg = False              # Does it have polyG tail?
        self.is_optical_dup = False         # Is it an optical duplicate?

    def to_dict(self):
        """Export as dictionary for JSON/TSV."""
        return {
            'read_id': self.read_id,
            'source_genome': self.source_genome,
            'taxonomy': self.taxonomy,
            'position': self.position,
            'read_type': self.read_type,
            'fragment_id': self.fragment_id,
            'vlp_enriched': self.vlp_enriched,
            'pcr_amplified': self.pcr_amplified,
            'has_polyg': self.has_polyg,
            'is_optical_duplicate': self.is_optical_dup
        }

# Usage:
metadata_db = []
for read in generated_reads:
    meta = ReadMetadata(
        read_id=read.id,
        source_genome=read.source,
        taxonomy=read.taxonomy,
        position=read.start_pos,
        read_type=read.contamination_type
    )
    metadata_db.append(meta)

# Export ground truth
import pandas as pd
df = pd.DataFrame([m.to_dict() for m in metadata_db])
df.to_csv('ground_truth_read_origins.tsv', sep='\t', index=False)
```

---

## Validation Strategy

### Statistical Comparison to Real Datasets

**Metrics to validate**:

1. **Alpha diversity metrics**
   - Shannon diversity
   - Simpson diversity
   - Species richness
   - Compare: Simulated gut virome vs real gut virome studies

2. **Taxonomic composition**
   - Family-level relative abundances
   - Genus-level composition
   - Presence of marker viruses (e.g., crAssphage in gut)
   - Statistical test: Jensen-Shannon divergence

3. **Read quality metrics**
   - Quality score distributions
   - Read length distributions
   - GC content distribution
   - Per-base quality scores

4. **Contamination patterns**
   - Host DNA percentages
   - rRNA percentages
   - PhiX detection rates
   - Compare to ViromeQC paper contamination survey

5. **Assembly statistics**
   - Contig length distributions
   - N50, L50 values
   - Genome recovery rates
   - Completeness estimates

### Empirical Validation Protocol

```python
# 1. Generate simulated gut virome
viroforge create --profile gut_virome_vlp --output sim_gut/

# 2. Download real gut virome data (e.g., from SRA)
# SRA accession: PRJNA1234567

# 3. Run identical analysis pipeline on both
snakemake --use-conda sim_gut/
snakemake --use-conda real_gut/

# 4. Compare statistics
python compare_datasets.py sim_gut/ real_gut/

# 5. Statistical tests
# - t-test for continuous metrics (diversity, coverage)
# - Chi-square for categorical (taxonomic composition)
# - Kolmogorov-Smirnov for distributions (read length, quality)

# 6. Criteria for validation:
# - No significant difference (p > 0.05) for most metrics
# - Effect sizes small (Cohen's d < 0.5)
# - Distributions overlap substantially
```

### Integration Testing with lab-virome-QC

**Test suite for QC pipeline validation**:

```bash
#!/bin/bash
# integration_test.sh

# Test 1: Clean VLP sample (should PASS all QC)
viroforge create --profile gut_virome_clean --output test1/
snakemake --use-conda --cores 8 test1/
assert_qc_pass test1/results/reports/sample_qc_flags.tsv

# Test 2: High host contamination (should FAIL host check)
viroforge create --profile gut_virome_vlp --host-dna 15.0 --output test2/
snakemake --use-conda --cores 8 test2/
assert_qc_fail_host test2/results/reports/sample_qc_flags.tsv

# Test 3: Failed VLP (should FAIL enrichment check)
viroforge create --profile gut_virome_bulk --output test3/
snakemake --use-conda --cores 8 test3/
assert_qc_fail_enrichment test3/results/reports/sample_qc_flags.tsv

# Test 4: Heavy polyG (should be removed by fastp)
viroforge create --profile novaseq_polyg_heavy --output test4/
snakemake --use-conda --cores 8 test4/
assert_polyg_removed test4/results/fastqc/trimmed/

# Test 5: Optical duplicates (should be removed by Clumpify)
viroforge create --profile optical_dup_high --output test5/
snakemake --use-conda --cores 8 test5/
assert_duplicates_removed test5/results/clumpify/

echo "All integration tests passed!"
```

---

## Risks and Mitigation Strategies

### Risk 1: Over-Complexity

**Risk**: Tool becomes too complex, users can't use it
**Probability**: Medium
**Impact**: High (limits adoption)

**Mitigation**:
- Start simple, add features incrementally
- Provide sensible defaults (80% of users use defaults)
- Pre-built scenarios (users don't set 100 parameters)
- Clear documentation with examples
- Tiered interface: simple CLI for beginners, config files for advanced users

### Risk 2: Poor Validation

**Risk**: Simulated data doesn't match real virome characteristics
**Probability**: Low (with proper validation)
**Impact**: High (tool not trusted)

**Mitigation**:
- Extensive validation against published datasets
- Statistical comparison (t-tests, KS-tests, effect sizes)
- User testing with real analysis pipelines
- Iterative refinement based on feedback
- Transparency: publish all validation results
- Open-source: community can verify

### Risk 3: Database Maintenance Burden

**Risk**: Viral databases constantly updating, tool becomes outdated
**Probability**: Medium
**Impact**: Medium

**Mitigation**:
- Design to use user-provided databases (not bundled)
- Clear documentation on updating databases
- Version tracking for databases
- Automated scripts to download/update databases
- Tool agnostic to database version

### Risk 4: Limited Adoption

**Risk**: Tool developed but nobody uses it
**Probability**: Low (solves real problem)
**Impact**: High

**Mitigation**:
- Early community engagement (pre-print, conferences)
- Solve real problems (validate QC pipelines, benchmark tools)
- Excellent documentation (quick start, tutorials, examples)
- Publication in high-impact journal
- Active maintenance and support
- Integration with popular tools
- Teaching/workshop demonstrations

### Risk 5: Scope Creep

**Risk**: Project never finishes, keeps adding features
**Probability**: Medium
**Impact**: High

**Mitigation**:
- Clear version 1.0 scope (defined in Phase 1-4)
- Parking lot for future features (v2.0, v3.0)
- Release early, iterate based on user feedback
- Monthly progress reviews
- Milestone-based development
- "Done is better than perfect" mentality for v1.0

### Risk 6: Lack of Ground Truth Validation

**Risk**: Users don't trust the simulated data
**Probability**: Low
**Impact**: High

**Mitigation**:
- Publish validation study (compare to physical synthetic communities)
- Open-source all code (transparency)
- Comprehensive testing
- User testimonials and case studies
- Benchmarking against known tools
- Reproducibility (same parameters → same output)

---

## Success Metrics

### Phase 1 Success Metrics

- [x] Repository created and structured
- [ ] Can sample 50 viral genomes from RefSeq
- [ ] Can add contamination at specified percentages
- [ ] Can generate 1M paired-end reads
- [ ] Ground truth metadata file produced
- [ ] 1 basic scenario works end-to-end
- [ ] Code passes basic tests

### Phase 2 Success Metrics

- [ ] VLP vs bulk comparison produces different results
- [ ] PCR bias model shows length-dependent abundances
- [ ] 10+ pre-built scenarios available
- [ ] NovaSeq artifacts (polyG, optical dupes) present in output
- [ ] Body-site specific profiles validate against literature

### Phase 3 Success Metrics

- [ ] Simulated data statistically similar to real data (p > 0.05)
- [ ] Integration tests pass with lab-virome-QC
- [ ] 5+ external users successfully use tool
- [ ] Documentation complete with tutorials
- [ ] Performance acceptable (<1 hour for 10M reads)

### Phase 4 Success Metrics

- [ ] Manuscript accepted for publication
- [ ] Tool cited by at least 3 external studies (within 1 year)
- [ ] 50+ GitHub stars
- [ ] Active community engagement (issues, PRs)
- [ ] Presentation at major conference

---

## Long-Term Vision (5 Years)

### Year 1
- ViroForge v1.0 released
- Published in high-impact journal
- Adopted by 10+ research groups
- Used in 3+ benchmarking studies

### Year 2
- ViroForge v2.0 with long-read support (PacBio, ONT)
- "CAVI Challenge" - first community benchmarking competition
- Integrated into teaching courses
- 100+ citations

### Year 3
- Standard reference datasets established (like HMP)
- Adopted by major virome analysis pipelines
- Expanded to other environments (marine, soil, wastewater)
- International collaboration on viral standards

### Year 4
- ViroForge v3.0 with advanced features
- Machine learning integration (train viral classifiers)
- Real-time simulation for experimental design
- 500+ citations

### Year 5
- Community-maintained resource
- Multiple derived tools and extensions
- Gold standard for virome simulation
- Lasting impact on virome research field

---

## Conclusion

ViroForge represents a unique opportunity to:

1. **Fill a genuine gap** in virome research tools
2. **Enable rigorous validation** of analysis pipelines
3. **Advance the field** through standardized benchmarking
4. **Train the lab** in collaborative software development
5. **Publish high-impact** methods paper
6. **Serve the community** with open-source tool

The project is:
- ✅ **Technically feasible** (leverages existing tools, well-defined algorithms)
- ✅ **Scientifically valuable** (addresses unmet need identified in literature)
- ✅ **Publishable** (novel contribution, high-impact potential)
- ✅ **Manageable** (clear phases, realistic timeline)
- ✅ **Impactful** (could become standard for virome simulation)

**Recommendation**: Proceed with development following the phased roadmap outlined above.

---

## References

### Key Papers Informing Design

1. **ViromeQC**: Zolfo et al. (2019). "Detecting contamination in viromes using ViromeQC." Nature Biotechnology 37:1408-1412.

2. **CAMISIM**: Fritz et al. (2019). "CAMISIM: simulating metagenomes and microbial communities." Microbiome 7:17.

3. **InSilicoSeq**: Gourlé et al. (2019). "Simulating Illumina metagenomic data with InSilicoSeq." Bioinformatics 35(3):521-522.

4. **Virome Benchmarking**: Multiple 2023-2024 papers on benchmarking virome analysis tools using synthetic communities.

5. **Contamination**: Salter et al. (2014). "Reagent and laboratory contamination can critically impact sequence-based microbiome analyses." BMC Biology 12:87.

6. **VLP Methods**: Conceição-Neto et al. (2015). "Modular approach to customise sample preparation procedures for viral metagenomics." BMC Genomics 16:617.

### Additional Resources

- RefSeq Viral: https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/
- IMG/VR: https://img.jgi.doe.gov/vr/
- Gut Virome Database: https://www.gutvirome.org/
- SILVA rRNA: https://www.arb-silva.de/
- InSilicoSeq GitHub: https://github.com/HadrienG/InSilicoSeq
- CAMISIM GitHub: https://github.com/CAMI-challenge/CAMISIM

---

**Document Version**: 1.0
**Last Updated**: 2025-01-30
**Next Review**: After Phase 1 completion

**Contact**: Scott Handley Lab
**Repository**: https://github.com/shandley/viroforge
