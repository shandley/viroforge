# Phase 10: Long-Read Sequencing Support - Research Summary

**Date**: 2025-11-10
**Status**: Research Complete
**Decision**: Use PBSIM3 for both PacBio HiFi and Nanopore simulation

---

## Executive Summary

After evaluating three major long-read simulators (PBSIM3, PBSIM2, NanoSim), **PBSIM3** is the clear choice for ViroForge. It supports both PacBio HiFi and Oxford Nanopore platforms in a single tool with state-of-the-art error models.

**Key Decision Points**:
- ✅ Single tool for both platforms (reduces complexity)
- ✅ Most recent and actively maintained (2022, NAR Genomics Bioinformatics)
- ✅ HMM-based error models (superior to quality score models)
- ✅ HiFi simulation support (>99.9% accuracy)
- ✅ Conda installable (easy user setup)
- ✅ FASTQ output (compatible with existing workflow)

---

## Simulator Comparison

### 1. PBSIM3 (2022) ⭐ RECOMMENDED

**Publication**: Ono et al. (2022) NAR Genomics and Bioinformatics
**GitHub**: https://github.com/yukiteruono/pbsim3

#### Platforms Supported
- PacBio RS II CLR
- PacBio Sequel CLR
- **PacBio Sequel HiFi** ✨
- Oxford Nanopore Technologies (ONT) reads

#### Key Features
- **Error Model**: Hidden Markov Model (HMM) using FIC-HMM methodology
  - Sequence-independent but context-aware
  - Homopolymer deletion bias parameter for ONT (`--hp-del-bias 6`)
- **Multi-pass Sequencing**: Simulates CLR with multiple passes, then uses PacBio's `ccs` software for consensus
- **HiFi Accuracy**: Real HiFi = 0.22-0.25% error, Simulated = comparable with 10 passes
- **Transcriptome Support**: Can simulate RNA-seq (compatible with ViroForge RNA workflow)
- **Read Length**: Gamma distribution for WGS/TS

#### Installation
```bash
# Via conda (recommended)
conda install bioconda::pbsim3

# From source
git clone https://github.com/yukiteruono/pbsim3.git
cd pbsim3
./configure
make
```

#### Dependencies
- SAMtools (for BAM output)
- gzip (for compression)
- PacBio ccs software (for HiFi consensus calling)

#### Command-Line Usage
```bash
# PacBio HiFi simulation (two-step process)
pbsim --strategy wgs --method errhmm --errhmm MODEL \
      --depth 10 --pass-num 10 --genome reference.fasta

# Oxford Nanopore simulation
pbsim --strategy wgs --method errhmm --errhmm ONT-MODEL \
      --hp-del-bias 6 --depth 10 --genome reference.fasta
```

#### Input/Output Formats
- **Input**: FASTA (genome/transcripts)
- **Output**:
  - `.fq.gz`: Simulated reads (FASTQ, gzip-compressed)
  - `.maf.gz`: Read-to-reference alignments (MAF)
  - `.bam`: Multi-pass sequencing output (for HiFi)
  - `.ref`: Reference sequence copy

#### Pros
- ✅ Unified solution for both PacBio and ONT
- ✅ Most modern error models (HMM-based)
- ✅ Active maintenance and community support
- ✅ Realistic HiFi simulation via actual ccs software
- ✅ Transcriptome sequencing support

#### Cons
- ⚠️ No Python bindings (command-line only, requires subprocess)
- ⚠️ HiFi requires two-step process (PBSIM3 → ccs)
- ⚠️ Requires external dependencies (SAMtools, ccs)

---

### 2. PBSIM2 (2020) ❌ NOT RECOMMENDED

**Publication**: Ono et al. (2020) Bioinformatics
**Status**: Superseded by PBSIM3

#### Platforms Supported
- PacBio RS II CLR (P4C2, P5C3, P6C4 chemistry)
- Oxford Nanopore Technologies (early models)

#### Key Limitations
- ❌ **No HiFi/CCS support** (critical missing feature)
- ❌ Quality score model doesn't distinguish error types well
- ❌ Only WGS (no transcriptome simulation)
- ❌ Outdated PacBio chemistry models

#### Verdict
Obsolete. PBSIM3 is a strict superset with better models.

---

### 3. NanoSim Family (2017-2023) ⚠️ SPECIALIZED

**Publication**: Yang et al. (2017) GigaScience
**GitHub**: https://github.com/bcgsc/NanoSim

#### Variants
- **NanoSim** (2017): DNA sequencing simulator
- **Trans-NanoSim** (2020): RNA-seq specific, mimics ONT transcriptome features
- **Meta-NanoSim** (2023): Metagenomic sequencing

#### Key Features
- Statistical characterization-based (trains on real data)
- Captures ONT-specific artifacts (homopolymer errors)
- Read length distributions from empirical data

#### Recent Developments (2024-2025)
- **Squigulator** (2024): Signal-level simulator for R10.4.1 chemistry (most recent pore)
- **AsaruSim** (2025): Single-cell RNA-seq applications

#### Pros
- ✅ Specialized for Nanopore (may be more accurate for ONT)
- ✅ Statistical models trained on real data
- ✅ Trans-NanoSim for RNA-seq applications

#### Cons
- ❌ **ONT only** (no PacBio support)
- ❌ Requires training data for accuracy
- ❌ Multiple tools needed (NanoSim + Trans-NanoSim)

#### Verdict
Excellent for ONT-only workflows, but PBSIM3 provides adequate ONT simulation with unified PacBio support. Consider NanoSim for future specialized Nanopore applications if needed.

---

## Recommendation: PBSIM3

### Rationale

1. **Platform Coverage**: Single tool supports both PacBio HiFi and Nanopore (primary targets)
2. **HiFi Quality**: Uses actual PacBio ccs software for consensus calling (maximum realism)
3. **Error Models**: HMM-based approach superior to quality score models
4. **Maintenance**: Recent publication (2022), actively maintained
5. **Integration Simplicity**: Command-line tool with FASTQ output (standard format)
6. **User Accessibility**: Conda installable (low barrier to entry)

### Integration Strategy

#### Architecture
```python
class LongReadSimulator:
    """Wrapper for PBSIM3 long-read simulation."""

    def simulate_pacbio_hifi(self, genome_fasta, depth, passes=10):
        """
        Two-step process:
        1. Generate CLR with PBSIM3 (--pass-num 10)
        2. Generate HiFi consensus with ccs software
        """

    def simulate_nanopore(self, genome_fasta, depth):
        """
        Single-step ONT simulation with PBSIM3
        Uses --hp-del-bias for homopolymer artifacts
        """
```

#### Workflow Integration
```
ViroForge Workflow (Existing)
    ↓
genome_fasta = extract_genomes_from_collection()
    ↓
    ├─→ [SHORT-READ] InSilicoSeq → FASTQ
    │
    └─→ [LONG-READ] PBSIM3 → FASTQ
            ↓
        [IF PACBIO HIFI]
            ↓
        PacBio ccs → HiFi FASTQ
```

#### VLP Modeling for Long Reads
**Key Difference**: Long reads are less affected by size bias than short reads
- **Short reads**: 150-300bp fragments → strong size selection during sequencing
- **Long reads**: 10-30kb reads → size bias primarily from VLP enrichment step

**Approach**:
- Apply VLP size filtration to genome selection (same as current)
- Reduce/eliminate read-length-dependent bias in sequencing step
- Model fragmentation differently (longer fragments, different shearing)

---

## Technical Implementation Plan

### Phase 10.1: PacBio HiFi Support (Week 1-2)

**Tasks**:
1. Add PBSIM3 to dependencies (conda environment)
2. Create `LongReadSimulator` class
3. Implement `simulate_pacbio_hifi()` method
   - Call PBSIM3 with `--pass-num` parameter
   - Call PacBio ccs for consensus generation
   - Handle intermediate files cleanup
4. Add `--platform pacbio-hifi` flag to `generate_fastq_dataset.py`
5. Update metadata to track long-read parameters
6. Create basic integration tests

### Phase 10.2: Nanopore Support (Week 2-3)

**Tasks**:
1. Implement `simulate_nanopore()` method
   - Use PBSIM3 ONT error models
   - Set `--hp-del-bias` for homopolymer errors
   - Model quality-length relationships
2. Add `--platform nanopore` flag
3. Update VLP modeling for long reads (reduced size bias)
4. Create Nanopore-specific tests

### Phase 10.3: Documentation & Testing (Week 3-4)

**Tasks**:
1. Create long-read tutorial: "Generating PacBio HiFi Datasets"
2. Create long-read tutorial: "Generating Nanopore Datasets"
3. Add assembly benchmarking examples
4. Comprehensive integration tests
5. Update README with long-read capabilities

---

## Dependencies to Add

### Required
```yaml
# environment.yml or requirements
dependencies:
  - pbsim3  # bioconda
  - samtools  # for BAM handling
  - pbccs  # PacBio ccs software (for HiFi)
```

### Optional (Future)
- NanoSim (if specialized ONT accuracy needed)
- Squigulator (for signal-level simulation, R10.4.1 chemistry)

---

## Future Enhancements

### Near-term (v1.0.0)
- Basic PacBio HiFi and Nanopore support
- Standard error models from PBSIM3
- Integration with existing VLP workflow

### Long-term (v1.1.0+)
- **Custom Error Models**: Train PBSIM3 on lab-specific data
- **Signal-Level Simulation**: Integrate Squigulator for basecalling benchmarks
- **Adaptive Sampling**: Model ONT's selective sequencing feature
- **Ultra-Long Reads**: Specialized handling for >100kb reads
- **Hybrid Assembly**: Generate matched short+long read datasets

---

## References

1. **PBSIM3**: Ono Y, Asai K, Hamada M. PBSIM3: a simulator for all types of PacBio and ONT long reads. *NAR Genomics and Bioinformatics* 2022;4(4):lqac092. https://doi.org/10.1093/nargab/lqac092

2. **PBSIM2**: Ono Y, Hamada M, Asai K. PBSIM2: a simulator for long-read sequencers with a novel generative model of quality scores. *Bioinformatics* 2021;37(5):589-595.

3. **NanoSim**: Yang C, Chu J, Warren RL, Birol I. NanoSim: nanopore sequence read simulator based on statistical characterization. *GigaScience* 2017;6(4):gix010.

4. **Trans-NanoSim**: Yang C, et al. Trans-NanoSim characterizes and simulates nanopore RNA-sequencing data. *GigaScience* 2020;9(6):giaa061.

5. **Squigulator**: Gamaarachchi H, et al. End-to-end simulation of nanopore sequencing signals with feed-forward transformers. *Bioinformatics* 2024;41(1):btae744.

---

## Conclusion

**Decision**: Implement PBSIM3 as the primary long-read simulator for ViroForge Phase 10.

**Timeline**: 3-4 weeks
- Week 1-2: PacBio HiFi implementation
- Week 2-3: Nanopore implementation
- Week 3-4: Testing, documentation, tutorials

**Expected Impact**:
- Enables complete genome assembly benchmarking
- Future-proofs ViroForge for long-read era
- Expands user base to assembly developers
- High impact with medium-high effort (realistic for 4-week timeline)

**Next Steps**: Design integration architecture and begin implementation.
