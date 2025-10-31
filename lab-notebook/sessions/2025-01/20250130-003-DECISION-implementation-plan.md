---
entry_id: 20250130-003-DECISION-implementation-plan
date: 2025-01-30
type: DECISION
status: complete
phase: 2 (planning)

author: Scott Handley + Claude

references:
  literature:
    - Shkoporov & Hill (2019) - VLP protocols
    - Zolfo et al. (2019) - ViromeQC enrichment metrics
    - Multiple papers on amplification methods
  prior_sessions:
    - 20250130-002
  related_topics:
    - vlp-enrichment
    - amplification-bias
    - platform-artifacts

tags:
  - phase2-planning
  - implementation-roadmap
  - architecture-design
  - 12-week-timeline

key_decisions:
  - Modular, composable architecture (flexible frameworks)
  - Lab-agnostic design (supports any protocol)
  - 12-week timeline for Phase 2
  - Pre-defined workflow templates for convenience

commits:
  - 3e96646  # Add community-focused Phase 2 implementation plan

raw_data: none
---

# Phase 2 Implementation Plan: Community-Focused, Lab-Agnostic Approach

**Date**: January 30, 2025
**Phase**: 2 (Planning)
**Status**: Complete

## Goals

Based on strategic review (Session 002), create detailed Phase 2 implementation plan that:
1. Implements core virome-specific features
2. Uses modular, composable architecture
3. Is lab-agnostic (not specific to one protocol)
4. Maintains focused timeline (~12 weeks)
5. Supports publication preparation

## Background

Strategic review revealed ViroForge only models 42% of real virome workflows. User wants to add virome-specific features BUT maintain community focus (not lab-specific).

**Key Challenge**: How to be both focused (fast implementation) AND flexible (lab-agnostic)?

**Solution**: Modular frameworks with pre-defined templates

## Architecture Decision: Modular, Composable Design

### Core Principle

**NOT**: Hard-coded protocols
```python
# BAD - Lab-specific, inflexible
def handley_lab_vlp_protocol():
    filter_0_2_um()
    dnase_treatment()
    # Only works for one lab
```

**INSTEAD**: Flexible frameworks
```python
# GOOD - Lab-agnostic, composable
workflow = ViromePipeline(
    enrichment=VLPEnrichment(
        filtration_cutoff_um=0.2,  # User configurable
        nuclease_efficiency=0.95,
        method='tangential_flow'
    ),
    amplification=RdABAmplification(
        cycles=40,
        length_bias='strong',
        gc_bias='moderate'
    ),
    library_prep=MechanicalShearing(
        target_size=300
    ),
    platform=NovaSeqPlatform(
        polyg_rate=0.2,
        optical_dup_rate=0.05
    )
)
composition = workflow.apply(base_composition)
```

### Benefits

**For Community**:
- Works for ANY lab's protocols
- Users configure parameters for their methods
- Extensible (add new methods easily)
- Supports method comparisons

**For Development**:
- Modular (test each component independently)
- Composable (mix and match)
- Maintainable (clean interfaces)
- Publication-ready (well-documented)

### Pre-Defined Templates

**Convenience wrappers** for common workflows:
```python
# Quick start - common protocol
from viroforge.workflows import gut_virome_vlp_rdab_novaseq

composition = gut_virome_vlp_rdab_novaseq(
    n_species=50,
    n_reads=10_000_000
)
```

Behind the scenes: Uses flexible framework with typical parameters

## Phase 2 Timeline: 12 Weeks

### Weeks 1-3: VLP Enrichment Framework ← **START HERE**

**Goal**: Flexible VLP enrichment supporting any protocol

**Components**:
1. `VLPEnrichment` base class
2. `FiltrationModel` (size-based viral enrichment)
3. `NucleaseModel` (host DNA/RNA removal)
4. Pre-defined protocols (0.2 μm standard, 0.45 μm alternative, etc.)

**Literature Support**:
- Shkoporov 2019: 0.1-0.45 μm filtration range
- Expected viral recovery: 60-95%
- Host removal: 95-99.9%

**Deliverables**:
- `viroforge/enrichment.py` (~500 lines)
- Unit tests (>90% coverage)
- Tutorial notebook
- Documentation

**Validation**:
- Viral recovery in expected range
- Host DNA removal matches literature
- VLP vs bulk comparison capability

### Weeks 4-6: Amplification Bias Framework

**Goal**: Support multiple amplification methods

**Components**:
1. `AmplificationBias` base class
2. `RdABAmplification` (random displacement, 40-cycle PCR)
3. `MDAAmplification` (multiple displacement, φ29)
4. `LinkerAmplification` (linker-based methods)
5. `NoAmplification` (baseline)

**Effects to Model**:
- **Length bias**: Short genomes amplify better (RdAB)
- **GC bias**: AT-rich amplify better (MDA extreme)
- **Coverage uniformity**: Varies by method

**Literature Support**:
- RdAB: Moderate length bias, 40 cycles typical
- MDA: Extreme GC bias (>60% GC underrepresented)
- Coverage: RdAB more uniform than MDA

**Deliverables**:
- `viroforge/amplification.py` (~600 lines)
- Bias models with literature-validated parameters
- Tests + documentation

### Weeks 7-8: Platform Artifacts & Library Prep

**Goal**: Realistic sequencing artifacts and library prep methods

**Platform Artifacts**:
1. `NovaSeqPlatform` (2-channel chemistry)
   - PolyG tails (G homopolymers in low-quality regions)
   - Optical duplicates (tile position-based)
2. `MiSeqPlatform` (4-channel, no polyG)
3. `HiSeqPlatform` (2-channel, older chemistry)

**Library Prep Methods**:
1. `MechanicalShearing` (Covaris, uniform)
2. `Tagmentation` (Nextera, GC bias)
3. `EnzymaticFragmentation` (dsDNA fragmentase)

**Literature Support**:
- NovaSeq polyG: 0.1-0.3% of bases (Chen 2021)
- Optical duplicates: 0.05-0.15% typical
- Tagmentation GC bias: Well documented

**Deliverables**:
- `viroforge/platforms.py` (~400 lines)
- `viroforge/library_prep.py` (~300 lines)
- Tests + documentation

### Weeks 9-10: Integration & Testing

**Goal**: Complete workflows work end-to-end

**Activities**:
1. Integration testing (all components together)
2. Pre-defined workflow templates
3. VLP vs bulk comparison validation
4. Performance optimization
5. User tutorials

**Validation**:
- Complete workflows generate realistic data
- Parameters within literature ranges
- Ground truth accurate throughout pipeline
- Comparison studies work correctly

### Weeks 11-12: Documentation & Publication Prep

**Goal**: Publication-ready documentation

**Deliverables**:
1. **User Guide**: Comprehensive tutorials
2. **API Reference**: Complete documentation
3. **Methods Section**: Publication-quality writeup
4. **Validation Report**: Literature comparison
5. **Example Workflows**: Common use cases

**Publication Materials**:
- Methods text (ready for manuscript)
- Validation figures
- Parameter tables with literature support
- Use case examples

## Key Architecture Principles

### 1. Modular
- Each component independently testable
- Clear interfaces between components
- No hidden dependencies

### 2. Composable
- Users mix and match components
- Pre-defined templates for convenience
- Custom workflows easy to create

### 3. Lab-Agnostic
- Configurable for any protocol
- Literature-based parameter ranges
- Not tied to specific lab methods

### 4. Literature-Validated
- Every parameter has citation
- Expected ranges from published studies
- Validation against real data

### 5. Ground Truth Throughout
- Metadata tracked through all transformations
- Complete tracking of:
  - Original composition
  - VLP enrichment effects
  - Amplification biases
  - Coverage patterns
  - Sequencing artifacts

## Implementation Order Rationale

**Why VLP first (Weeks 1-3)?**
- Most critical distinguishing feature of viromes
- Simplest conceptually (filtration + nuclease)
- Provides foundation for testing other features

**Why amplification second (Weeks 4-6)?**
- Builds on VLP-enriched composition
- More complex (multiple bias types)
- Critical for realistic diversity patterns

**Why platforms/library prep third (Weeks 7-8)?**
- Works on amplified composition
- Less critical for core biology
- Can be added modularly

**Why integration/docs last (Weeks 9-12)?**
- Requires all components complete
- Allows time for testing and refinement
- Publication prep needs complete system

## Success Criteria

**Technical**:
- ✅ All unit tests passing (>90% coverage)
- ✅ Integration tests working
- ✅ Performance acceptable (<2 hours for 10M reads)
- ✅ Ground truth accurate throughout

**Scientific**:
- ✅ Parameters within literature ranges
- ✅ VLP vs bulk comparison realistic
- ✅ Amplification biases match published patterns
- ✅ Platform artifacts realistic

**Community**:
- ✅ Lab-agnostic (works for any protocol)
- ✅ Well-documented (tutorials + API)
- ✅ Publication-ready (Nature Biotech level)
- ✅ Extensible (easy to add methods)

## Risks and Mitigation

### Risk 1: Scope Creep
**Mitigation**: Stick to 4 core features, defer others to Phase 3

### Risk 2: Literature Gaps
**Mitigation**: Document assumptions, make configurable

### Risk 3: Timeline Slippage
**Mitigation**: Weekly progress reviews, prioritize core functionality

### Risk 4: Biological Complexity
**Mitigation**: Simplify models while maintaining realism, validate key parameters

## Next Session

**Immediate Task**: Begin VLP enrichment framework (Week 1)

**Starting Point**:
1. Create `viroforge/enrichment.py`
2. Implement `VLPEnrichment` base class
3. Add `FiltrationModel` (size-based selection)
4. Add `NucleaseModel` (host DNA removal)
5. Write unit tests
6. Create tutorial

## Related Documents

- `docs/IMPLEMENTATION_PLAN.md` - Complete 997-line detailed plan (basis for this entry)
- `lab-notebook/sessions/2025-01/20250130-002-STRATEGIC-scope-review.md` - Strategic context
- `.claude.md` - Updated with Phase 2 start info

---

**Status**: Complete - Plan approved by user
**Impact**: Clear 12-week roadmap for Phase 2. Modular, composable architecture will produce lab-agnostic, community-focused tool while maintaining focused timeline. Publication-quality documentation planned throughout.

**User Approval**: "Yes I am Ok with this approach"
