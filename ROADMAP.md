# ViroForge Development Roadmap

**Version**: 0.4.0 → 1.0.0
**Timeline**: 3-6 months
**Goal**: Comprehensive virome simulation covering DNA/RNA, host/environmental, healthy/disease states

---

## Current Status (v0.6.0)

**Phase 7 Complete** - Critical Collections

✅ 14,423 RefSeq viral genomes with ICTV taxonomy
✅ 12 curated collections (9 host-associated, 3 environmental)
✅ 5 VLP enrichment protocols with size-based filtration
✅ Type-specific contamination reduction
✅ Amplification bias modeling (RdAB, MDA, Linker amplification)
✅ Wastewater virome for epidemiological surveillance
✅ Disease state collections (IBD, HIV+, CF)
✅ Progressive dysbiosis modeling (Healthy → IBD → HIV+)
✅ Platform-specific error models (NovaSeq, MiSeq, HiSeq)
✅ Complete ground truth tracking
✅ Comprehensive documentation with corrected citations

---

## Development Path B: Comprehensive Enhancement

### **PHASE 6: Enhanced Library Preparation**

**Timeline**: 1-2 weeks
**Status**: ✅ Complete

#### Objectives
- Integrate existing amplification bias code into user workflow
- Enable comparison of library prep methods
- Model realistic PCR artifacts

#### Tasks
- [x] Create development roadmap
- [ ] Add `--amplification` flag to `generate_fastq_dataset.py`
- [ ] Integrate RdAB, MDA, Linker amplification into workflow
- [ ] Add amplification to batch generation presets
- [ ] Create amplification comparison preset
- [ ] Update metadata to include amplification parameters
- [ ] Add integration tests for amplification
- [ ] Update documentation (QUICKSTART, PHASE4, tutorials)
- [ ] Create amplification comparison tutorial

#### Deliverables
- `--amplification {none,rdab,mda,linker}` flag in all scripts
- Amplification statistics in metadata.json
- Tutorial: "Comparing Library Preparation Methods"
- Tests: Integration tests for all 4 amplification methods

#### Impact
- **HIGH**: Affects all virome studies
- **Effort**: LOW (code exists, needs integration)
- **Users**: Pipeline developers, method comparison studies

---

### **PHASE 7: Critical Collections**

**Timeline**: 3-4 weeks
**Status**: ✅ Complete

#### 7.1 Wastewater Virome (Week 1)

**Rationale**: Highly timely for epidemiological surveillance

**Collection Design**:
- **Name**: Wastewater Virome - Urban Treatment Plant
- **Size**: 300-500 genomes
- **Composition**:
  - Human enteric viruses (40%): Norovirus, rotavirus, adenovirus, astrovirus, sapovirus
  - Bacteriophages (35%): Gut bacterial phages
  - Environmental viruses (15%): Plant viruses, insect viruses
  - Emerging pathogens (10%): SARS-CoV-2, mpox, polio, etc.
- **Literature**: Crits-Christoph et al. 2021, Crank et al. 2022

**Tasks**:
- [x] Literature review of wastewater virome composition
- [x] Curate genome list from RefSeq
- [x] Assign literature-validated abundances
- [x] Add to database as collection ID 17
- [x] Validate collection
- [ ] Create wastewater-specific contamination profile (future enhancement)
- [x] Document in Collection Implementation Guide

#### 7.2 Disease State Collections (Weeks 2-4)

**Collection 18: IBD Gut Virome** (Week 2)
- **Comparison**: vs healthy gut (collection 9)
- **Key differences**: Lower diversity, altered Caudovirales, increased temperate phages
- **Size**: 80-100 genomes
- **Literature**: Norman et al. 2015, Zuo et al. 2019

**Collection 19: HIV+ Gut Virome** (Week 3)
- **Comparison**: vs healthy gut (collection 9)
- **Key differences**: Dramatically reduced diversity, altered composition, increased eukaryotic viruses
- **Size**: 40-60 genomes
- **Literature**: Handley et al. 2012, Monaco et al. 2016

**Collection 20: Cystic Fibrosis Respiratory Virome** (Week 4)
- **Comparison**: vs healthy respiratory (collection 12)
- **Key differences**: Chronic viral infections, altered bacteriophage composition
- **Size**: 25-35 genomes
- **Literature**: Lim et al. 2014, Willner et al. 2009

#### Deliverables
- 4 new collections (1 wastewater, 3 disease states)
- Disease vs healthy comparison tutorials
- Dysbiosis modeling examples
- Updated Collection Implementation Guide

#### Impact
- **VERY HIGH**: Clinical validation, public health surveillance
- **Effort**: MEDIUM (1 week per collection)
- **Users**: Clinical researchers, epidemiologists, public health

---

### **PHASE 8: RNA Virome Workflow** (Current Phase)

**Timeline**: 3-4 weeks
**Status**: 🚧 In Progress

#### Objectives
- Support RNA virus sequencing workflows
- Model reverse transcription step
- Implement rRNA depletion for RNA viromes

#### 8.1 RNA Virus Collections (Weeks 1-2)

**Collection 21: Human Respiratory RNA Virome**
- **Viruses**: Influenza, RSV, HMPV, coronaviruses, rhinoviruses, enteroviruses
- **Size**: 50-70 genomes
- **Type**: ssRNA (negative and positive sense)

**Collection 22: Arbovirus Environmental (Mosquito)**
- **Viruses**: Flaviviruses, alphaviruses, bunyaviruses
- **Size**: 30-50 genomes
- **Type**: Mixed ssRNA, dsRNA

**Collection 23: Fecal RNA Virome**
- **Viruses**: Norovirus, rotavirus, astrovirus, picornaviruses
- **Size**: 40-60 genomes
- **Type**: Mixed ssRNA, dsRNA

#### 8.2 RNA Workflow Implementation (Weeks 3-4)

**Tasks**:
- [ ] Create RNAViromeWorkflow class
- [ ] Implement reverse transcription modeling
  - RT efficiency by virus type
  - Random vs oligo-dT vs specific primers
  - RT artifacts (template switching, truncation)
- [ ] Implement rRNA depletion
  - Different from DNA rRNA removal
  - Ribo-Zero or Ribominus modeling
  - Variable efficiency
- [ ] Add RNA degradation modeling
  - More fragile than DNA
  - RNase contamination
- [ ] Create `--molecule-type {dna,rna}` flag
- [ ] Update contamination profiles for RNA
- [ ] Add RNA-specific tests
- [ ] Create RNA virome tutorial

#### Deliverables
- RNA virome workflow
- 3 RNA virus collections
- RNA-specific documentation
- Tutorial: "Generating RNA Virome Datasets"

#### Impact
- **HIGH**: Opens entirely new user base
- **Effort**: MEDIUM-HIGH (new workflow components)
- **Users**: Respiratory researchers, clinical diagnostics, environmental surveillance

---

### **PHASE 9: Additional Host Niches**

**Timeline**: 2-3 weeks
**Status**: 📋 Planned

#### Collections to Add

**Collection 24: Vaginal Virome** (Week 1)
- **Size**: 20-30 genomes
- **Composition**: Papillomaviruses, herpesviruses, bacteriophages
- **Literature**: Wylie et al. 2014, Dols et al. 2016
- **Impact**: Women's health, pregnancy outcomes

**Collection 25: Blood/Plasma Virome** (Week 1)
- **Size**: 15-25 genomes
- **Composition**: Anelloviruses, herpesviruses, parvoviruses
- **Literature**: Moustafa et al. 2017, Young et al. 2016
- **Impact**: Viremia, transplant monitoring, diagnostics

**Collection 26: Ocular Surface Virome** (Week 2)
- **Size**: 10-20 genomes
- **Composition**: Adenoviruses, herpesviruses, bacteriophages
- **Literature**: Doan et al. 2016
- **Impact**: Ophthalmology, infectious keratitis

**Collection 27: Lower Respiratory (Lung) Virome** (Week 2)
- **Size**: 25-40 genomes
- **Composition**: Different from nasopharynx, includes lung-resident viruses
- **Literature**: Kitsios et al. 2018
- **Impact**: Pneumonia, COPD, lung transplant

**Collection 28: Urinary Virome** (Week 3)
- **Size**: 15-25 genomes
- **Composition**: BK virus, JC virus, adenoviruses
- **Literature**: Santiago-Rodriguez et al. 2015
- **Impact**: Transplant monitoring, UTI

#### Deliverables
- 5 new host-associated collections
- Coverage of major human body sites
- Updated collection comparison tools

#### Impact
- **MEDIUM**: Expands applicable research niches
- **Effort**: LOW-MEDIUM (similar to existing collections)
- **Users**: Clinical researchers across specialties

---

### **PHASE 10: Long-Read Sequencing Support**

**Timeline**: 3-4 weeks
**Status**: 📋 Planned (Future)

#### Objectives
- Support PacBio HiFi and Nanopore platforms
- Enable complete genome assembly benchmarking
- Model long-read specific artifacts

#### Tasks
- [ ] Research long-read simulators (pbsim3, NanoSim, PBSIM2)
- [ ] Integrate PacBio HiFi simulator
  - High accuracy (QV20+)
  - Different error profile than Illumina
  - Length distribution modeling
- [ ] Integrate Nanopore simulator
  - Homopolymer errors
  - Ultra-long reads
  - Quality-length relationships
- [ ] Add `--platform {novaseq,miseq,hiseq,pacbio-hifi,nanopore}` options
- [ ] Update VLP modeling for long reads (affects size bias differently)
- [ ] Create long-read specific tests
- [ ] Create long-read tutorial

#### Deliverables
- PacBio HiFi support
- Nanopore support
- Long-read assembly benchmarking capability
- Tutorial: "Long-Read Virome Benchmarking"

#### Impact
- **HIGH**: Future-proofing, better assemblies
- **Effort**: MEDIUM-HIGH (new external tools)
- **Users**: Assembly developers, complete genome studies

---

### **PHASE 11: Temporal Dynamics**

**Timeline**: 2-3 weeks
**Status**: 📋 Planned (Future)

#### Objectives
- Enable longitudinal study generation
- Model temporal variation in viral communities
- Support perturbation-recovery studies

#### Framework Design

**Temporal Models**:
1. **Infant Gut Development**: Birth → 2 years
2. **Antibiotic Perturbation**: Pre → during → recovery
3. **Seasonal Variation**: Environmental viromes over 12 months
4. **Infection Dynamics**: Acute infection → resolution
5. **Diet Change**: Baseline → intervention → new steady state

#### Tasks
- [ ] Create TemporalSeries class
- [ ] Implement turnover modeling
  - Virus appearance/disappearance
  - Colonization/extinction dynamics
- [ ] Implement abundance dynamics
  - Bloom-bust cycles
  - Gradual shifts
  - Stable vs transient viruses
- [ ] Add `--timeseries` flag with predefined trajectories
- [ ] Create visualization tools for temporal data
- [ ] Add temporal validation datasets
- [ ] Create temporal tutorial

#### Deliverables
- Temporal framework
- 5 predefined temporal trajectories
- Tutorial: "Generating Longitudinal Virome Studies"

#### Impact
- **MEDIUM-HIGH**: Enables longitudinal study design
- **Effort**: MEDIUM (complex modeling)
- **Users**: Microbiome researchers, intervention studies

---

### **PHASE 12: Additional Animal Models**

**Timeline**: 2-3 weeks
**Status**: 📋 Planned (Future)

#### Collections to Add

**Collection 29: Zebrafish Gut Virome**
- **Size**: 15-25 genomes
- **Literature**: Melancon et al. 2017
- **Impact**: Development biology, genetics

**Collection 30: Pig Gut Virome**
- **Size**: 60-80 genomes
- **Literature**: Shan et al. 2011
- **Impact**: Agricultural, translational research

**Collection 31: Chicken Gut Virome**
- **Size**: 40-60 genomes
- **Literature**: Day et al. 2010
- **Impact**: Poultry industry, food safety

**Collection 32: Non-Human Primate Gut Virome**
- **Size**: 80-100 genomes
- **Literature**: Handley et al. 2018
- **Impact**: Translational research, evolution

#### Deliverables
- 4 animal model collections
- Animal-specific contamination profiles
- Comparative viromics examples

#### Impact
- **MEDIUM**: Expands to animal research community
- **Effort**: LOW-MEDIUM (similar to existing)
- **Users**: Veterinary, agricultural, comparative biology

---

### **PHASE 13: Environmental Diversity**

**Timeline**: 2 weeks
**Status**: 📋 Planned (Future)

#### Collections to Add

**Collection 33: Hot Spring Virome** (Extreme environment)
- **Size**: 100-150 genomes
- **Literature**: Schoenfeld et al. 2008

**Collection 34: Hypersaline Virome** (Extreme environment)
- **Size**: 80-120 genomes
- **Literature**: Santos et al. 2010

**Collection 35: Hospital Surface Virome** (Built environment)
- **Size**: 40-60 genomes
- **Literature**: Lax et al. 2014

**Collection 36: Plant Phyllosphere Virome** (Agriculture)
- **Size**: 50-80 genomes
- **Literature**: Zablocki et al. 2016

#### Deliverables
- 4 environmental collections
- Extreme environment modeling
- Built environment applications

#### Impact
- **LOW-MEDIUM**: Niche applications
- **Effort**: LOW (standard curation)
- **Users**: Environmental microbiologists, ecology

---

## Version Milestones

### **v0.5.0 - Enhanced Library Prep** (2 weeks from now)
- Amplification bias integration
- Amplification comparison tools
- Updated documentation

### **v0.6.0 - Critical Collections** (6 weeks from now)
- Wastewater virome
- 3 disease state collections
- Disease vs healthy comparison framework

### **v0.7.0 - RNA Virome Support** (10 weeks from now)
- RNA workflow implementation
- 3 RNA virus collections
- RNA-specific documentation

### **v0.8.0 - Expanded Host Niches** (13 weeks from now)
- 5 additional host-associated collections
- Comprehensive human body site coverage

### **v1.0.0 - Feature Complete** (16-20 weeks from now)
- Long-read support (PacBio HiFi minimum)
- Temporal dynamics framework
- 30+ curated collections
- Comprehensive benchmarking platform

---

## Publication Strategy

### **Manuscript Timeline**

**Phase 1: Initial Submission** (After v0.7.0, ~10 weeks)
- Title: "ViroForge: A Comprehensive Synthetic Virome Data Generator for Pipeline Benchmarking"
- Content: DNA viromes, VLP modeling, amplification bias, disease states, RNA support
- Target: Bioinformatics, BMC Bioinformatics, or Microbiome

**Phase 2: Major Update** (After v1.0.0, ~20 weeks)
- Update with long-read support, temporal dynamics
- Expanded collection catalog
- Comprehensive validation studies

### **Preprint Strategy**

- Post to bioRxiv after v0.6.0 (wastewater + disease collections)
- Update preprint with v1.0.0 features
- Submit to peer review after v0.7.0 (RNA support)

---

## Success Metrics

### **Technical Metrics**
- [ ] 30+ curated virome collections
- [ ] DNA + RNA workflow support
- [ ] 5+ sequencing platforms (Illumina + PacBio + Nanopore)
- [ ] 4+ library prep methods
- [ ] Temporal dynamics framework
- [ ] 100% test coverage on critical paths

### **Community Metrics**
- [ ] 100+ GitHub stars
- [ ] 10+ citations in peer-reviewed literature
- [ ] Used in published benchmarking studies
- [ ] 50+ downloads per month
- [ ] Active community contributions

### **Scientific Impact**
- [ ] Enables standardized virome pipeline benchmarking
- [ ] Cited in clinical virome studies
- [ ] Used for method comparison papers
- [ ] Adopted by pipeline developers (Hecatomb, others)

---

## Resource Requirements

### **Computational**
- Database storage: ~500 MB (current) → ~2 GB (v1.0.0)
- Test datasets: ~10 GB
- CI/CD: GitHub Actions (free tier sufficient)

### **Time Commitment**
- **Phase 6**: 1-2 weeks (amplification)
- **Phase 7**: 3-4 weeks (critical collections)
- **Phase 8**: 3-4 weeks (RNA workflow)
- **Phase 9-13**: 8-10 weeks (additional collections, long-read, temporal)
- **Total**: 16-20 weeks for v1.0.0

### **Expertise Needed**
- Virome biology (curating collections)
- Bioinformatics (workflow implementation)
- Literature review (validation)
- Technical writing (documentation)

---

## Risk Assessment

### **Technical Risks**

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Long-read simulator integration issues | Medium | Medium | Start with PacBio HiFi only, well-documented tools |
| RNA workflow complexity | Medium | High | Leverage existing RT-PCR knowledge, literature |
| Database size growth | Low | Low | SQLite handles 2GB easily, can migrate to PostgreSQL |
| Test suite maintenance | Medium | Medium | Incremental testing, dry-run modes for CI/CD |

### **Scientific Risks**

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Literature gaps for collections | Low | Medium | Start with well-studied conditions, expert consultation |
| Parameter validation challenges | Low | High | Conservative estimates, sensitivity analysis |
| User adoption | Medium | High | Excellent documentation, tutorials, examples |

---

## Community Engagement

### **Documentation**
- Comprehensive tutorials for each phase
- Video tutorials (YouTube)
- Workshop materials
- Example workflows

### **Outreach**
- Present at conferences (ISME, ASM, Microbiome meetings)
- Twitter/X announcements for new features
- Engage with virome research community
- Collaborate with pipeline developers

### **Support**
- GitHub Discussions for Q&A
- Responsive issue tracking
- Monthly release notes
- Changelog documentation

---

## Conclusion

ViroForge v1.0.0 will be the **most comprehensive virome simulation platform** available, covering:

- **DNA + RNA viromes**
- **Host-associated + environmental samples**
- **Healthy + disease states**
- **Short-read + long-read sequencing**
- **Multiple library prep methods**
- **Temporal dynamics**
- **30+ curated collections**
- **Complete ground truth for rigorous benchmarking**

**Timeline**: 16-20 weeks to v1.0.0
**Impact**: Essential tool for virome research community
**Publication**: Strong methods paper in top-tier bioinformatics journal

---

**Let's build the future of virome benchmarking! 🦠🔨**
