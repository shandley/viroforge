# Body Site Virome Collections - Comprehensive Curation Plan

**Date**: November 1, 2025
**Status**: Planning Phase
**Target**: 8+ body site collections with literature-validated compositions

---

## Overview

This document provides comprehensive curation plans for creating realistic body site-specific virome collections in ViroForge. Each collection is based on recent literature (2023-2025) and designed to represent healthy adult humans unless otherwise specified.

**Total Collections Planned**: 8
**Total Genomes Across Collections**: ~2,200
**Completion Target**: Week 6-8

---

## Collection 1: Human Gut Virome (Adult, Healthy, Western) ✅

**Status**: Plan Complete (see `GUT_VIROME_CURATION.md`)
**Target Size**: 500 genomes
**Collection ID**: `gut_virome_adult_healthy_western`

**Quick Summary**:
- Crassvirales (crAssphage): 30% - Dominant, 90-95% in some individuals
- Caudoviricetes: 50% - Siphoviridae, Myoviridae, Podoviridae
- Microviridae: 10% - Stable colonizers
- Other phages: 6%
- Eukaryotic viruses: 4%

**Key References**: Dutilh et al. 2014, Guerin et al. 2018, Shkoporov et al. 2019

---

## Collection 2: Human Oral Virome (Saliva, Healthy Adult)

**Target Size**: 200 genomes
**Collection ID**: `oral_virome_saliva_healthy`
**Planned Completion**: Week 6, Day 7

### Literature-Based Composition

**Abundance**:
- Viral particles: ~10^8 VLPs/mL saliva
- Phage diversity: Hundreds to thousands of genotypes
- Individual variation: High - personalized virome

**Dominant Groups** (based on Pride et al. 2012, Lloyd-Price et al. 2017, recent 2024 studies):

**1. Siphoviridae - 80 genomes (40%)**
- Rationale: Most abundant in healthy oral cavity
- Hosts: Streptococcus, Actinomyces, Veillonella
- Characteristics: Temperate phages, long non-contractile tails
- Distribution:
  - Streptococcus phages: 40 genomes
  - Actinomyces phages: 20 genomes
  - Veillonella phages: 15 genomes
  - Other oral bacteria: 5 genomes

**2. Myoviridae - 40 genomes (20%)**
- Rationale: Increased in periodontal disease (baseline present in health)
- Hosts: Diverse oral bacteria
- Characteristics: Lytic phages, contractile tails
- Distribution:
  - Porphyromonas phages: 15 genomes
  - Fusobacterium phages: 10 genomes
  - Other oral bacteria: 15 genomes

**3. Podoviridae - 30 genomes (15%)**
- Hosts: Streptococcus, Lactobacillus
- Short-tailed phages
- Less abundant than Siphoviridae

**4. Microviridae - 20 genomes (10%)**
- ssDNA phages
- Stable colonizers
- Small genomes (4-8 kb)

**5. Inoviridae - 15 genomes (7.5%)**
- Filamentous phages
- Chronic infections
- Diverse oral bacteria hosts

**6. Eukaryotic Viruses - 15 genomes (7.5%)**
- Herpesviruses: 5 genomes (EBV, HSV-1, HHV-6, HHV-7)
- Papillomaviruses: 3 genomes (oral HPV types)
- Anelloviruses: 4 genomes (commensal)
- Adenoviruses: 3 genomes

### Anatomical Site Distribution

**Saliva** (this collection): Mixed from all oral sites
**Dental Plaque** (future): 10x higher phage abundance, different composition
**Tongue Dorsum** (future): Highest diversity, jumbo phages

### Abundance Model

**Tier 1** (Dominant, 5-20%): 10 genomes
- Major Streptococcus phages
- Persistent Actinomyces phages

**Tier 2** (Common, 1-5%): 40 genomes
- Diverse Siphoviridae
- Common Veillonella phages

**Tier 3** (Moderate, 0.1-1%): 80 genomes
- Myoviridae, Podoviridae
- Diverse Streptococcus variants

**Tier 4** (Rare, 0.01-0.1%): 50 genomes
- Eukaryotic viruses
- Rare phage families

**Tier 5** (Very rare, <0.01%): 20 genomes
- Novel/unclassified
- Transient viruses

### Selection Criteria

**Primary filters**:
```sql
-- Oral-associated families
WHERE family IN ('Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae', 'Inoviridae')
  AND (genome_name LIKE '%Streptococcus%'
    OR genome_name LIKE '%Actinomyces%'
    OR genome_name LIKE '%Veillonella%'
    OR genome_name LIKE '%oral%'
    OR genome_name LIKE '%saliva%')
```

**Host preferences** (from ecological metadata):
- Streptococcus spp. (most abundant oral bacteria)
- Actinomyces spp.
- Veillonella spp.
- Porphyromonas spp.
- Fusobacterium spp.

### Validation Metrics

- [ ] Siphoviridae dominant (40-50%)
- [ ] High individual variation modeled
- [ ] Streptococcus phages most abundant (20-30%)
- [ ] Herpesvirus presence (5-10% prevalence)
- [ ] Personal/temporal stability (for same seed)

### Key References

1. Pride et al. (2012). "Evidence of a robust resident bacteriophage population revealed through analysis of the human salivary virome." ISME J.
2. Abeles et al. (2014). "Altered oral viral ecology in association with periodontal disease." mBio.
3. Lloyd-Price et al. (2017). "The healthy human microbiome." Genome Med.
4. Recent 2024 review: "Oral bacteriophages: manifestations in the ecology of oral diseases."

---

## Collection 3: Human Skin Virome (Sebaceous Sites, Healthy)

**Target Size**: 150 genomes
**Collection ID**: `skin_virome_sebaceous_healthy`
**Planned Completion**: Week 7, Day 1

### Literature-Based Composition

**Abundance**:
- Lower viral load than gut/oral
- Site-specific variation (sebaceous vs. moist vs. dry)
- Dominated by Cutibacterium (Propionibacterium) phages

**Dominant Groups** (based on Hannigan et al. 2015, Oh et al. 2016, recent 2024 studies):

**1. Cutibacterium (Propionibacterium) Phages - 70 genomes (47%)**
- Rationale: Cutibacterium acnes is dominant sebaceous skin bacteria
- Type: Caudovirales (Siphoviridae-like)
- Distribution:
  - P1-like phages: 30 genomes
  - P2-like phages: 20 genomes
  - P3-like phages: 10 genomes
  - Novel Cutibacterium phages: 10 genomes

**2. Staphylococcus Phages - 40 genomes (27%)**
- Rationale: S. epidermidis, S. aureus common on skin
- Distribution:
  - S. epidermidis phages: 25 genomes (most abundant)
  - S. aureus phages: 15 genomes

**3. Corynebacterium Phages - 15 genomes (10%)**
- Skin-associated Corynebacterium species
- Less abundant but persistent

**4. Other Caudoviricetes - 10 genomes (7%)**
- Diverse skin bacteria hosts
- Myoviridae, Podoviridae

**5. Microviridae - 8 genomes (5%)**
- ssDNA phages
- Skin colonizers

**6. Eukaryotic Viruses - 7 genomes (4%)**
- Papillomaviruses: 3 genomes (skin HPV types)
- Polyomaviruses: 2 genomes (MCPyV, HPyV6, HPyV7)
- Herpesviruses: 2 genomes (HHV-6, HHV-7)

### Anatomical Site Variation

**Sebaceous** (this collection): Forehead, back
- High Cutibacterium phages
- Moderate Staphylococcus phages

**Moist** (future): Armpit, groin
- More Staphylococcus and Corynebacterium
- Different phage composition

**Dry** (future): Forearm, leg
- Lower viral abundance
- Different bacterial hosts

### Abundance Model

**Tier 1** (10-25%): 5-10 genomes
- Major Cutibacterium phage lineages
- Most prevalent S. epidermidis phages

**Tier 2** (1-10%): 30 genomes
- Diverse Cutibacterium variants
- Common Staphylococcus phages

**Tier 3** (0.1-1%): 60 genomes
- Other Caudoviricetes
- Corynebacterium phages

**Tier 4** (0.01-0.1%): 40 genomes
- Rare phages
- Eukaryotic viruses

**Tier 5** (<0.01%): 15 genomes
- Transient viruses
- Novel families

### Selection Criteria

```sql
WHERE (genome_name LIKE '%Propionibacterium%'
    OR genome_name LIKE '%Cutibacterium%'
    OR genome_name LIKE '%Staphylococcus%'
    OR genome_name LIKE '%skin%'
    OR genome_name LIKE '%sebum%')
  AND family IN ('Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae')
```

### Validation Metrics

- [ ] Cutibacterium phages dominant (40-50%)
- [ ] Staphylococcus phages common (25-30%)
- [ ] Auxiliary metabolic genes present (check annotations)
- [ ] Site-specific composition matches literature

### Key References

1. Hannigan et al. (2015). "The human skin double-stranded DNA virome." mBio.
2. Oh et al. (2016). "Temporal stability of the human skin microbiome." Cell.
3. Recent 2024 studies on skin virome and auxiliary genes.

---

## Collection 4: Human Respiratory Virome (Nasopharynx, Healthy)

**Target Size**: 200 genomes
**Collection ID**: `respiratory_virome_nasopharynx_healthy`
**Planned Completion**: Week 7, Day 2

### Literature-Based Composition

**Abundance**:
- Lower phage abundance than gut
- Dominated by eukaryotic viruses in some studies
- Bacteriophage diversity higher in healthy vs. disease

**Dominant Groups** (based on recent 2024 studies, Lim et al. 2019):

**1. Siphoviridae - 60 genomes (30%)**
- Rationale: Dominant in healthy respiratory tract
- Hosts: Streptococcus, Haemophilus, Moraxella
- Distribution:
  - Streptococcus phages: 30 genomes
  - Haemophilus phages: 15 genomes
  - Moraxella phages: 10 genomes
  - Other respiratory bacteria: 5 genomes

**2. Cutibacterium (Propionibacterium) Phages - 40 genomes (20%)**
- Rationale: Pahexavirus most prevalent in nasopharynx (90.6% of individuals)
- C. acnes colonizes nasopharynx
- Persistent colonizer

**3. Staphylococcus Phages - 30 genomes (15%)**
- Rationale: 60% of individuals in nasopharynx
- S. aureus and S. epidermidis common

**4. Myoviridae - 20 genomes (10%)**
- Lower in healthy, increased in disease
- Lytic phages

**5. Podoviridae - 15 genomes (7.5%)**
- Diverse respiratory bacteria hosts

**6. Microviridae - 10 genomes (5%)**
- ssDNA phages
- Small, persistent

**7. Eukaryotic Viruses - 25 genomes (12.5%)**
- Rationale: Eukaryotic viruses more prominent in respiratory
- Distribution:
  - Rhinoviruses: 5 genomes (common cold)
  - Coronaviruses: 3 genomes (endemic strains)
  - Adenoviruses: 4 genomes
  - Herpesviruses: 5 genomes (HSV, EBV, CMV, HHV-6)
  - Bocaviruses: 3 genomes
  - Anelloviruses: 5 genomes

### Anatomical Site Variation

**Nasopharynx** (this collection): Upper respiratory tract
- High Cutibacterium phages
- High eukaryotic virus prevalence

**Oropharynx** (future): Similar to oral
- More Streptococcus phages
- Overlap with oral virome

**Lung** (future): Lower respiratory
- Lower viral abundance
- Different composition

### Abundance Model

**Tier 1** (5-15%): 10 genomes
- Major Cutibacterium phages (Pahexavirus)
- Prevalent Streptococcus phages

**Tier 2** (1-5%): 40 genomes
- Diverse Siphoviridae
- Common Staphylococcus phages
- Some eukaryotic viruses

**Tier 3** (0.1-1%): 80 genomes
- Myoviridae, Podoviridae
- Diverse respiratory phages
- Eukaryotic viruses

**Tier 4** (0.01-0.1%): 50 genomes
- Rare phages
- Transient viruses

**Tier 5** (<0.01%): 20 genomes
- Novel families
- Very rare viruses

### Selection Criteria

```sql
WHERE (genome_name LIKE '%Streptococcus%'
    OR genome_name LIKE '%Haemophilus%'
    OR genome_name LIKE '%Propionibacterium%'
    OR genome_name LIKE '%Cutibacterium%'
    OR genome_name LIKE '%Staphylococcus%'
    OR genome_name LIKE '%respiratory%'
    OR genome_name LIKE '%nasopharyn%')
  AND (family IN ('Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae')
    OR family IN ('Picornaviridae', 'Coronaviridae', 'Adenoviridae', 'Herpesviridae'))
```

### Validation Metrics

- [ ] Cutibacterium phages highly prevalent (20-25%)
- [ ] Siphoviridae dominant among phages (30-40%)
- [ ] Eukaryotic viruses present (10-15%)
- [ ] Lower phage abundance than gut (modeled in abundance)

### Key References

1. Recent 2024: "Metagenomic profiling of nasopharyngeal samples from adults with acute respiratory infection."
2. Recent 2023: "Respiratory eukaryotic virome expansion and bacteriophage deficiency characterize childhood asthma." Sci Rep.
3. Lim et al. (2019). "Altered respiratory virome and serum cytokine profile associated with recurrent respiratory tract infections in children." Nat Commun.

---

## Collection 5: Marine Virome (Coastal Surface Water)

**Target Size**: 500 genomes
**Collection ID**: `marine_virome_coastal_surface`
**Planned Completion**: Week 7, Day 3-4

### Literature-Based Composition

**Abundance**:
- Highest viral abundance: ~10^7 VLPs/mL
- Virus-to-microbe ratio: ~10:1
- Dominated by bacteriophages (>95%)

**Dominant Groups** (based on Roux et al. 2016, Brum et al. 2015):

**1. Caudoviricetes - 350 genomes (70%)**
- Rationale: Dominant in marine environments
- Distribution:
  - Pelagiphages (SAR11 hosts): 100 genomes (most abundant bacteria)
  - Cyanophages (Prochlorococcus, Synechococcus): 100 genomes
  - Vibriophages: 50 genomes
  - Roseophages (Roseobacter): 50 genomes
  - Other marine bacteria: 50 genomes

**2. Microviridae - 80 genomes (16%)**
- ssDNA phages
- Very abundant in marine
- Diverse hosts

**3. Phycodnaviridae - 30 genomes (6%)**
- Large algal viruses
- Important in marine ecosystems

**4. Mimiviridae - 15 genomes (3%)**
- Giant viruses
- Infect protists

**5. Podoviridae - 15 genomes (3%)**
- Marine cyanophages
- T7-like phages

**6. Inoviridae - 10 genomes (2%)**
- Filamentous phages
- Diverse marine hosts

### Marine Bacteria Host Distribution

**SAR11 clade** (most abundant ocean bacteria) - 100 phages
**Cyanobacteria** (Prochlorococcus, Synechococcus) - 100 phages
**Vibrio** - 50 phages
**Roseobacter** - 50 phages
**Flavobacteria** - 40 phages
**SAR86** - 30 phages
**Alteromonas** - 30 phages
**Other marine bacteria** - 100 phages

### Abundance Model

**Very high turnover** - different model than human body sites

**Tier 1** (5-20%): 20 genomes
- Dominant pelagiphages
- Major cyanophages

**Tier 2** (1-5%): 80 genomes
- Common Caudoviricetes
- Abundant Microviridae

**Tier 3** (0.1-1%): 200 genomes
- Diverse phages
- Moderate abundance

**Tier 4** (0.01-0.1%): 150 genomes
- Rare phages
- Giant viruses

**Tier 5** (<0.01%): 50 genomes
- Very rare
- Novel families

### Selection Criteria

```sql
WHERE (genome_name LIKE '%marine%'
    OR genome_name LIKE '%ocean%'
    OR genome_name LIKE '%Pelagibacter%'
    OR genome_name LIKE '%SAR11%'
    OR genome_name LIKE '%Prochlorococcus%'
    OR genome_name LIKE '%Synechococcus%'
    OR genome_name LIKE '%cyano%'
    OR genome_name LIKE '%Vibrio%'
    OR genome_name LIKE '%Roseobacter%')
  AND (family IN ('Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae')
    OR family IN ('Phycodnaviridae', 'Mimiviridae', 'Inoviridae'))
```

### Validation Metrics

- [ ] Caudoviricetes dominant (70-80%)
- [ ] Pelagiphages highly abundant (20-30%)
- [ ] Cyanophages well-represented (20-25%)
- [ ] High virus-to-bacteria ratio modeled

### Key References

1. Roux et al. (2016). "Ecogenomics and potential biogeochemical impacts of globally abundant ocean viruses." Nature.
2. Brum et al. (2015). "Patterns and ecological drivers of ocean viral communities." Science.
3. Paez-Espino et al. (2016). "Uncovering Earth's virome." Nature.

---

## Collection 6: Soil Virome (Agricultural Soil)

**Target Size**: 300 genomes
**Collection ID**: `soil_virome_agricultural`
**Planned Completion**: Week 7, Day 5

### Literature-Based Composition

**Abundance**:
- ~10^8-10^9 VLPs/g soil
- High diversity
- Many novel families

**Dominant Groups**:

**1. Caudoviricetes - 200 genomes (67%)**
- Siphoviridae, Myoviridae, Podoviridae
- Hosts: Actinobacteria, Proteobacteria, Firmicutes

**2. Microviridae - 40 genomes (13%)**
- ssDNA phages

**3. Tectiviridae - 20 genomes (7%)**
- Membrane-containing
- Bacillus phages

**4. Inoviridae - 20 genomes (7%)**
- Filamentous

**5. Novel/Unclassified - 20 genomes (6%)**
- Soil-specific families

### Target Hosts

- Actinobacteria: 80 phages
- Proteobacteria: 70 phages
- Firmicutes (Bacillus): 50 phages
- Other soil bacteria: 100 phages

---

## Collection 7: Freshwater Virome (Lake Surface Water)

**Target Size**: 200 genomes
**Collection ID**: `freshwater_virome_lake_surface`
**Planned Completion**: Week 7, Day 6

### Literature-Based Composition

**Similar to marine but different hosts**:

**1. Caudoviricetes - 130 genomes (65%)**
- Actinophages: 40 genomes
- Cyanophages (freshwater): 40 genomes
- Bacteroidetes phages: 30 genomes
- Other: 20 genomes

**2. Microviridae - 30 genomes (15%)**

**3. Phycodnaviridae - 20 genomes (10%)**
- Freshwater algal viruses

**4. Inoviridae - 10 genomes (5%)**

**5. Other/Novel - 10 genomes (5%)**

---

## Collection 8: Mouse Gut Virome (Laboratory Mouse)

**Target Size**: 150 genomes
**Collection ID**: `mouse_gut_virome_lab`
**Planned Completion**: Week 8, Day 1

### Literature-Based Composition

**Similar to human gut but different bacteria**:

**1. Caudoviricetes - 100 genomes (67%)**
- Lactobacillus phages: 30 genomes
- Bacteroides phages: 25 genomes
- Clostridiales phages: 25 genomes
- Other: 20 genomes

**2. Microviridae - 25 genomes (17%)**

**3. Inoviridae - 10 genomes (6%)**

**4. Eukaryotic Viruses - 15 genomes (10%)**
- Murine norovirus: 5 genomes
- Mouse hepatitis virus: 5 genomes
- Other: 5 genomes

---

## Implementation Timeline

### Week 6 (Current)

- **Day 6-7**: Gut virome curation (500 genomes) ✅ Plan ready
- **Day 7**: Oral virome curation (200 genomes)

### Week 7

- **Day 1**: Skin virome curation (150 genomes)
- **Day 2**: Respiratory virome curation (200 genomes)
- **Day 3-4**: Marine virome curation (500 genomes)
- **Day 5**: Soil virome curation (300 genomes)
- **Day 6**: Freshwater virome curation (200 genomes)

### Week 8

- **Day 1**: Mouse gut virome curation (150 genomes)
- **Day 2-3**: Validation of all collections
- **Day 4-5**: Documentation and examples
- **Day 6**: Integration testing

---

## Summary Statistics

| Collection | Size | Phage % | Eukaryotic % | Top Family | Key Feature |
|------------|------|---------|--------------|------------|-------------|
| Gut (human) | 500 | 96% | 4% | Crassvirales | crAssphage dominant |
| Oral (human) | 200 | 92.5% | 7.5% | Siphoviridae | Streptococcus phages |
| Skin (human) | 150 | 96% | 4% | Cutibacterium | P. acnes phages |
| Respiratory | 200 | 87.5% | 12.5% | Siphoviridae | Mixed phage/eukaryotic |
| Marine | 500 | 97% | 3% | Caudoviricetes | Pelagiphages, cyanophages |
| Soil | 300 | 94% | 6% | Caudoviricetes | Actinophages |
| Freshwater | 200 | 95% | 5% | Caudoviricetes | Actinophages, cyanophages |
| Mouse gut | 150 | 90% | 10% | Caudoviricetes | Lactobacillus phages |

**Total Genomes**: 2,200
**Average Collection Size**: 275 genomes
**Total Families Represented**: 50+

---

## Validation Framework

### Cross-Collection Validation

**Phage abundance validation**:
- Human gut > Oral > Respiratory (phage dominance)
- Marine > Soil > Freshwater (total viral load)

**Diversity validation**:
- Soil > Freshwater > Marine (expected diversity)
- Gut > Oral > Skin > Respiratory (human body sites)

**Family distribution**:
- Siphoviridae most common across all
- Crassvirales unique to gut
- Pelagiphages unique to marine
- Giant viruses in aquatic only

### Literature Comparison

Each collection will be validated against at least 3 published studies.

**Statistical tests**:
- KS test vs published abundance distributions
- Family distribution chi-square tests
- Diversity metrics (Shannon, Simpson)

---

## Future Enhancements

### Disease States

- IBD gut virome
- Periodontal disease oral virome
- Acne skin virome
- Asthma respiratory virome

### Age Groups

- Infant gut virome (0-2 years)
- Elderly gut virome (65+ years)

### Geographical Variations

- Non-Western gut virome
- Different ocean regions
- Different soil types

---

## Documentation Deliverables

For each collection:

1. **Curation Report** - Methodology, selection criteria, validation
2. **Collection Metadata** - Database entries, abundance models
3. **Usage Examples** - Python code for creating communities
4. **Validation Report** - Comparison to literature
5. **Reference List** - All cited studies

---

**Last Updated**: November 1, 2025
**Status**: Planning Complete, Implementation Starting Week 6
**Total Collections**: 8 body sites, 2,200 genomes
