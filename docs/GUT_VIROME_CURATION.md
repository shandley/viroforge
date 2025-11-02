# Gut Virome Collection Curation Plan

**Date**: November 1, 2025
**Collection**: Human Gut Virome (Adult, Healthy, Western Diet)
**Target**: 500 genomes with literature-validated abundances
**Status**: Planning Phase

---

## Literature-Based Composition

### Key Findings from Recent Research (2023-2025)

**Overall Composition**:
- **Phages**: 97.7% of gut virome
- **Eukaryotic viruses**: 2.1%
- **Archaeal viruses**: 0.1%

**Dominant Viral Groups**:

1. **Crassvirales (crAssphage)** - MOST ABUNDANT
   - Prevalence: 98-100% of individuals
   - Abundance: Up to 90-95% of viral sequences in some individuals
   - Hosts: Bacteroidota (primarily Bacteroides species)
   - Status: Most abundant virus in human gut
   - Note: Multiple crAss-like phage families identified

2. **Caudoviricetes** (formerly Caudovirales) - DOMINANT CLASS
   - Comprises majority of gut phages
   - Traditional families (being reclassified):
     - **Siphoviridae**: Long non-contractile tails (most common, temperate)
     - **Myoviridae**: Contractile tails (complex structure)
     - **Podoviridae**: Short non-contractile tails
   - Note: Genome-based reclassification ongoing

3. **Microviridae** - STABLE COLONIZERS
   - Type: ssDNA bacteriophages
   - Characteristics: Stable, persistent gut colonizers
   - Size: Small genomes (~5-6 kb)

4. **Bacteroides-infecting phages**
   - Major target: Bacteroides (dominant gut bacteria)
   - Diversity: High, many novel families

### Age-Related Differences

**Adults** (Western diet):
- Dominated by Crassvirales (crAssphage)
- High Bacteroides-phage diversity
- Siphoviridae prevalent

**Infants**:
- Fewer crAssphages
- More Clostridiales-infecting phages
- More Bifidobacterium-infecting phages
- Different composition (excluded from this collection)

### Diversity and "Dark Matter"

- **Known families**: ~50-100 identified in gut
- **Novel diversity**: 232 previously unknown virus family-level clades (VFCs) identified in recent metagenomics
- **Total species**: ~10,000 viral species in comprehensive studies
- **Dark matter**: Many sequences don't match known viruses

---

## Curation Strategy

### Phase 1: Taxonomic Selection (Target: 500 genomes)

**1. Crassvirales - 150 genomes (30%)**
   - Rationale: Most abundant, 90-95% in some individuals
   - Distribution:
     - crAssphage prototypes: 20 genomes
     - crAss-like families: 130 genomes (diverse variants)
   - Host range: Bacteroidota (Bacteroides, Prevotella, etc.)
   - Selection criteria: Representatives from all major crAss-like clades

**2. Caudoviricetes (non-Crassvirales) - 250 genomes (50%)**
   - Rationale: Majority of remaining gut phages
   - Distribution by traditional family:
     - **Siphoviridae-like**: 140 genomes (28%)
       - Temperate phages (lysogenic)
       - Bacteroides, Clostridiales, Firmicutes hosts
       - Long non-contractile tails
     - **Myoviridae-like**: 70 genomes (14%)
       - Lytic phages
       - Complex contractile tails
       - Diverse bacterial hosts
     - **Podoviridae-like**: 40 genomes (8%)
       - Short-tailed phages
       - Various gut bacteria hosts

**3. Microviridae - 50 genomes (10%)**
   - Rationale: Stable colonizers, ssDNA
   - Distribution:
     - Gokushovirinae: 30 genomes
     - Bullavirinae: 20 genomes
   - Host range: Bacteroidetes, Firmicutes
   - Selection: Representatives of persistent strains

**4. Other Families - 30 genomes (6%)**
   - Inoviridae: 10 genomes (filamentous phages)
   - Tectiviridae: 5 genomes (membrane-containing)
   - Corticoviridae: 5 genomes
   - Novel/unclassified gut phages: 10 genomes

**5. Eukaryotic Viruses - 20 genomes (4%)**
   - Rationale: 2.1% of gut virome
   - Distribution:
     - Adenoviridae: 5 genomes (human pathogens)
     - Anelloviridae: 5 genomes (commensal)
     - Picobirnaviridae: 5 genomes (dsRNA)
     - Other eukaryotic: 5 genomes

### Phase 2: Host-Based Selection

**Target Bacterial Hosts** (reflecting gut microbiome composition):

1. **Bacteroidetes** (30-40% of gut bacteria) - 250 phages
   - Bacteroides: 150 phages
   - Prevotella: 50 phages
   - Parabacteroides: 30 phages
   - Other Bacteroidetes: 20 phages

2. **Firmicutes** (20-30% of gut bacteria) - 180 phages
   - Clostridiales: 80 phages
   - Faecalibacterium: 40 phages
   - Roseburia: 30 phages
   - Other Firmicutes: 30 phages

3. **Actinobacteria** (5-10%) - 40 phages
   - Bifidobacterium: 25 phages (important in some individuals)
   - Other Actinobacteria: 15 phages

4. **Proteobacteria** (<5%) - 20 phages
   - E. coli: 10 phages
   - Other Proteobacteria: 10 phages

5. **Other/Unknown hosts** - 10 phages

### Phase 3: Abundance Modeling

**Abundance Distribution Strategy**:

Based on literature (Manrique et al. 2016, Shkoporov et al. 2019, Guerin et al. 2018):

**Tier 1: Dominant (10-30% relative abundance each)** - 5-10 genomes
- crAssphage and major crAss-like variants
- Log-normal distribution with high mean
- Examples: crAssphage φ1, major Bacteroides phages

**Tier 2: Common (1-10% relative abundance)** - 50-100 genomes
- Common Siphoviridae
- Prevalent Microviridae
- Major Bacteroides/Firmicutes phages
- Log-normal distribution, moderate mean

**Tier 3: Moderate (0.1-1% relative abundance)** - 150-200 genomes
- Diverse Caudoviricetes
- Less common Bacteroides phages
- Clostridiales phages
- Log-normal distribution, lower mean

**Tier 4: Rare (0.01-0.1% relative abundance)** - 200-250 genomes
- Novel families
- Rare Myoviridae/Podoviridae
- Eukaryotic viruses
- Long-tail distribution

**Tier 5: Very Rare (<0.01% relative abundance)** - 50-100 genomes
- Archaeal viruses
- Rare eukaryotic viruses
- Transient phages
- Power-law tail

**Mathematical Model**:
```python
# Log-normal distribution for most abundant phages
# Tier 1-3: log-normal(μ=varies, σ=1.5)
# Tier 4-5: power-law tail (α=2.0)

# Relative abundance constraints:
# Sum(Tier 1) = 30-50% total virome
# Sum(Tier 2) = 20-30% total virome
# Sum(Tier 3) = 15-25% total virome
# Sum(Tier 4) = 5-10% total virome
# Sum(Tier 5) = 1-5% total virome
```

---

## Selection Criteria from Database

### Query Strategy

**Step 1: Identify Crassvirales**
```sql
SELECT genome_id, genome_name, family, genus, species
FROM genomes g
JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.order_name = 'Crassvirales'
   OR t.family LIKE '%crass%'
   OR g.genome_name LIKE '%crass%'
ORDER BY length DESC;
```

**Step 2: Identify Caudoviricetes**
```sql
SELECT genome_id, genome_name, family, genus, species, length, gc_content
FROM genomes g
JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.class = 'Caudoviricetes'
   AND t.order_name != 'Crassvirales'
ORDER BY family, length DESC;
```

**Step 3: Filter by Host (if available)**
```sql
SELECT g.genome_id, g.genome_name, h.host_species, h.host_taxid
FROM genomes g
JOIN host_associations h ON g.genome_id = h.genome_id
WHERE h.host_species LIKE '%Bacteroides%'
   OR h.host_species LIKE '%Clostrid%'
   OR h.host_species LIKE '%Firmicutes%';
```

**Step 4: Filter by Isolation Source**
```sql
SELECT g.genome_id, g.genome_name, e.body_site, e.environment
FROM genomes g
JOIN ecological_metadata e ON g.genome_id = e.genome_id
WHERE e.body_site LIKE '%gut%'
   OR e.body_site LIKE '%fec%'
   OR e.body_site LIKE '%intestin%'
   OR e.environment LIKE '%human%';
```

### Fallback Strategy (If Limited Metadata)

If host/ecological metadata is sparse, use taxonomic proxies:

**Proxy 1: Family-based**
- Families commonly found in gut: Siphoviridae, Myoviridae, Podoviridae, Microviridae, Inoviridae
- Exclude non-gut families: Marine viruses, plant viruses, etc.

**Proxy 2: Genome characteristics**
- Moderate GC content (30-50% - typical for gut bacteria)
- Genome size ranges:
  - Caudoviruses: 30-200 kb
  - Microviridae: 4-8 kb
  - Inoviridae: 6-12 kb

**Proxy 3: NCBI metadata mining**
- Parse genome names for "gut", "fecal", "intestinal", "human"
- Use NCBI BioSample information
- Cross-reference with IMGVR database

---

## Validation Against Literature

### Key Studies for Validation

1. **Guerin et al. (2018)** - "Biology and Taxonomy of crAss-like Bacteriophages"
   - Use for crAssphage abundance models
   - Validate dominance in collection

2. **Shkoporov et al. (2019)** - "The Human Gut Virome Is Highly Diverse, Stable, and Individual Specific"
   - Validate overall diversity
   - Check family-level distribution

3. **Manrique et al. (2016)** - "Healthy Human Gut Phageome"
   - Validate Caudovirales dominance
   - Compare Siphoviridae/Myoviridae/Podoviridae ratios

4. **Dutilh et al. (2014)** - "A highly abundant bacteriophage discovered in the unknown sequences of human faecal metagenomes"
   - Original crAssphage discovery
   - Validate crAssphage is most abundant

### Validation Metrics

**Composition Validation**:
- [ ] crAssphage represents 25-40% of collection
- [ ] Caudoviricetes represents 70-80% of collection
- [ ] Microviridae represents 5-15% of collection
- [ ] Eukaryotic viruses represent 2-5% of collection

**Diversity Validation**:
- [ ] 50+ viral families represented
- [ ] Bacteroides-phages most abundant (30-40%)
- [ ] Firmicutes-phages well represented (20-30%)
- [ ] Both temperate and lytic phages included

**Abundance Model Validation**:
- [ ] Log-normal distribution for main population
- [ ] Power-law tail for rare viruses
- [ ] Matches published virome abundance distributions

---

## Implementation Plan

### Week 6 Tasks

**Task 1: Database Query and Initial Selection** (Day 1)
1. Wait for full RefSeq download and database population
2. Query database for Crassvirales genomes
3. Query for Caudoviricetes genomes
4. Identify Microviridae genomes
5. Export candidate lists

**Task 2: Metadata Enrichment** (Day 2)
1. Parse NCBI BioSample information
2. Extract host information from genome descriptions
3. Mine isolation source data
4. Cross-reference with IMGVR database (if available)

**Task 3: Manual Curation** (Day 3)
1. Review top 100 candidates for accuracy
2. Ensure crAssphage is well-represented
3. Check for major gut phage families
4. Remove obvious non-gut viruses (marine, plant, etc.)

**Task 4: Abundance Model Creation** (Day 4)
1. Assign abundance tiers to genomes
2. Create log-normal + power-law model
3. Validate total abundances sum to 1.0
4. Compare to literature distributions

**Task 5: Collection Population** (Day 5)
1. Insert into body_site_collections table
2. Link genomes via collection_genomes table
3. Store abundance values
4. Validate collection integrity

**Task 6: Documentation and Validation** (Day 6)
1. Document curation methodology
2. Create validation report
3. Compare to published gut virome compositions
4. Prepare example usage code

---

## Expected Output

### Database Tables Updated

**body_site_collections**:
```sql
INSERT INTO body_site_collections VALUES (
    1,  -- collection_id
    'gut_virome_adult_healthy_western',  -- collection_name
    'Human gut virome - adult, healthy, Western diet',  -- description
    'gut',  -- body_site
    'Homo sapiens',  -- host_species
    9606,  -- host_taxid
    500,  -- total_genomes
    'healthy',  -- health_status
    'Western',  -- diet_type
    'adult',  -- age_group
    'Guerin et al. 2018; Shkoporov et al. 2019',  -- literature_references
    'log-normal + power-law',  -- abundance_model
    '2025-11-01',  -- created_date
    '2025-11-01'   -- last_updated
);
```

**collection_genomes** (500 rows):
```sql
INSERT INTO collection_genomes VALUES (
    1,  -- collection_id (gut_virome)
    'GCF_XXXXXXX.1',  -- genome_id
    0.15,  -- relative_abundance (e.g., crAssphage at 15%)
    0.98,  -- prevalence (98% of individuals)
    1,  -- is_core_member (always present)
    'Dominant crAssphage variant',  -- notes
    'Dutilh et al. 2014'  -- source_reference
);
```

### Usage Example

```python
from viroforge.core import create_body_site_profile

# Use pre-curated gut virome collection
gut_virome = create_body_site_profile(
    body_site='gut',
    collection='gut_virome_adult_healthy_western',
    random_seed=42
)

# Results:
# - 500 genomes with literature-based abundances
# - crAssphage dominant (25-40%)
# - Diverse Caudoviricetes
# - Realistic abundance distribution
# - Validated against published studies
```

---

## Quality Metrics

### Completeness Metrics

- [ ] **Taxonomic Coverage**: 50+ families represented
- [ ] **Abundance Range**: 5 orders of magnitude (10% to 0.001%)
- [ ] **Host Diversity**: Bacteroidetes, Firmicutes, Actinobacteria, Proteobacteria
- [ ] **Genome Size Range**: 4 kb (Microviridae) to 200 kb (large Myoviridae)
- [ ] **GC Diversity**: 30-55% (reflecting host GC content)

### Accuracy Metrics

- [ ] **Literature Alignment**: Family distributions match published studies (±10%)
- [ ] **Dominance**: crAssphage is most abundant viral group
- [ ] **Realism**: Abundance distribution matches log-normal + power-law
- [ ] **Composition**: 97-98% phages, 2-3% eukaryotic viruses

### Validation Metrics

- [ ] **Cross-Study Comparison**: Matches at least 3 independent gut virome studies
- [ ] **Statistical Tests**: KS-test p>0.05 vs published abundances
- [ ] **Expert Review**: Curated by virome researcher
- [ ] **Reproducibility**: Same seed produces same composition

---

## Future Enhancements

### Additional Gut Virome Collections

1. **Gut Virome - Infant** (0-2 years)
   - More Bifidobacterium phages
   - Fewer crAssphages
   - Different Clostridiales phages

2. **Gut Virome - Non-Western** (Non-industrialized)
   - Different Bacteroides phage composition
   - More Prevotella phages
   - Regional diversity

3. **Gut Virome - IBD** (Inflammatory Bowel Disease)
   - Dysbiotic composition
   - Altered Caudoviricetes diversity
   - Reduced crAssphage in some studies

4. **Gut Virome - Post-Antibiotic**
   - Reduced diversity
   - Siphoviridae depletion
   - Recovery dynamics

5. **Gut Virome - Elderly**
   - Age-related changes
   - Different stability patterns

---

## References

1. Dutilh et al. (2014). "A highly abundant bacteriophage discovered in the unknown sequences of human faecal metagenomes." PNAS.

2. Guerin et al. (2018). "Biology and Taxonomy of crAss-like Bacteriophages, the Most Abundant Virus in the Human Gut." Cell Host & Microbe.

3. Shkoporov et al. (2019). "The Human Gut Virome Is Highly Diverse, Stable, and Individual Specific." Cell Host & Microbe.

4. Manrique et al. (2016). "Healthy Human Gut Phageome." PNAS.

5. Clooney et al. (2019). "Whole-Virome Analysis Sheds Light on Viral Dark Matter in Inflammatory Bowel Disease." Cell Host & Microbe.

6. Camarillo-Guerrero et al. (2021). "Massive expansion of human gut bacteriophage diversity." Cell.

7. Benler et al. (2021). "Thousands of previously unknown phages discovered in whole-community human gut metagenomes." Microbiome.

---

## Status

**Current Phase**: Planning Complete ✅
**Next Phase**: Database Query and Selection (awaiting full RefSeq database)
**Estimated Completion**: Week 6, Day 5
**Target Collection Size**: 500 genomes
**Expected Completion Date**: November 6, 2025

---

**Last Updated**: November 1, 2025
**Author**: Scott Handley + Claude Code
**Document Version**: 1.0
