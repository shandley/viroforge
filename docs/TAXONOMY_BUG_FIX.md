# Taxonomy Bug Fix Documentation

**Date**: 2025-11-09
**Phase**: 7-8 Collection Review
**Impact**: Critical - 4 collections regenerated

## Problem Discovery

During Phase 8 (RNA Virome Collections), we discovered that **46% of the database (6,651/14,423 genomes) had `family='Unknown'`**. This was caused by mismatches between:
- **RefSeq naming**: Strain-specific names (e.g., "Influenza A virus (A/California/07/2009(H1N1))")
- **ICTV taxonomy**: General species names (e.g., "influenza A virus")

### Initial Impact

**Collection 21 (Human Respiratory RNA)**:
- Only 31 genomes obtained (vs target 60-70)
- **ZERO influenza viruses** despite being critical respiratory pathogen
- Missing seasonal flu strains essential for benchmarking

This prompted a systematic review of all previous collections.

---

## Taxonomy Bug Review - Collections Affected

### Collection 19: HIV+ Gut Virome - **CRITICAL** ❌

**Status**: Scientifically invalid without fix

**Problem**:
- 0 Herpesviridae (family was 'Unknown' for all herpesviruses)
- HIV+ patients characteristically show herpesvirus reactivation:
  - **CMV (Human herpesvirus 5)**: Most common reactivation
  - **EBV (Human herpesvirus 4)**: Associated with lymphomas
  - **KSHV/HHV-8**: Kaposi's sarcoma (AIDS-defining illness)

**Before Fix**: 49 genomes, 0 Herpesviridae
**After Fix**: 55 genomes, 6 human Herpesviridae (EBV, KSHV, HSV-1, VZV, HHV-6B, HHV-7)

**Fix Applied**:
```python
# Updated herpesvirus query to:
WHERE t.family IN ('Herpesviridae', 'Orthoherpesviridae', 'Alloherpesviridae')
  AND (g.genome_name LIKE 'Human herpesvirus%'
   OR g.genome_name LIKE 'Human betaherpesvirus%'
   OR g.genome_name LIKE 'Human gammaherpesvirus%')
```

---

### Collection 17: Wastewater Virome - **MISSING KEY PATHOGENS** ⚠️

**Problem**: Missing rotavirus (major sewage pathogen)

**Before Fix**: 321 genomes, 0 Sedoreoviridae (rotavirus)
**After Fix**: 352 genomes (+9.7%), 8 rotavirus + 32 norovirus

**Fix Applied**: Added Sedoreoviridae query to enteric virus selection

---

### Collection 23: Fecal RNA Virome - **SEVERELY INCOMPLETE** ⚠️

**Problem**: Missing primary enteric virus targets

**Before Fix**: 32 genomes, 1 norovirus, 0 rotavirus
**After Fix**: 58 genomes (+81%!), 15 norovirus + 12 rotavirus

**Fix Applied**:
```python
# Updated reovirus query from:
WHERE t.family = 'Reoviridae'
# To:
WHERE (t.family = 'Reoviridae' OR t.family = 'Sedoreoviridae')
```

---

### Collection 20: CF Respiratory Virome - **INCOMPLETE** ⚠️

**Problem**: Missing 40% of available influenza viruses

**Before Fix**: 77 genomes, 6 influenza (60% of available)
**After Fix**: 81 genomes (+5.2%), 10 influenza (100% coverage)

**Fix Applied**: Increased influenza limit from 6 to 10

---

### Collection 18: IBD Gut Virome - ✓ NO CHANGES NEEDED

**Status**: Working as designed (phage-dominated dysbiotic gut)
**Reason**: Enteric viruses (rotavirus/norovirus) not primary focus for chronic IBD

---

## Taxonomy Fix Solution

### Phase 1: Enhanced Fuzzy Matching (443 genomes fixed)

**Added pattern-based family matching**:
```python
family_patterns = {
    'herpesvirus': ['Orthoherpesviridae', 'Alloherpesviridae'],
    'cytomegalovirus': ['Orthoherpesviridae'],
    'papillomavirus': ['Papillomaviridae'],
    'polyomavirus': ['Polyomaviridae'],
    'adenovirus': ['Adenoviridae'],
    'parvovirus': ['Parvoviridae'],
    'retrovirus': ['Retroviridae'],
    'lentivirus': ['Retroviridae'],
    'nucleopolyhedrovirus': ['Baculoviridae'],
    'respirovirus': ['Paramyxoviridae'],
    'morbillivirus': ['Paramyxoviridae'],
    # ... 20+ patterns total
}
```

**Enhanced normalization**:
- Remove "type" keyword: "papillomavirus type 60" → "papillomavirus 60"
- Remove trailing numbers: "virus - 4" → "virus"
- Remove subtype information: "subtype 1" → ""
- Case-insensitive matching

**Results**:
- Fixed 443 genomes (6.7% of 6,651 unmatched)
- Success examples: 10 influenza, herpesviruses, papillomaviruses, retroviruses

### Phase 2: Rotavirus/Norovirus Fix (30 genomes fixed, 96.8% success)

**Problem**: Critical enteric viruses still unmatched

**Added patterns**:
```python
'rotavirus': ['Sedoreoviridae'],  # Rotavirus A, B, C, etc.
'norovirus': ['Caliciviridae'],   # Norovirus GI, GII, etc.
```

**Why needed**:
- ICTV: "rotavirus A", "Norwalk virus"
- RefSeq: "Rotavirus B strain Bang373", "Norovirus GII.17"

**Results**:
- Caliciviridae (norovirus): 0 → 20 genomes
- Sedoreoviridae (rotavirus): 4 → 10 genomes
- Only 1 genome remains unmatched (96.8% success)

---

## Final Statistics

### Database Coverage
- **Before fixes**: 7,772/14,423 genomes assigned (53.9%)
- **After fixes**: 8,241/14,423 genomes assigned (57.1%)
- **Total fixed**: 469 genomes (3.2% improvement)
- **Still Unknown**: 6,182 genomes (mostly phages, plant viruses)

### Collection Impact Summary

| Collection | Before | After | Change | Status |
|------------|--------|-------|--------|--------|
| 19 (HIV+) | 49 genomes, 0 herpes | 55 genomes, 6 herpes | +12% | **CRITICAL FIX** |
| 17 (Wastewater) | 321, 0 rotavirus | 352, 8 rotavirus | +10% | Fixed |
| 23 (Fecal RNA) | 32, 1 norovirus | 58, 15 norovirus | +81% | **MAJOR FIX** |
| 20 (CF Respiratory) | 77, 6 influenza | 81, 10 influenza | +5% | Fixed |
| 18 (IBD) | 90 genomes | 90 genomes | 0% | ✓ No change needed |

---

## Lessons Learned

### 1. Always Review Collection Sizes Against Targets
- Collection 21 had only 31/60 genomes (52% of target) - RED FLAG
- Collection 23 had only 32/40 genomes (80% of target) - Should have investigated

### 2. Missing Key Pathogens = Database Issue
- No influenza in respiratory collection → investigate database
- No rotavirus in fecal collection → investigate database
- Don't assume RefSeq doesn't have the data

### 3. Taxonomy Assignment is Critical
- 46% unmatched is a severe problem
- Affects collection quality and scientific validity
- Pattern-based matching essential for strain-specific nomenclature

### 4. Systematic Review is Essential
When a bug is discovered:
1. Review ALL previous collections
2. Identify which are affected
3. Prioritize by scientific impact (Collection 19 was CRITICAL)
4. Regenerate affected collections
5. Document thoroughly

---

## Prevention for Future Work

### Quality Checks After Curation

```python
# Always check collection size vs target
if len(collection) < target * 0.8:
    logger.warning(f"Collection below 80% of target: {len(collection)}/{target}")

# Always verify key families are present
expected_families = ['Orthomyxoviridae', 'Caliciviridae', 'Sedoreoviridae']
for family in expected_families:
    count = len([g for g in collection if g['family'] == family])
    if count == 0:
        logger.warning(f"Expected family {family} has 0 genomes!")
```

### When Adding New Collections

1. Check database coverage for target families FIRST
2. If key families missing, investigate taxonomy assignment
3. Don't proceed with incomplete collections
4. Document any known limitations

---

## Scripts Modified

1. **scripts/fix_taxonomy_unmatched.py**
   - Enhanced fuzzy matching with 20+ virus family patterns
   - Improved normalization (type, subtype, trailing numbers)
   - Added rotavirus/norovirus patterns

2. **scripts/curate_hiv_gut_collection.py**
   - Updated to query Orthoherpesviridae (not Herpesviridae)
   - Filter for HUMAN herpesviruses specifically

3. **scripts/curate_wastewater_collection.py**
   - Added Sedoreoviridae (rotavirus) query
   - Adjusted enteric virus percentages

4. **scripts/curate_fecal_rna_collection.py**
   - Updated to include Sedoreoviridae OR Reoviridae
   - Fixed rotavirus query

5. **scripts/curate_cf_respiratory_collection.py**
   - Increased influenza limit to capture all 10 genomes

---

## Git Commits

1. `df2344a` - feat: enhance taxonomy fuzzy matching for unmatched genomes
2. `e78f344` - fix: CRITICAL - Regenerate Collection 19 (HIV+) with human herpesviruses
3. `79a4fd4` - fix: Add rotavirus/norovirus patterns to taxonomy fuzzy matching
4. `452c636` - fix: Regenerate Collections 17 & 23 with rotavirus/norovirus
5. `a71a74f` - fix: Regenerate Collection 20 (CF) with all 10 influenza viruses

---

## Future Considerations

### Remaining Unknown Genomes (6,182)

**Why still unmatched**:
- Plant viruses with complex nomenclature
- Bacteriophages with genus-specific naming
- Novel/unclassified viruses
- Possible ICTV vs RefSeq nomenclature divergence

**Not urgent** because:
- Mostly affect specialized collections (plant, environmental phages)
- All critical human pathogens now assigned
- Further enhancement requires deeper ICTV analysis

### Potential Future Enhancements

1. **Genus-level matching**: Use genus information when species fails
2. **Synonym handling**: Build ICTV synonym database
3. **Machine learning**: Train classifier on matched genomes
4. **Manual curation**: For remaining important viruses

---

## Contact

**Issue discovered by**: ViroForge Development Team (Phase 8)
**Fixed by**: Claude Code assistant
**Documentation**: 2025-11-09

**For questions**: Review git commits df2344a through a71a74f
