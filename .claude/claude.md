# ViroForge - Claude Code Context

## Project Overview

**ViroForge** is a synthetic viral metagenomics dataset generator for benchmarking virome analysis pipelines. It generates realistic FASTQ sequencing data with known ground truth for validation.

**Key Features**:
- 23 curated viral genome collections (body sites, disease states, environmental)
- Realistic abundance distributions (log-normal, power-law)
- Multi-step library preparation modeling (VLP extraction, amplification bias, sequencing)
- RNA and DNA virome workflows
- 14,423 viral genomes from RefSeq with ICTV taxonomy

**Target Users**: Bioinformatics researchers developing/benchmarking virome analysis tools

---

## Current Status (v0.6.0)

### Phases Completed

- ✅ **Phase 1-5**: Core infrastructure, 8 original collections, FASTQ generation
- ✅ **Phase 6**: Amplification bias modeling (MDA, SISPA, TruSeq)
- ✅ **Phase 7**: Critical collections (wastewater, IBD, HIV+, CF respiratory)
- ✅ **Phase 8.1**: RNA virome collections (respiratory, arbovirus, fecal RNA)
- ✅ **Taxonomy Bug Fix**: Enhanced fuzzy matching, fixed 469 genomes

### Active Development

**Phase 8.2**: RNA Workflow Components (PENDING)
- Reverse transcription modeling
- rRNA depletion (Ribo-Zero)
- RNA degradation modeling
- `--molecule-type {dna,rna}` flag

---

## CRITICAL: Taxonomy Bug (2025-11-09)

**⚠️ IMPORTANT**: We discovered and fixed a major taxonomy assignment bug affecting 46% of the database.

### What Happened

**Problem**: RefSeq uses strain-specific names, ICTV uses general species names
- RefSeq: "Influenza A virus (A/California/07/2009(H1N1))"
- ICTV: "influenza A virus"
- Result: 6,651/14,423 genomes (46%) had `family='Unknown'`

### Collections Affected

**CRITICAL** - Collection 19 (HIV+ Gut):
- Had ZERO herpesviruses (scientifically invalid - HIV+ patients show herpesvirus reactivation)
- Fixed: 0 → 6 human herpesviruses (EBV, KSHV, HSV-1, VZV, HHV-6B, HHV-7)

**MAJOR** - Collection 23 (Fecal RNA):
- Only 32 genomes, missing rotavirus/norovirus (primary enteric pathogens)
- Fixed: 32 → 58 genomes (+81%), now has 15 norovirus + 12 rotavirus

**Collections 17 (Wastewater) & 20 (CF Respiratory)**: Also fixed

### Solution Applied

Enhanced `scripts/fix_taxonomy_unmatched.py` with:
1. Pattern-based family matching (20+ virus families)
2. Improved normalization (remove "type", trailing numbers, etc.)
3. Fuzzy matching for strain-specific names

**Results**: Fixed 469 genomes (7.1% of unmatched)

### Full Documentation

📄 **See**: `docs/TAXONOMY_BUG_FIX.md` for complete details, lessons learned, and prevention strategies

**Git commits**: df2344a, e78f344, 79a4fd4, 452c636, a71a74f

---

## Database Structure

### Location
```
viroforge/data/viral_genomes.db
```

### Key Tables

**genomes**: 14,423 viral genomes
- genome_id (RefSeq accession)
- genome_name, length, gc_content, sequences

**taxonomy**: ICTV taxonomy (NOW 57.1% assigned after fix)
- genome_id → realm, kingdom, phylum, class, order_name, **family**, genus, species
- **NOTE**: `family` has NOT NULL constraint - use 'Unknown' for unassigned
- **IMPORTANT**: After taxonomy fix, check family != 'Unknown' for critical viruses

**body_site_collections**: 23 curated collections
- collection_id, collection_name, description, n_genomes
- literature_references (verify citations!)

**collection_genomes**: Genome-collection associations
- collection_id, genome_id, relative_abundance, prevalence, abundance_rank

---

## Collections Overview (23 Total)

### Original Collections (1-8)
1. Healthy Human Gut Virome (134 genomes)
2. Healthy Human Skin Virome (42 genomes)
3. Healthy Human Oral Virome (67 genomes)
4. Healthy Human Urogenital Virome (31 genomes)
5. Healthy Human Respiratory Virome (58 genomes)
6. Marine Virome (78 genomes)
7. Soil Virome (82 genomes)
8. Freshwater Virome (71 genomes)

### VLP Comparison Collections (9-15)
9. Healthy Gut (Comparison baseline) (134 genomes)
10. High Bacterial Lysis (142 genomes)
11. Low Bacterial Lysis (128 genomes)
12. High Prophage Induction (149 genomes)
13. Low Prophage Induction (127 genomes)
14. High Eukaryotic Virus Shedding (139 genomes)
15. Low Eukaryotic Virus Shedding (127 genomes)

### Amplification Comparison (16)
16. Pre-Amplification Control (100 genomes)

### Critical Disease/Environmental Collections (17-20) - **FIXED AFTER TAXONOMY BUG**
17. Wastewater Virome (352 genomes) - Now includes rotavirus
18. IBD Gut Virome (90 genomes)
19. HIV+ Gut Virome (55 genomes) - **CRITICAL FIX**: Now has herpesviruses
20. CF Respiratory Virome (81 genomes) - Now has all 10 influenza

### RNA Virome Collections (21-23) - **FIXED AFTER TAXONOMY BUG**
21. Human Respiratory RNA Virome (56 genomes) - Now includes influenza
22. Arbovirus Environmental (39 genomes)
23. Fecal RNA Virome (58 genomes) - **MAJOR FIX**: +81% size, now has rotavirus/norovirus

---

## Common Workflows

### Creating a New Collection

```python
# 1. Check database coverage FIRST
sqlite3 viroforge/data/viral_genomes.db "
SELECT t.family, COUNT(*)
FROM genomes g
JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.family = 'TargetFamily'
GROUP BY t.family"

# 2. If count is 0 or low, investigate taxonomy
# Don't assume RefSeq doesn't have the data - check for family='Unknown'

# 3. Create curation script in scripts/
# Follow pattern: curate_<collection_name>_collection.py

# 4. Always verify collection size vs target
if len(collection) < target * 0.8:
    logger.warning(f"Collection below target - investigate!")

# 5. Verify critical families are present
for family in expected_families:
    count = len([g for g in collection if g['family'] == family])
    if count == 0:
        logger.warning(f"Missing expected family: {family}")
```

### Running Taxonomy Fix

```bash
# Enhance fuzzy matching (if adding new patterns)
# Edit: scripts/fix_taxonomy_unmatched.py
python3 scripts/fix_taxonomy_unmatched.py

# Verify results
sqlite3 viroforge/data/viral_genomes.db "
SELECT COUNT(*) FROM taxonomy WHERE family = 'Unknown'"
```

### Generating FASTQ Datasets

```bash
# Test with dry-run first
python scripts/generate_fastq_dataset.py --collection-id 21 --dry-run

# Generate full dataset
python scripts/generate_fastq_dataset.py \
  --collection-id 21 \
  --output-dir datasets/respiratory_rna \
  --num-reads 1000000 \
  --read-length 150 \
  --sequencing-platform novaseq \
  --vlp-protocol standard \
  --amplification-method MDA
```

---

## Important Files

### Core Scripts
- `scripts/fix_taxonomy_unmatched.py` - **CRITICAL**: Fuzzy matching for taxonomy assignment
- `scripts/generate_fastq_dataset.py` - Main FASTQ generation
- `scripts/batch_generate_fastq.py` - Generate multiple datasets
- `scripts/curate_*_collection.py` - Collection curation scripts (23 total)

### Documentation
- `docs/TAXONOMY_BUG_FIX.md` - **READ FIRST** if revisiting taxonomy issues
- `docs/COLLECTION_IMPLEMENTATION_GUIDE.md` - All 23 collections documented
- `docs/PHASE4_FASTQ_GENERATION.md` - FASTQ generation guide
- `ROADMAP.md` - Development roadmap (current: Phase 8)

### Database
- `viroforge/data/viral_genomes.db` - Main database (14,423 genomes)
- `data/ictv/ictv_taxonomy.json` - ICTV reference (17,925 records)

---

## Development Guidelines

### Citations
**CRITICAL**: Always verify citations with web search
- Include DOI and PMID when possible
- Use proper journal names (not abbreviations)
- Check publication years
- See `docs/CITATION_CORRECTIONS.md` for past errors

### Code Quality
- Follow existing patterns in curation scripts
- Add comprehensive docstrings
- Log progress at INFO level
- Include example outputs in docstrings

### Testing
```bash
# Always test collections with dry-run
python scripts/generate_fastq_dataset.py --collection-id <ID> --dry-run

# Run integration tests
pytest tests/test_fastq_integration.py -v
```

---

## Known Issues

### Taxonomy Assignment
- **6,182 genomes still have `family='Unknown'` (42.8%)**
  - Mostly: Plant viruses, bacteriophages, novel viruses
  - Not critical: All major human pathogens now assigned
  - See `data/ictv/still_unmatched_after_fix.tsv`

### Collection Gaps
- No long-read sequencing support (PacBio/Nanopore) - planned for Phase 10
- No temporal dynamics - planned for Phase 11
- Limited to human and environmental viromes currently

---

## Quick Reference

### Check Collection Status
```bash
sqlite3 viroforge/data/viral_genomes.db "
SELECT c.collection_id, c.collection_name, c.n_genomes
FROM body_site_collections c
ORDER BY c.collection_id"
```

### Check Taxonomy Coverage
```bash
sqlite3 viroforge/data/viral_genomes.db "
SELECT
  CASE WHEN family = 'Unknown' THEN 'Unassigned' ELSE 'Assigned' END as status,
  COUNT(*) as count,
  ROUND(COUNT(*) * 100.0 / (SELECT COUNT(*) FROM taxonomy), 1) as percent
FROM taxonomy
GROUP BY status"
```

### Find Genomes by Family
```bash
sqlite3 viroforge/data/viral_genomes.db "
SELECT g.genome_name, t.family, t.genus
FROM genomes g
JOIN taxonomy t ON g.genome_id = t.genome_id
WHERE t.family = 'Orthomyxoviridae'
LIMIT 10"
```

---

## Next Steps (from ROADMAP.md)

**Immediate (Phase 8.2)**:
- Implement RNA workflow components (RT, rRNA depletion)
- Add `--molecule-type {dna,rna}` flag
- Test RNA collections (21-23)

**Short-term (Phase 9)**:
- 5 additional host-associated collections
- Blood/plasma, ocular, lung, urinary viromes
- Clinical research applications

**Medium-term (Phase 10-11)**:
- Long-read sequencing support (PacBio, Nanopore)
- Temporal dynamics modeling
- Longitudinal study generation

---

## Contact & History

**Project**: ViroForge
**Repository**: hecatomb/viroforge
**Current Version**: 0.6.0
**Last Major Update**: 2025-11-09 (Taxonomy bug fix)

**Development Team**: ViroForge Development Team
**Assistant**: Claude Code (Anthropic)

**Major Milestones**:
- 2025-11-09: Taxonomy bug discovered and fixed (Phase 7-8)
- Earlier: Phases 1-6 completed (core functionality)
