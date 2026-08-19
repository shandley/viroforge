# ViroForge - Claude Code Context

## Project Overview

**ViroForge** is a synthetic viral metagenomics dataset generator for benchmarking virome analysis pipelines. It generates realistic FASTQ sequencing data with known ground truth for validation.

**Key Features**:
- 20 curated viral genome collections (body sites, disease states, environmental)
- Realistic abundance distributions (log-normal, power-law)
- Multi-step library preparation modeling (VLP extraction, amplification bias, sequencing)
- RNA and DNA virome workflows
- 14,423 viral genomes from RefSeq with ICTV taxonomy
- Real reference contamination (rRNA, host DNA, PhiX, Illumina adapters)
- Sequencing artifact injection (adapters, low-complexity, PCR duplicates, ERVs)
- Per-read ground truth labels for exact QC validation metrics

**Target Users**: Bioinformatics researchers developing/benchmarking virome analysis tools

---

## Current Status (v0.20.0)

### Phases Completed

- ✅ **Phase 1-5**: Core infrastructure, 8 original collections, FASTQ generation
- ✅ **Phase 6**: Amplification bias modeling (MDA, SISPA, TruSeq)
- ✅ **Phase 7**: Critical collections (wastewater, IBD, HIV+, CF respiratory)
- ✅ **Phase 8**: RNA virome workflow (RT, rRNA depletion, degradation)
- ✅ **Phase 9**: Additional host niches (vaginal, blood, ocular, lung, urinary)
- ✅ **Phase 10**: Long-read sequencing (PacBio HiFi, Oxford Nanopore)
- ✅ **Phase 11**: Hybrid assembly support (matched short + long reads)
- ✅ **Phase 12**: CLI enhancements & web interface
- ✅ **Taxonomy Bug Fix**: Enhanced fuzzy matching, fixed 469 genomes
- ✅ **QC Validation Toolkit** (2026-03-22):
  - Real reference contamination (rRNA, host DNA, PhiX, adapters)
  - Per-read source labels in FASTQ headers
  - Adapter read-through injection with normal insert size distribution
  - Low-complexity artifact injection with controlled entropy spectrum
  - PCR duplicate injection with geometric copy distribution
  - Retroviral read injection (endogenous HERV + exogenous Retroviridae)
  - PhiX174 curation fix (removed lab phage from 8 collections)
  - Insert-size-driven adapter contamination (emergent rate from insert size distribution)
  - GC/length-biased PCR duplicate template selection
  - Chimeric adapter support (internal adapter from ligation events)
  - ERV injection integrated into contamination profile factory
  - virome-qc profile batch configs (4 profiles, 5 datasets each)
  - `viroforge summary` command for expected range derivation

### Project Status: Production Ready

**20 curated collections** covering:
- Host-associated viromes (gut, oral, skin, respiratory, vaginal, blood, etc.)
- Environmental viromes (marine, soil, freshwater, wastewater)
- Disease states (IBD, HIV+, CF)
- RNA viromes (respiratory, arbovirus, fecal)

**5 sequencing platforms**: NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore

**Complete CLI**: browse, generate, batch, report, compare, presets, web, setup-db,
summary, benchmark

### Benchmarking Framework (Phase 13) - Modules 1/2/4 built

**Status**: Phase 13A metadata complete; Modules 1 (QC), 2 (Assembly), 4 (Taxonomy)
implemented in `viroforge/benchmarking/` with the `viroforge benchmark` command.

Implemented (JSON + markdown reports, independent-oracle tests):
- **Module 1 QC** (`benchmark qc`): contamination removal + viral retention
  against per-read `source=` labels; read-name match-rate gate; dedup scored
  separately. See `docs/PHASE13B_QC_BENCHMARK.md`.
- **Module 2 Assembly** (`benchmark assembly`): genome recovery/completeness/
  identity, chimeras, N50/L50, observed-vs-expected completeness, abundance
  accuracy. minimap2 via mappy (`benchmarking/align.py`, the `benchmark` extra).
  See `docs/PHASE13B_ASSEMBLY_BENCHMARK.md`.
- **Module 4 Taxonomy** (`benchmark taxonomy`): read- and contig-based; taxid-exact
  + genus/family (NCBI taxdump via `--taxdump-dir`); abundance profile;
  known/dark-matter stratification; Kraken2/Centrifuge/DIAMOND/MMseqs2/generic
  formats. Metadata exports `benchmarking.taxonomy`. See `docs/PHASE13C_TAXONOMY_BENCHMARK.md`.

Data quality was independently evaluated before QC (`docs/DATA_QUALITY_EVALUATION.md`,
`scripts/evaluate_dataset.py`): read-level ground truth is accurate and reproducible.

**Next**: Module 5 (completeness across coverage), HTML reports + visualizations,
length-weighted contig taxonomy metrics, remaining modules (binning, annotation,
host, discovery, end-to-end). See `ROADMAP.md` and `docs/PHASE13_BENCHMARKING_FRAMEWORK.md`.

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

**taxonomy**: ICTV/NCBI taxonomy (71.7% family, 88.9% class assigned after the
2026-07-16 NCBI rank enrichment)
- genome_id → realm, kingdom, phylum, class, order_name, **family**, genus, species
- **NOTE**: `family` has NOT NULL constraint - use 'Unknown' for unassigned
- **IMPORTANT**: After taxonomy fix, check family != 'Unknown' for critical viruses

**body_site_collections**: 20 curated collections
- collection_id, collection_name, description, n_genomes
- literature_references (verify citations!)

**collection_genomes**: Genome-collection associations
- collection_id, genome_id, relative_abundance, prevalence, abundance_rank

---

## Collections Overview (20 Total)

Collections use contiguous IDs 1-20 (renumbered from the legacy 9-28 scheme in
PR #6; see scripts/migrate_renumber_collections.py). There are no separate
"VLP comparison" collections: VLP enrichment is applied to any collection via
`--vlp-protocol`. Genome counts below reflect the live database after the
animal-virus (#41) and phage-host (#50) cleanups. Verify with
`viroforge browse` or a query against `body_site_collections`.

### Body-site collections (1-8)
1. Gut Virome - Adult Healthy (Western Diet) (134 genomes)
2. Oral Virome - Saliva (Healthy) (47 genomes)
3. Skin Virome - Sebaceous Sites (Healthy) (15 genomes)
4. Respiratory Virome - Nasopharynx (Healthy) (43 genomes)
5. Marine Virome - Coastal Surface Water (445 genomes)
6. Soil Virome - Agricultural (290 genomes)
7. Freshwater Virome - Lake Surface Water (200 genomes)
8. Mouse Gut Virome - Laboratory (C57BL/6) (22 genomes)

### Disease and environmental collections (9-12)
9. Wastewater Virome - Urban Treatment Plant (352 genomes)
10. IBD Gut Virome (Inflammatory Bowel Disease) (90 genomes)
11. HIV+ Gut Virome (55 genomes) - includes human herpesviruses
12. Cystic Fibrosis (CF) Respiratory Virome (79 genomes)

### RNA virome collections (13-15)
13. Human Respiratory RNA Virome (53 genomes)
14. Arbovirus Environmental (Mosquito Virome) (39 genomes)
15. Fecal RNA Virome (46 genomes) - includes rotavirus/norovirus

### Additional host niches (16-20)
16. Vaginal Virome (Healthy) (26 genomes) - proxy Lactobacillus phages (see notes)
17. Blood/Plasma Virome (Healthy) (21 genomes)
18. Ocular Surface Virome (Healthy) (17 genomes)
19. Lower Respiratory (Lung) Virome (Healthy) (31 genomes)
20. Urinary Virome (Healthy) (20 genomes)

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
python scripts/generate_fastq_dataset.py --collection-id 13 --dry-run

# Generate full dataset
python scripts/generate_fastq_dataset.py \
  --collection-id 13 \
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
- `scripts/curate_*_collection.py` - Collection curation scripts (20 total)

### Documentation
- `docs/TAXONOMY_BUG_FIX.md` - **READ FIRST** if revisiting taxonomy issues
- `docs/COLLECTION_IMPLEMENTATION_GUIDE.md` - All 20 collections documented
- `docs/PHASE4_FASTQ_GENERATION.md` - FASTQ generation guide
- `ROADMAP.md` - Development roadmap (current: Phase 8)

### Database
- `viroforge/data/viral_genomes.db` - Main database (14,423 genomes)
- `data/ictv/ictv_taxonomy.json` - ICTV reference (17,925 records)

---

## Development Guidelines

### Releases and reproducibility

ViroForge's value is reproducible ground truth, so any change that alters
generator output for a fixed seed is treated as breaking even when the API is
untouched. Bump the minor version and record it under "Reproducibility" in
`CHANGELOG.md`. `viroforge/__init__.py` holds `__version__`; `setup.py` reads it
and `CITATION.cff` must be updated to match. There are no git tags and nothing
is published to PyPI, so versions are repo state only.

The lab-notebook workflow is retired. There is no longer a pre-commit hook
requiring an entry, and the `.claude/hooks/` scripts have been removed.

### Citations
**CRITICAL**: Run `/verify-references` before committing anything containing a
PMID, DOI or accession. Never write an identifier from memory.

Two separate failures to guard against, because the first check does not catch
the second:
1. The identifier resolves to a different paper. Both PMIDs submitted in PRs
   during 2026-07 were wrong (one was TreeDyn cited as a phi29 study, one was
   Fast UniFrac cited as a VLP protocol paper).
2. The identifier is right but the paper does not contain the number attributed
   to it. PR #65 justified a parameter change with "~2-10x bias" that appears
   nowhere in the cited work, whose actual conclusion pointed the other way.

So verify the accession, then confirm the paper really says what you are citing
it for. Include DOI and PMID, use full journal names, check publication years.
See `docs/CITATION_CORRECTIONS.md` for past errors.

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
- **4,080 genomes have `family='Unknown'` (28.3%)** after the 2026-07-16 NCBI
  enrichment (was 6,182 / 42.8%)
  - Mostly: modern bacteriophages ICTV assigns to a genus but no family (families
    abolished 2021) - correct taxonomy, not a gap. Only 1,481 have no class at all.
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

## Next Steps

**Open on GitHub**:
- Issue #37: bacterial/fungal background, so `--no-vlp` bulk metagenomes are
  realistic. Today they come out ~79% viral where a real bulk stool metagenome
  is 1-5%. Build the reference FASTA on demand through `resolver.py` and verify
  every accession; PR #39 stalled on 4.4M lines of FASTA committed to git.
- Issue #30: replace the batch YAML text box with a visual form builder.
- PR #51: `host_associations` table, deferred pending a real `body_site` column.

**Known modelling gaps**:
- MDA and RdAB produce nearly identical GC efficiency curves at their defaults.
  The docstring claim that MDA bias is 2-5x stronger than PCR no longer holds
  and needs its own measurement.
- The pore-less and column branches of VLP enrichment favour larger particles,
  opposite to the filtration curves. Pelleting and column capture are not
  sieving so this may be right, but it has not been checked against literature.
- `--contamination-level failed` exists in `LEVEL_MULTIPLIERS` but neither CLI
  parser accepts it.
- Amplification (MDA/RdAB) adjusts relative abundances between genomes but
  does not model intra-genome coverage unevenness. Read Simulation samples
  positions uniformly across the full reference. Real MDA produces uneven
  coverage (over-represented near primer binding sites, gaps elsewhere),
  which breaks contigs during assembly. This makes ViroForge's assembly
  benchmarks optimistic relative to real amplified data.

**Benchmarking framework**: Module 5 (completeness across coverage), HTML reports
and visualizations, then Modules 3 and 6-9. See `ROADMAP.md`.

---

## Contact & History

**Project**: ViroForge
**Repository**: shandley/viroforge
**Current Version**: 0.20.0
**Last Major Update**: 2026-07-30 (microbial background; issue #37 closed)

**Development Team**: ViroForge Development Team
**Assistant**: Claude Code (Anthropic)

**Major Milestones**:
- 2026-07-30: Bacterial/fungal/archaeal background, issue #37 closed; amplification
  GC bias was compounding over cycles (fixed); every preset pointed at the wrong
  collection (fixed); v0.17.0-v0.20.0
- 2026-07-29: MDA GC-bias optimum moved 40% -> 50% GC; v0.16.0
- 2026-07-28: PR backlog cleared (10 merged); VLP filtration direction corrected;
  CHANGELOG started; v0.15.0
- 2026-07-16: Benchmarking framework (QC, assembly, taxonomy); v0.14.0
- 2026-07-14: Collections renumbered to 1-20; v0.13.0 realistic default output
- 2025-11-09: Taxonomy bug discovered and fixed (Phase 7-8)
- Earlier: Phases 1-6 completed (core functionality)
