# ViroForge Database Setup Guide

**Complete guide for building the ViroForge viral genome database from scratch**

---

## Overview

The ViroForge database contains 14,568 RefSeq viral genomes with rich metadata including ICTV taxonomy, genome characteristics, and curated body site collections. The database is not included in the repository due to its size (~2 GB) and must be built locally.

## Requirements

### System Requirements
- **Disk Space:** ~2-2.5 GB
  - RefSeq downloads (compressed): ~198 MB
  - Parsed genome data: ~300 MB
  - SQLite database: ~1.7 GB
  - ICTV VMR file: ~3.5 MB
- **RAM:** 2 GB minimum (4 GB recommended)
- **Time:** ~3 hours for full database with ICTV taxonomy
- **Internet:** Stable connection for NCBI downloads and ICTV VMR

### Software Requirements
- Python 3.9+
- Biopython
- NumPy
- Pandas
- openpyxl (for ICTV Excel parsing)

```bash
pip install biopython numpy pandas openpyxl
```

---

## Quick Setup (Recommended)

### Step 1: Create Database Schema

Creates the SQLite database with tables for genomes, taxonomy, collections, etc.

```bash
python viroforge/data/database_schema.py viroforge/data/viral_genomes.db
```

**Output:**
```
Creating database: viroforge/data/viral_genomes.db
✓ Schema verification passed
Database statistics:
  total_genomes: 0
  schema_version: 1.0.0
```

**Time:** ~5 seconds

---

### Step 2: Download RefSeq Viral Genomes

Downloads complete viral genomes from NCBI RefSeq FTP server.

```bash
python scripts/download_refseq.py --output data/refseq --all
```

**What it does:**
1. Downloads assembly summary from NCBI (~15,000 genome records)
2. Filters for complete genomes only
3. Downloads each genome's FASTA file (gzipped)
4. Progress updates every 10 genomes
5. Automatic retry on failures (3 attempts)
6. Resume capability (skips existing files)

**Output:**
```
Step 1: Downloading assembly summary...
✓ Downloaded assembly summary: data/refseq/assembly_summary.txt

Step 2: Parsing assembly summary...
✓ Parsed 14568 genomes matching criteria

Step 3: Downloading genome sequences...
Progress: 100/14568 (0.7%) | Success: 100 | Failed: 0 | Skipped: 0
Progress: 200/14568 (1.4%) | Success: 200 | Failed: 0 | Skipped: 0
...
```

**Time:** ~90-120 minutes (depends on internet speed)

**Troubleshooting:**
- **Slow download?** Normal. NCBI throttles requests. Script downloads 2-3 genomes/second.
- **Download failed?** Script auto-retries 3 times. If persistent failures, check internet connection.
- **Want to resume?** Just re-run the command. It skips existing files automatically.

---

### Step 3: Parse Downloaded Genomes

Extracts sequences and calculates genome statistics (length, GC content, genome type).

```bash
python scripts/parse_genomes.py --input data/refseq --output data/parsed
```

**What it does:**
1. Reads gzipped FASTA files
2. Extracts sequences and metadata
3. Calculates GC content, length
4. Detects genome type (dsDNA, ssDNA, dsRNA, ssRNA)
5. Saves JSON metadata + FASTA sequences

**Output:**
```
Parsing genomes from: data/refseq
✓ Parsed 14568/14568 genomes (100.0%)
✓ Saved to: data/parsed/
```

**Time:** ~5 minutes

---

### Step 4: Populate Database

Inserts genomes into database with quality filtering.

```bash
python scripts/populate_database.py \
    --input data/parsed \
    --database viroforge/data/viral_genomes.db
```

**What it does:**
1. Loads parsed genome metadata
2. Applies quality filters (see below)
3. Inserts genomes in batches (100 at a time)
4. Updates database metadata

**Quality Filters:**
- Minimum length: 1,000 bp
- Maximum length: 500,000 bp
- GC content: 15-75%
- Max ambiguous bases: 5%

**Output:**
```
Loading genomes from: data/parsed
✓ Loaded 14568 genomes

Applying quality filters...
✓ Passed: 14423 genomes (99.0%)
✗ Failed: 145 genomes (1.0%)

Populating database...
Progress: 100/14423 | Success: 100
Progress: 200/14423 | Success: 200
...
✓ Inserted 14423 genomes successfully
```

**Time:** ~3 seconds

---

### Step 5: Download and Integrate ICTV Taxonomy

Downloads the official ICTV Virus Metadata Resource (VMR) and maps taxonomy to genomes.

```bash
# Download ICTV VMR (3.5 MB)
mkdir -p data/ictv
curl -L -o data/ictv/VMR_MSL40_v2.xlsx \
  "https://ictv.global/sites/default/files/VMR/VMR_MSL40.v2.20251013.xlsx"

# Parse and integrate taxonomy
python scripts/parse_ictv_taxonomy.py \
  --vmr data/ictv/VMR_MSL40_v2.xlsx \
  --output data/ictv \
  --database viroforge/data/viral_genomes.db \
  --update-db
```

**What it does:**
1. Downloads ICTV Master Species List 40 (17,925 virus records)
2. Parses ICTV taxonomy hierarchy (realm → species)
3. Maps ICTV species names to RefSeq genomes (~69% match rate)
4. Updates database taxonomy table with proper viral families

**Output:**
```
Parsed 17925 virus taxonomy records
Unique families: 368
Mapping Statistics:
  Total genomes: 14,423
  Matched: 9,962 (69.1%)
  With ICTV family: 7,773 (53.9%)
✓ Database taxonomy updated
```

**Time:** ~1 minute

---

### Step 6: Curate Body Site Collections

Creates 8 literature-validated virome collections for different body sites.

```bash
python scripts/curate_body_site_collections.py --all
```

**What it does:**
1. Selects genomes for each body site based on taxonomic composition
2. Assigns realistic abundance distributions (log-normal, power-law)
3. Creates collection metadata and composition tables
4. Saves collections to database

**Collections Created:**
- Gut Virome - 134 genomes
- Oral Virome - 47 genomes
- Skin Virome - 15 genomes
- Respiratory Virome - 41 genomes
- Marine Virome - 448 genomes (90% of target!)
- Soil Virome - 291 genomes (97% of target!)
- Freshwater Virome - 200 genomes (100% of target!)
- Mouse Gut Virome - 22 genomes

**Output:**
```
✓ Collection curated: 134 genomes (Gut)
✓ Collection curated: 47 genomes (Oral)
✓ Collection curated: 448 genomes (Marine)
...
Curation complete!
```

**Time:** ~2 seconds

---

### Step 7: Verify Setup

Check that collections are available:

```bash
python scripts/generate_fastq_dataset.py --list-collections
```

**Expected Output:**
```
Available Collections:
================================================================================
ID: 9
  Name: Gut Virome - Adult Healthy (Western Diet)
  Genomes: 134
  Description: Human gut virome from healthy adults...

ID: 10
  Name: Oral Virome - Saliva (Healthy)
  Genomes: 47
  Description: Oral virome from healthy saliva...

[... 6 more collections ...]
```

If you see 8 collections with genome counts, **setup is complete!** 🎉

---

## Advanced Options

### Quick Test Setup (100 genomes)

For testing or exploring ViroForge without waiting 2 hours:

```bash
# Step 1: Create schema
python viroforge/data/database_schema.py viroforge/data/viral_genomes.db

# Step 2: Download only 100 genomes (~2 minutes)
python scripts/download_refseq.py --output data/refseq --limit 100

# Step 3: Parse
python scripts/parse_genomes.py --input data/refseq --output data/parsed

# Step 4: Populate
python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db
```

**Time:** ~5 minutes total

**Note:** This won't include all body site collections (not enough diversity), but you can test FASTQ generation with custom genome sets.

---

### Resuming Failed Downloads

If download was interrupted:

```bash
# Just re-run - automatically skips existing files
python scripts/download_refseq.py --output data/refseq --all
```

Progress will show:
```
Progress: 500/14568 | Success: 0 | Failed: 0 | Skipped: 500
```

---

### Custom Quality Filters

Adjust filters for specific needs:

```bash
python scripts/populate_database.py \
    --input data/parsed \
    --database viroforge/data/viral_genomes.db \
    --min-length 5000 \      # Longer genomes only
    --max-length 300000 \    # Exclude large genomes
    --min-gc 0.20 \          # More stringent GC
    --max-gc 0.70
```

---

## Troubleshooting

### Database file not found

**Error:**
```
FileNotFoundError: Database not found: viroforge/data/viral_genomes.db
```

**Solution:** Run Step 1 to create the database schema first.

---

### No collections found

**Error:** `--list-collections` shows empty output

**Cause:** Database doesn't have collections populated yet

**Solution:** Collections are automatically included when you populate the database in Step 4. If they're missing:

```bash
# Check database stats
python scripts/viroforge_db.py --database viroforge/data/viral_genomes.db --stats
```

If `total_collections: 0`, the database wasn't fully populated. Re-run Step 4.

---

### Download failures

**Error:**
```
✗ Failed to download: GCF_XXXXXX (attempt 3/3)
```

**Possible causes:**
- Transient network issues (script auto-retries)
- NCBI server temporarily unavailable
- Genome removed from RefSeq

**Solution:**
1. Wait a few minutes and re-run the download command
2. Script will skip successful downloads and retry failures
3. A few failures (<1%) are normal and won't affect functionality

---

### Insufficient disk space

**Error:**
```
OSError: [Errno 28] No space left on device
```

**Solution:**
1. Ensure you have 2.5 GB free before starting
2. You can delete `data/refseq/` and `data/parsed/` after Step 4 completes (saves ~550 MB)
3. The database file (`viroforge/data/viral_genomes.db`) is the only required file

**To clean up after setup:**
```bash
# Keep only the database, remove intermediate files
rm -rf data/refseq/
rm -rf data/parsed/

# Database is still fully functional
python scripts/generate_fastq_dataset.py --list-collections
```

---

### Parse failures

**Error:**
```
ERROR - Failed to parse: GCF_XXXXXX.1
```

**Causes:**
- Corrupted download
- Invalid FASTA format
- Unexpected file structure

**Solution:**
```bash
# Delete the specific genome file
rm data/refseq/GCF_XXXXXX.1_genomic.fna.gz

# Re-download that genome
python scripts/download_refseq.py --output data/refseq --all

# Re-run parse
python scripts/parse_genomes.py --input data/refseq --output data/parsed
```

---

## Database Statistics

After setup completes, check statistics:

```bash
python scripts/viroforge_qc.py --database viroforge/data/viral_genomes.db
```

**Expected output:**
```
ViroForge Database Quality Report
================================================================================

Overall Statistics
------------------
Total genomes: 14,423
Total families: 145
Total collections: 8
Schema version: 1.0.0

Genome Characteristics
----------------------
Length range: 1,028 - 499,997 bp
Mean length: 52,431 bp
GC content: 15.0% - 74.8%
Mean GC: 45.3%

Genome Types
------------
dsDNA: 13,891 (96.3%)
ssDNA: 387 (2.7%)
dsRNA: 98 (0.7%)
ssRNA: 47 (0.3%)

ICTV Taxonomy Coverage
----------------------
With ICTV taxonomy: 10,756 (74.6%)
Missing taxonomy: 3,667 (25.4%)

Body Site Collections
---------------------
Gut Virome: 134 genomes
Oral Virome: 47 genomes
Skin Virome: 15 genomes
Respiratory Virome: 41 genomes
Marine Virome: 448 genomes
Soil Virome: 291 genomes
Freshwater Virome: 200 genomes
Mouse Gut Virome: 22 genomes
```

---

## Database Schema

The SQLite database contains 8 tables:

1. **genomes** - Core genome data (14,423 rows)
   - genome_id, accession, species, length, gc_content, sequence, etc.

2. **taxonomy** - ICTV taxonomy (10,756 genomes with taxonomy)
   - realm, kingdom, phylum, class, order, family, genus, species

3. **host_associations** - Virus-host relationships
   - genome_id, host_name, host_taxid

4. **ecological_metadata** - Environmental context
   - genome_id, environment, body_site, sample_type

5. **genome_annotations** - Gene-level data
   - genome_id, gene_name, product, sequence

6. **body_site_collections** - Pre-curated virome collections (8 collections)
   - collection_id, collection_name, n_genomes, description

7. **collection_genomes** - Collection membership
   - collection_id, genome_id, abundance, rank

8. **database_metadata** - Version tracking
   - version, created_date, last_updated

**Explore schema:**
```bash
sqlite3 viroforge/data/viral_genomes.db .schema
```

---

## Next Steps

Once database setup is complete:

1. **Generate your first dataset:**
   ```bash
   python scripts/generate_fastq_dataset.py \
       --collection-id 9 \
       --output test_gut_virome \
       --coverage 5 \
       --platform novaseq
   ```

2. **Explore the database:**
   ```bash
   # Interactive database browser
   python scripts/viroforge_db.py --database viroforge/data/viral_genomes.db

   # Search for specific viruses
   python scripts/viroforge_search.py --database viroforge/data/viral_genomes.db --family Siphoviridae

   # Compare collections
   python scripts/viroforge_compare.py --database viroforge/data/viral_genomes.db --collections 9 10
   ```

3. **Read the documentation:**
   - [FASTQ Generation Guide](PHASE4_FASTQ_GENERATION.md)
   - [Collection Implementation Guide](COLLECTION_IMPLEMENTATION_GUIDE.md)
   - [User Guide](USER_GUIDE.md)

---

## Frequently Asked Questions

### Do I need to rebuild the database?

**No.** Once built, the database is permanent. Only rebuild if:
- You want to update to newer RefSeq releases
- Database file is corrupted
- You want different quality filters

### Can I share the database?

**Yes!** The database is portable. You can:
- Copy `viroforge/data/viral_genomes.db` to another machine
- Share it with collaborators
- Back it up for archival

**Note:** It's 1.7 GB, so use cloud storage or external drives.

### How do I update to new RefSeq releases?

```bash
# Delete old downloads
rm -rf data/refseq/ data/parsed/

# Re-run setup pipeline
python scripts/download_refseq.py --output data/refseq --all
python scripts/parse_genomes.py --input data/refseq --output data/parsed
python scripts/populate_database.py --input data/parsed --database viroforge/data/viral_genomes.db --force
```

RefSeq is updated monthly, but ViroForge doesn't require the absolute latest version.

### Can I add custom genomes?

**Yes!** See `scripts/README_HELPER_UTILITIES.md` for instructions on adding custom genomes to existing collections or creating new collections.

---

## Support

If you encounter issues not covered here:

1. Check [GitHub Issues](https://github.com/shandley/viroforge/issues)
2. Search documentation in `docs/`
3. Open a new issue with:
   - Setup step where error occurred
   - Full error message
   - Output of `python --version` and `pip list`

---

**Last Updated:** 2025-11-06
**Version:** 1.0
