# ViroForge Database Exploration Tools

**Created**: 2025-11-01
**Phase**: Phase 3, Week 6

## Overview

ViroForge includes three command-line tools for exploring and analyzing the viral genome database:

| Tool | Purpose | Key Features |
|------|---------|--------------|
| `viroforge_db.py` | Database inspector | Stats, search, taxonomy browsing, export |
| `viroforge_compare.py` | Collection comparison | Compare taxonomic composition, diversity metrics |
| `viroforge_search.py` | Advanced search | Complex queries, batch export (planned) |

## 1. Database Inspector (`viroforge_db.py`)

### Overview

Interactive CLI tool for exploring the ViroForge database. Provides quick access to statistics, search functionality, taxonomy browsing, and genome export.

### Commands

#### `stats` - Database Statistics

Get comprehensive database statistics including genome counts, taxonomy coverage, and collection summaries.

**Usage**:
```bash
python scripts/viroforge_db.py stats
python scripts/viroforge_db.py stats --verbose  # More detail
```

**Output Example**:
```
======================================================================
ViroForge Database Statistics
======================================================================

📊 Genome Statistics:
  Total genomes:      13,847
  Mean length:        35,915 bp
  Length range:       1,023 - 2,473,870 bp
  Mean GC content:    45.1%

  Genome types:
    dsDNA             11,982 (86.5%)
    ssRNA              1,752 (12.7%)
    dsRNA                 98 ( 0.7%)
    ssDNA                 15 ( 0.1%)

🧬 Taxonomy Coverage:
  With ICTV realm:    10,234 (73.9%)
  With species:       13,847 (100.0%)

  Top 5 Realms:
    Riboviria                         7,123 (51.4%)
    Monodnaviria                      1,543 (11.1%)
    Varidnaviria                        891 ( 6.4%)
    ...

📚 Collections:
  Total collections:  8

  Available collections:
    Gut Virome - Adult Healthy (Western Diet)          500 genomes
    Oral Virome - Saliva (Healthy)                     200 genomes
    ...
```

---

#### `search` - Search Genomes

Search for genomes using multiple filters. Supports taxonomy, host, genome characteristics, and more.

**Usage**:
```bash
# Search by taxonomy
python scripts/viroforge_db.py search --family Siphoviridae

# Search by host
python scripts/viroforge_db.py search --host "Escherichia coli"

# Complex search
python scripts/viroforge_db.py search \
  --family Crassvirales \
  --length-min 100000 \
  --gc-min 0.4 \
  --limit 20

# Search with different output formats
python scripts/viroforge_db.py search --family Potyviridae --format tsv
python scripts/viroforge_db.py search --realm Riboviria --format json
```

**Search Filters**:
| Filter | Description | Example |
|--------|-------------|---------|
| `--realm` | ICTV realm | `--realm Riboviria` |
| `--family` | ICTV family | `--family Siphoviridae` |
| `--genus` | ICTV genus | `--genus Orthopoxvirus` |
| `--species` | ICTV species or genome name | `--species "cowpox virus"` |
| `--host` | Host organism | `--host Streptococcus` |
| `--genome-type` | Genome type | `--genome-type dsDNA` |
| `--length-min` | Minimum length (bp) | `--length-min 50000` |
| `--length-max` | Maximum length (bp) | `--length-max 150000` |
| `--gc-min` | Minimum GC content | `--gc-min 0.4` |
| `--gc-max` | Maximum GC content | `--gc-max 0.6` |
| `--limit` | Max results (default: 100) | `--limit 500` |

**Output Formats**:
- `table` (default) - Pretty-printed table
- `tsv` - Tab-separated values
- `json` - JSON format

**Output Example** (table format):
```
Found 56 genomes

Genome ID          Species                                  Type         Length Family
--------------------------------------------------------------------------------------------------------------
GCF_000851785.1    Yam mosaic virus                         dsDNA         9,608 Potyviridae
GCF_000854025.1    Sugarcane mosaic virus                   dsDNA         9,596 Potyviridae
...
```

---

#### `collection` - Collection Details

View detailed information about a specific body site collection.

**Usage**:
```bash
# Full collection ID
python scripts/viroforge_db.py collection gut_virome_adult_healthy_western

# Partial match supported
python scripts/viroforge_db.py collection gut

# Show more genomes
python scripts/viroforge_db.py collection gut --top 50
```

**Output Example**:
```
======================================================================
Collection: Gut Virome - Adult Healthy (Western Diet)
======================================================================

ID:               gut_virome_adult_healthy_western
Body Site:        gut
Sample Type:      stool
Health Status:    healthy
Genome Count:     500
Phage Fraction:   96.0%
Curation Date:    2025-11-01

Description:
  Human gut virome from healthy adults with Western diet.
  crAssphage-dominated composition typical of industrialized populations.

Top 20 Most Abundant Genomes:
Rank   Abundance  Genome ID          Species                                  Family
------------------------------------------------------------------------------------------------------
1         28.45%  NC_024711          Bacteroides phage crAssphage            Crassvirales
2         15.67%  NC_019724          Escherichia phage T7                    Podoviridae
...

Taxonomic Composition (Top 10 Families):
  Crassvirales                                    150 (30.0%)
  Siphoviridae                                    140 (28.0%)
  Myoviridae                                       70 (14.0%)
  ...
```

---

#### `taxonomy` - Browse Taxonomy

Browse ICTV taxonomy at any rank, showing genome counts.

**Usage**:
```bash
# Browse families (default)
python scripts/viroforge_db.py taxonomy

# Browse other ranks
python scripts/viroforge_db.py taxonomy --rank realm
python scripts/viroforge_db.py taxonomy --rank genus --limit 100
```

**Available Ranks**:
- `realm`
- `kingdom`
- `phylum`
- `class`
- `order`
- `family` (default)
- `genus`
- `species`

**Output Example**:
```
Taxonomy Browser: Family
======================================================================

Rank Family                                             Count
----------------------------------------------------------------------
1    Potyviridae                                          834
2    Flaviviridae                                         721
3    Geminiviridae                                        658
4    Siphoviridae                                         542
...

Showing top 50 of 234 familys
```

---

#### `export` - Export Genomes

Export genomes matching search criteria to FASTA, TSV, or JSON format.

**Usage**:
```bash
# Export to FASTA
python scripts/viroforge_db.py export \
  --family Crassvirales \
  --output crassvirales.fasta

# Export to TSV
python scripts/viroforge_db.py export \
  --realm Riboviria \
  --genome-type ssRNA \
  --format tsv \
  --output rna_viruses.tsv

# Export large set
python scripts/viroforge_db.py export \
  --family Siphoviridae \
  --length-min 30000 \
  --limit 1000 \
  --output siphoviridae_large.fasta
```

**Export Formats**:
- `fasta` (default) - FASTA format with sequences
- `tsv` - Tab-separated metadata (no sequences)
- `json` - JSON format with metadata

**Output**:
```
Exporting genomes to crassvirales.fasta...
✓ Exported 237 genomes to crassvirales.fasta
```

---

### Common Workflows

#### Workflow 1: Explore Database

```bash
# 1. Get overview
python scripts/viroforge_db.py stats

# 2. Browse taxonomy
python scripts/viroforge_db.py taxonomy --rank family

# 3. Search for genomes of interest
python scripts/viroforge_db.py search --family Siphoviridae --host Streptococcus

# 4. View collection
python scripts/viroforge_db.py collection oral
```

#### Workflow 2: Find and Export Phages

```bash
# 1. Search for large Siphoviridae phages
python scripts/viroforge_db.py search \
  --family Siphoviridae \
  --length-min 40000 \
  --format tsv > candidates.tsv

# 2. Review candidates
cat candidates.tsv

# 3. Export sequences
python scripts/viroforge_db.py export \
  --family Siphoviridae \
  --length-min 40000 \
  --output large_siphoviridae.fasta
```

#### Workflow 3: Collection Analysis

```bash
# 1. List all collections
python scripts/viroforge_db.py stats | grep -A 20 "Collections:"

# 2. View specific collection
python scripts/viroforge_db.py collection gut --top 50

# 3. Export collection genomes
# (Use collection_genomes table directly or build custom query)
```

---

## 2. Collection Comparison Tool (`viroforge_compare.py`)

### Overview

Compare taxonomic composition, abundance distributions, and diversity metrics between body site collections.

### Usage

```bash
# Compare two collections
python scripts/viroforge_compare.py gut oral

# Compare multiple collections
python scripts/viroforge_compare.py gut oral skin respiratory

# Save to file
python scripts/viroforge_compare.py gut oral --report comparison_report.txt
```

### Output

**Collection Overview**:
```
Collection Overview:
--------------------------------------------------------------------------------
Collection                               Genomes  Body Site
--------------------------------------------------------------------------------
Gut Virome - Adult Healthy (Western D        500  gut
Oral Virome - Saliva (Healthy)               200  oral_cavity
```

**Diversity Metrics**:
```
Diversity Metrics:
--------------------------------------------------------------------------------
Collection                                  Shannon    Simpson  Evenness
--------------------------------------------------------------------------------
Gut Virome - Adult Healthy (Western D         5.23     0.9856     0.841
Oral Virome - Saliva (Healthy)                 4.67     0.9782     0.880
```

**Genome Characteristics**:
```
Genome Characteristics:
--------------------------------------------------------------------------------
Collection                                Mean Length   Mean GC
--------------------------------------------------------------------------------
Gut Virome - Adult Healthy (Western D   45,234 bp      43.2%
Oral Virome - Saliva (Healthy)          38,921 bp      45.8%
```

**Taxonomic Composition** (Top 10 families per collection)

**Genome Overlap** (for 2-collection comparisons):
```
Genome Overlap:
--------------------------------------------------------------------------------

Gut Virome vs Oral Virome:
  Shared genomes:             23
  Unique to first:           477
  Unique to second:          177
  Jaccard index:           0.033
```

### Use Cases

#### Compare Body Sites

```bash
# Compare human body sites
python scripts/viroforge_compare.py gut oral skin respiratory

# Compare aquatic environments
python scripts/viroforge_compare.py marine freshwater

# Compare human vs mouse gut
python scripts/viroforge_compare.py gut mouse_gut
```

#### Assess Collection Similarity

```bash
# Calculate Jaccard index for two collections
python scripts/viroforge_compare.py gut oral

# Look for shared genomes between sites
```

#### Validate Curation

```bash
# Check diversity metrics match expectations
python scripts/viroforge_compare.py gut --report gut_validation.txt

# Compare taxonomic composition to literature
```

---

## 3. Advanced Search Tool (`viroforge_search.py`)

**Status**: Placeholder created, implementation planned

**Planned Features**:
- Complex boolean queries
- Batch export with metadata
- Statistical summaries
- Comparison between search results

---

## Installation

These tools are included with ViroForge and require no additional installation beyond the main ViroForge dependencies.

```bash
# From ViroForge root directory
source venv/bin/activate

# Make scripts executable (optional)
chmod +x scripts/viroforge_db.py
chmod +x scripts/viroforge_compare.py
```

---

## Requirements

- Python 3.8+
- SQLite3
- NumPy
- ViroForge database (`viroforge/data/viral_genomes.db`)

All dependencies are installed with ViroForge:
```bash
pip install -e .
```

---

## Database Schema Compatibility

These tools work with both old and new database schemas:

**Old Schema** (Phase 1-2):
- `collection_name` column
- `n_genomes` column
- `genome_name` column in genomes table

**New Schema** (Phase 3+):
- `name` column
- `genome_count` column
- `species_name` column in genomes table
- Additional metadata: `body_site`, `sample_type`, `health_status`, `phage_fraction`

The tools automatically detect and adapt to the schema version.

---

## Tips and Best Practices

### Performance

- **Large exports**: Use `--limit` to control output size
- **Complex searches**: Be specific with filters to reduce results
- **Multiple queries**: Run searches in parallel for different criteria

### Search Strategy

1. **Start broad**: Begin with taxonomy searches at family level
2. **Refine**: Add genome characteristic filters (length, GC)
3. **Validate**: Check result counts before exporting
4. **Export**: Use appropriate format for downstream analysis

### Output Formats

- **`table`**: Human-readable, terminal viewing
- **`tsv`**: Spreadsheet import, scripting
- **`json`**: Programmatic access, web apps
- **`fasta`**: Sequence analysis, BLAST, alignment

### Common Patterns

```bash
# Find all phages for a specific host
python scripts/viroforge_db.py search --host "Escherichia coli" --format tsv

# Get large dsDNA viruses
python scripts/viroforge_db.py search --genome-type dsDNA --length-min 100000

# Export all Crassvirales
python scripts/viroforge_db.py export --family Crassvirales --output crass.fasta

# Compare related families
# (Use compare tool or multiple search commands)
```

---

## Troubleshooting

### Error: Database not found

```
❌ Error: Database not found at viroforge/data/viral_genomes.db
```

**Solution**: Check database path or specify with `--database` flag:
```bash
python scripts/viroforge_db.py --database path/to/viral_genomes.db stats
```

### Error: Collection not found

```
❌ Error: Collection 'xyz' not found
```

**Solution**: List available collections:
```bash
python scripts/viroforge_db.py stats | grep -A 10 "Collections:"
```

### Error: No such column

```
sqlite3.OperationalError: no such column: species_name
```

**Solution**: Database schema mismatch. Tools should auto-detect schema, but if error persists:
1. Check database was created with latest pipeline
2. Report issue with database version info

### No results returned

```
⚠️  No genomes match the search criteria.
```

**Causes**:
- Filters too restrictive
- Taxonomy not in database
- Host associations not populated

**Solutions**:
- Try broader search (e.g., family instead of genus)
- Check taxonomy coverage: `python scripts/viroforge_db.py stats`
- Verify host_associations table exists

---

## Examples Gallery

### Example 1: Find Streptococcus Phages for Oral Virome

```bash
# Search for candidates
python scripts/viroforge_db.py search \
  --host Streptococcus \
  --family Siphoviridae \
  --length-min 30000 \
  --length-max 60000 \
  --format tsv \
  > strep_phages.tsv

# Count results
wc -l strep_phages.tsv

# Export sequences
python scripts/viroforge_db.py export \
  --host Streptococcus \
  --family Siphoviridae \
  --length-min 30000 \
  --length-max 60000 \
  --output strep_phages.fasta
```

### Example 2: Compare All Human Body Site Collections

```bash
# Full comparison
python scripts/viroforge_compare.py \
  gut oral skin respiratory \
  --report human_body_sites_comparison.txt

# View report
cat human_body_sites_comparison.txt
```

### Example 3: Export All crAssphage Variants

```bash
# Find all Crassvirales
python scripts/viroforge_db.py taxonomy --rank order | grep -i crass

# Search
python scripts/viroforge_db.py search --order Crassvirales

# Export
python scripts/viroforge_db.py export \
  --order Crassvirales \
  --output crassvirales_all.fasta \
  --format fasta
```

### Example 4: Analyze Genome Size Distribution by Family

```bash
# Export family data
python scripts/viroforge_db.py search --family Siphoviridae --format tsv > sipho.tsv
python scripts/viroforge_db.py search --family Myoviridae --format tsv > myo.tsv
python scripts/viroforge_db.py search --family Podoviridae --format tsv > podo.tsv

# Analyze in R/Python
# (Calculate statistics, plot distributions)
```

---

## Future Enhancements

Planned features for future versions:

- **Visualization**: Generate plots (taxonomy pie charts, length distributions)
- **Batch operations**: Process multiple searches in single command
- **Advanced statistics**: Beta diversity, PCoA, clustering
- **Interactive mode**: REPL-style interface for exploratory analysis
- **Web interface**: Browser-based database explorer
- **Export to other formats**: GenBank, GFF, Newick trees

---

## Related Documentation

- **[Database Schema Reference](../docs/DATABASE_SCHEMA.md)** - Complete schema documentation (planned)
- **[Body Site Collections](../docs/BODY_SITE_COLLECTIONS.md)** - Collection specifications
- **[Curation Workflow](../docs/BODY_SITE_CURATION_WORKFLOW.md)** - How collections are curated
- **[API Reference](../docs/API.md)** - Python API for programmatic access

---

## Support

- **Issues**: Open a GitHub issue
- **Questions**: Check documentation or ask in discussions
- **Feature requests**: Submit via GitHub issues

---

**Document Version**: 1.0
**Last Updated**: 2025-11-01
**Author**: ViroForge Development Team
