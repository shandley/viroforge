# Genome Database Design Document

**Phase 3: Genome Database Expansion**

**Version**: 1.0
**Date**: October 31, 2025
**Status**: Planning & Design

---

## Executive Summary

This document details the design for ViroForge's comprehensive viral genome database, which will serve as the foundation for realistic virome simulations. The database will contain 10,000+ viral genomes with rich metadata including ICTV taxonomy, host range, ecological niches, and body site associations.

---

## Goals & Requirements

### Primary Goals

1. **Comprehensive Coverage**: 10,000+ high-quality viral genomes
2. **Rich Metadata**: ICTV taxonomy, host, ecology, isolation source
3. **Fast Queries**: <100ms for typical queries
4. **Easy Updates**: Automated RefSeq synchronization
5. **Body Site Collections**: Pre-curated collections for common use cases
6. **Version Control**: Track changes, updates, additions

### Requirements

**Functional Requirements**:
- Query genomes by taxonomy (family, genus, species)
- Filter by body site, host, ecological niche
- Select representative genomes for body sites
- Export genome sequences and metadata
- Support custom genome additions
- Track provenance and quality

**Non-Functional Requirements**:
- Query performance: <100ms for filtering
- Database size: <5GB including sequences
- Update frequency: Monthly (automated)
- Data quality: Only verified, complete genomes
- Reproducibility: Version tracking for citations

---

## Data Sources

### 1. NCBI RefSeq Viral Genomes

**URL**: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/

**What We Get**:
- Complete viral genomes (verified, curated)
- ~15,000 viral genome assemblies
- NCBI taxonomy IDs
- GenBank annotations
- Sequence data (FASTA)

**Download Strategy**:
```bash
# Download RefSeq viral genome catalog
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt

# Download genomes in batches
# Parse assembly_summary.txt for FTP paths
# Download genomic.fna.gz files
```

**Quality Filtering**:
- Only "Complete Genome" assembly level
- Exclude draft/partial genomes
- Exclude metagenome-assembled genomes (MAGs)
- Minimum length: 1,000 bp
- Maximum length: 500,000 bp (exclude giant viruses for now)

### 2. ICTV Taxonomy Database

**URL**: https://ictv.global/taxonomy

**What We Get**:
- Official virus taxonomy (realm → species)
- Virus names and classifications
- Regular updates (annual)

**Format**: VMR (Virus Metadata Resource) - Excel spreadsheet

**Fields We Need**:
- Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species
- Virus name(s)
- Genome composition (DNA/RNA, ss/ds)
- Host source

### 3. Literature-Derived Metadata

**Sources**:
- Body site-specific papers (gut, oral, skin virome studies)
- Environmental virome papers (marine, soil, etc.)
- Host-virus association databases
- Ecological metadata from publications

**Curation Strategy**:
- Manual curation for major body sites (gut, oral, skin)
- Automated extraction from papers (where possible)
- Community contributions (future)

---

## Database Schema

### Technology Choice: SQLite

**Why SQLite**:
- ✅ Serverless (embedded in Python)
- ✅ Fast for read-heavy workloads
- ✅ Single file (easy distribution)
- ✅ ACID compliant
- ✅ Well-supported in Python
- ✅ No server setup required

**Schema Design**:

### Table 1: `genomes`

Core genome information and sequences.

```sql
CREATE TABLE genomes (
    -- Primary identifiers
    genome_id TEXT PRIMARY KEY,           -- NC_001416 (RefSeq accession)
    genome_name TEXT NOT NULL,            -- Enterobacteria phage T7

    -- Sequence data
    sequence TEXT NOT NULL,               -- Full genome sequence (FASTA)
    length INTEGER NOT NULL,              -- Genome length in bp
    gc_content REAL NOT NULL,             -- GC content (0-1)

    -- Genome characteristics
    genome_type TEXT NOT NULL,            -- dsDNA, ssDNA, dsRNA, ssRNA, ssRNA-RT, dsDNA-RT
    genome_structure TEXT,                -- linear, circular, segmented
    n_segments INTEGER DEFAULT 1,        -- Number of segments

    -- Quality metrics
    assembly_level TEXT,                  -- Complete Genome, Scaffold, Contig
    quality_score REAL,                   -- Custom quality metric

    -- Provenance
    source_database TEXT NOT NULL,        -- RefSeq, custom, etc.
    refseq_category TEXT,                 -- reference, representative
    genbank_accession TEXT,               -- GenBank accession if different

    -- Timestamps
    date_added TEXT NOT NULL,             -- ISO 8601 date
    date_modified TEXT,                   -- Last modification date

    -- Version control
    version INTEGER DEFAULT 1,            -- Version number

    -- Indexes
    CHECK (gc_content >= 0 AND gc_content <= 1),
    CHECK (length > 0)
);

CREATE INDEX idx_genome_name ON genomes(genome_name);
CREATE INDEX idx_genome_type ON genomes(genome_type);
CREATE INDEX idx_length ON genomes(length);
CREATE INDEX idx_gc_content ON genomes(gc_content);
```

### Table 2: `taxonomy`

ICTV taxonomy information.

```sql
CREATE TABLE taxonomy (
    genome_id TEXT PRIMARY KEY,

    -- ICTV taxonomy hierarchy
    realm TEXT,                           -- Riboviria, Duplodnaviria, etc.
    kingdom TEXT,
    phylum TEXT,
    class TEXT,
    order_name TEXT,                      -- 'order' is SQL keyword
    family TEXT NOT NULL,                 -- Minimum required
    subfamily TEXT,
    genus TEXT,
    species TEXT,

    -- NCBI taxonomy
    ncbi_taxid INTEGER,                   -- NCBI Taxonomy ID

    -- Alternative names
    common_names TEXT,                    -- JSON array of common names
    synonyms TEXT,                        -- JSON array of synonyms

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);

CREATE INDEX idx_family ON taxonomy(family);
CREATE INDEX idx_genus ON taxonomy(genus);
CREATE INDEX idx_species ON taxonomy(species);
CREATE INDEX idx_order ON taxonomy(order_name);
```

### Table 3: `host_associations`

Host range and virus-host relationships.

```sql
CREATE TABLE host_associations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    genome_id TEXT NOT NULL,

    -- Host taxonomy
    host_domain TEXT,                     -- Bacteria, Archaea, Eukaryota
    host_phylum TEXT,
    host_class TEXT,
    host_order TEXT,
    host_family TEXT,
    host_genus TEXT,
    host_species TEXT,

    -- Host name
    host_name TEXT NOT NULL,              -- Escherichia coli, Homo sapiens

    -- Association type
    association_type TEXT,                -- experimental, predicted, literature
    evidence TEXT,                        -- PubMed ID, database reference

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);

CREATE INDEX idx_host_domain ON host_associations(host_domain);
CREATE INDEX idx_host_species ON host_associations(host_species);
CREATE INDEX idx_genome_host ON host_associations(genome_id);
```

### Table 4: `ecological_metadata`

Ecological and environmental metadata.

```sql
CREATE TABLE ecological_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    genome_id TEXT NOT NULL,

    -- Body sites (human/animal)
    body_site TEXT,                       -- gut, oral, skin, respiratory, etc.
    body_site_specific TEXT,              -- colon, small_intestine, saliva, etc.

    -- Environmental niches
    environment_type TEXT,                -- marine, freshwater, soil, etc.
    environment_subtype TEXT,             -- coastal, deep_sea, agricultural_soil

    -- Sample metadata
    isolation_source TEXT,                -- human_feces, seawater, soil_sample
    geographic_origin TEXT,               -- Country/region
    latitude REAL,
    longitude REAL,

    -- Collection metadata
    collection_date TEXT,                 -- ISO 8601 date
    collected_by TEXT,                    -- Institution/person

    -- Evidence
    evidence_type TEXT,                   -- isolated, sequenced, MAG, literature
    pubmed_id TEXT,                       -- PubMed reference

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);

CREATE INDEX idx_body_site ON ecological_metadata(body_site);
CREATE INDEX idx_environment ON ecological_metadata(environment_type);
CREATE INDEX idx_isolation ON ecological_metadata(isolation_source);
```

### Table 5: `genome_annotations`

Functional annotations and gene information.

```sql
CREATE TABLE genome_annotations (
    genome_id TEXT PRIMARY KEY,

    -- Gene counts
    n_genes INTEGER,                      -- Total number of genes
    n_coding_genes INTEGER,              -- Protein-coding genes
    n_rna_genes INTEGER,                 -- RNA genes (tRNA, rRNA)

    -- Genome features
    coding_density REAL,                  -- Proportion of genome that codes

    -- Functional categories (counts)
    n_structural_genes INTEGER,          -- Head, tail, capsid proteins
    n_replication_genes INTEGER,         -- DNA/RNA polymerase, etc.
    n_regulatory_genes INTEGER,          -- Transcription factors, etc.

    -- Gene information (JSON)
    genes_json TEXT,                      -- Detailed gene information

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);
```

### Table 6: `body_site_collections`

Pre-curated genome collections for body sites.

```sql
CREATE TABLE body_site_collections (
    collection_id INTEGER PRIMARY KEY AUTOINCREMENT,

    -- Collection metadata
    collection_name TEXT UNIQUE NOT NULL, -- gut_virome, oral_virome, etc.
    description TEXT,

    -- Collection properties
    n_genomes INTEGER NOT NULL,           -- Target number of genomes
    selection_criteria TEXT,              -- How genomes were selected

    -- Curation
    curated_by TEXT,
    curation_date TEXT,
    literature_references TEXT,           -- JSON array of references

    -- Version
    version INTEGER DEFAULT 1
);

CREATE TABLE collection_genomes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    collection_id INTEGER NOT NULL,
    genome_id TEXT NOT NULL,

    -- Relative abundance in this body site
    relative_abundance REAL,              -- Expected abundance (0-1)
    prevalence REAL,                      -- Prevalence across samples (0-1)

    -- Rank within collection
    abundance_rank INTEGER,               -- 1 = most abundant

    FOREIGN KEY (collection_id) REFERENCES body_site_collections(collection_id) ON DELETE CASCADE,
    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE,

    UNIQUE(collection_id, genome_id)
);

CREATE INDEX idx_collection ON collection_genomes(collection_id);
```

### Table 7: `database_metadata`

Database version and update information.

```sql
CREATE TABLE database_metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);

-- Initial values:
-- version: 1.0.0
-- created_date: 2025-10-31
-- last_updated: 2025-10-31
-- refseq_release: 223
-- ictv_release: 2024
-- total_genomes: 15000
-- total_families: 150
```

---

## Query Interface Design

### Python API

```python
from viroforge.data.genome_database import GenomeDatabase

# Initialize database
db = GenomeDatabase('viroforge/data/viral_genomes.db')

# Query by taxonomy
gut_caudovirales = db.query(
    family='Caudovirales',
    body_site='gut',
    limit=100
)

# Query by genome properties
large_gc_rich = db.query(
    length_min=50000,
    gc_min=0.6,
    gc_max=0.7
)

# Query by host
e_coli_phages = db.query(
    host_species='Escherichia coli',
    genome_type='dsDNA'
)

# Get body site collection
gut_virome = db.get_collection('gut_virome')
oral_virome = db.get_collection('oral_virome')

# Get specific genome
t7_phage = db.get_genome('NC_001416')

# Random sampling
random_genomes = db.random_sample(
    n=50,
    family='Siphoviridae',
    body_site='gut'
)

# Statistics
stats = db.get_statistics()
# Returns: {
#     'total_genomes': 15234,
#     'families': 127,
#     'genera': 892,
#     'genome_types': {'dsDNA': 12000, 'ssDNA': 800, ...},
#     'body_sites': {'gut': 2300, 'oral': 800, ...}
# }
```

### GenomeDatabase Class Design

```python
class GenomeDatabase:
    """Interface to ViroForge viral genome database."""

    def __init__(self, db_path: str):
        """Initialize connection to SQLite database."""

    def query(
        self,
        genome_id: Optional[str] = None,
        family: Optional[str] = None,
        genus: Optional[str] = None,
        species: Optional[str] = None,
        genome_type: Optional[str] = None,
        host_species: Optional[str] = None,
        host_domain: Optional[str] = None,
        body_site: Optional[str] = None,
        environment: Optional[str] = None,
        length_min: Optional[int] = None,
        length_max: Optional[int] = None,
        gc_min: Optional[float] = None,
        gc_max: Optional[float] = None,
        limit: Optional[int] = None,
        random_seed: Optional[int] = None
    ) -> List[ViralGenome]:
        """Query database with filters."""

    def get_genome(self, genome_id: str) -> ViralGenome:
        """Get specific genome by ID."""

    def get_collection(self, collection_name: str) -> ViralCollection:
        """Get pre-curated genome collection."""

    def random_sample(
        self,
        n: int,
        **filters
    ) -> List[ViralGenome]:
        """Randomly sample genomes with filters."""

    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics."""

    def add_genome(
        self,
        genome: ViralGenome,
        taxonomy: Dict,
        metadata: Dict
    ) -> None:
        """Add custom genome to database."""

    def update_from_refseq(self) -> Dict[str, int]:
        """Update database from RefSeq (monthly)."""
```

---

## Body Site Collection Specifications

### Human Body Sites (Priority 1)

#### 1. Gut Virome (Target: 500 genomes)

**Selection Criteria**:
- Family distribution based on literature:
  - Caudovirales (60%): Siphoviridae, Podoviridae, Myoviridae
  - Microviridae (20%): ssDNA phages
  - Other families (20%): Inoviridae, etc.
- Host range: Primarily bacterial hosts (gut bacteria)
- Isolation source: Human feces, gut samples
- Literature references: Shkoporov et al. 2019, Manrique et al. 2016

**Abundance Model**:
- Log-normal distribution
- Few dominant families, many rare genomes
- Matches real gut virome structure

#### 2. Oral Virome (Target: 200 genomes)

**Selection Criteria**:
- Family distribution:
  - Siphoviridae (50%)
  - Myoviridae (30%)
  - Other (20%)
- Isolation source: Saliva, oral swabs, dental plaque
- Host range: Oral bacteria (Streptococcus, etc.)
- Literature: Pride et al. 2012, Abeles et al. 2014

#### 3. Skin Virome (Target: 150 genomes)

**Selection Criteria**:
- Family distribution:
  - Podoviridae (40%)
  - Siphoviridae (30%)
  - Staphylococcus phages (high representation)
- Isolation source: Skin swabs
- Literature: Oh et al. 2014, Hannigan et al. 2015

#### 4. Respiratory Virome (Target: 200 genomes)

**Selection Criteria**:
- Mixed: Bacterial phages + eukaryotic viruses
- Isolation source: Nasal swabs, lung samples
- Include common respiratory viruses
- Literature: Wylie et al. 2014

### Environmental Viromes (Priority 2)

#### 5. Marine Virome (Target: 500 genomes)

**Selection Criteria**:
- Marine cyanophages (abundant)
- Heterotrophic bacterial phages
- Isolation source: Seawater
- Include diverse marine viral families
- Literature: Brum et al. 2015, Roux et al. 2016

#### 6. Soil Virome (Target: 300 genomes)

**Selection Criteria**:
- Soil bacterial phages
- Plant viruses (some)
- Isolation source: Agricultural soil, forest soil
- Literature: Emerson et al. 2018

#### 7. Freshwater Virome (Target: 200 genomes)

**Selection Criteria**:
- Freshwater phages
- Cyanophages
- Isolation source: Lakes, rivers

### Animal Models (Priority 3)

#### 8. Mouse Gut Virome (Target: 150 genomes)

**Selection Criteria**:
- Mouse gut microbiome phages
- Host: Mouse gut bacteria
- Literature: Reyes et al. 2013

---

## Implementation Plan

### Phase 3A: Database Foundation (Weeks 1-2)

**Tasks**:
1. Set up SQLite database with schema
2. Create GenomeDatabase class
3. Implement basic query methods
4. Write unit tests for database operations

**Deliverables**:
- `viroforge/data/genome_database.py`
- `tests/test_genome_database.py`
- Empty database with schema

### Phase 3B: RefSeq Data Acquisition (Weeks 3-4)

**Tasks**:
1. Download RefSeq viral genome catalog
2. Parse assembly_summary.txt
3. Download genome FASTA files (batched)
4. Extract ICTV taxonomy from VMR
5. Quality filtering

**Deliverables**:
- `scripts/download_refseq.py`
- `scripts/parse_ictv_taxonomy.py`
- Raw data in `data/refseq/`

### Phase 3C: Database Population (Week 5)

**Tasks**:
1. Parse genome sequences
2. Calculate genome statistics (length, GC)
3. Populate genomes table
4. Populate taxonomy table
5. Basic host associations (from RefSeq metadata)

**Deliverables**:
- `scripts/populate_database.py`
- Populated database: `viroforge/data/viral_genomes.db`
- ~10,000+ genomes loaded

### Phase 3D: Body Site Curation (Week 6)

**Tasks**:
1. Literature review for body site compositions
2. Curate gut virome collection (500 genomes)
3. Curate oral virome collection (200 genomes)
4. Curate skin virome collection (150 genomes)
5. Populate collection tables

**Deliverables**:
- `data/body_site_collections/gut_virome.tsv`
- `data/body_site_collections/oral_virome.tsv`
- Pre-curated collections in database

### Phase 3E: Integration (Week 6)

**Tasks**:
1. Update `create_body_site_profile()` to use database
2. Implement collection-based loading
3. Update tests
4. Documentation

**Deliverables**:
- Updated `viroforge/core/community.py`
- All tests passing
- Documentation updated

---

## Data Quality Standards

### Genome Inclusion Criteria

**Must Have**:
- ✅ Complete genome sequence
- ✅ RefSeq accession OR verified source
- ✅ Minimum length: 1,000 bp
- ✅ Maximum length: 500,000 bp
- ✅ Valid ICTV taxonomy (at least family)

**Should Have**:
- ✅ Host association information
- ✅ Isolation source
- ✅ Gene annotations
- ✅ Literature reference

**Quality Metrics**:
- No ambiguous bases (N) > 1%
- GC content reasonable (0.2 - 0.8)
- Assembly level: Complete Genome preferred

---

## Performance Targets

### Query Performance

- Simple queries (by family): <10ms
- Complex queries (multi-filter): <100ms
- Random sampling (n=100): <50ms
- Collection loading: <200ms

### Database Size

- Genome sequences: ~2GB (compressed)
- Metadata: <500MB
- Indexes: <200MB
- Total: <3GB

### Update Frequency

- RefSeq synchronization: Monthly (automated)
- ICTV updates: Annually (manual)
- Body site collections: Quarterly (curated)

---

## Testing Strategy

### Unit Tests

- Database connection and initialization
- Query methods (all filter combinations)
- Genome retrieval
- Collection loading
- Statistics calculation
- Data validation

### Integration Tests

- End-to-end: Download → Parse → Populate → Query
- Collection-based community creation
- Performance benchmarks
- Data integrity checks

### Validation Tests

- Genome sequence validity
- Taxonomy consistency
- Metadata completeness
- Cross-reference integrity

---

## Documentation Requirements

### User Documentation

- **Database Schema**: Complete table descriptions
- **Query Guide**: How to query the database
- **Collection Guide**: Using pre-curated collections
- **API Reference**: GenomeDatabase class methods

### Developer Documentation

- **Data Sources**: Where data comes from
- **Update Procedures**: How to update the database
- **Curation Guidelines**: How to curate body site collections
- **Quality Standards**: Data quality requirements

---

## Version Control & Releases

### Versioning Scheme

Format: `MAJOR.MINOR.PATCH`

- MAJOR: Database schema changes
- MINOR: Significant genome additions (>1000 genomes)
- PATCH: Bug fixes, small additions

Example: `1.0.0` → `1.1.0` (added 2000 genomes)

### Release Process

1. Update RefSeq genomes
2. Update ICTV taxonomy
3. Quality check (automated tests)
4. Update version number
5. Tag release in git
6. Update documentation
7. Announce to users

---

## Success Metrics

### Database Completeness

- ✅ 10,000+ genomes
- ✅ 100+ viral families
- ✅ 80% genomes with host information
- ✅ 60% genomes with body site associations
- ✅ 8+ pre-curated body site collections

### Functionality

- ✅ Fast queries (<100ms)
- ✅ All API methods working
- ✅ Comprehensive test coverage (>90%)
- ✅ Documentation complete

### Integration

- ✅ Seamless integration with existing ViroForge
- ✅ All existing tests passing
- ✅ No performance regression

---

## Risks & Mitigation

### Risk 1: RefSeq Data Quality Issues

**Mitigation**: Strict quality filtering, validation tests

### Risk 2: Database Size Too Large

**Mitigation**: Sequence compression, optional full download

### Risk 3: ICTV Taxonomy Changes

**Mitigation**: Version tracking, update procedures

### Risk 4: Query Performance

**Mitigation**: Proper indexing, query optimization

---

## Next Steps

1. **Review this design** with stakeholders
2. **Create database schema** (SQLite)
3. **Implement GenomeDatabase class** (basic)
4. **Download RefSeq sample** (100 genomes for testing)
5. **Test population pipeline**
6. **Proceed with full data acquisition**

---

**Status**: Ready for implementation ✅

**Next Document**: RefSeq Data Acquisition Plan

**Last Updated**: October 31, 2025
