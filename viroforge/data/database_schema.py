"""
Database schema for ViroForge viral genome database.

This module defines the SQLite database schema for storing viral genomes,
taxonomy, host associations, ecological metadata, and body site collections.
"""

import sqlite3
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


# SQL schema definitions
SCHEMA_VERSION = "1.0.0"

CREATE_GENOMES_TABLE = """
CREATE TABLE IF NOT EXISTS genomes (
    -- Primary identifiers
    genome_id TEXT PRIMARY KEY,
    genome_name TEXT NOT NULL,

    -- Sequence data
    sequence TEXT NOT NULL,
    length INTEGER NOT NULL,
    gc_content REAL NOT NULL,

    -- Genome characteristics
    genome_type TEXT NOT NULL,
    genome_structure TEXT,
    n_segments INTEGER DEFAULT 1,

    -- Quality metrics
    assembly_level TEXT,
    quality_score REAL,

    -- Provenance
    source_database TEXT NOT NULL,
    refseq_category TEXT,
    genbank_accession TEXT,

    -- Timestamps
    date_added TEXT NOT NULL,
    date_modified TEXT,

    -- Version control
    version INTEGER DEFAULT 1,

    -- Constraints
    CHECK (gc_content >= 0 AND gc_content <= 1),
    CHECK (length > 0)
);
"""

CREATE_GENOMES_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_genome_name ON genomes(genome_name);
CREATE INDEX IF NOT EXISTS idx_genome_type ON genomes(genome_type);
CREATE INDEX IF NOT EXISTS idx_length ON genomes(length);
CREATE INDEX IF NOT EXISTS idx_gc_content ON genomes(gc_content);
CREATE INDEX IF NOT EXISTS idx_source ON genomes(source_database);
"""

CREATE_TAXONOMY_TABLE = """
CREATE TABLE IF NOT EXISTS taxonomy (
    genome_id TEXT PRIMARY KEY,

    -- ICTV taxonomy hierarchy
    realm TEXT,
    kingdom TEXT,
    phylum TEXT,
    class TEXT,
    order_name TEXT,
    family TEXT NOT NULL,
    subfamily TEXT,
    genus TEXT,
    species TEXT,

    -- NCBI taxonomy
    ncbi_taxid INTEGER,

    -- Alternative names (JSON arrays)
    common_names TEXT,
    synonyms TEXT,

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);
"""

CREATE_TAXONOMY_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_family ON taxonomy(family);
CREATE INDEX IF NOT EXISTS idx_genus ON taxonomy(genus);
CREATE INDEX IF NOT EXISTS idx_species ON taxonomy(species);
CREATE INDEX IF NOT EXISTS idx_order ON taxonomy(order_name);
CREATE INDEX IF NOT EXISTS idx_ncbi_taxid ON taxonomy(ncbi_taxid);
"""

CREATE_HOST_ASSOCIATIONS_TABLE = """
CREATE TABLE IF NOT EXISTS host_associations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    genome_id TEXT NOT NULL,

    -- Host taxonomy
    host_domain TEXT,
    host_phylum TEXT,
    host_class TEXT,
    host_order TEXT,
    host_family TEXT,
    host_genus TEXT,
    host_species TEXT,

    -- Host name
    host_name TEXT NOT NULL,

    -- Association type
    association_type TEXT,
    evidence TEXT,

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);
"""

CREATE_HOST_ASSOCIATIONS_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_host_domain ON host_associations(host_domain);
CREATE INDEX IF NOT EXISTS idx_host_species ON host_associations(host_species);
CREATE INDEX IF NOT EXISTS idx_genome_host ON host_associations(genome_id);
"""

CREATE_ECOLOGICAL_METADATA_TABLE = """
CREATE TABLE IF NOT EXISTS ecological_metadata (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    genome_id TEXT NOT NULL,

    -- Body sites (human/animal)
    body_site TEXT,
    body_site_specific TEXT,

    -- Environmental niches
    environment_type TEXT,
    environment_subtype TEXT,

    -- Sample metadata
    isolation_source TEXT,
    geographic_origin TEXT,
    latitude REAL,
    longitude REAL,

    -- Collection metadata
    collection_date TEXT,
    collected_by TEXT,

    -- Evidence
    evidence_type TEXT,
    pubmed_id TEXT,

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);
"""

CREATE_ECOLOGICAL_METADATA_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_body_site ON ecological_metadata(body_site);
CREATE INDEX IF NOT EXISTS idx_environment ON ecological_metadata(environment_type);
CREATE INDEX IF NOT EXISTS idx_isolation ON ecological_metadata(isolation_source);
"""

CREATE_GENOME_ANNOTATIONS_TABLE = """
CREATE TABLE IF NOT EXISTS genome_annotations (
    genome_id TEXT PRIMARY KEY,

    -- Gene counts
    n_genes INTEGER,
    n_coding_genes INTEGER,
    n_rna_genes INTEGER,

    -- Genome features
    coding_density REAL,

    -- Functional categories (counts)
    n_structural_genes INTEGER,
    n_replication_genes INTEGER,
    n_regulatory_genes INTEGER,

    -- Gene information (JSON)
    genes_json TEXT,

    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE
);
"""

CREATE_BODY_SITE_COLLECTIONS_TABLE = """
CREATE TABLE IF NOT EXISTS body_site_collections (
    collection_id INTEGER PRIMARY KEY AUTOINCREMENT,

    -- Collection metadata
    collection_name TEXT UNIQUE NOT NULL,
    description TEXT,

    -- Collection properties
    n_genomes INTEGER NOT NULL,
    selection_criteria TEXT,

    -- Curation
    curated_by TEXT,
    curation_date TEXT,
    literature_references TEXT,

    -- Version
    version INTEGER DEFAULT 1
);
"""

CREATE_COLLECTION_GENOMES_TABLE = """
CREATE TABLE IF NOT EXISTS collection_genomes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    collection_id INTEGER NOT NULL,
    genome_id TEXT NOT NULL,

    -- Relative abundance in this body site
    relative_abundance REAL,
    prevalence REAL,

    -- Rank within collection
    abundance_rank INTEGER,

    FOREIGN KEY (collection_id) REFERENCES body_site_collections(collection_id) ON DELETE CASCADE,
    FOREIGN KEY (genome_id) REFERENCES genomes(genome_id) ON DELETE CASCADE,

    UNIQUE(collection_id, genome_id)
);
"""

CREATE_COLLECTION_GENOMES_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_collection ON collection_genomes(collection_id);
CREATE INDEX IF NOT EXISTS idx_collection_genome ON collection_genomes(genome_id);
"""

CREATE_DATABASE_METADATA_TABLE = """
CREATE TABLE IF NOT EXISTS database_metadata (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


def create_database(db_path: str) -> None:
    """
    Create a new ViroForge genome database with schema.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database file to create

    Examples
    --------
    >>> create_database('viroforge/data/viral_genomes.db')
    """
    # Create parent directory if needed
    db_file = Path(db_path)
    db_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Creating database at: {db_path}")

    # Connect to database (creates file if doesn't exist)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Create tables
        logger.info("Creating genomes table...")
        cursor.execute(CREATE_GENOMES_TABLE)
        for statement in CREATE_GENOMES_INDEXES.strip().split(';'):
            if statement.strip():
                cursor.execute(statement)

        logger.info("Creating taxonomy table...")
        cursor.execute(CREATE_TAXONOMY_TABLE)
        for statement in CREATE_TAXONOMY_INDEXES.strip().split(';'):
            if statement.strip():
                cursor.execute(statement)

        logger.info("Creating host associations table...")
        cursor.execute(CREATE_HOST_ASSOCIATIONS_TABLE)
        for statement in CREATE_HOST_ASSOCIATIONS_INDEXES.strip().split(';'):
            if statement.strip():
                cursor.execute(statement)

        logger.info("Creating ecological metadata table...")
        cursor.execute(CREATE_ECOLOGICAL_METADATA_TABLE)
        for statement in CREATE_ECOLOGICAL_METADATA_INDEXES.strip().split(';'):
            if statement.strip():
                cursor.execute(statement)

        logger.info("Creating genome annotations table...")
        cursor.execute(CREATE_GENOME_ANNOTATIONS_TABLE)

        logger.info("Creating body site collections tables...")
        cursor.execute(CREATE_BODY_SITE_COLLECTIONS_TABLE)
        cursor.execute(CREATE_COLLECTION_GENOMES_TABLE)
        for statement in CREATE_COLLECTION_GENOMES_INDEXES.strip().split(';'):
            if statement.strip():
                cursor.execute(statement)

        logger.info("Creating database metadata table...")
        cursor.execute(CREATE_DATABASE_METADATA_TABLE)

        # Initialize metadata
        from datetime import datetime
        now = datetime.now().isoformat()

        metadata = [
            ('schema_version', SCHEMA_VERSION),
            ('created_date', now),
            ('last_updated', now),
            ('total_genomes', '0'),
            ('total_families', '0'),
            ('total_collections', '0')
        ]

        cursor.executemany(
            "INSERT INTO database_metadata (key, value) VALUES (?, ?)",
            metadata
        )

        # Commit changes
        conn.commit()
        logger.info(f"Database created successfully at: {db_path}")
        logger.info(f"Schema version: {SCHEMA_VERSION}")

    except Exception as e:
        conn.rollback()
        logger.error(f"Error creating database: {e}")
        raise

    finally:
        conn.close()


def verify_schema(db_path: str) -> bool:
    """
    Verify that database has correct schema.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database file

    Returns
    -------
    bool
        True if schema is valid, False otherwise
    """
    expected_tables = [
        'genomes',
        'taxonomy',
        'host_associations',
        'ecological_metadata',
        'genome_annotations',
        'body_site_collections',
        'collection_genomes',
        'database_metadata'
    ]

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Get list of tables
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]

        # Check all expected tables exist
        for table in expected_tables:
            if table not in tables:
                logger.error(f"Missing table: {table}")
                return False

        # Check schema version
        cursor.execute(
            "SELECT value FROM database_metadata WHERE key = 'schema_version'"
        )
        result = cursor.fetchone()
        if result is None:
            logger.error("Missing schema_version in metadata")
            return False

        version = result[0]
        if version != SCHEMA_VERSION:
            logger.warning(
                f"Schema version mismatch: expected {SCHEMA_VERSION}, got {version}"
            )

        conn.close()
        logger.info("Schema verification passed")
        return True

    except Exception as e:
        logger.error(f"Error verifying schema: {e}")
        return False


def get_database_stats(db_path: str) -> dict:
    """
    Get database statistics.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database file

    Returns
    -------
    dict
        Dictionary of database statistics
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    stats = {}

    # Get genome count
    cursor.execute("SELECT COUNT(*) FROM genomes")
    stats['total_genomes'] = cursor.fetchone()[0]

    # Get family count
    cursor.execute("SELECT COUNT(DISTINCT family) FROM taxonomy WHERE family IS NOT NULL")
    stats['total_families'] = cursor.fetchone()[0]

    # Get collection count
    cursor.execute("SELECT COUNT(*) FROM body_site_collections")
    stats['total_collections'] = cursor.fetchone()[0]

    # Get metadata
    cursor.execute("SELECT key, value FROM database_metadata")
    for key, value in cursor.fetchall():
        stats[key] = value

    # Get genome type distribution
    cursor.execute("SELECT genome_type, COUNT(*) FROM genomes GROUP BY genome_type")
    genome_types = {}
    for genome_type, count in cursor.fetchall():
        genome_types[genome_type] = count
    stats['genome_types'] = genome_types

    # Get body site distribution
    cursor.execute(
        "SELECT body_site, COUNT(DISTINCT genome_id) FROM ecological_metadata "
        "WHERE body_site IS NOT NULL GROUP BY body_site"
    )
    body_sites = {}
    for site, count in cursor.fetchall():
        body_sites[site] = count
    stats['body_sites'] = body_sites

    conn.close()
    return stats


if __name__ == "__main__":
    # Example usage: create test database
    import sys

    if len(sys.argv) > 1:
        db_path = sys.argv[1]
    else:
        db_path = "viroforge/data/viral_genomes.db"

    print(f"Creating database: {db_path}")
    create_database(db_path)

    print("\nVerifying schema...")
    if verify_schema(db_path):
        print("✓ Schema verification passed")
    else:
        print("✗ Schema verification failed")
        sys.exit(1)

    print("\nDatabase statistics:")
    stats = get_database_stats(db_path)
    for key, value in stats.items():
        print(f"  {key}: {value}")
