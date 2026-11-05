"""Tests for CollectionLoader backward compatibility (issue #78)."""

import sqlite3
from datetime import datetime

import pytest

from viroforge.core.collection import CollectionLoader


class TestHostPctColumnRename:
    """Old databases use default_host_dna_pct; new ones use default_host_pct."""

    @pytest.fixture
    def db_old_column_name(self, tmp_path):
        db_path = str(tmp_path / "old_host.db")
        conn = sqlite3.connect(db_path)
        conn.execute("""
            CREATE TABLE body_site_collections (
                collection_id INTEGER PRIMARY KEY,
                collection_name TEXT NOT NULL,
                description TEXT,
                n_genomes INTEGER DEFAULT 0,
                default_host_dna_pct REAL
            )
        """)
        conn.execute("""
            INSERT INTO body_site_collections
                (collection_id, collection_name, n_genomes, default_host_dna_pct)
            VALUES (1, 'Blood Virome', 1, 40.0)
        """)
        conn.execute("""
            CREATE TABLE genomes (
                genome_id TEXT PRIMARY KEY,
                genome_name TEXT NOT NULL,
                sequence TEXT NOT NULL,
                length INTEGER NOT NULL,
                gc_content REAL NOT NULL,
                genome_type TEXT NOT NULL,
                genome_structure TEXT,
                n_segments INTEGER DEFAULT 1,
                assembly_level TEXT,
                quality_score REAL,
                source_database TEXT,
                refseq_category TEXT,
                genbank_accession TEXT,
                date_added TEXT,
                date_modified TEXT,
                version INTEGER DEFAULT 1
            )
        """)
        conn.execute("""
            CREATE TABLE taxonomy (
                genome_id TEXT PRIMARY KEY,
                realm TEXT, kingdom TEXT, phylum TEXT, class TEXT,
                order_name TEXT, family TEXT NOT NULL DEFAULT 'Unknown',
                subfamily TEXT, genus TEXT, species TEXT,
                ncbi_taxid INTEGER, common_names TEXT, synonyms TEXT
            )
        """)
        conn.execute("""
            CREATE TABLE collection_genomes (
                collection_id INTEGER,
                genome_id TEXT,
                relative_abundance REAL DEFAULT 0.0,
                abundance_rank INTEGER,
                PRIMARY KEY (collection_id, genome_id)
            )
        """)
        now = datetime.now().isoformat()
        conn.execute("""
            INSERT INTO genomes (genome_id, genome_name, sequence, length,
                gc_content, genome_type, date_added, date_modified)
            VALUES ('NC_001', 'TTV', 'ATGC', 4, 0.5, 'ssDNA', ?, ?)
        """, (now, now))
        conn.execute("INSERT INTO taxonomy (genome_id, family) VALUES ('NC_001', 'Anelloviridae')")
        conn.execute("""
            INSERT INTO collection_genomes (collection_id, genome_id, relative_abundance, abundance_rank)
            VALUES (1, 'NC_001', 1.0, 1)
        """)
        conn.commit()
        conn.close()
        return db_path

    def test_old_column_normalized_to_new_name(self, db_old_column_name):
        loader = CollectionLoader(db_old_column_name)
        meta, _ = loader.load_collection(1)
        assert meta.get('default_host_pct') == 40.0
        assert 'default_host_dna_pct' not in meta

    def test_generator_check_passes_with_old_db(self, db_old_column_name):
        """The check in generator.py line 1864 should find the value."""
        loader = CollectionLoader(db_old_column_name)
        meta, _ = loader.load_collection(1)
        collection_defaults = meta if meta.get('default_host_pct') is not None else None
        assert collection_defaults is not None
        assert collection_defaults['default_host_pct'] == 40.0
