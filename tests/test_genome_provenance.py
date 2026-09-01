"""Tests for genome provenance detection and propagation."""

import sqlite3
import tempfile
from pathlib import Path
from datetime import datetime

import pytest

from scripts.populate_database import detect_genome_provenance


class TestDetectGenomeProvenance:

    def test_prophage(self):
        assert detect_genome_provenance("Bacillus phage SPbeta prophage") == "prophage"
        assert detect_genome_provenance("Stx2-converting prophage 1717") == "prophage"

    def test_proviral(self):
        assert detect_genome_provenance("HIV-1 proviral DNA") == "provirus"
        assert detect_genome_provenance("Provirus sequence from host") == "provirus"

    def test_provirus(self):
        assert detect_genome_provenance("HTLV-1 provirus sequence") == "provirus"

    def test_endogenous(self):
        assert detect_genome_provenance("Endogenous retrovirus K113") == "endogenous"

    def test_satellite(self):
        assert detect_genome_provenance("Tobacco necrosis satellite virus") == "satellite"

    def test_viroid(self):
        assert detect_genome_provenance("Potato spindle tuber viroid") == "viroid"

    def test_isolate_default(self):
        assert detect_genome_provenance("Enterobacteria phage T4") == "isolate"
        assert detect_genome_provenance("Human herpesvirus 5") == "isolate"
        assert detect_genome_provenance("Pandoravirus salinus") == "isolate"

    def test_case_insensitive(self):
        assert detect_genome_provenance("PROPHAGE element") == "prophage"
        assert detect_genome_provenance("Endogenous Retrovirus") == "endogenous"

    def test_empty_string(self):
        assert detect_genome_provenance("") == "isolate"

    def test_priority_endogenous_over_provirus(self):
        # "endogenous" checked before "provirus" keywords
        assert detect_genome_provenance("Endogenous provirus element") == "endogenous"


class TestProvenanceInSchema:

    def test_genome_provenance_column_exists(self):
        from viroforge.data.database_schema import CREATE_GENOMES_TABLE
        assert "genome_provenance" in CREATE_GENOMES_TABLE

    def test_provenance_default_is_isolate(self):
        from viroforge.data.database_schema import CREATE_GENOMES_TABLE
        assert "DEFAULT 'isolate'" in CREATE_GENOMES_TABLE

    def test_provenance_index_exists(self):
        from viroforge.data.database_schema import CREATE_GENOMES_INDEXES
        assert "idx_provenance" in CREATE_GENOMES_INDEXES


class TestProvenanceInDatabase:

    @pytest.fixture
    def db_with_provenance(self, tmp_path):
        """Create a minimal database with genomes of different provenance."""
        from viroforge.data.database_schema import (
            CREATE_GENOMES_TABLE, CREATE_GENOMES_INDEXES,
            CREATE_TAXONOMY_TABLE, CREATE_TAXONOMY_INDEXES,
            CREATE_BODY_SITE_COLLECTIONS_TABLE,
            CREATE_COLLECTION_GENOMES_TABLE, CREATE_COLLECTION_GENOMES_INDEXES,
        )
        db_path = str(tmp_path / "test.db")
        conn = sqlite3.connect(db_path)
        conn.execute(CREATE_GENOMES_TABLE)
        for stmt in CREATE_GENOMES_INDEXES.strip().split(';'):
            if stmt.strip():
                conn.execute(stmt)
        conn.execute(CREATE_TAXONOMY_TABLE)
        for stmt in CREATE_TAXONOMY_INDEXES.strip().split(';'):
            if stmt.strip():
                conn.execute(stmt)
        conn.execute(CREATE_BODY_SITE_COLLECTIONS_TABLE)
        conn.execute(CREATE_COLLECTION_GENOMES_TABLE)
        for stmt in CREATE_COLLECTION_GENOMES_INDEXES.strip().split(';'):
            if stmt.strip():
                conn.execute(stmt)

        now = datetime.now().isoformat()
        genomes = [
            ("NC_001", "Enterobacteria phage T4", "ATGC" * 1000, 4000, 0.5,
             "dsDNA", "linear", 1, "isolate", "Complete Genome", None,
             "RefSeq", None, None, now, now, 1),
            ("NC_002", "Bacillus phage SPbeta prophage", "GCTA" * 500, 2000, 0.45,
             "dsDNA", "linear", 1, "prophage", "Complete Genome", None,
             "RefSeq", None, None, now, now, 1),
            ("NC_003", "Endogenous retrovirus K113", "ATAT" * 750, 3000, 0.42,
             "ssRNA(RT)", "linear", 1, "endogenous", "Complete Genome", None,
             "RefSeq", None, None, now, now, 1),
        ]
        conn.executemany("""
            INSERT INTO genomes (genome_id, genome_name, sequence, length,
                gc_content, genome_type, genome_structure, n_segments,
                genome_provenance, assembly_level, quality_score,
                source_database, refseq_category, genbank_accession,
                date_added, date_modified, version)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, genomes)

        taxonomy = [
            ("NC_001", None, None, None, None, None, "Straboviridae", None,
             "Tequatrovirus", "Tequatrovirus T4", None, None, None),
            ("NC_002", None, None, None, None, None, "Unknown", None,
             None, None, None, None, None),
            ("NC_003", None, None, None, None, None, "Retroviridae", None,
             "Betaretrovirus", None, None, None, None),
        ]
        conn.executemany("""
            INSERT INTO taxonomy (genome_id, realm, kingdom, phylum, class,
                order_name, family, subfamily, genus, species, ncbi_taxid,
                common_names, synonyms)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, taxonomy)

        conn.execute("""
            INSERT INTO body_site_collections
                (collection_id, collection_name, description, n_genomes)
            VALUES (1, 'Test Collection', 'Test', 3)
        """)
        for i, gid in enumerate(["NC_001", "NC_002", "NC_003"]):
            conn.execute("""
                INSERT INTO collection_genomes
                    (collection_id, genome_id, relative_abundance, abundance_rank)
                VALUES (1, ?, ?, ?)
            """, (gid, 1.0 / 3, i + 1))

        conn.commit()
        conn.close()
        return db_path

    def test_load_collection_includes_provenance(self, db_with_provenance):
        from viroforge.core.collection import CollectionLoader
        loader = CollectionLoader(db_with_provenance)
        _, genomes = loader.load_collection(1)

        provenance_by_id = {g['genome_id']: g['genome_provenance'] for g in genomes}
        assert provenance_by_id['NC_001'] == 'isolate'
        assert provenance_by_id['NC_002'] == 'prophage'
        assert provenance_by_id['NC_003'] == 'endogenous'

    def test_source_label_uses_provenance(self, db_with_provenance):
        from viroforge.core.collection import CollectionLoader
        loader = CollectionLoader(db_with_provenance)
        _, genomes = loader.load_collection(1)

        provenance_map = {g['genome_id']: g.get('genome_provenance', 'isolate') for g in genomes}

        # Simulate the source label logic from generator.py
        labels = {}
        for g in genomes:
            provenance = provenance_map.get(g['genome_id'], 'isolate')
            if provenance != 'isolate':
                labels[g['genome_id']] = f"viral_{provenance}"
            else:
                labels[g['genome_id']] = "viral"

        assert labels['NC_001'] == 'viral'
        assert labels['NC_002'] == 'viral_prophage'
        assert labels['NC_003'] == 'viral_endogenous'

    def test_provenance_labels_backward_compatible(self, db_with_provenance):
        """All provenance-aware labels start with 'viral' for backward compat."""
        from viroforge.core.collection import CollectionLoader
        loader = CollectionLoader(db_with_provenance)
        _, genomes = loader.load_collection(1)

        provenance_map = {g['genome_id']: g.get('genome_provenance', 'isolate') for g in genomes}
        for gid, prov in provenance_map.items():
            if prov != 'isolate':
                label = f"viral_{prov}"
            else:
                label = "viral"
            assert label.startswith("viral")
