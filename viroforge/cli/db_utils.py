#!/usr/bin/env python3
"""
Database utility functions for CLI

Provides convenient functions for accessing ViroForge genome database
from CLI tools.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sqlite3
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json


def get_database_path() -> Path:
    """Get path to ViroForge genome database."""
    # Try multiple possible locations
    possible_paths = [
        Path(__file__).parent.parent / "data" / "viral_genomes.db",
        Path(__file__).parent.parent.parent / "viroforge" / "data" / "viral_genomes.db",
        Path.cwd() / "viroforge" / "data" / "viral_genomes.db",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    raise FileNotFoundError(
        "Could not find viral_genomes.db. "
        "Make sure ViroForge is properly installed."
    )


def load_all_collections() -> List[Dict]:
    """
    Load all collections from database.

    Returns
    -------
    List[Dict]
        List of collection dictionaries with metadata
    """
    db_path = get_database_path()
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Query collections with genome counts
    query = """
    SELECT
        c.collection_id,
        c.collection_name,
        c.description,
        c.n_genomes,
        c.selection_criteria,
        c.literature_references,
        c.curation_date,
        COUNT(DISTINCT cg.genome_id) as actual_genome_count
    FROM body_site_collections c
    LEFT JOIN collection_genomes cg ON c.collection_id = cg.collection_id
    GROUP BY c.collection_id
    ORDER BY c.collection_id
    """

    cursor.execute(query)
    collections = []

    for row in cursor.fetchall():
        collection = {
            'id': row['collection_id'],
            'name': row['collection_name'],
            'description': row['description'],
            'n_genomes': row['actual_genome_count'],
            'selection_criteria': row['selection_criteria'],
            'literature_references': row['literature_references'],
            'curation_date': row['curation_date'],
        }
        collections.append(collection)

    conn.close()
    return collections


def load_collection_details(collection_id: int) -> Dict:
    """
    Load detailed information about a collection.

    Parameters
    ----------
    collection_id : int
        Collection ID

    Returns
    -------
    Dict
        Detailed collection information including genome list
    """
    db_path = get_database_path()
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    # Get collection metadata
    cursor.execute(
        """
        SELECT * FROM body_site_collections
        WHERE collection_id = ?
        """,
        (collection_id,)
    )
    collection_row = cursor.fetchone()

    if not collection_row:
        conn.close()
        raise ValueError(f"Collection {collection_id} not found")

    collection = dict(collection_row)

    # Get genomes in this collection
    cursor.execute(
        """
        SELECT
            cg.genome_id,
            cg.relative_abundance,
            cg.abundance_rank,
            g.genome_name,
            g.genome_type,
            g.length,
            g.gc_content,
            t.family,
            t.genus,
            t.species
        FROM collection_genomes cg
        JOIN genomes g ON cg.genome_id = g.genome_id
        LEFT JOIN taxonomy t ON cg.genome_id = t.genome_id
        WHERE cg.collection_id = ?
        ORDER BY cg.abundance_rank
        """,
        (collection_id,)
    )

    genomes = []
    genome_types = {}
    total_abundance = 0.0

    for row in cursor.fetchall():
        genome = dict(row)
        genomes.append(genome)

        # Track genome type distribution
        genome_type = genome['genome_type']
        genome_types[genome_type] = genome_types.get(genome_type, 0) + 1

        # Track abundance
        if genome['relative_abundance']:
            total_abundance += genome['relative_abundance']

    collection['genomes'] = genomes
    collection['genome_types'] = genome_types
    collection['actual_genome_count'] = len(genomes)

    # Get family distribution
    cursor.execute(
        """
        SELECT t.family, COUNT(*) as count
        FROM collection_genomes cg
        JOIN taxonomy t ON cg.genome_id = t.genome_id
        WHERE cg.collection_id = ? AND t.family IS NOT NULL
        GROUP BY t.family
        ORDER BY count DESC
        LIMIT 10
        """,
        (collection_id,)
    )

    families = []
    for row in cursor.fetchall():
        families.append({
            'family': row[0],
            'count': row[1]
        })

    collection['top_families'] = families

    conn.close()
    return collection


def search_collections(query: str) -> List[Dict]:
    """
    Search collections by name or description.

    Parameters
    ----------
    query : str
        Search query

    Returns
    -------
    List[Dict]
        List of matching collections
    """
    collections = load_all_collections()

    if not query:
        return collections

    query_lower = query.lower()
    filtered = []

    for collection in collections:
        name_match = query_lower in collection['name'].lower()
        desc_match = query_lower in (collection['description'] or '').lower()

        if name_match or desc_match:
            filtered.append(collection)

    return filtered


def get_collection_categories() -> Dict[str, List[Dict]]:
    """
    Get collections organized by category.

    Returns
    -------
    Dict[str, List[Dict]]
        Collections grouped by category
    """
    collections = load_all_collections()

    categories = {
        'Host-Associated': [],
        'Environmental': [],
        'Disease-Associated': [],
        'RNA Viromes': [],
        'Other': []
    }

    for collection in collections:
        name = collection['name']
        desc = collection['description'] or ''

        # Categorize based on name/description
        if any(keyword in name.lower() for keyword in ['gut', 'oral', 'skin', 'vaginal', 'urinary', 'blood', 'ocular', 'lung', 'respiratory', 'fecal']):
            if any(keyword in name.lower() for keyword in ['ibd', 'hiv', 'cf', 'disease']):
                categories['Disease-Associated'].append(collection)
            else:
                categories['Host-Associated'].append(collection)
        elif 'rna' in name.lower() or 'rna' in desc.lower():
            categories['RNA Viromes'].append(collection)
        elif any(keyword in name.lower() for keyword in ['soil', 'marine', 'water', 'environment']):
            categories['Environmental'].append(collection)
        else:
            categories['Other'].append(collection)

    # Remove empty categories
    categories = {k: v for k, v in categories.items() if v}

    return categories


def get_database_stats() -> Dict:
    """
    Get database statistics.

    Returns
    -------
    Dict
        Database statistics
    """
    db_path = get_database_path()
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    stats = {}

    # Total genomes
    cursor.execute("SELECT COUNT(*) FROM genomes")
    stats['total_genomes'] = cursor.fetchone()[0]

    # Total collections
    cursor.execute("SELECT COUNT(*) FROM body_site_collections")
    stats['total_collections'] = cursor.fetchone()[0]

    # Total families
    cursor.execute("SELECT COUNT(DISTINCT family) FROM taxonomy WHERE family IS NOT NULL")
    stats['total_families'] = cursor.fetchone()[0]

    # Genome types
    cursor.execute("SELECT genome_type, COUNT(*) FROM genomes GROUP BY genome_type")
    stats['genome_types'] = dict(cursor.fetchall())

    conn.close()
    return stats
