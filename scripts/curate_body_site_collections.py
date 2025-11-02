#!/usr/bin/env python3
"""
Body Site Virome Collection Curation Script

Creates realistic body site-specific virome collections based on literature-validated
compositions and abundance distributions.

Usage:
    python scripts/curate_body_site_collections.py --collection gut --output-name gut_virome_adult_healthy_western
    python scripts/curate_body_site_collections.py --all  # Curate all collections

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import random
import numpy as np
from dataclasses import dataclass


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class TaxonTarget:
    """Target composition for a taxonomic group."""
    taxon_name: str
    target_count: int
    rank: str  # 'family', 'genus', 'species', etc.
    host_filter: Optional[str] = None
    genome_type: Optional[str] = None  # 'dsDNA', 'ssRNA', etc.


@dataclass
class CollectionSpec:
    """Specification for a body site collection."""
    collection_id: str
    name: str
    description: str
    body_site: str
    sample_type: str
    health_status: str
    target_size: int
    phage_fraction: float
    composition: List[TaxonTarget]
    abundance_model: str  # 'log_normal', 'power_law', 'mixed'
    dominant_tier_count: int  # Number of dominant genomes (10-30%)
    common_tier_count: int    # Number of common genomes (1-10%)
    moderate_tier_count: int  # Number of moderate genomes (0.1-1%)
    rare_tier_count: int      # Number of rare genomes (0.01-0.1%)
    very_rare_tier_count: int # Number of very rare genomes (<0.01%)


# Collection specifications based on docs/BODY_SITE_COLLECTIONS.md
COLLECTION_SPECS = {
    'gut': CollectionSpec(
        collection_id='gut_virome_adult_healthy_western',
        name='Gut Virome - Adult Healthy (Western Diet)',
        description='Human gut virome from healthy adults with Western diet. crAssphage-dominated composition typical of industrialized populations.',
        body_site='gut',
        sample_type='stool',
        health_status='healthy',
        target_size=500,
        phage_fraction=0.96,
        composition=[
            TaxonTarget('Crassvirales', 150, 'order'),  # 30%
            TaxonTarget('Siphoviridae', 140, 'family'),  # 28%
            TaxonTarget('Myoviridae', 70, 'family'),     # 14%
            TaxonTarget('Podoviridae', 40, 'family'),    # 8%
            TaxonTarget('Microviridae', 50, 'family'),   # 10%
            TaxonTarget('Inoviridae', 15, 'family'),     # 3%
            TaxonTarget('Adenoviridae', 10, 'family', genome_type='dsDNA'),  # 2%
            TaxonTarget('Anelloviridae', 10, 'family', genome_type='ssDNA'), # 2%
            TaxonTarget('Other', 15, 'any'),  # 3%
        ],
        abundance_model='mixed',
        dominant_tier_count=8,    # 10-30% each
        common_tier_count=75,     # 1-10% each
        moderate_tier_count=175,  # 0.1-1% each
        rare_tier_count=217,      # 0.01-0.1% each
        very_rare_tier_count=25   # <0.01% each
    ),

    'oral': CollectionSpec(
        collection_id='oral_virome_saliva_healthy',
        name='Oral Virome - Saliva (Healthy)',
        description='Oral virome from healthy saliva. Streptococcus phage-dominated with high Siphoviridae prevalence.',
        body_site='oral_cavity',
        sample_type='saliva',
        health_status='healthy',
        target_size=200,
        phage_fraction=0.925,
        composition=[
            TaxonTarget('Siphoviridae', 80, 'family'),  # 40%
            TaxonTarget('Myoviridae', 40, 'family'),    # 20%
            TaxonTarget('Podoviridae', 30, 'family'),   # 15%
            TaxonTarget('Microviridae', 20, 'family'),  # 10%
            TaxonTarget('Inoviridae', 15, 'family'),    # 7.5%
            TaxonTarget('Herpesviridae', 8, 'family', genome_type='dsDNA'),  # 4%
            TaxonTarget('Papillomaviridae', 4, 'family', genome_type='dsDNA'), # 2%
            TaxonTarget('Anelloviridae', 3, 'family', genome_type='ssDNA'),    # 1.5%
        ],
        abundance_model='log_normal',
        dominant_tier_count=5,
        common_tier_count=30,
        moderate_tier_count=70,
        rare_tier_count=80,
        very_rare_tier_count=15
    ),

    'skin': CollectionSpec(
        collection_id='skin_virome_sebaceous_healthy',
        name='Skin Virome - Sebaceous Sites (Healthy)',
        description='Skin virome from sebaceous sites (forehead, back). Cutibacterium (P. acnes) phage-dominated.',
        body_site='skin',
        sample_type='swab',
        health_status='healthy',
        target_size=150,
        phage_fraction=0.96,
        composition=[
            TaxonTarget('Cutibacterium', 70, 'genus', host_filter='Cutibacterium'),  # 47%
            TaxonTarget('Staphylococcus', 40, 'genus', host_filter='Staphylococcus'), # 27%
            TaxonTarget('Corynebacterium', 15, 'genus', host_filter='Corynebacterium'), # 10%
            TaxonTarget('Other_Caudoviricetes', 10, 'class'),  # 7%
            TaxonTarget('Microviridae', 8, 'family'),  # 5%
            TaxonTarget('Papillomaviridae', 4, 'family', genome_type='dsDNA'), # 3%
            TaxonTarget('Polyomaviridae', 3, 'family', genome_type='dsDNA'),  # 1%
        ],
        abundance_model='power_law',
        dominant_tier_count=3,
        common_tier_count=20,
        moderate_tier_count=50,
        rare_tier_count=60,
        very_rare_tier_count=17
    ),

    'respiratory': CollectionSpec(
        collection_id='respiratory_virome_nasopharynx_healthy',
        name='Respiratory Virome - Nasopharynx (Healthy)',
        description='Respiratory virome from nasopharyngeal swabs. Mixed phage and eukaryotic viral composition.',
        body_site='respiratory',
        sample_type='nasopharyngeal_swab',
        health_status='healthy',
        target_size=200,
        phage_fraction=0.875,
        composition=[
            TaxonTarget('Cutibacterium', 40, 'genus', host_filter='Cutibacterium'),  # 20%
            TaxonTarget('Siphoviridae', 50, 'family'),  # 25%
            TaxonTarget('Myoviridae', 30, 'family'),    # 15%
            TaxonTarget('Podoviridae', 25, 'family'),   # 12.5%
            TaxonTarget('Microviridae', 20, 'family'),  # 10%
            TaxonTarget('Adenoviridae', 10, 'family', genome_type='dsDNA'),  # 5%
            TaxonTarget('Herpesviridae', 8, 'family', genome_type='dsDNA'),  # 4%
            TaxonTarget('Orthomyxoviridae', 7, 'family', genome_type='ssRNA'), # 3.5%
            TaxonTarget('Anelloviridae', 5, 'family', genome_type='ssDNA'),    # 2.5%
            TaxonTarget('Coronaviridae', 5, 'family', genome_type='ssRNA'),    # 2.5%
        ],
        abundance_model='log_normal',
        dominant_tier_count=4,
        common_tier_count=25,
        moderate_tier_count=70,
        rare_tier_count=85,
        very_rare_tier_count=16
    ),

    'marine': CollectionSpec(
        collection_id='marine_virome_coastal_surface',
        name='Marine Virome - Coastal Surface Water',
        description='Marine virome from coastal surface waters. Pelagiphage and cyanophage-dominated composition.',
        body_site='marine',
        sample_type='seawater',
        health_status='not_applicable',
        target_size=500,
        phage_fraction=0.97,
        composition=[
            TaxonTarget('Caudoviricetes', 400, 'class'),  # 80% - Pelagiphages, Cyanophages
            TaxonTarget('Microviridae', 50, 'family'),    # 10%
            TaxonTarget('Phycodnaviridae', 30, 'family', genome_type='dsDNA'), # 6%
            TaxonTarget('Other', 20, 'any'),  # 4%
        ],
        abundance_model='power_law',
        dominant_tier_count=10,
        common_tier_count=80,
        moderate_tier_count=200,
        rare_tier_count=180,
        very_rare_tier_count=30
    ),

    'soil': CollectionSpec(
        collection_id='soil_virome_agricultural',
        name='Soil Virome - Agricultural',
        description='Soil virome from agricultural fields. Actinophage-dominated with high diversity.',
        body_site='soil',
        sample_type='soil_extract',
        health_status='not_applicable',
        target_size=300,
        phage_fraction=0.94,
        composition=[
            TaxonTarget('Caudoviricetes', 240, 'class'),  # 80% - Actinophages, diverse hosts
            TaxonTarget('Microviridae', 30, 'family'),    # 10%
            TaxonTarget('Inoviridae', 12, 'family'),      # 4%
            TaxonTarget('Other', 18, 'any'),  # 6%
        ],
        abundance_model='log_normal',
        dominant_tier_count=5,
        common_tier_count=50,
        moderate_tier_count=120,
        rare_tier_count=105,
        very_rare_tier_count=20
    ),

    'freshwater': CollectionSpec(
        collection_id='freshwater_virome_lake_surface',
        name='Freshwater Virome - Lake Surface Water',
        description='Freshwater virome from lake surface waters. Actinophage and cyanophage-dominated.',
        body_site='freshwater',
        sample_type='lake_water',
        health_status='not_applicable',
        target_size=200,
        phage_fraction=0.95,
        composition=[
            TaxonTarget('Caudoviricetes', 160, 'class'),  # 80%
            TaxonTarget('Microviridae', 20, 'family'),    # 10%
            TaxonTarget('Inoviridae', 10, 'family'),      # 5%
            TaxonTarget('Other', 10, 'any'),  # 5%
        ],
        abundance_model='log_normal',
        dominant_tier_count=4,
        common_tier_count=30,
        moderate_tier_count=80,
        rare_tier_count=70,
        very_rare_tier_count=16
    ),

    'mouse_gut': CollectionSpec(
        collection_id='mouse_gut_virome_lab',
        name='Mouse Gut Virome - Laboratory (C57BL/6)',
        description='Mouse gut virome from laboratory C57BL/6 mice. Lactobacillus phage-dominated.',
        body_site='gut',
        sample_type='feces',
        health_status='healthy',
        target_size=150,
        phage_fraction=0.90,
        composition=[
            TaxonTarget('Lactobacillus', 45, 'genus', host_filter='Lactobacillus'),  # 30%
            TaxonTarget('Siphoviridae', 45, 'family'),  # 30%
            TaxonTarget('Myoviridae', 23, 'family'),    # 15%
            TaxonTarget('Podoviridae', 15, 'family'),   # 10%
            TaxonTarget('Microviridae', 7, 'family'),   # 5%
            TaxonTarget('Murine_viruses', 15, 'any', host_filter='Mus'),  # 10%
        ],
        abundance_model='log_normal',
        dominant_tier_count=3,
        common_tier_count=20,
        moderate_tier_count=50,
        rare_tier_count=60,
        very_rare_tier_count=17
    ),
}


class BodySiteCurator:
    """Curates body site-specific virome collections from RefSeq database."""

    def __init__(self, db_path: str, random_seed: int = 42):
        self.db_path = Path(db_path)
        self.random_seed = random_seed
        random.seed(random_seed)
        np.random.seed(random_seed)

        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

    def curate_collection(self, spec: CollectionSpec) -> Dict:
        """
        Curate a body site collection based on specification.

        Returns:
            Dictionary with collection statistics and genome IDs
        """
        logger.info(f"Curating collection: {spec.name}")
        logger.info(f"  Target size: {spec.target_size} genomes")
        logger.info(f"  Body site: {spec.body_site}")
        logger.info(f"  Phage fraction: {spec.phage_fraction:.1%}")

        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        selected_genomes = []

        # Process each taxon target
        for target in spec.composition:
            logger.info(f"  Selecting {target.target_count} genomes for {target.taxon_name}...")

            genomes = self._select_genomes_for_taxon(
                conn, target, spec
            )

            # Randomly select target count
            if len(genomes) > target.target_count:
                genomes = random.sample(genomes, target.target_count)

            logger.info(f"    Found {len(genomes)} genomes")
            selected_genomes.extend(genomes)

        conn.close()

        # Generate abundance distribution
        abundances = self._generate_abundances(spec)

        # Truncate abundances to match actual genome count (in case we didn't reach target)
        if len(selected_genomes) < len(abundances):
            abundances = abundances[:len(selected_genomes)]
            # Renormalize
            abundances = abundances / abundances.sum()

        # Combine genomes with abundances
        collection = {
            'spec': spec,
            'genomes': selected_genomes,
            'abundances': abundances,
            'total_count': len(selected_genomes)
        }

        logger.info(f"✓ Collection curated: {len(selected_genomes)} genomes")

        return collection

    def _select_genomes_for_taxon(
        self,
        conn: sqlite3.Connection,
        target: TaxonTarget,
        spec: CollectionSpec
    ) -> List[str]:
        """Select genomes matching taxon criteria."""

        if target.taxon_name == 'Other' or target.rank == 'any':
            # Select random genomes not already selected
            query = """
                SELECT genome_id FROM genomes
                WHERE genome_id NOT IN (
                    SELECT genome_id FROM collection_genomes
                    WHERE collection_id = ?
                )
                ORDER BY RANDOM()
                LIMIT ?
            """
            cursor = conn.execute(query, (spec.collection_id, target.target_count))
            return [row['genome_id'] for row in cursor.fetchall()]

        # Build query based on rank
        if target.rank == 'order':
            taxonomy_field = 'order_name'
        elif target.rank == 'family':
            taxonomy_field = 'family'
        elif target.rank == 'genus':
            taxonomy_field = 'genus'
        elif target.rank == 'species':
            taxonomy_field = 'species'
        elif target.rank == 'class':
            taxonomy_field = 'class'
        else:
            raise ValueError(f"Unknown rank: {target.rank}")

        query = f"""
            SELECT g.genome_id, g.genome_name, t.{taxonomy_field}, g.genome_type
            FROM genomes g
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE t.{taxonomy_field} LIKE ?
        """
        params = [f"%{target.taxon_name}%"]

        # Add genome type filter if specified
        if target.genome_type:
            query += " AND g.genome_type = ?"
            params.append(target.genome_type)

        # Add host filter if specified
        if target.host_filter:
            query += """
                AND EXISTS (
                    SELECT 1 FROM host_associations h
                    WHERE h.genome_id = g.genome_id
                    AND h.host_name LIKE ?
                )
            """
            params.append(f"%{target.host_filter}%")

        query += " ORDER BY RANDOM()"

        cursor = conn.execute(query, params)
        genomes = [row['genome_id'] for row in cursor.fetchall()]

        # If not enough genomes found, try broader search
        if len(genomes) < target.target_count:
            logger.warning(f"    Only found {len(genomes)}/{target.target_count} for {target.taxon_name}")
            logger.warning(f"    Trying broader search...")

            # Try without host filter
            if target.host_filter:
                query_broad = f"""
                    SELECT g.genome_id
                    FROM genomes g
                    LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
                    WHERE t.{taxonomy_field} LIKE ?
                """
                params_broad = [f"%{target.taxon_name}%"]

                if target.genome_type:
                    query_broad += " AND g.genome_type = ?"
                    params_broad.append(target.genome_type)

                query_broad += " ORDER BY RANDOM()"

                cursor = conn.execute(query_broad, params_broad)
                genomes = [row['genome_id'] for row in cursor.fetchall()]

        return genomes

    def _generate_abundances(self, spec: CollectionSpec) -> np.ndarray:
        """
        Generate realistic abundance distribution for collection.

        Uses tiered abundance model:
        - Tier 1 (Dominant): 10-30% each
        - Tier 2 (Common): 1-10% each
        - Tier 3 (Moderate): 0.1-1% each
        - Tier 4 (Rare): 0.01-0.1% each
        - Tier 5 (Very rare): <0.01% each
        """
        abundances = []

        # Tier 1: Dominant (10-30%)
        for _ in range(spec.dominant_tier_count):
            abundances.append(np.random.uniform(0.10, 0.30))

        # Tier 2: Common (1-10%)
        for _ in range(spec.common_tier_count):
            abundances.append(np.random.uniform(0.01, 0.10))

        # Tier 3: Moderate (0.1-1%)
        for _ in range(spec.moderate_tier_count):
            abundances.append(np.random.uniform(0.001, 0.01))

        # Tier 4: Rare (0.01-0.1%)
        for _ in range(spec.rare_tier_count):
            abundances.append(np.random.uniform(0.0001, 0.001))

        # Tier 5: Very rare (<0.01%)
        for _ in range(spec.very_rare_tier_count):
            abundances.append(np.random.uniform(0.00001, 0.0001))

        # Normalize to sum to 1.0
        abundances = np.array(abundances)
        abundances = abundances / abundances.sum()

        # Shuffle to randomize order
        np.random.shuffle(abundances)

        return abundances

    def save_collection(self, collection: Dict, output_dir: Path):
        """Save curated collection to database and output files."""
        spec = collection['spec']
        genomes = collection['genomes']
        abundances = collection['abundances']

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        conn = sqlite3.connect(self.db_path)

        # Insert collection metadata
        logger.info(f"Saving collection to database...")

        # Build selection criteria description
        selection_criteria = f"Body site: {spec.body_site}, Sample type: {spec.sample_type}, Health status: {spec.health_status}, Phage fraction: {spec.phage_fraction:.1%}"

        cursor = conn.execute("""
            INSERT OR REPLACE INTO body_site_collections
            (collection_name, description, n_genomes, selection_criteria,
             curated_by, curation_date)
            VALUES (?, ?, ?, ?, ?, datetime('now'))
        """, (
            spec.name,
            spec.description,
            len(genomes),
            selection_criteria,
            'ViroForge automated curation'
        ))

        # Get the collection_id of the inserted collection
        collection_id = cursor.lastrowid

        # Insert genome associations with abundances
        for rank, (genome_id, abundance) in enumerate(zip(genomes, abundances), start=1):
            conn.execute("""
                INSERT OR REPLACE INTO collection_genomes
                (collection_id, genome_id, relative_abundance, abundance_rank)
                VALUES (?, ?, ?, ?)
            """, (
                collection_id,
                genome_id,
                float(abundance),
                rank
            ))

        conn.commit()

        # Export to TSV
        tsv_path = output_dir / f"{spec.collection_id}_composition.tsv"
        logger.info(f"Exporting to: {tsv_path}")

        with open(tsv_path, 'w') as f:
            f.write("genome_id\tabundance\trank\n")

            # Sort by abundance
            sorted_indices = np.argsort(abundances)[::-1]

            for rank, idx in enumerate(sorted_indices, 1):
                genome_id = genomes[idx]
                abundance = abundances[idx]
                f.write(f"{genome_id}\t{abundance:.8f}\t{rank}\n")

        # Generate summary statistics
        summary_path = output_dir / f"{spec.collection_id}_summary.txt"
        self._write_summary(collection, summary_path)

        conn.close()

        logger.info(f"✓ Collection saved successfully")

    def _write_summary(self, collection: Dict, output_path: Path):
        """Write collection summary statistics."""
        spec = collection['spec']
        abundances = collection['abundances']

        with open(output_path, 'w') as f:
            f.write(f"Body Site Virome Collection Summary\n")
            f.write(f"=" * 70 + "\n\n")
            f.write(f"Collection: {spec.name}\n")
            f.write(f"ID: {spec.collection_id}\n")
            f.write(f"Body Site: {spec.body_site}\n")
            f.write(f"Sample Type: {spec.sample_type}\n")
            f.write(f"Health Status: {spec.health_status}\n\n")

            f.write(f"Composition:\n")
            f.write(f"  Total genomes: {spec.target_size}\n")
            f.write(f"  Phage fraction: {spec.phage_fraction:.1%}\n\n")

            f.write(f"Taxonomic Targets:\n")
            for target in spec.composition:
                f.write(f"  - {target.taxon_name}: {target.target_count} genomes ")
                f.write(f"({target.target_count/spec.target_size:.1%})\n")

            f.write(f"\nAbundance Distribution:\n")
            f.write(f"  Model: {spec.abundance_model}\n")
            f.write(f"  Dominant tier (10-30%): {spec.dominant_tier_count} genomes\n")
            f.write(f"  Common tier (1-10%): {spec.common_tier_count} genomes\n")
            f.write(f"  Moderate tier (0.1-1%): {spec.moderate_tier_count} genomes\n")
            f.write(f"  Rare tier (0.01-0.1%): {spec.rare_tier_count} genomes\n")
            f.write(f"  Very rare tier (<0.01%): {spec.very_rare_tier_count} genomes\n\n")

            f.write(f"Statistics:\n")
            f.write(f"  Max abundance: {abundances.max():.4%}\n")
            f.write(f"  Mean abundance: {abundances.mean():.4%}\n")
            f.write(f"  Median abundance: {np.median(abundances):.4%}\n")
            f.write(f"  Min abundance: {abundances.min():.6%}\n")
            f.write(f"  Shannon diversity: {self._calculate_shannon(abundances):.2f}\n")

    def _calculate_shannon(self, abundances: np.ndarray) -> float:
        """Calculate Shannon diversity index."""
        # Remove zeros
        abundances = abundances[abundances > 0]
        return -np.sum(abundances * np.log(abundances))


def main():
    parser = argparse.ArgumentParser(
        description='Curate body site-specific virome collections'
    )
    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to viral genomes database'
    )
    parser.add_argument(
        '--collection',
        choices=list(COLLECTION_SPECS.keys()),
        help='Collection to curate (gut, oral, skin, respiratory, marine, soil, freshwater, mouse_gut)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Curate all collections'
    )
    parser.add_argument(
        '--output',
        default='data/body_site_collections',
        help='Output directory for collection files'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )

    args = parser.parse_args()

    if not args.collection and not args.all:
        parser.error("Must specify either --collection or --all")

    # Create curator
    curator = BodySiteCurator(args.database, random_seed=args.seed)

    # Determine which collections to curate
    if args.all:
        collections_to_curate = list(COLLECTION_SPECS.keys())
    else:
        collections_to_curate = [args.collection]

    logger.info(f"Curating {len(collections_to_curate)} collection(s)")
    logger.info(f"Output directory: {args.output}")

    # Curate each collection
    for collection_name in collections_to_curate:
        spec = COLLECTION_SPECS[collection_name]

        try:
            collection = curator.curate_collection(spec)
            curator.save_collection(collection, Path(args.output))
            logger.info(f"✓ Successfully curated: {collection_name}\n")

        except Exception as e:
            logger.error(f"✗ Failed to curate {collection_name}: {e}")
            continue

    logger.info(f"Curation complete!")


if __name__ == '__main__':
    main()
