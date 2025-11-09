#!/usr/bin/env python3
"""
Curate Collection 18: IBD Gut Virome (Inflammatory Bowel Disease)

Comparison to healthy gut (Collection 9):
- Lower diversity (80-100 genomes vs 134)
- Altered Caudovirales composition
- Increased temperate phages
- Reduced crAssphage diversity
- Different abundance patterns (dysbiosis)

Literature basis:
- Norman et al. 2015 (Cell): IBD gut virome shows reduced diversity
- Zuo et al. 2019 (Gut): Altered phage-bacteria dynamics in IBD
- Clooney et al. 2019 (Cell Host Microbe): Crohn's disease virome expansion

Target size: 80-100 genomes

Author: ViroForge Development Team
Date: 2025-11-09
"""

import sqlite3
import numpy as np
from pathlib import Path
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class IBDGutCurator:
    """Curate IBD gut virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_healthy_gut_genomes(self) -> List[Dict]:
        """Get genomes from healthy gut collection (ID 9) for reference."""
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, cg.relative_abundance
        FROM collection_genomes cg
        JOIN genomes g ON cg.genome_id = g.genome_id
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE cg.collection_id = 9
        ORDER BY cg.relative_abundance DESC
        """
        return [dict(row) for row in self.conn.execute(query)]

    def get_reduced_crassphage(self, n_target: int = 20) -> List[Dict]:
        """
        Get reduced crAssphage diversity (IBD characteristic).

        Healthy has ~36 crAssphage genomes, IBD should have ~20 (44% reduction)
        reflecting loss of diversity documented in Norman et al. 2015.
        """
        logger.info("Selecting reduced crAssphage diversity...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Intestiviridae', 'Suoliviridae', 'Steigviridae', 'Crevaviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """

        crassphage = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  crAssphage-like: {len(crassphage)} (vs 36 in healthy)")
        return crassphage

    def get_temperate_phages(self, n_target: int = 35) -> List[Dict]:
        """
        Get temperate phages (increased in IBD).

        Siphoviridae and related families are often temperate (lysogenic).
        IBD shows increased temperate phage activity (Clooney et al. 2019).
        """
        logger.info("Selecting temperate phages (increased in IBD)...")

        # Get various gut phage families, biased toward temperate types
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Ackermannviridae', 'Drexlerviridae', 'Demerecviridae')
           OR (t.genus LIKE '%phage%' AND g.genome_name LIKE '%Lactobacillus%')
           OR (t.genus LIKE '%phage%' AND g.genome_name LIKE '%Enterococcus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        temperate = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Temperate phages: {len(temperate)}")
        return temperate

    def get_microviridae(self, n_target: int = 15) -> List[Dict]:
        """
        Get Microviridae (coliphage).

        Maintained in IBD but potentially altered composition.
        """
        logger.info("Selecting Microviridae...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Microviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """

        microvir = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Microviridae: {len(microvir)}")
        return microvir

    def get_inoviridae(self, n_target: int = 12) -> List[Dict]:
        """
        Get Inoviridae (filamentous phages).

        Present but reduced diversity in IBD.
        """
        logger.info("Selecting Inoviridae...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Inoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """

        inovir = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Inoviridae: {len(inovir)}")
        return inovir

    def get_eukaryotic_viruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get eukaryotic viruses (potentially increased in IBD).

        Some studies show increased eukaryotic viral load in IBD,
        possibly due to altered immune function or mucosal damage.
        """
        logger.info("Selecting eukaryotic viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Adenoviridae', 'Picornaviridae', 'Anelloviridae', 'Parvoviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """

        eukaryotic = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Eukaryotic viruses: {len(eukaryotic)}")
        return eukaryotic

    def assign_dysbiotic_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign abundances reflecting dysbiotic state.

        IBD gut virome characteristics:
        - Less even distribution (higher dominance)
        - Some viral blooms (few very abundant species)
        - Overall lower diversity

        Uses more skewed log-normal distribution than healthy gut.
        """
        n = len(genomes)

        # More skewed than healthy gut
        # Healthy uses tiered random, IBD uses log-normal with stronger skew
        mu = -2.5  # Lower = more skewed
        sigma = 2.5  # Higher = more variable

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete IBD gut virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING IBD GUT VIROME COLLECTION")
        logger.info("=" * 80)

        # First, show healthy gut for comparison
        healthy = self.get_healthy_gut_genomes()
        logger.info(f"\nHealthy gut (Collection 9): {len(healthy)} genomes")

        # IBD characteristics: reduced diversity across all groups
        crassphage = self.get_reduced_crassphage(n_target=20)
        temperate = self.get_temperate_phages(n_target=35)
        microviridae = self.get_microviridae(n_target=15)
        inoviridae = self.get_inoviridae(n_target=12)
        eukaryotic = self.get_eukaryotic_viruses(n_target=8)

        # Combine all
        collection = crassphage + temperate + microviridae + inoviridae + eukaryotic

        # Remove any duplicates (shouldn't be any, but safety check)
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign dysbiotic abundances
        logger.info("\nAssigning dysbiotic abundance pattern...")
        collection = self.assign_dysbiotic_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in IBD collection: {len(collection)}")
        logger.info(f"Diversity reduction: {(1 - len(collection)/len(healthy)) * 100:.1f}% fewer genomes than healthy")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        # Show top 10
        logger.info("\nTop 10 most abundant genomes:")
        for i, genome in enumerate(collection[:10], 1):
            family = genome.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {genome['genome_name'][:50]:50s} {genome['relative_abundance']:.6f} ({family})")

        # Compare diversity to healthy
        logger.info("\n" + "=" * 80)
        logger.info("COMPARISON TO HEALTHY GUT (Collection 9)")
        logger.info("=" * 80)
        logger.info(f"Healthy:  {len(healthy)} genomes")
        logger.info(f"IBD:      {len(collection)} genomes ({len(collection)/len(healthy)*100:.1f}% of healthy)")
        logger.info(f"Reduction: {len(healthy) - len(collection)} genomes lost")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 18,
            'collection_name': 'IBD Gut Virome (Inflammatory Bowel Disease)',
            'description': (
                'Gut virome from patients with inflammatory bowel disease (IBD). '
                'Characterized by reduced viral diversity, altered Caudovirales composition, '
                'increased temperate phage activity, and dysbiotic abundance patterns. '
                'Compare to Collection 9 (healthy gut) to study disease effects. '
                'Based on Norman et al. 2015, Zuo et al. 2019, and Clooney et al. 2019.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'Reduced diversity vs healthy (80-100 genomes vs 134). '
                'Reduced crAssphage (20 vs 36). Increased temperate phages. '
                'Altered Microviridae and Inoviridae. Dysbiotic abundance distribution.'
            ),
            'curated_by': 'ViroForge Phase 7',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Norman et al. 2015 (Cell), Zuo et al. 2019 (Gut), '
                'Clooney et al. 2019 (Cell Host Microbe)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 18")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 18 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 18")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 18")

        # Insert collection
        cursor.execute("""
            INSERT INTO body_site_collections
            (collection_id, collection_name, description, n_genomes, selection_criteria,
             curated_by, curation_date, literature_references, version)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            collection_meta['collection_id'],
            collection_meta['collection_name'],
            collection_meta['description'],
            collection_meta['n_genomes'],
            collection_meta['selection_criteria'],
            collection_meta['curated_by'],
            collection_meta['curation_date'],
            collection_meta['literature_references'],
            collection_meta['version']
        ))

        logger.info(f"✓ Inserted collection metadata")

        # Insert genomes
        for genome in collection:
            cursor.execute("""
                INSERT INTO collection_genomes
                (collection_id, genome_id, relative_abundance, prevalence, abundance_rank)
                VALUES (?, ?, ?, ?, ?)
            """, (
                18,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 18 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = IBDGutCurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ IBD GUT VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 18 --dry-run")
        logger.info("  2. Compare to healthy: Generate both collections and compare compositions")
        logger.info("  3. Update documentation")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
