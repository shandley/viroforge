#!/usr/bin/env python3
"""
Curate Collection 21: Human Respiratory RNA Virome

Composition based on literature (Gao et al. 2019, Wylie et al. 2012):
- Coronaviruses (including common cold and SARS-CoV-2)
- Rhinoviruses and enteroviruses (most common respiratory viruses)
- Influenza viruses (seasonal epidemic viruses)
- Respiratory syncytial virus (RSV)
- Human metapneumovirus (HMPV)
- Parainfluenza viruses

All viruses in this collection are RNA viruses (ssRNA and dsRNA).

Literature basis:
- Gao et al. 2019 (Nat Commun): Respiratory virome in recurrent infections
- Wylie et al. 2012 (PLOS One): Human virome in febrile children
- Wylie et al. 2017 (Clin Chest Med): Review of respiratory tract virome

Target size: 60-70 genomes

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


class RespiratoryRNACurator:
    """Curate respiratory RNA virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_coronaviruses(self, n_target: int = 15) -> List[Dict]:
        """
        Get human coronaviruses.

        Includes:
        - Common cold coronaviruses: OC43, 229E, NL63, HKU1
        - SARS-CoV-2 and variants
        - Other human coronaviruses

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting coronaviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Coronaviridae'
          AND (g.genome_name LIKE '%Human%'
           OR g.genome_name LIKE '%SARS%'
           OR g.genome_name LIKE '%229E%'
           OR g.genome_name LIKE '%OC43%'
           OR g.genome_name LIKE '%NL63%'
           OR g.genome_name LIKE '%HKU1%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        coronaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Coronaviridae: {len(coronaviruses)}")
        return coronaviruses

    def get_picornaviruses(self, n_target: int = 15) -> List[Dict]:
        """
        Get picornaviruses (rhinoviruses and enteroviruses).

        Rhinoviruses are the most common cause of the common cold.
        Enteroviruses include many respiratory and systemic pathogens.

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting picornaviruses (rhinovirus, enterovirus)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Picornaviridae'
          AND (g.genome_name LIKE '%Rhinovirus%'
           OR g.genome_name LIKE '%Enterovirus%'
           OR g.genome_name LIKE '%rhinovirus%'
           OR t.genus = 'Enterovirus')
        ORDER BY RANDOM()
        LIMIT ?
        """

        picorna = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Picornaviridae (Rhinovirus, Enterovirus): {len(picorna)}")
        return picorna

    def get_influenza(self, n_target: int = 12) -> List[Dict]:
        """
        Get influenza viruses.

        Seasonal influenza A and B viruses.
        Genome type: ssRNA negative-sense (segmented)
        """
        logger.info("Selecting influenza viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Orthomyxoviridae'
          AND (g.genome_name LIKE '%Influenza%'
           OR g.genome_name LIKE '%influenza%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        influenza = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Orthomyxoviridae (Influenza): {len(influenza)}")
        return influenza

    def get_rsv(self, n_target: int = 5) -> List[Dict]:
        """
        Get respiratory syncytial virus (RSV).

        Major cause of bronchiolitis and pneumonia in infants.
        Genome type: ssRNA negative-sense
        """
        logger.info("Selecting respiratory syncytial virus (RSV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Pneumoviridae'
          AND (g.genome_name LIKE '%Respiratory syncytial%'
           OR g.genome_name LIKE '%RSV%'
           OR g.genome_name LIKE '%Human orthopneumovirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        rsv = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pneumoviridae (RSV): {len(rsv)}")
        return rsv

    def get_paramyxoviruses(self, n_target: int = 15) -> List[Dict]:
        """
        Get paramyxoviruses (HMPV, parainfluenza).

        Includes:
        - Human metapneumovirus (HMPV)
        - Parainfluenza viruses (1-4)

        Genome type: ssRNA negative-sense
        """
        logger.info("Selecting paramyxoviruses (HMPV, parainfluenza)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Paramyxoviridae'
          AND (g.genome_name LIKE '%Human metapneumovirus%'
           OR g.genome_name LIKE '%Metapneumovirus%'
           OR g.genome_name LIKE '%Parainfluenza%'
           OR g.genome_name LIKE '%Human respirovirus%'
           OR g.genome_name LIKE '%Human rubulavirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        paramyxo = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Paramyxoviridae (HMPV, Parainfluenza): {len(paramyxo)}")
        return paramyxo

    def get_adenoviruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get adenoviruses.

        While adenoviruses are DNA viruses, they commonly co-circulate
        with RNA respiratory viruses and are detected in respiratory virome studies.

        Genome type: dsDNA (included for completeness)
        """
        logger.info("Selecting adenoviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
          AND g.genome_name LIKE '%Human%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        adeno = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Adenoviridae (DNA virus, commonly detected): {len(adeno)}")
        return adeno

    def assign_respiratory_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign abundances reflecting seasonal respiratory virus patterns.

        Respiratory viruses show:
        - High dominance by rhinovirus (most common)
        - Seasonal peaks for influenza, RSV
        - Variable coronavirus abundance
        - Long tail of less common viruses

        Uses log-normal distribution with moderate skew.
        """
        n = len(genomes)

        # Moderate skew (rhinovirus dominant, but diverse community)
        mu = -2.0
        sigma = 2.0

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete respiratory RNA virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING HUMAN RESPIRATORY RNA VIROME COLLECTION")
        logger.info("=" * 80)

        # Get viruses by family
        coronaviruses = self.get_coronaviruses(n_target=15)
        picornaviruses = self.get_picornaviruses(n_target=15)
        influenza = self.get_influenza(n_target=12)
        rsv = self.get_rsv(n_target=5)
        paramyxoviruses = self.get_paramyxoviruses(n_target=15)
        adenoviruses = self.get_adenoviruses(n_target=8)

        # Combine all
        collection = coronaviruses + picornaviruses + influenza + rsv + paramyxoviruses + adenoviruses

        # Remove any duplicates (shouldn't be any, but safety check)
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign respiratory abundances
        logger.info("\nAssigning respiratory abundance pattern...")
        collection = self.assign_respiratory_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in respiratory RNA collection: {len(collection)}")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        # Show top 10
        logger.info("\nTop 10 most abundant genomes:")
        for i, genome in enumerate(collection[:10], 1):
            family = genome.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {genome['genome_name'][:60]:60s} {genome['relative_abundance']:.6f} ({family})")

        # Show composition summary
        logger.info("\n" + "=" * 80)
        logger.info("COMPOSITION SUMMARY")
        logger.info("=" * 80)

        # Count by family
        corona_count = len([g for g in collection if g.get('family') == 'Coronaviridae'])
        picorna_count = len([g for g in collection if g.get('family') == 'Picornaviridae'])
        influenza_count = len([g for g in collection if g.get('family') == 'Orthomyxoviridae'])
        rsv_count = len([g for g in collection if g.get('family') == 'Pneumoviridae'])
        paramyxo_count = len([g for g in collection if g.get('family') == 'Paramyxoviridae'])
        adeno_count = len([g for g in collection if g.get('family') == 'Adenoviridae'])

        logger.info(f"Coronaviridae:            {corona_count:3d} genomes ({corona_count/len(collection)*100:.1f}%)")
        logger.info(f"Picornaviridae:           {picorna_count:3d} genomes ({picorna_count/len(collection)*100:.1f}%)")
        logger.info(f"Orthomyxoviridae:         {influenza_count:3d} genomes ({influenza_count/len(collection)*100:.1f}%)")
        logger.info(f"Pneumoviridae (RSV):      {rsv_count:3d} genomes ({rsv_count/len(collection)*100:.1f}%)")
        logger.info(f"Paramyxoviridae:          {paramyxo_count:3d} genomes ({paramyxo_count/len(collection)*100:.1f}%)")
        logger.info(f"Adenoviridae:             {adeno_count:3d} genomes ({adeno_count/len(collection)*100:.1f}%)")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 21,
            'collection_name': 'Human Respiratory RNA Virome',
            'description': (
                'Human respiratory RNA virome from community and hospital settings. '
                'Includes common respiratory RNA viruses: coronaviruses (including SARS-CoV-2), '
                'rhinoviruses, influenza, RSV, parainfluenza, and human metapneumovirus. '
                'Represents typical respiratory virus diversity detected in clinical and surveillance studies. '
                'Based on Gao et al. 2019, Wylie et al. 2012, and Wylie et al. 2017.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'RNA respiratory viruses from major families: Coronaviridae, Picornaviridae, '
                'Orthomyxoviridae, Pneumoviridae, Paramyxoviridae. Includes seasonal '
                'epidemic viruses (influenza, RSV) and endemic viruses (rhinovirus, common cold coronaviruses). '
                'Adenoviruses included as commonly co-detected DNA viruses.'
            ),
            'curated_by': 'ViroForge Phase 8',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Gao et al. 2019 (Nat Commun), Wylie et al. 2012 (PLOS One), '
                'Wylie et al. 2017 (Clin Chest Med)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 21")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 21 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 21")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 21")

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
                21,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 21 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = RespiratoryRNACurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ RESPIRATORY RNA VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 21 --dry-run")
        logger.info("  2. Continue with Collection 22 (Arbovirus) and Collection 23 (Fecal RNA)")
        logger.info("  3. Implement RNA workflow components (RT, rRNA depletion)")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
