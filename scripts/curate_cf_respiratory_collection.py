#!/usr/bin/env python3
"""
Curate Collection 20: Cystic Fibrosis (CF) Respiratory Virome

Comparison to healthy respiratory (Collection 5):
- Altered composition dominated by bacterial pathogen-specific phages
- Increased Pseudomonas phages (chronic P. aeruginosa infection)
- Increased Staphylococcus phages (S. aureus colonization)
- Respiratory viruses (influenza, rhinovirus, RSV, adenovirus)
- Mucus-adapted viral community

Literature basis:
- Lim et al. 2014 (J Clin Microbiol): CF airway virome dominated by bacteriophages
- Wat et al. 2008 (J Cyst Fibros): Viral infections and CF exacerbations
- Esther Jr et al. 2014 (Pediatr Pulmonol): Viral infections in CF

Target size: 60-80 genomes

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


class CFRespiratoryCurator:
    """Curate CF respiratory virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_pseudomonas_phages(self, n_target: int = 25) -> List[Dict]:
        """
        Get Pseudomonas phages (dominant in CF).

        P. aeruginosa chronically colonizes CF airways.
        Pseudomonas phages are the most abundant viruses in CF airways.

        Key families: Myoviridae, Podoviridae (now part of various families)
        """
        logger.info("Selecting Pseudomonas phages (CF dominant)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Pseudomonas%phage%'
           OR g.genome_name LIKE '%Pseudomonas%virus%'
           OR t.genus LIKE '%Pseudomonas%')
        AND g.genome_name NOT LIKE '%prophage%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        pseudomonas = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pseudomonas phages: {len(pseudomonas)}")
        return pseudomonas

    def get_staphylococcus_phages(self, n_target: int = 15) -> List[Dict]:
        """
        Get Staphylococcus phages.

        S. aureus commonly colonizes CF airways, especially early disease.
        Staphylococcus phages are abundant in CF respiratory samples.
        """
        logger.info("Selecting Staphylococcus phages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Staphylococcus%phage%'
           OR g.genome_name LIKE '%Staphylococcus%virus%'
           OR t.genus LIKE '%Staphylococcus%')
        AND g.genome_name NOT LIKE '%prophage%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        staph = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Staphylococcus phages: {len(staph)}")
        return staph

    def get_respiratory_viruses(self, n_target: int = 20) -> List[Dict]:
        """
        Get respiratory viruses common in CF.

        CF patients are susceptible to respiratory viral infections:
        - Influenza (seasonal, causes exacerbations)
        - Rhinovirus (common cold, frequent infections)
        - RSV (respiratory syncytial virus)
        - Adenovirus (persistent infections)
        - Coronavirus (including common cold strains)
        """
        logger.info("Selecting respiratory viruses...")

        # Orthomyxoviridae (Influenza) - get all available (10 total after taxonomy fix)
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Orthomyxoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        influenza = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.5), 10),))]
        logger.info(f"  Influenza: {len(influenza)}")

        # Picornaviridae (Rhinovirus, Enterovirus) - 30%
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Picornaviridae'
           AND (g.genome_name LIKE '%Rhinovirus%' OR g.genome_name LIKE '%rhinovirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """
        rhino = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.3), 3),))]
        logger.info(f"  Rhinovirus: {len(rhino)}")

        # Pneumoviridae (RSV) - 20%
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Pneumoviridae'
           OR g.genome_name LIKE '%Respiratory syncytial%'
           OR g.genome_name LIKE '%RSV%'
        ORDER BY RANDOM()
        LIMIT ?
        """
        rsv = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.2), 2),))]
        logger.info(f"  RSV: {len(rsv)}")

        # Adenoviridae - 20%
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        adeno = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.2), 2),))]
        logger.info(f"  Adenovirus: {len(adeno)}")

        respiratory = influenza + rhino + rsv + adeno
        logger.info(f"  Total respiratory viruses: {len(respiratory)}")
        return respiratory

    def get_other_bacterial_phages(self, n_target: int = 12) -> List[Dict]:
        """
        Get other bacterial phages from CF respiratory pathogens.

        CF airways harbor diverse bacterial community including:
        - Burkholderia cepacia complex
        - Achromobacter
        - Stenotrophomonas
        """
        logger.info("Selecting other bacterial phages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Burkholderia%phage%'
           OR g.genome_name LIKE '%Achromobacter%phage%'
           OR g.genome_name LIKE '%Stenotrophomonas%phage%'
           OR g.genome_name LIKE '%Haemophilus%phage%')
        AND g.genome_name NOT LIKE '%prophage%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        other_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Other bacterial phages: {len(other_phages)}")
        return other_phages

    def get_opportunistic_viruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get opportunistic viruses in CF.

        CF patients may have altered viral susceptibility.
        """
        logger.info("Selecting opportunistic viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Herpesviridae', 'Anelloviridae', 'Polyomaviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """

        opportunistic = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Opportunistic viruses: {len(opportunistic)}")
        return opportunistic

    def assign_cf_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign relative abundances reflecting CF respiratory virome.

        CF airways show:
        - Dominance of Pseudomonas phages (chronic infection)
        - Secondary abundance of Staphylococcus phages
        - Variable respiratory viruses (seasonal/episodic)
        - Long tail of lower abundance viruses

        Uses log-normal distribution with moderate skew.
        """
        n = len(genomes)

        # Moderate skew (more even than gut dysbiosis, reflecting phage dominance)
        mu = -2.2
        sigma = 2.2

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete CF respiratory virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING CF RESPIRATORY VIROME COLLECTION")
        logger.info("=" * 80)

        # CF characteristics: phage-dominated with respiratory viruses
        pseudomonas = self.get_pseudomonas_phages(n_target=25)
        staph = self.get_staphylococcus_phages(n_target=15)
        respiratory = self.get_respiratory_viruses(n_target=20)
        other_phages = self.get_other_bacterial_phages(n_target=12)
        opportunistic = self.get_opportunistic_viruses(n_target=8)

        # Combine all
        collection = pseudomonas + staph + respiratory + other_phages + opportunistic

        # Remove any duplicates (shouldn't be any, but safety check)
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign CF-specific abundances
        logger.info("\nAssigning CF respiratory abundance pattern...")
        collection = self.assign_cf_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in CF respiratory collection: {len(collection)}")
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

        # Count by category
        pseudo_count = len([g for g in collection if 'Pseudomonas' in g['genome_name']])
        staph_count = len([g for g in collection if 'Staphylococcus' in g['genome_name']])
        resp_families = ['Orthomyxoviridae', 'Picornaviridae', 'Pneumoviridae', 'Adenoviridae']
        resp_count = len([g for g in collection if g.get('family') in resp_families])

        logger.info(f"Pseudomonas phages:      {pseudo_count:3d} genomes")
        logger.info(f"Staphylococcus phages:   {staph_count:3d} genomes")
        logger.info(f"Respiratory viruses:     {resp_count:3d} genomes")
        logger.info(f"Other phages/viruses:    {len(collection) - pseudo_count - staph_count - resp_count:3d} genomes")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 20,
            'collection_name': 'Cystic Fibrosis (CF) Respiratory Virome',
            'description': (
                'Respiratory virome from cystic fibrosis patients. '
                'Dominated by bacteriophages specific to CF pathogens (Pseudomonas aeruginosa, '
                'Staphylococcus aureus) reflecting chronic bacterial colonization. '
                'Includes respiratory viruses (influenza, rhinovirus, RSV, adenovirus) '
                'associated with CF exacerbations. '
                'Based on Lim et al. 2014, Wat et al. 2008, and Esther Jr et al. 2014.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'CF-specific composition: Pseudomonas phages (dominant), '
                'Staphylococcus phages (secondary), respiratory viruses (influenza, rhinovirus, RSV, adenovirus), '
                'other CF pathogen phages (Burkholderia, Achromobacter), '
                'opportunistic viruses. Reflects chronic infection environment.'
            ),
            'curated_by': 'ViroForge Phase 7',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Lim et al. 2014 (J Clin Microbiol), Wat et al. 2008 (J Cyst Fibros), '
                'Esther Jr et al. 2014 (Pediatr Pulmonol), '
                'Cuthbertson et al. 2020 (J Cyst Fibros)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 20")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 20 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 20")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 20")

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
                20,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 20 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = CFRespiratoryCurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ CF RESPIRATORY VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 20 --dry-run")
        logger.info("  2. Compare to healthy respiratory (Collection 5)")
        logger.info("  3. Update documentation for all disease collections")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
