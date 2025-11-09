#!/usr/bin/env python3
"""
Curate Collection 17: Wastewater Virome - Urban Treatment Plant

Composition based on literature (Crits-Christoph et al. 2021, Crank et al. 2022):
- Human enteric viruses (40%): Norovirus, rotavirus, adenovirus, astrovirus, sapovirus
- Bacteriophages (35%): Gut bacterial phages (crAssphage, coliphage, etc.)
- Environmental viruses (15%): Plant viruses, insect viruses
- Emerging pathogens (10%): SARS-CoV-2, mpox, polio, etc.

Target size: ~400 genomes

Author: ViroForge Development Team
Date: 2025-11-09
"""

import sqlite3
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class WastewaterCurator:
    """Curate wastewater virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_enteric_viruses(self, n_target: int = 160) -> List[Dict]:
        """
        Get human enteric viruses (40% of collection).

        Families: Caliciviridae, Adenoviridae, Astrovir

idae, Picornaviridae

        Common in wastewater:
        - Norovirus GII (dominant)
        - Adenovirus (highly abundant)
        - Rotavirus (seasonal)
        - Astrovirus (moderate)
        - Enterovirus/Poliovirus (surveillance)
        """
        logger.info("Selecting human enteric viruses...")

        # Caliciviridae (Norovirus, Sapovirus) - ~30% of enteric
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Caliciviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        calici = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.3),))]
        logger.info(f"  Caliciviridae (Norovirus, Sapovirus): {len(calici)}")

        # Adenoviridae - ~35% of enteric (very stable in wastewater)
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        adeno = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.35),))]
        logger.info(f"  Adenoviridae: {len(adeno)}")

        # Astroviridae - ~10% of enteric
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Astroviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        astro = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.10), 5),))]
        logger.info(f"  Astroviridae: {len(astro)}")

        # Picornaviridae (Enterovirus, Poliovirus) - ~20% of enteric
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Picornaviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        picorna = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.20),))]
        logger.info(f"  Picornaviridae (Enterovirus, Polio): {len(picorna)}")

        # Sedoreoviridae (Rotavirus) - ~5% of enteric (seasonal, newly available after taxonomy fix)
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Sedoreoviridae'
           OR g.genome_name LIKE '%Rotavirus%'
           OR g.genome_name LIKE '%rotavirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """
        rotavirus = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.05), 5),))]
        logger.info(f"  Sedoreoviridae (Rotavirus): {len(rotavirus)}")

        enteric = calici + adeno + astro + picorna + rotavirus
        logger.info(f"  Total enteric viruses: {len(enteric)}")
        return enteric

    def get_bacteriophages(self, n_target: int = 140) -> List[Dict]:
        """
        Get gut-associated bacteriophages (35% of collection).

        Dominated by:
        - crAssphage (highly abundant in human gut, stable in wastewater)
        - Microviridae (coliphage)
        - Inoviridae (filamentous phages)
        """
        logger.info("Selecting bacteriophages...")

        # crAssphage variants - ~50% of phages (dominant in wastewater)
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%crAss%'
           OR t.genus LIKE '%crAss%'
           OR t.family IN ('Intestiviridae', 'Suoliviridae', 'Steigviridae', 'Crevaviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """
        crassphage = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.5),))]
        logger.info(f"  crAssphage-like: {len(crassphage)}")

        # Microviridae (coliphage) - ~30% of phages
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Microviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        microvir = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.3),))]
        logger.info(f"  Microviridae (coliphage): {len(microvir)}")

        # Inoviridae (filamentous phages) - ~20% of phages
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Inoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        inovir = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.2),))]
        logger.info(f"  Inoviridae: {len(inovir)}")

        phages = crassphage + microvir + inovir
        logger.info(f"  Total bacteriophages: {len(phages)}")
        return phages

    def get_environmental_viruses(self, n_target: int = 60) -> List[Dict]:
        """
        Get environmental viruses (15% of collection).

        Plant and insect viruses washed in from environment:
        - Plant viruses from agricultural runoff
        - Insect viruses from urban environment
        """
        logger.info("Selecting environmental viruses...")

        # Insect viruses - ~50%
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Baculoviridae', 'Iflaviridae', 'Dicistroviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """
        insect = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.5),))]
        logger.info(f"  Insect viruses: {len(insect)}")

        # Plant viruses - ~50%
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Virgaviridae', 'Tombusviridae', 'Bromoviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """
        plant = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.5),))]
        logger.info(f"  Plant viruses: {len(plant)}")

        env = insect + plant
        logger.info(f"  Total environmental viruses: {len(env)}")
        return env

    def get_emerging_pathogens(self, n_target: int = 40) -> List[Dict]:
        """
        Get emerging pathogens (10% of collection).

        Public health surveillance targets:
        - SARS-CoV-2 (COVID-19 surveillance)
        - Mpox virus (monkeypox)
        - Poliovirus (included in Picornaviridae above, but add specific)
        """
        logger.info("Selecting emerging pathogens...")

        # Coronaviridae (SARS-CoV-2) - ~60%
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Coronaviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        corona = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.6),))]
        logger.info(f"  Coronaviridae (SARS-CoV-2): {len(corona)}")

        # Poxviridae (Mpox) - ~40%
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Poxviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        pox = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.4),))]
        logger.info(f"  Poxviridae (Mpox): {len(pox)}")

        emerging = corona + pox
        logger.info(f"  Total emerging pathogens: {len(emerging)}")
        return emerging

    def assign_abundances(self, genomes: List[Dict], group_name: str, group_fraction: float) -> List[Dict]:
        """
        Assign relative abundances within a group using realistic distributions.

        Uses log-normal distribution to simulate:
        - Few dominant species (crAssphage, adenovirus, norovirus)
        - Many rare species (long tail)
        """
        n = len(genomes)

        # Generate log-normal abundances
        # Higher mu = more even, lower mu = more skewed
        mu = -2.0  # Creates strong dominance pattern
        sigma = 2.0

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to group fraction
        normalized = raw_abundances / raw_abundances.sum() * group_fraction

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)
            genome['group'] = group_name

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete wastewater virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING WASTEWATER VIROME COLLECTION")
        logger.info("=" * 80)

        # Get genomes for each category
        enteric = self.get_enteric_viruses(n_target=160)
        phages = self.get_bacteriophages(n_target=140)
        environmental = self.get_environmental_viruses(n_target=60)
        emerging = self.get_emerging_pathogens(n_target=40)

        # Assign abundances
        logger.info("\nAssigning relative abundances...")
        enteric = self.assign_abundances(enteric, 'enteric', 0.40)
        phages = self.assign_abundances(phages, 'phage', 0.35)
        environmental = self.assign_abundances(environmental, 'environmental', 0.15)
        emerging = self.assign_abundances(emerging, 'emerging', 0.10)

        # Combine all
        collection = enteric + phages + environmental + emerging

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in collection: {len(collection)}")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        # Show top 10
        logger.info("\nTop 10 most abundant genomes:")
        for i, genome in enumerate(collection[:10], 1):
            logger.info(f"  {i:2d}. {genome['genome_name'][:60]:60s} {genome['relative_abundance']:.6f} ({genome['group']})")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        # Insert collection metadata
        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 17,
            'collection_name': 'Wastewater Virome - Urban Treatment Plant',
            'description': (
                'Municipal wastewater viral community from urban treatment plant. '
                'Composition based on Crits-Christoph et al. 2021 and Crank et al. 2022. '
                'Includes human enteric viruses (40%), gut bacteriophages (35%), '
                'environmental viruses (15%), and emerging pathogens (10%) for '
                'epidemiological surveillance applications.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'Literature-based composition: Caliciviridae, Adenoviridae, Astroviridae, '
                'Picornaviridae (enteric); crAssphage, Microviridae, Inoviridae (phages); '
                'Plant and insect viruses (environmental); SARS-CoV-2, Mpox (emerging)'
            ),
            'curated_by': 'ViroForge Phase 7',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Crits-Christoph et al. 2021 (mBio), Crank et al. 2022 (Sci Total Environ), '
                'Symonds et al. 2019 (Curr Opin Virol)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 17")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 17 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 17")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 17")

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
                17,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present (can adjust for prevalence patterns)
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 17 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = WastewaterCurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ WASTEWATER COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 17 --dry-run")
        logger.info("  2. Create wastewater-specific contamination profile")
        logger.info("  3. Update documentation")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
