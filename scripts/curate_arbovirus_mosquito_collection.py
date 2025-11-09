#!/usr/bin/env python3
"""
Curate Collection 22: Arbovirus Environmental (Mosquito Virome)

Composition based on literature (Shi et al. 2016, Bolling et al. 2012):
- Flaviviruses (dengue, Zika, West Nile, yellow fever, Japanese encephalitis)
- Alphaviruses (chikungunya, equine encephalitis viruses)
- Bunyaviruses (La Crosse, Rift Valley fever)
- Mosquito-specific insect viruses (densovirus, iflavirus, negevirus)
- Environmental RNA viruses from vector pools

All viruses are RNA viruses (mostly ssRNA, some dsRNA).

Literature basis:
- Shi et al. 2016 (Nature): Redefining the invertebrate RNA virome
- Bolling et al. 2012 (Vector Borne Zoonotic Dis): Mosquito viromes
- Coffey et al. 2014 (Trends Microbiol): Arbovirus-vector interactions

Target size: 30-50 genomes

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


class ArbovirusCurator:
    """Curate arbovirus environmental collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_flaviviruses(self, n_target: int = 12) -> List[Dict]:
        """
        Get flaviviruses (arboviruses).

        Includes:
        - Dengue virus
        - Zika virus
        - West Nile virus
        - Yellow fever virus
        - Japanese encephalitis virus
        - Tick-borne encephalitis virus

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting flaviviruses (arboviruses)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Flaviviridae'
          AND (g.genome_name LIKE '%Dengue%'
           OR g.genome_name LIKE '%Zika%'
           OR g.genome_name LIKE '%West Nile%'
           OR g.genome_name LIKE '%Yellow fever%'
           OR g.genome_name LIKE '%Japanese encephalitis%'
           OR g.genome_name LIKE '%Tick-borne encephalitis%'
           OR g.genome_name LIKE '%encephalitis virus%'
           OR t.genus = 'Flavivirus')
        ORDER BY RANDOM()
        LIMIT ?
        """

        flaviviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Flaviviridae (arboviruses): {len(flaviviruses)}")
        return flaviviruses

    def get_alphaviruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get alphaviruses (Togaviridae).

        Includes:
        - Chikungunya virus
        - Eastern equine encephalitis virus
        - Western equine encephalitis virus
        - Venezuelan equine encephalitis virus
        - Sindbis virus
        - Ross River virus

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting alphaviruses (Togaviridae)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Togaviridae'
          AND (g.genome_name LIKE '%Chikungunya%'
           OR g.genome_name LIKE '%equine encephalitis%'
           OR g.genome_name LIKE '%Sindbis%'
           OR g.genome_name LIKE '%Ross River%'
           OR g.genome_name LIKE '%O%nyong-nyong%'
           OR g.genome_name LIKE '%Mayaro%'
           OR t.genus = 'Alphavirus')
        ORDER BY RANDOM()
        LIMIT ?
        """

        alphaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Togaviridae (Alphavirus): {len(alphaviruses)}")
        return alphaviruses

    def get_bunyaviruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get bunyaviruses (Peribunyaviridae, Phenuiviridae, Nairoviridae).

        Includes:
        - La Crosse virus (California encephalitis)
        - Rift Valley fever virus
        - Crimean-Congo hemorrhagic fever virus
        - Hantaviruses

        Genome type: ssRNA negative-sense (segmented)
        """
        logger.info("Selecting bunyaviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (t.family = 'Peribunyaviridae'
           OR t.family = 'Phenuiviridae'
           OR t.family = 'Nairoviridae')
          AND (g.genome_name LIKE '%La Crosse%'
           OR g.genome_name LIKE '%Rift Valley%'
           OR g.genome_name LIKE '%Crimean-Congo%'
           OR g.genome_name LIKE '%California encephalitis%'
           OR g.genome_name LIKE '%Hantavirus%'
           OR g.genome_name LIKE '%Bunyamwera%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        bunyaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Bunyaviruses (Peribunyaviridae, Phenuiviridae, Nairoviridae): {len(bunyaviruses)}")
        return bunyaviruses

    def get_insect_specific_viruses(self, n_target: int = 12) -> List[Dict]:
        """
        Get mosquito-specific insect viruses.

        Includes:
        - Insect-specific flaviviruses (ISF)
        - Mosquito densoviruses (Parvoviridae)
        - Mosquito iflaviruses (Iflaviridae)
        - Mosquito negeviruses (Negeviruses)

        These viruses replicate only in mosquitoes, not vertebrates.
        Genome type: mixed RNA
        """
        logger.info("Selecting mosquito-specific insect viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (
            (t.family = 'Flaviviridae' AND (
                g.genome_name LIKE '%insect%'
                OR g.genome_name LIKE '%Culex%'
                OR g.genome_name LIKE '%Aedes%'
                OR g.genome_name LIKE '%Anopheles%'
                OR g.genome_name LIKE '%cell fusing agent%'
            ))
            OR (t.family = 'Parvoviridae' AND (
                g.genome_name LIKE '%Aedes%'
                OR g.genome_name LIKE '%Culex%'
                OR g.genome_name LIKE '%mosquito%'
            ))
            OR t.family = 'Iflaviridae'
            OR t.family = 'Mesoniviridae'
        )
        ORDER BY RANDOM()
        LIMIT ?
        """

        insect_viruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Mosquito-specific viruses (ISF, densovirus, iflavirus): {len(insect_viruses)}")
        return insect_viruses

    def assign_arbovirus_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign abundances reflecting mosquito pool virome patterns.

        Mosquito viromes show:
        - Dominated by insect-specific viruses
        - Variable arbovirus abundance (depends on infection prevalence)
        - Uneven distribution (some pools have high virus, many have low)

        Uses log-normal distribution with high skew.
        """
        n = len(genomes)

        # High skew (insect-specific viruses dominant, arboviruses rarer)
        mu = -2.5
        sigma = 2.5

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete arbovirus environmental collection."""
        logger.info("=" * 80)
        logger.info("CURATING ARBOVIRUS ENVIRONMENTAL (MOSQUITO) VIROME COLLECTION")
        logger.info("=" * 80)

        # Get viruses by family
        flaviviruses = self.get_flaviviruses(n_target=12)
        alphaviruses = self.get_alphaviruses(n_target=8)
        bunyaviruses = self.get_bunyaviruses(n_target=8)
        insect_specific = self.get_insect_specific_viruses(n_target=12)

        # Combine all
        collection = flaviviruses + alphaviruses + bunyaviruses + insect_specific

        # Remove any duplicates
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign mosquito virome abundances
        logger.info("\nAssigning mosquito virome abundance pattern...")
        collection = self.assign_arbovirus_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in arbovirus mosquito collection: {len(collection)}")
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
        flavi_count = len([g for g in collection if g.get('family') == 'Flaviviridae'])
        toga_count = len([g for g in collection if g.get('family') == 'Togaviridae'])
        bunya_count = len([g for g in collection if g.get('family') in ['Peribunyaviridae', 'Phenuiviridae', 'Nairoviridae']])
        other_count = len(collection) - flavi_count - toga_count - bunya_count

        logger.info(f"Flaviviridae:             {flavi_count:3d} genomes ({flavi_count/len(collection)*100:.1f}%)")
        logger.info(f"Togaviridae (Alphavirus): {toga_count:3d} genomes ({toga_count/len(collection)*100:.1f}%)")
        logger.info(f"Bunyaviruses:             {bunya_count:3d} genomes ({bunya_count/len(collection)*100:.1f}%)")
        logger.info(f"Insect-specific/Other:    {other_count:3d} genomes ({other_count/len(collection)*100:.1f}%)")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 22,
            'collection_name': 'Arbovirus Environmental (Mosquito Virome)',
            'description': (
                'Arbovirus-enriched environmental virome from mosquito vector pools. '
                'Includes human-pathogenic arboviruses (flaviviruses, alphaviruses, bunyaviruses) '
                'and mosquito-specific insect viruses. Represents typical mosquito pool composition '
                'from arbovirus surveillance programs. Based on Shi et al. 2016, Bolling et al. 2012.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'RNA arboviruses from major families: Flaviviridae (dengue, Zika, West Nile, yellow fever), '
                'Togaviridae Alphavirus (chikungunya, equine encephalitis), '
                'Bunyaviruses (La Crosse, Rift Valley fever). '
                'Includes mosquito-specific insect viruses (insect-specific flaviviruses, densoviruses) '
                'commonly found in vector surveillance.'
            ),
            'curated_by': 'ViroForge Phase 8',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Shi et al. 2016 (Nature 540:539-543), '
                'Bolling et al. 2012 (Vector Borne Zoonotic Dis 12:657-664), '
                'Coffey et al. 2014 (Trends Microbiol 22:119-127)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 22")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 22 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 22")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 22")

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
                22,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 22 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = ArbovirusCurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ ARBOVIRUS MOSQUITO VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 22 --dry-run")
        logger.info("  2. Continue with Collection 23 (Fecal RNA Virome)")
        logger.info("  3. Implement RNA workflow components (RT, rRNA depletion)")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
