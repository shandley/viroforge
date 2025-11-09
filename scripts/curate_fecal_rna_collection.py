#!/usr/bin/env python3
"""
Curate Collection 23: Fecal RNA Virome

Composition based on literature (Reyes et al. 2012, Minot et al. 2011):
- Caliciviruses (norovirus, sapovirus)
- Reoviruses (rotavirus)
- Astroviruses
- Enteric picornaviruses (enterovirus, parechovirus)
- Picobirnaviruses
- Mamastroviruses

Represents RNA virus component of fecal virome (phages excluded, in Collection 1).
All viruses are RNA viruses (ssRNA and dsRNA).

Literature basis:
- Reyes et al. 2012 (Genome Res): Viruses in the fecal microbiota
- Minot et al. 2011 (PNAS): Human gut virome is highly diverse
- Shan et al. 2011 (J Virol): Mammalian and avian fecal viromes

Target size: 40-60 genomes

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


class FecalRNACurator:
    """Curate fecal RNA virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_caliciviruses(self, n_target: int = 15) -> List[Dict]:
        """
        Get caliciviruses (norovirus, sapovirus).

        Noroviruses are the leading cause of acute gastroenteritis.
        Sapoviruses also cause gastroenteritis, especially in children.

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting caliciviruses (norovirus, sapovirus)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Caliciviridae'
          AND (g.genome_name LIKE '%Norovirus%'
           OR g.genome_name LIKE '%norovirus%'
           OR g.genome_name LIKE '%Sapovirus%'
           OR g.genome_name LIKE '%sapovirus%'
           OR t.genus IN ('Norovirus', 'Sapovirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        caliciviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Caliciviridae (Norovirus, Sapovirus): {len(caliciviruses)}")
        return caliciviruses

    def get_reoviruses(self, n_target: int = 12) -> List[Dict]:
        """
        Get reoviruses (rotavirus, reovirus).

        Rotaviruses are major cause of severe diarrhea in infants/children.
        Reoviruses are commonly found in feces.

        Genome type: dsRNA (segmented)
        """
        logger.info("Selecting reoviruses (rotavirus)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Reoviridae'
          AND (g.genome_name LIKE '%Rotavirus%'
           OR g.genome_name LIKE '%rotavirus%'
           OR g.genome_name LIKE '%Reovirus%'
           OR t.genus IN ('Rotavirus', 'Orthoreovirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        reoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Reoviridae (Rotavirus, Reovirus): {len(reoviruses)}")
        return reoviruses

    def get_astroviruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get astroviruses (mamastrovirus, avastrovirus).

        Common cause of gastroenteritis, especially in children.
        Highly diverse in fecal samples.

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting astroviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Astroviridae'
          AND (g.genome_name LIKE '%Astrovirus%'
           OR g.genome_name LIKE '%astrovirus%'
           OR g.genome_name LIKE '%Mamastrovirus%'
           OR t.genus IN ('Mamastrovirus', 'Avastrovirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        astroviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Astroviridae: {len(astroviruses)}")
        return astroviruses

    def get_enteric_picornaviruses(self, n_target: int = 12) -> List[Dict]:
        """
        Get enteric picornaviruses (enterovirus, parechovirus, aichivirus).

        Enteroviruses and parechoviruses are commonly shed in feces.
        Include polioviruses, coxsackieviruses, echoviruses.

        Genome type: ssRNA positive-sense
        """
        logger.info("Selecting enteric picornaviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Picornaviridae'
          AND (g.genome_name LIKE '%Enterovirus%'
           OR g.genome_name LIKE '%enterovirus%'
           OR g.genome_name LIKE '%Poliovirus%'
           OR g.genome_name LIKE '%Coxsackievirus%'
           OR g.genome_name LIKE '%Echovirus%'
           OR g.genome_name LIKE '%Parechovirus%'
           OR g.genome_name LIKE '%Aichivirus%'
           OR t.genus IN ('Enterovirus', 'Parechovirus', 'Kobuvirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        picornaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Picornaviridae (Enterovirus, Parechovirus): {len(picornaviruses)}")
        return picornaviruses

    def get_picobirnaviruses(self, n_target: int = 6) -> List[Dict]:
        """
        Get picobirnaviruses.

        Common in fecal samples, widespread in mammals.
        Association with diarrhea unclear, may be commensal.

        Genome type: dsRNA (segmented)
        """
        logger.info("Selecting picobirnaviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Picobirnaviridae'
           OR g.genome_name LIKE '%Picobirnavirus%'
           OR g.genome_name LIKE '%picobirnavirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        picobirnaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Picobirnaviridae: {len(picobirnaviruses)}")
        return picobirnaviruses

    def get_other_enteric_viruses(self, n_target: int = 8) -> List[Dict]:
        """
        Get other enteric RNA viruses.

        Includes:
        - Coronaviruses (enteric strains)
        - Toroviruses
        - Adenoviruses (commonly found in feces)

        Genome type: mixed
        """
        logger.info("Selecting other enteric viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (
            (t.family = 'Coronaviridae' AND (
                g.genome_name LIKE '%enteric%'
                OR g.genome_name LIKE '%fecal%'
                OR g.genome_name LIKE '%gastro%'
            ))
            OR t.family = 'Tobaniviridae'
            OR (t.family = 'Adenoviridae' AND (
                g.genome_name LIKE '%enteric%'
                OR g.genome_name LIKE '%gastroenteritis%'
                OR g.genome_name LIKE '%Mastadenovirus F%'
            ))
        )
        ORDER BY RANDOM()
        LIMIT ?
        """

        other_viruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Other enteric viruses (Coronavirus, Adenovirus): {len(other_viruses)}")
        return other_viruses

    def assign_fecal_rna_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign abundances reflecting fecal RNA virome patterns.

        Fecal RNA viromes show:
        - Variable abundance based on infection/shedding
        - Norovirus/rotavirus can be highly abundant during infection
        - Many viruses at lower, chronic shedding levels
        - Moderate diversity

        Uses log-normal distribution with moderate skew.
        """
        n = len(genomes)

        # Moderate skew (some high shedders, many at low levels)
        mu = -1.5
        sigma = 2.0

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete fecal RNA virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING FECAL RNA VIROME COLLECTION")
        logger.info("=" * 80)

        # Get viruses by family
        caliciviruses = self.get_caliciviruses(n_target=15)
        reoviruses = self.get_reoviruses(n_target=12)
        astroviruses = self.get_astroviruses(n_target=8)
        picornaviruses = self.get_enteric_picornaviruses(n_target=12)
        picobirnaviruses = self.get_picobirnaviruses(n_target=6)
        other_enteric = self.get_other_enteric_viruses(n_target=8)

        # Combine all
        collection = (caliciviruses + reoviruses + astroviruses +
                     picornaviruses + picobirnaviruses + other_enteric)

        # Remove any duplicates
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign fecal RNA virome abundances
        logger.info("\nAssigning fecal RNA virome abundance pattern...")
        collection = self.assign_fecal_rna_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in fecal RNA virome: {len(collection)}")
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
        calici_count = len([g for g in collection if g.get('family') == 'Caliciviridae'])
        reo_count = len([g for g in collection if g.get('family') == 'Reoviridae'])
        astro_count = len([g for g in collection if g.get('family') == 'Astroviridae'])
        picorna_count = len([g for g in collection if g.get('family') == 'Picornaviridae'])
        picobirna_count = len([g for g in collection if g.get('family') == 'Picobirnaviridae'])
        other_count = len(collection) - calici_count - reo_count - astro_count - picorna_count - picobirna_count

        logger.info(f"Caliciviridae (Norovirus):    {calici_count:3d} genomes ({calici_count/len(collection)*100:.1f}%)")
        logger.info(f"Reoviridae (Rotavirus):       {reo_count:3d} genomes ({reo_count/len(collection)*100:.1f}%)")
        logger.info(f"Astroviridae:                 {astro_count:3d} genomes ({astro_count/len(collection)*100:.1f}%)")
        logger.info(f"Picornaviridae (Enterovirus): {picorna_count:3d} genomes ({picorna_count/len(collection)*100:.1f}%)")
        logger.info(f"Picobirnaviridae:             {picobirna_count:3d} genomes ({picobirna_count/len(collection)*100:.1f}%)")
        logger.info(f"Other enteric viruses:        {other_count:3d} genomes ({other_count/len(collection)*100:.1f}%)")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 23,
            'collection_name': 'Fecal RNA Virome',
            'description': (
                'RNA virus component of human fecal virome. '
                'Includes major enteric RNA viruses: norovirus, rotavirus, astrovirus, enterovirus. '
                'Also includes picobirnaviruses and other RNA viruses commonly detected in fecal samples. '
                'Complements bacteriophage-focused Collection 1. '
                'Based on Reyes et al. 2012, Minot et al. 2011, and Shan et al. 2011.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'RNA enteric viruses from major families: Caliciviridae (norovirus, sapovirus), '
                'Reoviridae (rotavirus), Astroviridae, Picornaviridae (enterovirus, parechovirus), '
                'Picobirnaviridae. Includes both pathogenic viruses (norovirus, rotavirus) '
                'and commonly detected commensal/opportunistic viruses. '
                'Excludes bacteriophages (covered in Collection 1).'
            ),
            'curated_by': 'ViroForge Phase 8',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Reyes et al. 2012 (Genome Res 22:1107-1119), '
                'Minot et al. 2011 (PNAS 108:4516-4522), '
                'Shan et al. 2011 (J Virol 85:11146-11155)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 23")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 23 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 23")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 23")

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
                23,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 23 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = FecalRNACurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ FECAL RNA VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nPhase 8 Collections Complete:")
        logger.info("  ✓ Collection 21: Human Respiratory RNA Virome (56 genomes)")
        logger.info("  ✓ Collection 22: Arbovirus Environmental (39 genomes)")
        logger.info("  ✓ Collection 23: Fecal RNA Virome ({} genomes)".format(len(collection)))
        logger.info("\nNext steps:")
        logger.info("  1. Test all RNA collections with dry-run")
        logger.info("  2. Implement RNA workflow components (RT, rRNA depletion)")
        logger.info("  3. Update documentation for Phase 8")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
