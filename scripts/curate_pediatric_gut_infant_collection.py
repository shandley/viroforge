#!/usr/bin/env python3
"""
Curate Collection 21: Pediatric Gut Virome - Healthy Infant (0-3 years)

The infant gut virome is distinct from adults: dominated by eukaryotic
viruses (Anelloviridae, Genomoviridae, CRESS-DNA) in early life, with
bacteriophages colonising alongside the developing microbiome. crAssphage
emerges around 6-12 months and rises to dominance by year 2-3.

Key differences from adult gut (Collection 1):
- Anelloviridae dominant (30-60% of virome), rare in adults
- CRESS-DNA viruses (Genomoviridae, Circoviridae) highly prevalent
- crAssphage-like phages present but lower diversity than adults
- Pioneer bacteriophages: Bifidobacterium, Lactobacillus, E. coli phages
- Lower overall phage diversity reflecting simpler microbiome
- Adenovirus and astrovirus common from environmental/dietary exposure

Literature basis:
- Liang et al. 2020 (Nature Medicine): Stepwise assembly of the neonatal virome
- Shkoporov et al. 2019 (Cell Host Microbe): crAssphage temporal dynamics
- Lim et al. 2015 (Nature): Early-life gut virome and immunity
- Reyes et al. 2015 (Nature Medicine): Viral diversity in malnourished infants

Target size: 60-80 genomes

Author: ViroForge Development Team
Date: 2026-09-01
"""

import sqlite3
import random
import numpy as np
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

COLLECTION_ID = 21


class PediatricGutInfantCurator:
    """Curate healthy infant gut virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.create_function("seeded_rand", 0, random.Random(42).random)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_anelloviridae(self, n_target: int = 15) -> List[Dict]:
        """
        Anelloviridae (TTV, TTMV, TTMDV) - dominant in infant gut.

        Vertical transmission and environmental acquisition; prevalence
        reaches 90%+ in infants by 12 months (Lim et al. 2015).
        """
        logger.info("Selecting Anelloviridae (dominant in infants)...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Anelloviridae'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Anelloviridae: {len(results)}")
        return results

    def get_cress_dna_viruses(self, n_target: int = 10) -> List[Dict]:
        """
        CRESS-DNA viruses (Genomoviridae, Circoviridae) - highly prevalent
        in infant gut, likely from dietary and environmental exposure.
        """
        logger.info("Selecting CRESS-DNA viruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Genomoviridae', 'Circoviridae')
          AND g.genome_name NOT LIKE '%Porcine%'
          AND g.genome_name NOT LIKE '%Beak and feather%'
          AND g.genome_name NOT LIKE '%Pigeon%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  CRESS-DNA (Genomoviridae + Circoviridae): {len(results)}")
        return results

    def get_crassphage(self, n_target: int = 8) -> List[Dict]:
        """
        crAssphage-like phages: present from ~6 months, lower diversity
        than adults. Emerges with Bacteroides colonisation.
        """
        logger.info("Selecting crAssphage-like phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Intestiviridae', 'Suoliviridae',
                           'Steigviridae', 'Crevaviridae')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  crAssphage-like: {len(results)}")
        return results

    def get_pioneer_phages(self, n_target: int = 12) -> List[Dict]:
        """
        Pioneer bacteriophages infecting early gut colonisers:
        Bifidobacterium, Lactobacillus, Escherichia coli.

        These track the developing microbiome (Liang et al. 2020).
        """
        logger.info("Selecting pioneer bacteriophages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Bifidobacterium phage%'
           OR g.genome_name LIKE 'Lactobacillus phage%'
           OR g.genome_name LIKE 'Lactococcus phage%'
           OR (g.genome_name LIKE 'Escherichia phage%'
               AND g.genome_name NOT LIKE '%phiX174%')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pioneer phages: {len(results)}")
        return results

    def get_microviridae(self, n_target: int = 5) -> List[Dict]:
        """Small ssDNA coliphages - present early in infant gut."""
        logger.info("Selecting Microviridae...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Microviridae'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Microviridae: {len(results)}")
        return results

    def get_enteric_eukaryotic(self, n_target: int = 10) -> List[Dict]:
        """
        Enteric eukaryotic viruses commonly detected in infants:
        adenovirus, astrovirus, calicivirus (norovirus/sapovirus),
        parechovirus. Often asymptomatic in healthy infants.
        """
        logger.info("Selecting enteric eukaryotic viruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (t.family = 'Adenoviridae' AND g.genome_name LIKE 'Human adenovirus%')
           OR (t.family = 'Astroviridae' AND g.genome_name LIKE 'Human astrovirus%')
           OR (t.family = 'Caliciviridae'
               AND (g.genome_name LIKE 'Norovirus%' OR g.genome_name LIKE 'Sapovirus%'))
           OR (t.family = 'Parvoviridae' AND g.genome_name LIKE 'Human bocavirus%')
           OR g.genome_name LIKE 'Human parechovirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Enteric eukaryotic: {len(results)}")
        return results

    def get_parvoviridae(self, n_target: int = 5) -> List[Dict]:
        """
        Parvoviridae (non-bocavirus): dependoparvoviruses (AAV) and
        erythroparvoviruses detected in paediatric stool.
        """
        logger.info("Selecting Parvoviridae...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Parvoviridae'
          AND g.genome_name NOT LIKE '%bocavirus%'
          AND g.genome_name NOT LIKE '%Porcine%'
          AND g.genome_name NOT LIKE '%Canine%'
          AND g.genome_name NOT LIKE '%Feline%'
          AND g.genome_name NOT LIKE '%Bovine%'
          AND g.genome_name NOT LIKE '%Rat%'
          AND g.genome_name NOT LIKE '%Mouse%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Parvoviridae: {len(results)}")
        return results

    def assign_infant_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Infant gut abundance: Anelloviridae-dominated with moderate
        evenness among phages. Less skewed than adult gut but with a
        strong eukaryotic virus component.
        """
        n = len(genomes)
        mu = -2.0
        sigma = 2.0
        raw = np.random.lognormal(mu, sigma, n)

        # Boost Anelloviridae abundances (dominant in infants)
        for i, g in enumerate(genomes):
            if g.get('family') == 'Anelloviridae':
                raw[i] *= 3.0

        normalized = raw / raw.sum()
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)
        return genomes

    def create_collection(self) -> List[Dict]:
        logger.info("=" * 80)
        logger.info("CURATING PEDIATRIC GUT VIROME - HEALTHY INFANT (0-3 years)")
        logger.info("=" * 80)

        anello = self.get_anelloviridae(n_target=15)
        cress = self.get_cress_dna_viruses(n_target=10)
        crass = self.get_crassphage(n_target=8)
        pioneer = self.get_pioneer_phages(n_target=12)
        micro = self.get_microviridae(n_target=5)
        enteric = self.get_enteric_eukaryotic(n_target=10)
        parvo = self.get_parvoviridae(n_target=5)

        collection = anello + cress + crass + pioneer + micro + enteric + parvo

        seen_ids = set()
        unique = []
        for g in collection:
            if g['genome_id'] not in seen_ids:
                unique.append(g)
                seen_ids.add(g['genome_id'])
        collection = unique

        logger.info(f"\nPre-abundance total: {len(collection)} genomes")
        collection = self.assign_infant_abundances(collection)
        collection.sort(key=lambda x: x['relative_abundance'], reverse=True)
        for i, g in enumerate(collection, 1):
            g['abundance_rank'] = i

        logger.info(f"Total genomes in infant gut collection: {len(collection)}")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        logger.info("\nTop 10 most abundant genomes:")
        for i, g in enumerate(collection[:10], 1):
            fam = g.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {g['genome_name'][:50]:50s} {g['relative_abundance']:.6f} ({fam})")

        return collection

    def insert_collection(self, collection: List[Dict]):
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        cursor.execute(
            "SELECT collection_id FROM body_site_collections WHERE collection_id = ?",
            (COLLECTION_ID,))
        if cursor.fetchone():
            logger.info(f"Collection {COLLECTION_ID} already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = ?",
                           (COLLECTION_ID,))
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = ?",
                           (COLLECTION_ID,))

        cursor.execute("""
            INSERT INTO body_site_collections
            (collection_id, collection_name, description, n_genomes, selection_criteria,
             curated_by, curation_date, literature_references, host_organism, version)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            COLLECTION_ID,
            'Pediatric Gut Virome - Healthy Infant (0-3 years)',
            'Gut virome of healthy infants aged 0-3 years. Dominated by '
            'eukaryotic viruses (Anelloviridae, Genomoviridae) and pioneer '
            'bacteriophages (Bifidobacterium, Lactobacillus phages). crAssphage '
            'diversity lower than adults, reflecting the developing microbiome. '
            'Based on Liang et al. 2020, Lim et al. 2015, Shkoporov et al. 2019.',
            len(collection),
            'Anelloviridae-dominated (30-60%). CRESS-DNA viruses prevalent. '
            'Pioneer phages from early colonisers. crAssphage present but reduced. '
            'Enteric eukaryotic viruses (adeno, astro, calici) common.',
            'ViroForge Pediatric',
            '2026-09-01',
            'Liang et al. 2020 (Nature Medicine), Shkoporov et al. 2019 '
            '(Cell Host Microbe), Lim et al. 2015 (Nature), '
            'Reyes et al. 2015 (Nature Medicine)',
            'Homo sapiens (infant)',
            1,
        ))
        logger.info("Inserted collection metadata")

        for g in collection:
            cursor.execute("""
                INSERT INTO collection_genomes
                (collection_id, genome_id, relative_abundance, prevalence, abundance_rank)
                VALUES (?, ?, ?, ?, ?)
            """, (COLLECTION_ID, g['genome_id'], g['relative_abundance'],
                  1.0, g['abundance_rank']))

        self.conn.commit()
        logger.info(f"Inserted {len(collection)} genome associations")
        logger.info(f"\nCollection {COLLECTION_ID} successfully created!")

    def close(self):
        self.conn.close()


def main():
    curator = PediatricGutInfantCurator()
    try:
        collection = curator.create_collection()
        curator.insert_collection(collection)
        logger.info("\n" + "=" * 80)
        logger.info("PEDIATRIC GUT INFANT VIROME CURATION COMPLETE!")
        logger.info("=" * 80)
    finally:
        curator.close()


if __name__ == '__main__':
    main()
