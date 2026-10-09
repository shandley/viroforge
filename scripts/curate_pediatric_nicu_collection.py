#!/usr/bin/env python3
"""
Curate Collection 22: Pediatric Gut Virome - Preterm/NICU

The NICU infant gut virome is dominated by nosocomial bacteriophages
infecting hospital-acquired bacteria (Staphylococcus, Klebsiella,
Enterococcus, E. coli). Eukaryotic viruses are less diverse than
in healthy term infants — the immature immune system and hospital
environment shape a distinct viral community.

Key differences from healthy infant gut (Collection 21):
- Hospital-associated phages dominate (60-70% of virome)
- Lower Anelloviridae prevalence (delayed environmental acquisition)
- Enteroviruses and parechoviruses clinically significant
- Lower overall diversity reflecting restricted exposure
- Phage composition tracks antibiotic-resistant bacterial colonisation

Literature basis:
- Liang et al. 2020 (Nature Medicine): Preterm vs term virome assembly
- Mukhopadhya et al. 2019 (Gut Microbes): Preterm infant virome
- Lim et al. 2015 (Nature): Early-life immune-virome dynamics
- Garmaeva et al. 2019 (Cell Host Microbe): Mother-to-infant virome

Target size: 50-70 genomes

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

COLLECTION_ID = 22


class PediatricNICUCurator:
    """Curate preterm/NICU infant gut virome collection."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.create_function("seeded_rand", 0, random.Random(42).random)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_staphylococcus_phages(self, n_target: int = 8) -> List[Dict]:
        """
        Staphylococcus phages — nosocomial. S. epidermidis and S. aureus
        are among the earliest NICU colonisers.
        """
        logger.info("Selecting Staphylococcus phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Staphylococcus phage%'
           OR g.genome_name LIKE 'Staphylococcus virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Staphylococcus phages: {len(results)}")
        return results

    def get_klebsiella_phages(self, n_target: int = 7) -> List[Dict]:
        """
        Klebsiella phages — nosocomial. K. pneumoniae is a major NICU
        pathogen associated with necrotising enterocolitis (NEC).
        """
        logger.info("Selecting Klebsiella phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Klebsiella phage%'
           OR g.genome_name LIKE 'Klebsiella virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Klebsiella phages: {len(results)}")
        return results

    def get_enterococcus_phages(self, n_target: int = 6) -> List[Dict]:
        """
        Enterococcus phages — nosocomial. VRE (vancomycin-resistant
        Enterococcus) colonisation is common in NICU infants.
        """
        logger.info("Selecting Enterococcus phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Enterococcus phage%'
           OR g.genome_name LIKE 'Enterococcus virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Enterococcus phages: {len(results)}")
        return results

    def get_ecoli_phages(self, n_target: int = 6) -> List[Dict]:
        """
        E. coli phages — both commensal and pathogenic strains colonise
        the NICU infant gut.
        """
        logger.info("Selecting Escherichia phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE 'Escherichia phage%'
               OR g.genome_name LIKE 'Escherichia virus%')
          AND g.genome_name NOT LIKE '%phiX174%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Escherichia phages: {len(results)}")
        return results

    def get_enterobacter_phages(self, n_target: int = 4) -> List[Dict]:
        """
        Enterobacter phages — Enterobacter cloacae is a common NICU
        organism associated with late-onset sepsis.
        """
        logger.info("Selecting Enterobacter phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Enterobacter%phage%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Enterobacter phages: {len(results)}")
        return results

    def get_pseudomonas_phages(self, n_target: int = 3) -> List[Dict]:
        """Pseudomonas phages — P. aeruginosa is a NICU concern."""
        logger.info("Selecting Pseudomonas phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Pseudomonas phage%'
           OR g.genome_name LIKE 'Pseudomonas virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pseudomonas phages: {len(results)}")
        return results

    def get_anelloviridae(self, n_target: int = 5) -> List[Dict]:
        """
        Anelloviridae — lower prevalence than term infants due to
        limited environmental exposure in hospital setting.
        """
        logger.info("Selecting Anelloviridae (reduced in NICU)...")
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

    def get_clinical_enteric(self, n_target: int = 8) -> List[Dict]:
        """
        Clinically significant enteric viruses: parechovirus,
        enterovirus, rotavirus. These can cause severe disease
        in preterm infants.
        """
        logger.info("Selecting clinically significant enteric viruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human parechovirus%'
           OR g.genome_name LIKE 'Enterovirus%'
           OR g.genome_name LIKE 'Human enterovirus%'
           OR (t.family = 'Sedoreoviridae' AND g.genome_name LIKE 'Rotavirus%')
           OR (t.family = 'Caliciviridae' AND g.genome_name LIKE 'Norovirus%')
           OR (t.family = 'Astroviridae' AND g.genome_name LIKE 'Human astrovirus%')
           OR (t.family = 'Polyomaviridae' AND g.genome_name LIKE 'Human polyomavirus%')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Clinical enteric viruses: {len(results)}")
        return results

    def get_microviridae(self, n_target: int = 3) -> List[Dict]:
        """Microviridae (ssDNA coliphages) — present in NICU infants."""
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

    def assign_nicu_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        NICU abundance: hospital phage-dominated with high dominance
        by a few taxa (low evenness). Reflects limited microbial
        diversity in the hospital environment.
        """
        n = len(genomes)
        mu = -2.5
        sigma = 2.5
        raw = np.random.lognormal(mu, sigma, n)

        # Boost hospital-associated phages
        hospital_hosts = ('Staphylococcus', 'Klebsiella', 'Enterococcus',
                          'Escherichia', 'Enterobacter', 'Pseudomonas')
        for i, g in enumerate(genomes):
            name = g.get('genome_name', '')
            if any(name.startswith(h) for h in hospital_hosts):
                raw[i] *= 2.5

        normalized = raw / raw.sum()
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)
        return genomes

    def create_collection(self) -> List[Dict]:
        logger.info("=" * 80)
        logger.info("CURATING PEDIATRIC GUT VIROME - PRETERM/NICU")
        logger.info("=" * 80)

        staph = self.get_staphylococcus_phages(n_target=8)
        kleb = self.get_klebsiella_phages(n_target=7)
        entero = self.get_enterococcus_phages(n_target=6)
        ecoli = self.get_ecoli_phages(n_target=6)
        enterob = self.get_enterobacter_phages(n_target=4)
        pseudo = self.get_pseudomonas_phages(n_target=3)
        anello = self.get_anelloviridae(n_target=5)
        clinical = self.get_clinical_enteric(n_target=8)
        micro = self.get_microviridae(n_target=3)

        collection = (staph + kleb + entero + ecoli + enterob +
                      pseudo + anello + clinical + micro)

        seen_ids = set()
        unique = []
        for g in collection:
            if g['genome_id'] not in seen_ids:
                unique.append(g)
                seen_ids.add(g['genome_id'])
        collection = unique

        logger.info(f"\nPre-abundance total: {len(collection)} genomes")
        collection = self.assign_nicu_abundances(collection)
        collection.sort(key=lambda x: x['relative_abundance'], reverse=True)
        for i, g in enumerate(collection, 1):
            g['abundance_rank'] = i

        logger.info(f"Total genomes in NICU collection: {len(collection)}")
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
            'Pediatric Gut Virome - Preterm/NICU',
            'Gut virome of preterm infants in NICU settings. Dominated by '
            'nosocomial bacteriophages (Staphylococcus, Klebsiella, '
            'Enterococcus, E. coli phages) reflecting hospital-acquired '
            'bacterial colonisation. Lower eukaryotic virus diversity than '
            'term infants. Clinically relevant enteroviruses and '
            'parechoviruses included. '
            'Based on Liang et al. 2020, Mukhopadhya et al. 2019.',
            len(collection),
            'Hospital phage-dominated (60-70%). Staphylococcus + Klebsiella '
            'phages most abundant. Reduced Anelloviridae. Parechovirus and '
            'enterovirus as clinical markers. Low crAssphage (immature microbiome).',
            'ViroForge Pediatric',
            '2026-09-01',
            'Liang et al. 2020 (Nature Medicine), '
            'Mukhopadhya et al. 2019 (Gut Microbes), '
            'Lim et al. 2015 (Nature), '
            'Garmaeva et al. 2019 (Cell Host Microbe)',
            'Homo sapiens (preterm infant)',
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
    curator = PediatricNICUCurator()
    try:
        collection = curator.create_collection()
        curator.insert_collection(collection)
        logger.info("\n" + "=" * 80)
        logger.info("PEDIATRIC NICU VIROME CURATION COMPLETE!")
        logger.info("=" * 80)
    finally:
        curator.close()


if __name__ == '__main__':
    main()
