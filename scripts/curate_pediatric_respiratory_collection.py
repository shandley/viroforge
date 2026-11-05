#!/usr/bin/env python3
"""
Curate Collection 23: Pediatric Respiratory Virome

The pediatric respiratory virome is richer and more dynamic than the
adult respiratory virome (Collection 4/13). Children under 5 experience
6-8 respiratory infections per year, driving higher prevalence and
diversity of respiratory viruses.

Key differences from adult respiratory (Collections 4, 13):
- Higher rhinovirus diversity (>100 serotypes circulate)
- RSV and HMPV more prevalent and pathogenic in <2 year olds
- Parainfluenza viruses more common than in adults
- Bocavirus (HBoV1) nearly universal in children <5
- Adenoviruses more diverse (species A, B, C, E)
- Anelloviridae present in upper airways
- Higher co-detection rates (2+ viruses simultaneously common)

Literature basis:
- Jartti et al. 2012 (Journal of Allergy and Clinical Immunology):
  Respiratory viral epidemiology in children
- Wylie et al. 2018 (mBio): Human airway virome across development
- Arden et al. 2010 (PLoS One): Bocavirus prevalence in children
- Van der Hoek et al. 2005 (Nature Medicine): Novel coronaviruses in children

Target size: 40-55 genomes

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

COLLECTION_ID = 23


class PediatricRespiratoryCurator:
    """Curate pediatric respiratory virome collection."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.create_function("seeded_rand", 0, random.Random(42).random)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_rhinoviruses(self, n_target: int = 6) -> List[Dict]:
        """
        Rhinoviruses — most frequent respiratory virus in children.
        >100 serotypes; species A and C drive most paediatric illness.
        """
        logger.info("Selecting rhinoviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Rhinovirus%'
           OR g.genome_name LIKE 'Human rhinovirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Rhinoviruses: {len(results)}")
        return results

    def get_rsv(self, n_target: int = 2) -> List[Dict]:
        """
        RSV — leading cause of bronchiolitis and pneumonia in infants.
        Limited RefSeq diversity but critical inclusion.
        """
        logger.info("Selecting RSV...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human respiratory syncytial virus%'
           OR g.genome_name LIKE 'Human orthopneumovirus%'
           OR g.genome_name LIKE 'Respiratory syncytial virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  RSV: {len(results)}")
        return results

    def get_parainfluenza(self, n_target: int = 4) -> List[Dict]:
        """
        Parainfluenza viruses (types 1-4) — cause croup and
        bronchiolitis; more prevalent in children than adults.
        """
        logger.info("Selecting parainfluenza viruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human parainfluenza%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Parainfluenza: {len(results)}")
        return results

    def get_hmpv(self, n_target: int = 1) -> List[Dict]:
        """
        Human metapneumovirus — second leading cause of bronchiolitis
        after RSV. Only 1 genome in RefSeq.
        """
        logger.info("Selecting HMPV...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human metapneumovirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  HMPV: {len(results)}")
        return results

    def get_adenoviruses(self, n_target: int = 6) -> List[Dict]:
        """
        Human adenoviruses — species A-E cause respiratory, GI, and
        ocular infections in children. More diverse than in adults.
        """
        logger.info("Selecting adenoviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
          AND g.genome_name LIKE 'Human adenovirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Adenoviruses: {len(results)}")
        return results

    def get_bocavirus(self, n_target: int = 3) -> List[Dict]:
        """
        Human bocavirus (HBoV1-4) — nearly universal seroconversion by
        age 5. HBoV1 is respiratory; HBoV2-4 are enteric.
        """
        logger.info("Selecting bocaviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human bocavirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Bocaviruses: {len(results)}")
        return results

    def get_coronaviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Coronaviruses: seasonal HCoV (229E, OC43, NL63, HKU1) plus
        SARS-CoV-2. NL63 was discovered through paediatric surveillance.
        """
        logger.info("Selecting coronaviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Coronaviridae'
          AND (g.genome_name LIKE 'Human coronavirus%'
               OR g.genome_name LIKE 'Severe acute respiratory%')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Coronaviruses: {len(results)}")
        return results

    def get_influenza(self, n_target: int = 4) -> List[Dict]:
        """
        Influenza A and B — children are major drivers of
        community transmission; broader serotype exposure.
        """
        logger.info("Selecting influenza viruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Influenza A virus%'
           OR g.genome_name LIKE 'Influenza B virus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Influenza: {len(results)}")
        return results

    def get_enteroviruses(self, n_target: int = 4) -> List[Dict]:
        """
        Enteroviruses (EV-D68, coxsackieviruses) — respiratory
        tropism in some species; paediatric outbreaks documented.
        """
        logger.info("Selecting enteroviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Enterovirus%'
           OR g.genome_name LIKE 'Coxsackievirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Enteroviruses: {len(results)}")
        return results

    def get_anelloviridae(self, n_target: int = 5) -> List[Dict]:
        """
        Anelloviridae — persistent in upper airways; prevalence
        correlates inversely with immune maturation (Wylie 2018).
        """
        logger.info("Selecting Anelloviridae...")
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

    def get_polyomaviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Respiratory polyomaviruses (KIPyV, WUPyV) — discovered in
        paediatric respiratory samples. Often co-detected with other
        respiratory viruses.
        """
        logger.info("Selecting polyomaviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Polyomaviridae'
          AND (g.genome_name LIKE '%KI polyomavirus%'
               OR g.genome_name LIKE '%WU polyomavirus%'
               OR g.genome_name LIKE 'Human polyomavirus%')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Polyomaviruses: {len(results)}")
        return results

    def assign_respiratory_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Pediatric respiratory abundance: rhinovirus-dominated with
        moderate diversity. RSV and HMPV get seasonal boosts when
        present. More even than adult respiratory due to higher
        co-detection rates.
        """
        n = len(genomes)
        mu = -1.5
        sigma = 1.8
        raw = np.random.lognormal(mu, sigma, n)

        # Boost the most clinically important respiratory viruses
        for i, g in enumerate(genomes):
            name = g.get('genome_name', '')
            if 'Rhinovirus' in name or 'rhinovirus' in name:
                raw[i] *= 3.0
            elif 'respiratory syncytial' in name.lower() or 'orthopneumovirus' in name.lower():
                raw[i] *= 2.5
            elif 'metapneumovirus' in name.lower():
                raw[i] *= 2.0
            elif 'parainfluenza' in name.lower():
                raw[i] *= 1.5

        normalized = raw / raw.sum()
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)
        return genomes

    def create_collection(self) -> List[Dict]:
        logger.info("=" * 80)
        logger.info("CURATING PEDIATRIC RESPIRATORY VIROME")
        logger.info("=" * 80)

        rhino = self.get_rhinoviruses(n_target=6)
        rsv = self.get_rsv(n_target=2)
        piv = self.get_parainfluenza(n_target=4)
        hmpv = self.get_hmpv(n_target=1)
        adeno = self.get_adenoviruses(n_target=6)
        boca = self.get_bocavirus(n_target=3)
        corona = self.get_coronaviruses(n_target=5)
        flu = self.get_influenza(n_target=4)
        entero = self.get_enteroviruses(n_target=4)
        anello = self.get_anelloviridae(n_target=5)
        polyoma = self.get_polyomaviruses(n_target=3)

        collection = (rhino + rsv + piv + hmpv + adeno + boca +
                      corona + flu + entero + anello + polyoma)

        seen_ids = set()
        unique = []
        for g in collection:
            if g['genome_id'] not in seen_ids:
                unique.append(g)
                seen_ids.add(g['genome_id'])
        collection = unique

        logger.info(f"\nPre-abundance total: {len(collection)} genomes")
        collection = self.assign_respiratory_abundances(collection)
        collection.sort(key=lambda x: x['relative_abundance'], reverse=True)
        for i, g in enumerate(collection, 1):
            g['abundance_rank'] = i

        logger.info(f"Total genomes in pediatric respiratory collection: {len(collection)}")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        logger.info("\nTop 10 most abundant genomes:")
        for i, g in enumerate(collection[:10], 1):
            fam = g.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {g['genome_name'][:50]:50s} {g['relative_abundance']:.6f} ({fam})")

        # Show family breakdown
        families = {}
        for g in collection:
            fam = g.get('family', 'Unknown')
            families[fam] = families.get(fam, 0) + 1
        logger.info("\nFamily breakdown:")
        for fam, count in sorted(families.items(), key=lambda x: -x[1]):
            logger.info(f"  {fam}: {count}")

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
            'Pediatric Respiratory Virome',
            'Respiratory virome of children under 5 years. Higher diversity '
            'than adult respiratory virome with rhinovirus, RSV, parainfluenza, '
            'HMPV, adenovirus, bocavirus, and coronaviruses. Anelloviridae '
            'persistent in upper airways. Higher co-detection rates than adults. '
            'Based on Jartti et al. 2012, Wylie et al. 2018.',
            len(collection),
            'Rhinovirus-dominated. RSV and HMPV critical for <2 years. '
            'Parainfluenza types 1-4. Bocavirus nearly universal <5 years. '
            'Adenovirus species A-E. Seasonal coronaviruses. '
            'Anelloviridae and respiratory polyomaviruses.',
            'ViroForge Pediatric',
            '2026-09-01',
            'Jartti et al. 2012 (Journal of Allergy and Clinical Immunology), '
            'Wylie et al. 2018 (mBio), '
            'Arden et al. 2010 (PLoS One), '
            'Van der Hoek et al. 2005 (Nature Medicine)',
            'Homo sapiens (child, <5 years)',
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
    curator = PediatricRespiratoryCurator()
    try:
        collection = curator.create_collection()
        curator.insert_collection(collection)
        logger.info("\n" + "=" * 80)
        logger.info("PEDIATRIC RESPIRATORY VIROME CURATION COMPLETE!")
        logger.info("=" * 80)
    finally:
        curator.close()


if __name__ == '__main__':
    main()
