#!/usr/bin/env python3
"""
Curate Collection 24: HPV-Infected Vaginal Virome

The HPV-positive vaginal virome differs from the healthy vaginal virome
(Collection 16) in several key ways:
- HPV dominates the eukaryotic viral fraction (high-risk types 16, 18, etc.)
- Lactobacillus-dominated CST shifts toward CST IV (BV-like) in persistent
  HPV infection, altering the phageome
- Co-infections with multiple HPV types are common (20-40% of HPV+ women)
- HSV reactivation and other sexually transmitted viruses co-occur
- Altered Anelloviridae prevalence
- Increased overall viral diversity compared to healthy

Key differences from healthy vaginal (Collection 16):
- HPV genomes present at high abundance (absent/rare in healthy)
- Higher phage diversity reflecting dysbiotic bacterial community
- HSV-1 as STI co-infection marker
- Broader enteric virus shedding
- More diverse Anelloviridae (immune perturbation)

Literature basis:
- Wylie et al. 2014 (PLOS Pathogens): Cervicovaginal virome and HPV
- Virtanen et al. 2016 (International Journal of Cancer): HPV type
  distribution in cervical infections
- Brusselaers et al. 2019 (BMC Medicine): Vaginal dysbiosis and HPV
- Laniewski et al. 2020 (mBio): Vaginal microbiome and HPV persistence

Target size: 40-55 genomes

Author: ViroForge Development Team
Date: 2026-09-01
"""

import re
import sqlite3
import random
import numpy as np
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

COLLECTION_ID = 24


class HPVVaginalCurator:
    """Curate HPV-infected vaginal virome collection."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.create_function("seeded_rand", 0, random.Random(42).random)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def _fetch_all_hpv(self) -> List[Dict]:
        """Fetch all HPV genomes once, extract type number in Python."""
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human papillomavirus type%'
        ORDER BY seeded_rand()
        """
        results = []
        for row in self.conn.execute(query):
            g = dict(row)
            m = re.search(r'type (\d+)', g['genome_name'])
            if m:
                g['_hpv_type'] = int(m.group(1))
                results.append(g)
        return results

    def get_hpv_high_risk(self, all_hpv: List[Dict], n_target: int = 10) -> List[Dict]:
        """
        High-risk HPV types (IARC Group 1/2A carcinogens).

        HPV 16 and 18 cause ~70% of cervical cancers. Types 31, 33,
        35, 39, 45, 51, 52, 56, 58, 59, 68 are also high-risk.
        Multiple type co-infection is common.
        """
        logger.info("Selecting high-risk HPV types...")
        hr_types = {16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 68}
        results = [g for g in all_hpv if g['_hpv_type'] in hr_types][:n_target]
        logger.info(f"  High-risk HPV: {len(results)} (available types: "
                    f"{sorted(set(g['_hpv_type'] for g in results))})")
        return results

    def get_hpv_low_risk(self, all_hpv: List[Dict], n_target: int = 4) -> List[Dict]:
        """
        Low-risk HPV types (condylomata, low-grade lesions).
        Types 6 and 11 cause >90% of genital warts.
        """
        logger.info("Selecting low-risk HPV types...")
        lr_types = {6, 11, 42, 43, 44}
        results = [g for g in all_hpv if g['_hpv_type'] in lr_types][:n_target]
        logger.info(f"  Low-risk HPV: {len(results)} (available types: "
                    f"{sorted(set(g['_hpv_type'] for g in results))})")
        return results

    def get_hpv_other(self, all_hpv: List[Dict], selected_ids: set,
                      n_target: int = 4) -> List[Dict]:
        """
        Other HPV types — cutaneous and mucosal types found in the
        vaginal tract at lower prevalence.
        """
        logger.info("Selecting other HPV types...")
        results = [g for g in all_hpv
                   if g['genome_id'] not in selected_ids][:n_target]
        logger.info(f"  Other HPV: {len(results)}")
        return results

    def get_lactobacillus_phages(self, n_target: int = 6) -> List[Dict]:
        """
        Lactobacillus phages — still present but composition shifts
        in HPV+ women as Lactobacillus dominance decreases (CST IV
        transition in persistent infection).
        """
        logger.info("Selecting Lactobacillus phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Lactobacillus phage%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Lactobacillus phages: {len(results)}")
        return results

    def get_bv_associated_phages(self, n_target: int = 5) -> List[Dict]:
        """
        Phages associated with BV-related bacteria. HPV persistence
        correlates with vaginal dysbiosis (Brusselaers 2019). No
        Gardnerella phages in RefSeq, so we use Prevotella, E. coli,
        and Enterococcus phages as proxies for the dysbiotic phageome.
        """
        logger.info("Selecting BV-associated phages...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Prevotella%phage%'
           OR g.genome_name LIKE 'Enterococcus phage%'
           OR (g.genome_name LIKE 'Escherichia phage%'
               AND g.genome_name NOT LIKE '%phiX174%')
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  BV-associated phages: {len(results)}")
        return results

    def get_hsv(self) -> List[Dict]:
        """
        HSV-1 — sexually transmitted co-infection. HSV-2 is more
        common vaginally but absent from RefSeq; HSV-1 genital
        infections are increasingly prevalent.
        """
        logger.info("Selecting HSV...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human herpesvirus 1%'
        LIMIT 1
        """
        results = [dict(row) for row in self.conn.execute(query)]
        logger.info(f"  HSV: {len(results)}")
        return results

    def get_other_herpesviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Other herpesviruses detected in the vaginal tract: CMV (HHV-5),
        EBV (HHV-4), HHV-6. Reactivation more common in HPV+ women.
        """
        logger.info("Selecting other herpesviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human herpesvirus 5%'
           OR g.genome_name LIKE 'Human gammaherpesvirus 4%'
           OR g.genome_name LIKE 'Human herpesvirus 6B%'
           OR g.genome_name LIKE 'Human betaherpesvirus 6A%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Other herpesviruses: {len(results)}")
        return results

    def get_anelloviridae(self, n_target: int = 4) -> List[Dict]:
        """
        Anelloviridae (TTV) — altered prevalence in HPV+ women,
        potentially reflecting immune modulation.
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

    def get_polyomaviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Human polyomaviruses — MCPyV and HPyV detected in vaginal
        samples, co-infection with HPV reported.
        """
        logger.info("Selecting polyomaviruses...")
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus,
               t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human polyomavirus%'
           OR g.genome_name LIKE 'Merkel cell polyomavirus%'
        ORDER BY seeded_rand()
        LIMIT ?
        """
        results = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Polyomaviruses: {len(results)}")
        return results

    def assign_hpv_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        HPV-infected vaginal abundance: HPV types dominate the
        eukaryotic viral fraction. HPV 16 typically the most abundant
        in high-risk infections. Phage community more diverse but
        each at lower abundance than in healthy.
        """
        n = len(genomes)
        mu = -2.0
        sigma = 2.0
        raw = np.random.lognormal(mu, sigma, n)

        max_raw = raw.max()
        for i, g in enumerate(genomes):
            name = g.get('genome_name', '')
            m = re.search(r'papillomavirus type (\d+)', name)
            if m:
                hpv_type = int(m.group(1))
                if hpv_type == 16:
                    raw[i] = max_raw * 10.0
                elif hpv_type == 18:
                    raw[i] = max_raw * 6.0
                elif hpv_type in {31, 33, 35, 45, 52, 58}:
                    raw[i] *= 3.0
                else:
                    raw[i] *= 2.0

        normalized = raw / raw.sum()
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)
        return genomes

    def create_collection(self) -> List[Dict]:
        logger.info("=" * 80)
        logger.info("CURATING HPV-INFECTED VAGINAL VIROME")
        logger.info("=" * 80)

        all_hpv = self._fetch_all_hpv()
        hpv_hr = self.get_hpv_high_risk(all_hpv, n_target=10)
        hpv_lr = self.get_hpv_low_risk(all_hpv, n_target=4)
        selected_ids = {g['genome_id'] for g in hpv_hr + hpv_lr}
        hpv_other = self.get_hpv_other(all_hpv, selected_ids, n_target=4)
        lacto = self.get_lactobacillus_phages(n_target=6)
        bv_phages = self.get_bv_associated_phages(n_target=5)
        hsv = self.get_hsv()
        herpes = self.get_other_herpesviruses(n_target=3)
        anello = self.get_anelloviridae(n_target=4)
        polyoma = self.get_polyomaviruses(n_target=2)

        collection = (hpv_hr + hpv_lr + hpv_other + lacto + bv_phages +
                      hsv + herpes + anello + polyoma)

        seen_ids = set()
        unique = []
        for g in collection:
            if g['genome_id'] not in seen_ids:
                unique.append(g)
                seen_ids.add(g['genome_id'])
        collection = unique

        logger.info(f"\nPre-abundance total: {len(collection)} genomes")
        collection = self.assign_hpv_abundances(collection)
        collection.sort(key=lambda x: x['relative_abundance'], reverse=True)
        for i, g in enumerate(collection, 1):
            g['abundance_rank'] = i

        logger.info(f"Total genomes in HPV-infected vaginal collection: {len(collection)}")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        # Verify HPV 16 is present and dominant
        hpv16 = [g for g in collection if 'type 16' in g.get('genome_name', '')]
        if hpv16:
            logger.info(f"\nHPV 16 abundance: {hpv16[0]['relative_abundance']:.4f} (rank {hpv16[0]['abundance_rank']})")
        else:
            logger.warning("HPV 16 NOT FOUND — scientifically invalid!")

        logger.info("\nTop 10 most abundant genomes:")
        for i, g in enumerate(collection[:10], 1):
            fam = g.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {g['genome_name'][:50]:50s} {g['relative_abundance']:.6f} ({fam})")

        # Family breakdown
        families = {}
        for g in collection:
            fam = g.get('family', 'Unknown')
            families[fam] = families.get(fam, 0) + 1
        logger.info("\nFamily breakdown:")
        for fam, count in sorted(families.items(), key=lambda x: -x[1]):
            logger.info(f"  {fam}: {count}")

        # Compare to healthy
        healthy_count = self.conn.execute(
            "SELECT n_genomes FROM body_site_collections WHERE collection_id = 16"
        ).fetchone()
        if healthy_count:
            logger.info(f"\nComparison: healthy vaginal has {healthy_count[0]} genomes, "
                        f"HPV-infected has {len(collection)}")

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
            'HPV-Infected Vaginal Virome',
            'Vaginal virome of women with HPV infection. HPV dominates the '
            'eukaryotic viral fraction with high-risk types (16, 18, 31, 33, '
            '45, 52, 58) and low-risk types (6, 11). Phageome shifted toward '
            'BV-associated bacteria reflecting dysbiosis linked to HPV '
            'persistence. HSV and herpesvirus co-infections included. '
            'Compare to Collection 16 (healthy vaginal) to study HPV impact. '
            'Based on Wylie et al. 2014, Brusselaers et al. 2019.',
            len(collection),
            'HPV-dominated (high-risk types 16/18 most abundant). '
            'Multi-type co-infection modeled (12+ HPV types). '
            'BV-associated phages (dysbiotic shift). HSV co-infection. '
            'Altered Anelloviridae and polyomavirus prevalence.',
            'ViroForge Pediatric/Gynecologic',
            '2026-09-01',
            'Wylie et al. 2014 (PLOS Pathogens), '
            'Virtanen et al. 2016 (International Journal of Cancer), '
            'Brusselaers et al. 2019 (BMC Medicine), '
            'Laniewski et al. 2020 (mBio)',
            'Homo sapiens (female)',
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
    curator = HPVVaginalCurator()
    try:
        collection = curator.create_collection()
        curator.insert_collection(collection)
        logger.info("\n" + "=" * 80)
        logger.info("HPV-INFECTED VAGINAL VIROME CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nKnown limitations:")
        logger.info("  - No HSV-2 in RefSeq (HSV-1 used as proxy)")
        logger.info("  - No Gardnerella phages in RefSeq (Prevotella/Enterococcus phages as proxy)")
        logger.info("  - Lactobacillus phages are from dairy species (see Collection 16 notes)")
    finally:
        curator.close()


if __name__ == '__main__':
    main()
