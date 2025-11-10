#!/usr/bin/env python3
"""
Curate Collection 24: Vaginal Virome

Composition based on literature (Wylie et al. 2014, Dols et al. 2016,
Nature Microbiology 2024 VMGC study):
- Papillomaviruses (HPV) - 78.3% prevalence, most common eukaryotic virus
- Anelloviruses - 69.6% prevalence, ubiquitous commensal viruses
- Herpesviruses (HSV-1, HSV-2, CMV, EBV)
- Bacteriophages (Lactobacillus phages - Siphoviridae, Myoviridae, Podoviridae)
- Other: Polyomaviruses (BK, JC), Adenoviruses

Represents DNA virome of healthy vaginal microbiome (cervicovaginal samples).
Bacteriophages dominate (~80% of reads), Lactobacillus phages reflect dominance
of vaginal Lactobacillus (L. crispatus, L. iners, L. jensenii, L. gasseri).

Literature basis:
- Wylie et al. 2014 (BMC Biol): Metagenomic analysis of dsDNA viruses
- Dols et al. 2016: Vaginal virome and bacterial vaginosis
- Nature Microbiology 2024: Vaginal Microbial Genome Collection (VMGC)
- Frontiers 2025: Vaginal virome diversity and vaginitis
- Virology Journal 2020: Comparison of viromes in vaginal secretion

Target size: 20-30 genomes (as per roadmap)

Women's health impact:
- HPV: Cervical cancer screening, vaccine effectiveness
- Bacterial vaginosis: Altered phageome, reduced Lactobacillus phages
- Pregnancy outcomes: Virome dysbiosis associations
- Transplant monitoring: Anellovirus immune status marker

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


class VaginalViromeCurator:
    """Curate vaginal virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_papillomaviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Get human papillomaviruses (HPV).

        HPV is the most prevalent eukaryotic virus in the vaginal virome
        (78.3% of women). Includes high-risk types (HPV16, 18) and low-risk
        types (HPV6, 11). Critical for cervical cancer screening.

        Genome type: dsDNA circular
        """
        logger.info("Selecting papillomaviruses (HPV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Papillomaviridae'
          AND (g.genome_name LIKE '%Human papillomavirus%'
           OR g.genome_name LIKE '%human papillomavirus%'
           OR g.genome_name LIKE '%HPV%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        papillomaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Papillomaviridae (HPV): {len(papillomaviruses)}")
        return papillomaviruses

    def get_anelloviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get anelloviruses (TTV, TTMV, TTMDV).

        Anelloviruses are highly prevalent (69.6%) ubiquitous commensal viruses.
        Small circular ssDNA viruses, not pathogenic. Marker of immune status.
        Extensive diversity in vaginal samples.

        Genome type: ssDNA circular
        """
        logger.info("Selecting anelloviruses (TTV family)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Anelloviridae'
          AND (g.genome_name LIKE '%Torque teno%'
           OR g.genome_name LIKE '%torque teno%'
           OR g.genome_name LIKE '%TTV%'
           OR g.genome_name LIKE '%TTMV%'
           OR t.genus IN ('Alphatorquevirus', 'Betatorquevirus', 'Gammatorquevirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        anelloviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Anelloviridae (TTV): {len(anelloviruses)}")
        return anelloviruses

    def get_herpesviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get human herpesviruses.

        Includes HSV-1, HSV-2 (genital herpes), CMV, EBV, VZV.
        HSV-2 is sexually transmitted, common in reproductive tract.
        HSV-1 also increasingly detected in genital infections.

        Genome type: dsDNA linear
        """
        logger.info("Selecting herpesviruses (HSV-1/2, CMV, EBV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Herpesviridae'
          AND (g.genome_name LIKE '%Human alphaherpesvirus 1%'  -- HSV-1
           OR g.genome_name LIKE '%Human alphaherpesvirus 2%'   -- HSV-2
           OR g.genome_name LIKE '%herpes simplex%'
           OR g.genome_name LIKE '%Human betaherpesvirus 5%'    -- CMV
           OR g.genome_name LIKE '%cytomegalovirus%'
           OR g.genome_name LIKE '%Human gammaherpesvirus 4%'   -- EBV
           OR g.genome_name LIKE '%Epstein-Barr%'
           OR t.genus IN ('Simplexvirus', 'Cytomegalovirus', 'Lymphocryptovirus'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        herpesviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Herpesviridae: {len(herpesviruses)}")
        return herpesviruses

    def get_lactobacillus_siphoviruses(self, n_target: int = 4) -> List[Dict]:
        """
        Get Lactobacillus Siphoviridae phages.

        Siphoviridae is the dominant phage family in healthy vaginal microbiome.
        Targets Lactobacillus (L. crispatus, L. iners, L. jensenii, L. gasseri).
        Reduced in bacterial vaginosis.

        Genome type: dsDNA
        """
        logger.info("Selecting Lactobacillus Siphoviridae phages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Siphoviridae'
          AND (g.genome_name LIKE '%Lactobacillus%'
           OR g.genome_name LIKE '%lactobacillus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        siphoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Siphoviridae (Lactobacillus): {len(siphoviruses)}")

        # If no specific Lactobacillus phages, get general Siphoviridae
        if len(siphoviruses) < n_target:
            additional = n_target - len(siphoviruses)
            logger.info(f"  Finding {additional} additional Siphoviridae phages...")

            query = """
            SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
            FROM genomes g
            JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE t.family = 'Siphoviridae'
            ORDER BY RANDOM()
            LIMIT ?
            """

            additional_phages = [dict(row) for row in self.conn.execute(query, (additional,))]
            siphoviruses.extend(additional_phages)
            logger.info(f"  Added {len(additional_phages)} general Siphoviridae")

        return siphoviruses

    def get_lactobacillus_myoviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get Lactobacillus Myoviridae phages.

        Myoviridae is another major phage family in vaginal virome.
        Also targets Lactobacillus species.

        Genome type: dsDNA
        """
        logger.info("Selecting Lactobacillus Myoviridae phages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Myoviridae'
          AND (g.genome_name LIKE '%Lactobacillus%'
           OR g.genome_name LIKE '%lactobacillus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        myoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Myoviridae (Lactobacillus): {len(myoviruses)}")

        # If no specific Lactobacillus phages, get general Myoviridae
        if len(myoviruses) < n_target:
            additional = n_target - len(myoviruses)
            logger.info(f"  Finding {additional} additional Myoviridae phages...")

            query = """
            SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
            FROM genomes g
            JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE t.family = 'Myoviridae'
            ORDER BY RANDOM()
            LIMIT ?
            """

            additional_phages = [dict(row) for row in self.conn.execute(query, (additional,))]
            myoviruses.extend(additional_phages)
            logger.info(f"  Added {len(additional_phages)} general Myoviridae")

        return myoviruses

    def get_microviridae_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Microviridae phages.

        Microviridae is also a leading bacteriophage family in vaginal virome
        (alongside Siphoviridae). Small ssDNA phages.

        Genome type: ssDNA circular
        """
        logger.info("Selecting Microviridae phages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Microviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """

        microviridae = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Microviridae: {len(microviridae)}")
        return microviridae

    def get_polyomaviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human polyomaviruses (BK, JC).

        Polyomaviruses are part of the eukaryotic virome.
        BK and JC viruses are ubiquitous, often detected in healthy individuals.

        Genome type: dsDNA circular
        """
        logger.info("Selecting polyomaviruses (BK, JC)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Polyomaviridae'
          AND (g.genome_name LIKE '%Human polyomavirus%'
           OR g.genome_name LIKE '%BK polyomavirus%'
           OR g.genome_name LIKE '%JC polyomavirus%'
           OR g.genome_name LIKE '%BKV%'
           OR g.genome_name LIKE '%JCV%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        polyomaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Polyomaviridae: {len(polyomaviruses)}")
        return polyomaviruses

    def get_adenoviruses(self, n_target: int = 1) -> List[Dict]:
        """
        Get human adenoviruses.

        Adenoviruses are occasionally detected in vaginal samples.
        Part of the eukaryotic virome.

        Genome type: dsDNA linear
        """
        logger.info("Selecting adenoviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
          AND (g.genome_name LIKE '%Human adenovirus%'
           OR g.genome_name LIKE '%human adenovirus%'
           OR g.genome_name LIKE '%HAdV%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        adenoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Adenoviridae: {len(adenoviruses)}")
        return adenoviruses

    def assign_abundances(self, genomes: List[Dict], categories: Dict) -> List[Dict]:
        """
        Assign relative abundances to genomes based on literature.

        Vaginal virome composition (healthy):
        - Bacteriophages: ~80% (Siphoviridae, Microviridae dominant)
        - Papillomaviruses: High prevalence (78.3%), moderate abundance
        - Anelloviruses: High prevalence (69.6%), moderate abundance
        - Herpesviruses: Low abundance (episodic shedding)
        - Polyoma/Adenoviruses: Low abundance

        Abundance model: Moderate skew (phages dominant)
        """
        logger.info("Assigning relative abundances...")

        n_total = len(genomes)

        # Calculate category abundances (literature-based)
        bacteriophages_abundance = 0.80  # 80% phages (Siphoviridae + Microviridae + Myoviridae)
        papilloma_abundance = 0.10       # 10% HPV (high prevalence, moderate abundance)
        anello_abundance = 0.05          # 5% Anelloviruses
        herpes_abundance = 0.03          # 3% Herpesviruses (episodic)
        other_abundance = 0.02           # 2% Other (polyoma, adeno)

        # Distribute within categories
        abundances = []

        for genome in genomes:
            family = genome.get('family', '')

            if family in ['Siphoviridae', 'Myoviridae', 'Podoviridae', 'Microviridae']:
                # Bacteriophages: 80% total, distributed among phages
                n_phages = categories['bacteriophages']
                base_abundance = bacteriophages_abundance / n_phages
                # Add variation: log-normal distribution
                abundance = base_abundance * np.random.lognormal(0, 0.5)

            elif family == 'Papillomaviridae':
                # Papillomaviruses: 10% total
                n_hpv = categories['papillomaviruses']
                base_abundance = papilloma_abundance / n_hpv
                abundance = base_abundance * np.random.lognormal(0, 0.3)

            elif family == 'Anelloviridae':
                # Anelloviruses: 5% total
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello
                abundance = base_abundance * np.random.lognormal(0, 0.3)

            elif family == 'Herpesviridae':
                # Herpesviruses: 3% total (episodic shedding)
                n_herpes = categories['herpesviruses']
                base_abundance = herpes_abundance / n_herpes
                abundance = base_abundance * np.random.lognormal(0, 0.4)

            else:
                # Other viruses: 2% total
                n_other = categories['other']
                base_abundance = other_abundance / n_other
                abundance = base_abundance * np.random.lognormal(0, 0.3)

            abundances.append(abundance)

        # Normalize to sum to 1.0
        total = sum(abundances)
        abundances = [a / total for a in abundances]

        # Add abundances to genomes
        for genome, abundance in zip(genomes, abundances):
            genome['relative_abundance'] = abundance

        # Rank by abundance
        genomes_sorted = sorted(genomes, key=lambda x: x['relative_abundance'], reverse=True)
        for rank, genome in enumerate(genomes_sorted, 1):
            genome['abundance_rank'] = rank

        logger.info(f"  Total abundance: {sum(abundances):.6f}")
        logger.info(f"  Top 3 by abundance:")
        for genome in genomes_sorted[:3]:
            logger.info(f"    {genome['genome_name'][:60]:60s} {genome['relative_abundance']:.4f}")

        return genomes

    def create_collection(self):
        """Create Collection 24: Vaginal Virome."""
        logger.info("=" * 70)
        logger.info("Creating Collection 24: Vaginal Virome")
        logger.info("=" * 70)

        # Collect genomes by category
        papillomaviruses = self.get_papillomaviruses(n_target=5)
        anelloviruses = self.get_anelloviruses(n_target=3)
        herpesviruses = self.get_herpesviruses(n_target=3)
        siphoviruses = self.get_lactobacillus_siphoviruses(n_target=4)
        myoviruses = self.get_lactobacillus_myoviruses(n_target=2)
        microviridae = self.get_microviridae_phages(n_target=2)
        polyomaviruses = self.get_polyomaviruses(n_target=2)
        adenoviruses = self.get_adenoviruses(n_target=1)

        # Combine all genomes
        all_genomes = (
            papillomaviruses +
            anelloviruses +
            herpesviruses +
            siphoviruses +
            myoviruses +
            microviridae +
            polyomaviruses +
            adenoviruses
        )

        logger.info(f"\nTotal genomes collected: {len(all_genomes)}")

        if len(all_genomes) == 0:
            logger.error("No genomes found! Check database content.")
            return

        # Track category sizes for abundance assignment
        categories = {
            'papillomaviruses': len(papillomaviruses),
            'anelloviruses': len(anelloviruses),
            'herpesviruses': len(herpesviruses),
            'bacteriophages': len(siphoviruses) + len(myoviruses) + len(microviridae),
            'other': len(polyomaviruses) + len(adenoviruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Insert collection
        collection_meta = {
            'collection_id': 24,
            'collection_name': 'Vaginal Virome (Healthy)',
            'description': 'Healthy human vaginal/cervicovaginal virome with papillomaviruses, anelloviruses, herpesviruses, and Lactobacillus bacteriophages. Bacteriophages ~80% (Lactobacillus-targeting), HPV 78.3% prevalence, Anelloviruses 69.6% prevalence. Host: Homo sapiens, Body site: Vagina. Women\'s health applications: cervical cancer screening, bacterial vaginosis, pregnancy outcomes.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Wylie et al. 2014 (BMC Biol), Dols et al. 2016, Nature Microbiology 2024 VMGC study, Frontiers 2025. Includes HPV (high prevalence 78.3%), anelloviruses (69.6%), herpesviruses (HSV-1/2), and Lactobacillus bacteriophages.',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-09',
            'literature_references': 'Wylie et al. 2014 BMC Biol: Metagenomic analysis of dsDNA viruses; Dols et al. 2016: Vaginal virome and BV; Nature Microbiology 2024: Vaginal Microbial Genome Collection; Frontiers 2025: Vaginal virome diversity',
            'version': 1
        }

        # Insert into database
        self._insert_collection(collection_meta, all_genomes)

        logger.info("\n" + "=" * 70)
        logger.info("Collection 24: Vaginal Virome - COMPLETE")
        logger.info("=" * 70)

    def _insert_collection(self, collection_meta: Dict, genomes: List[Dict]):
        """Insert collection and genomes into database."""
        logger.info("\nInserting collection into database...")

        cursor = self.conn.cursor()

        # Insert collection metadata
        cursor.execute("""
            INSERT INTO body_site_collections (
                collection_id, collection_name, description, n_genomes,
                selection_criteria, curated_by, curation_date,
                literature_references, version
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
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

        # Insert collection genomes
        for genome in genomes:
            cursor.execute("""
                INSERT INTO collection_genomes (
                    collection_id, genome_id, relative_abundance, abundance_rank
                ) VALUES (?, ?, ?, ?)
            """, (
                collection_meta['collection_id'],
                genome['genome_id'],
                genome['relative_abundance'],
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"  Collection metadata inserted")
        logger.info(f"  {len(genomes)} genome associations inserted")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main function."""
    curator = VaginalViromeCurator()

    try:
        curator.create_collection()
    finally:
        curator.close()

    logger.info("\nDone!")


if __name__ == '__main__':
    main()
