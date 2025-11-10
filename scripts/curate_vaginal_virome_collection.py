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

    def get_herpesviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Get human herpesviruses.

        Includes HSV-1 (HHV-1), VZV (HHV-3), EBV (HHV-4), CMV (HHV-5),
        HHV-6A/6B, HHV-7, KSHV (HHV-8).
        NOTE: HSV-2 (HHV-2) is not in RefSeq database.

        Family: Orthoherpesviridae (modern ICTV name, not "Herpesviridae")
        Genome type: dsDNA linear
        """
        logger.info("Selecting human herpesviruses (HSV-1, VZV, EBV, CMV, KSHV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Human%herpesvirus%'
           OR g.genome_name LIKE '%Human%herpes%'
           OR g.genome_name LIKE '%Human betaherpes%'
           OR g.genome_name LIKE '%Human gammaherpes%'
           OR g.genome_name LIKE '%Human alphaherpes%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        herpesviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Orthoherpesviridae (Human): {len(herpesviruses)}")

        # Log what we found
        for hv in herpesviruses:
            logger.info(f"    - {hv['genome_name'][:60]}")

        return herpesviruses

    def get_lactobacillus_phages(self, n_target: int = 8) -> List[Dict]:
        """
        Get Lactobacillus bacteriophages.

        NOTE: Old morphology-based families (Siphoviridae, Myoviridae, Podoviridae)
        were abolished by ICTV in 2022 and replaced with genome-based families.
        Lactobacillus phages now classified in families like Herelleviridae or
        have "Unknown" family (not yet reclassified).

        64 Lactobacillus phages available in database:
        - 12 in Herelleviridae (genome-based family)
        - 52 with Unknown family (awaiting reclassification)

        Targets Lactobacillus (L. crispatus, L. iners, L. jensenii, L. gasseri).
        Dominant in healthy vaginal microbiome, reduced in bacterial vaginosis.

        Genome type: dsDNA
        """
        logger.info("Selecting Lactobacillus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Lactobacillus phage%'
           OR g.genome_name LIKE '%Lactobacillus prophage%'
           OR g.genome_name LIKE '%lactobacillus phage%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        lactobacillus_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Lactobacillus phages: {len(lactobacillus_phages)}")

        # Log family distribution
        from collections import Counter
        families = Counter(p.get('family') or 'Unknown' for p in lactobacillus_phages)
        for family, count in families.most_common():
            logger.info(f"    - {family}: {count} phages")

        return lactobacillus_phages

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
            genome_name = genome.get('genome_name', '')

            # Determine category by content, not just family name
            # (since many phages have Unknown family)
            is_phage = ('Lactobacillus phage' in genome_name or
                        'phage' in genome_name.lower() or
                        family in ['Herelleviridae', 'Microviridae'])

            if is_phage:
                # Bacteriophages: 80% total, distributed among phages
                # Note: Many Lactobacillus phages have "Unknown" family
                n_phages = categories['bacteriophages']
                base_abundance = bacteriophages_abundance / n_phages if n_phages > 0 else 0
                # Add variation: log-normal distribution
                abundance = base_abundance * np.random.lognormal(0, 0.5) if n_phages > 0 else 0

            elif family == 'Papillomaviridae':
                # Papillomaviruses: 10% total
                n_hpv = categories['papillomaviruses']
                base_abundance = papilloma_abundance / n_hpv if n_hpv > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_hpv > 0 else 0

            elif family == 'Anelloviridae':
                # Anelloviruses: 5% total
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello if n_anello > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_anello > 0 else 0

            elif family == 'Orthoherpesviridae' or 'herpesvirus' in genome_name.lower():
                # Herpesviruses: 3% total (episodic shedding)
                # Note: Modern ICTV name is "Orthoherpesviridae"
                n_herpes = categories['herpesviruses']
                base_abundance = herpes_abundance / n_herpes if n_herpes > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.4) if n_herpes > 0 else 0

            else:
                # Other viruses: 2% total (polyoma, adeno)
                n_other = categories['other']
                base_abundance = other_abundance / n_other if n_other > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_other > 0 else 0

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
        herpesviruses = self.get_herpesviruses(n_target=5)  # Increased from 3
        lactobacillus_phages = self.get_lactobacillus_phages(n_target=8)  # NEW: Replaces sipho/myo
        microviridae = self.get_microviridae_phages(n_target=2)
        polyomaviruses = self.get_polyomaviruses(n_target=2)
        adenoviruses = self.get_adenoviruses(n_target=1)

        # Combine all genomes
        all_genomes = (
            papillomaviruses +
            anelloviruses +
            herpesviruses +
            lactobacillus_phages +  # NEW
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
            'bacteriophages': len(lactobacillus_phages) + len(microviridae),  # Updated
            'other': len(polyomaviruses) + len(adenoviruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Insert collection
        collection_meta = {
            'collection_id': 24,
            'collection_name': 'Vaginal Virome (Healthy)',
            'description': f'Healthy human vaginal/cervicovaginal virome with {len(all_genomes)} genomes. Includes papillomaviruses (HPV), human herpesviruses (HSV-1, VZV, EBV, CMV, KSHV), anelloviruses, Lactobacillus bacteriophages, polyomaviruses, and adenoviruses. Bacteriophages ~80% abundance (literature-consistent), HPV 78.3% prevalence, Anelloviruses 69.6% prevalence. Host: Homo sapiens, Body site: Vagina. Women\'s health applications: cervical cancer screening, bacterial vaginosis, pregnancy outcomes, anellovirus immune markers.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Wylie et al. 2014 (BMC Biol), Dols et al. 2016, Nature Microbiology 2024 VMGC study (4,263 viral OTUs, 85.8% vaginal-specific), Frontiers 2025 (267 women). Includes HPV (5 types), human herpesviruses (HSV-1, VZV, EBV, CMV, HHV-6A/6B, HHV-7, KSHV), anelloviruses (3 types), and Lactobacillus bacteriophages (64 available in RefSeq). NOTE: Family names updated for ICTV 2022 taxonomy (Orthoherpesviridae, Herelleviridae).',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-09',
            'literature_references': 'Wylie et al. 2014 BMC Biol; Dols et al. 2016; Fu et al. 2024 Nat Microbiol (VMGC); Zhang et al. 2025 Front Cell Infect Microbiol; ICTV 2022 taxonomy update (bacteriophage families)',
            'version': 2
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
