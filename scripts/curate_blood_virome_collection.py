#!/usr/bin/env python3
"""
Curate Collection 25: Blood/Plasma Virome (Healthy)

Literature basis:
- Moustafa et al. 2017 (PLOS Pathogens) - Blood DNA virome in 8,000 humans
- Young et al. 2016 - Polyomavirus monitoring in transplant patients
- PMC 2023 - Atlas of blood virome in healthy individuals

Composition (prevalence in healthy blood):
- Anelloviruses (TTV, TLMV): 8.91% prevalence (most abundant)
- Herpesviruses: 14-20% prevalence
  - HHV-7: 20.37% (most common)
  - EBV (HHV-4): 14.45%
  - HHV-6B: 4.80%
  - HHV-6A: 1.47%
  - CMV (HHV-5): 0.35%
  - HSV-1 (HHV-1): 0.12%
  - KSHV (HHV-8): 0.04%
- Polyomaviruses: Important for transplant monitoring
  - BK virus: 26.4% in transplant patients
  - JC virus: 18% in blood donors
  - Merkel cell polyomavirus: 0.59%
- Parvovirus B19: 0.12% (can have very high viral loads)
- HPV: 0.19% (7 types detected)
- Adenoviruses: Low prevalence
- Hepatitis viruses (HBV, HCV): <0.02%

Target: 15-25 genomes
Applications: Viremia detection, transplant monitoring, blood safety screening
"""

import sqlite3
import logging
from pathlib import Path
from typing import List, Dict
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BloodViromeCurator:
    """Curator for blood/plasma virome collection."""

    def __init__(self):
        """Initialize connection to viral genomes database."""
        db_path = Path(__file__).parent.parent / "viroforge" / "data" / "viral_genomes.db"
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        logger.info(f"Connected to database: {db_path}")

    def get_anelloviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Get anelloviruses (Torque teno viruses).

        Most abundant virus family in blood (8.91% prevalence).
        Includes TTV (Torque teno virus) and TLMV (Torque teno midi virus).

        Family: Anelloviridae
        Genome type: Circular ssDNA
        """
        logger.info("Selecting anelloviruses (TTV, TLMV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Anelloviridae'
           OR g.genome_name LIKE '%Torque teno%'
           OR g.genome_name LIKE '%TTV%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        anelloviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Anelloviridae: {len(anelloviruses)}")

        for av in anelloviruses:
            logger.info(f"    - {av['genome_name'][:60]}")

        return anelloviruses

    def get_herpesviruses(self, n_target: int = 6) -> List[Dict]:
        """
        Get human herpesviruses found in blood.

        Includes HHV-7 (20.37%), EBV/HHV-4 (14.45%), HHV-6B (4.80%),
        HHV-6A (1.47%), CMV/HHV-5 (0.35%), HSV-1/HHV-1 (0.12%), KSHV/HHV-8 (0.04%).

        Family: Orthoherpesviridae (modern ICTV name, not "Herpesviridae")
        Genome type: dsDNA linear
        """
        logger.info("Selecting human herpesviruses (HHV-7, EBV, CMV, HHV-6)...")

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

        for hv in herpesviruses:
            logger.info(f"    - {hv['genome_name'][:60]}")

        return herpesviruses

    def get_polyomaviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get human polyomaviruses.

        Critical for transplant monitoring.
        Includes BK virus (26.4% in transplant patients), JC virus (18% in donors),
        Merkel cell polyomavirus (0.59%).

        Family: Polyomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting polyomaviruses (BK, JC, Merkel cell)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Polyomaviridae'
           AND (g.genome_name LIKE '%BK polyomavirus%'
                OR g.genome_name LIKE '%JC polyomavirus%'
                OR g.genome_name LIKE '%Merkel cell polyomavirus%'
                OR g.genome_name LIKE '%Human polyomavirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        polyomaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Polyomaviridae (Human): {len(polyomaviruses)}")

        for pv in polyomaviruses:
            logger.info(f"    - {pv['genome_name'][:60]}")

        return polyomaviruses

    def get_parvoviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human parvoviruses.

        Parvovirus B19: 0.12% prevalence but can reach very high viral loads
        (up to 302 million copies per 100,000 cells).

        Family: Parvoviridae
        Genome type: ssDNA linear
        """
        logger.info("Selecting parvoviruses (Parvovirus B19)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Parvoviridae'
           AND (g.genome_name LIKE '%Human parvovirus B19%'
                OR g.genome_name LIKE '%Human parvovirus 4%'
                OR g.genome_name LIKE '%Primate erythroparvovirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        parvoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Parvoviridae (Human): {len(parvoviruses)}")

        for pv in parvoviruses:
            logger.info(f"    - {pv['genome_name'][:60]}")

        return parvoviruses

    def get_papillomaviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human papillomaviruses.

        7 HPV types detected in blood (0.19% prevalence).
        Lower prevalence than in other body sites.

        Family: Papillomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting papillomaviruses (HPV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Papillomaviridae'
           AND g.genome_name LIKE '%Human papillomavirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        papillomaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Papillomaviridae (Human): {len(papillomaviruses)}")

        for pv in papillomaviruses:
            logger.info(f"    - {pv['genome_name'][:60]}")

        return papillomaviruses

    def get_adenoviruses(self, n_target: int = 1) -> List[Dict]:
        """
        Get human adenoviruses.

        Low prevalence in healthy blood, but important for immunocompromised.

        Family: Adenoviridae
        Genome type: dsDNA linear
        """
        logger.info("Selecting adenoviruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
           AND g.genome_name LIKE '%Human adenovirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        adenoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Adenoviridae (Human): {len(adenoviruses)}")

        for av in adenoviruses:
            logger.info(f"    - {av['genome_name'][:60]}")

        return adenoviruses

    def get_bloodborne_viruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get bloodborne viruses (HCV, HIV).

        HCV: Low prevalence in healthy populations, important for blood safety
        HIV: 0.06% prevalence (Moustafa et al. 2017)

        NOTE: Human HBV not available in RefSeq database (only animal HBV)
        NOTE: HCV has family="Unknown" in database, not "Flaviviridae"

        Families: Unknown (HCV), Retroviridae (HIV)
        """
        logger.info("Selecting bloodborne viruses (HCV, HIV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE 'Hepatitis C virus%')
           OR (g.genome_name LIKE 'Human immunodeficiency virus%' OR g.genome_name LIKE 'HIV-1%' OR g.genome_name LIKE 'HIV-2%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        bloodborne = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Bloodborne viruses: {len(bloodborne)}")

        for bv in bloodborne:
            logger.info(f"    - {bv['genome_name'][:60]}")

        return bloodborne

    def assign_abundances(self, genomes: List[Dict], categories: Dict[str, int]) -> List[Dict]:
        """
        Assign relative abundances based on blood virome composition.

        Distribution (literature-based):
        - Anelloviruses: 40% (most abundant, 8.91% prevalence)
        - Herpesviruses: 35% (HHV-7 20.37%, EBV 14.45%)
        - Polyomaviruses: 15% (important for transplant monitoring)
        - Others: 10% (parvoviruses, HPV, adenoviruses, hepatitis)
        """
        logger.info("\nAssigning relative abundances...")

        # Target abundances (will be perturbed with log-normal)
        anello_abundance = 0.40
        herpes_abundance = 0.35
        polyoma_abundance = 0.15
        other_abundance = 0.10

        for genome in genomes:
            family = genome.get('family') or 'Unknown'
            genome_name = genome.get('genome_name', '')

            # Categorize by family
            if family == 'Anelloviridae' or 'Torque teno' in genome_name:
                # Anelloviruses: Most abundant (40%)
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello if n_anello > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_anello > 0 else 0

            elif family == 'Orthoherpesviridae' or 'herpesvirus' in genome_name.lower():
                # Herpesviruses: High abundance (35%)
                # HHV-7 and EBV should be most abundant
                n_herpes = categories['herpesviruses']
                base_abundance = herpes_abundance / n_herpes if n_herpes > 0 else 0

                # Boost HHV-7 and EBV
                if 'herpesvirus 7' in genome_name.lower() or 'HHV-7' in genome_name:
                    abundance = base_abundance * np.random.lognormal(0.5, 0.3)  # Higher
                elif 'gammaherpesvirus 4' in genome_name or 'EBV' in genome_name:
                    abundance = base_abundance * np.random.lognormal(0.3, 0.3)  # High
                else:
                    abundance = base_abundance * np.random.lognormal(0, 0.3)

            elif family == 'Polyomaviridae':
                # Polyomaviruses: Moderate (15%)
                n_polyoma = categories['polyomaviruses']
                base_abundance = polyoma_abundance / n_polyoma if n_polyoma > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_polyoma > 0 else 0

            else:
                # Others (parvoviruses, HPV, adenoviruses, hepatitis): Low (10%)
                n_other = categories['other']
                base_abundance = other_abundance / n_other if n_other > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.4) if n_other > 0 else 0

            genome['relative_abundance'] = abundance

        # Normalize to sum to 1.0
        total = sum(g['relative_abundance'] for g in genomes)
        for genome in genomes:
            genome['relative_abundance'] /= total

        # Sort by abundance and assign ranks
        genomes.sort(key=lambda x: x['relative_abundance'], reverse=True)
        for rank, genome in enumerate(genomes, 1):
            genome['abundance_rank'] = rank

        # Log abundance distribution
        logger.info("\nAbundance distribution:")
        for genome in genomes[:10]:  # Top 10
            logger.info(f"  {genome['abundance_rank']:2d}. {genome['genome_name'][:50]:50s} {genome['relative_abundance']:7.1%}")

        if len(genomes) > 10:
            logger.info(f"  ... and {len(genomes) - 10} more genomes")

        return genomes

    def create_collection(self):
        """Create Collection 25: Blood/Plasma Virome."""
        logger.info("=" * 70)
        logger.info("Creating Collection 25: Blood/Plasma Virome")
        logger.info("=" * 70)

        # Collect genomes by category
        anelloviruses = self.get_anelloviruses(n_target=5)
        herpesviruses = self.get_herpesviruses(n_target=6)
        polyomaviruses = self.get_polyomaviruses(n_target=3)
        parvoviruses = self.get_parvoviruses(n_target=2)
        papillomaviruses = self.get_papillomaviruses(n_target=2)
        adenoviruses = self.get_adenoviruses(n_target=1)
        bloodborne_viruses = self.get_bloodborne_viruses(n_target=2)

        # Combine all genomes
        all_genomes = (
            anelloviruses +
            herpesviruses +
            polyomaviruses +
            parvoviruses +
            papillomaviruses +
            adenoviruses +
            bloodborne_viruses
        )

        logger.info(f"\nTotal genomes collected: {len(all_genomes)}")

        if len(all_genomes) == 0:
            logger.error("No genomes found! Check database content.")
            return

        # Track category sizes for abundance assignment
        categories = {
            'anelloviruses': len(anelloviruses),
            'herpesviruses': len(herpesviruses),
            'polyomaviruses': len(polyomaviruses),
            'other': len(parvoviruses) + len(papillomaviruses) + len(adenoviruses) + len(bloodborne_viruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Delete existing Collection 25 if present
        cursor = self.conn.cursor()
        cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 25")
        cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 25")
        self.conn.commit()
        logger.info("\nDeleted existing Collection 25 (if present)")

        # Insert collection
        collection_meta = {
            'collection_id': 25,
            'collection_name': 'Blood/Plasma Virome (Healthy)',
            'description': f'Healthy human blood/plasma virome with {len(all_genomes)} genomes. Includes anelloviruses (TTV, TLMV - most abundant at 8.91% prevalence), human herpesviruses (HHV-7, EBV, CMV, HHV-6 - 14-20% prevalence), polyomaviruses (BK, JC, Merkel cell - critical for transplant monitoring), parvovirus B19, bloodborne viruses (HCV, HIV), HPV, and adenoviruses. Host: Homo sapiens, Body site: Blood/Plasma. Applications: viremia detection, transplant monitoring, blood safety screening, immunocompromised patient diagnostics. NOTE: Human HBV not in RefSeq database.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Moustafa et al. 2017 (PLOS Pathogens - 8,000 individuals), Young et al. 2016 (polyomavirus transplant monitoring), PMC 2023 (Atlas of blood virome). Anelloviruses most abundant (40% - TTV/TLMV 8.91% prevalence), herpesviruses high (35% - HHV-7 20.37%, EBV 14.45%), polyomaviruses moderate (15% - BK 26.4% in transplant, JC 18% in donors), others low (10% - parvovirus B19, hepatitis, HPV, adenoviruses). Modern ICTV taxonomy used (Orthoherpesviridae).',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-09',
            'literature_references': 'Moustafa et al. 2017 PLOS Pathog (blood DNA virome 8,000 humans); Young et al. 2016 (polyomavirus monitoring); PMC 2023 (atlas of blood virome); ICTV modern taxonomy',
            'version': 1
        }

        # Insert into database
        self._insert_collection(collection_meta, all_genomes)

        logger.info("\n" + "=" * 70)
        logger.info("Collection 25: Blood/Plasma Virome - COMPLETE")
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
    curator = BloodViromeCurator()

    try:
        curator.create_collection()
    finally:
        curator.close()

    logger.info("\nDone!")


if __name__ == '__main__':
    main()
