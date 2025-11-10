#!/usr/bin/env python3
"""
Curate Collection 28: Urinary Virome (Healthy)

Literature basis:
- Santiago-Rodriguez et al. 2015 (Frontiers Microbiol) - Human urine virome in UTI
- Frontiers Medicine 2018 - Urinary virome perturbations in kidney transplant
- Scientific Reports 2016 - Diverse virome in kidney transplant patients

Composition (healthy urine and transplant monitoring):
- Papillomaviruses (HPV): 95% prevalence in healthy subjects (DOMINANT)
- Bacteriophages: Most identifiable viruses in urine
- BK virus: 60-80% in transplant complications, critical biomarker
- JC virus: Important for transplant monitoring
- Anelloviruses: 53% in BKV-negative transplant patients
- Adenoviruses: Present in urine
- Herpesviruses: Present in urine

Transplant applications:
- BK virus prevalence: 60-70% acute rejection, 70-80% chronic nephropathy, 100% BK nephritis
- Virome diversity higher in transplant patients (20-32 unique viruses)

Target: 15-25 genomes
Applications: Transplant monitoring (BK/JC viremia), UTI diagnosis, hemorrhagic cystitis
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


class UrinaryViromeCurator:
    """Curator for urinary virome collection."""

    def __init__(self):
        """Initialize connection to viral genomes database."""
        db_path = Path(__file__).parent.parent / "viroforge" / "data" / "viral_genomes.db"
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        logger.info(f"Connected to database: {db_path}")

    def get_papillomaviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Get human papillomaviruses.

        DOMINANT in healthy urine (95% prevalence).
        Both low-risk genotypes and novel strains detected.

        Family: Papillomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting papillomaviruses (HPV - DOMINANT 95% prevalence)...")

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

    def get_bk_virus(self, n_target: int = 1) -> List[Dict]:
        """
        Get BK polyomavirus.

        CRITICAL for transplant monitoring.
        Prevalence: 60-80% in transplant complications, 100% in BK nephritis.

        Family: Polyomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting BK polyomavirus (CRITICAL transplant biomarker)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'BK polyomavirus%'
           OR g.genome_name LIKE 'BK virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        bk_virus = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  BK polyomavirus: {len(bk_virus)}")

        for bv in bk_virus:
            logger.info(f"    - {bv['genome_name'][:60]}")

        return bk_virus

    def get_jc_virus(self, n_target: int = 1) -> List[Dict]:
        """
        Get JC polyomavirus.

        Important for transplant monitoring.
        Shares 75% nucleotide sequence identity with BK virus.

        Family: Polyomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting JC polyomavirus (transplant monitoring)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'JC polyomavirus%'
           OR g.genome_name LIKE 'JC virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        jc_virus = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  JC polyomavirus: {len(jc_virus)}")

        for jv in jc_virus:
            logger.info(f"    - {jv['genome_name'][:60]}")

        return jc_virus

    def get_anelloviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get anelloviruses (Torque teno viruses).

        53% prevalence in BKV-negative transplant patients.
        Lower in patients with BK viremia.

        Family: Anelloviridae
        Genome type: Circular ssDNA
        """
        logger.info("Selecting anelloviruses (TTV - 53% in transplant patients)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Torque teno virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        anelloviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Anelloviridae (TTV): {len(anelloviruses)}")

        for av in anelloviruses:
            logger.info(f"    - {av['genome_name'][:60]}")

        return anelloviruses

    def get_adenoviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human adenoviruses.

        Present in urine, pathogenic in immunocompromised.
        Can cause hemorrhagic cystitis in transplant patients.

        Family: Adenoviridae
        Genome type: dsDNA linear
        """
        logger.info("Selecting adenoviruses (hemorrhagic cystitis)...")

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

    def get_herpesviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human herpesviruses.

        Present in urinary virome.
        CMV particularly relevant for transplant patients.

        Family: Orthoherpesviridae (modern ICTV name)
        Genome type: dsDNA linear
        """
        logger.info("Selecting herpesviruses (CMV for transplant)...")

        # Get CMV (HHV-5) - CRITICAL for transplant
        query_cmv = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Human herpesvirus 5%'
           OR g.genome_name LIKE '%Human betaherpesvirus 5%')
        ORDER BY RANDOM()
        LIMIT 1
        """

        # Get EBV (HHV-4)
        query_ebv = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Human herpesvirus 4%'
           OR g.genome_name LIKE '%Human gammaherpesvirus 4%')
        ORDER BY RANDOM()
        LIMIT 1
        """

        herpesviruses = []
        herpesviruses.extend([dict(row) for row in self.conn.execute(query_cmv)])
        herpesviruses.extend([dict(row) for row in self.conn.execute(query_ebv)])

        logger.info(f"  Orthoherpesviridae (CMV, EBV): {len(herpesviruses)}")

        for hv in herpesviruses:
            logger.info(f"    - {hv['genome_name'][:60]}")

        return herpesviruses

    def get_escherichia_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Escherichia (E. coli) bacteriophages.

        E. coli is major UTI pathogen (80-90% of UTIs).
        Bacteriophages are most identifiable viruses in urine.

        Genome type: dsDNA
        """
        logger.info("Selecting Escherichia (E. coli) bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Escherichia phage%'
           OR g.genome_name LIKE '%Escherichia virus%'
           OR g.genome_name LIKE '%Enterobacteria phage%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        escherichia_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Escherichia/Enterobacteria phages: {len(escherichia_phages)}")

        for ep in escherichia_phages:
            logger.info(f"    - {ep['genome_name'][:60]}")

        return escherichia_phages

    def get_enterococcus_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Enterococcus bacteriophages.

        Enterococcus is common UTI pathogen (10-20% of UTIs).

        Genome type: dsDNA
        """
        logger.info("Selecting Enterococcus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Enterococcus phage%'
           OR g.genome_name LIKE '%Enterococcus virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        enterococcus_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Enterococcus phages: {len(enterococcus_phages)}")

        for ep in enterococcus_phages:
            logger.info(f"    - {ep['genome_name'][:60]}")

        return enterococcus_phages

    def get_staphylococcus_phages(self, n_target: int = 1) -> List[Dict]:
        """
        Get Staphylococcus bacteriophages.

        Staphylococcus saprophyticus causes UTIs (5-10%).

        Genome type: dsDNA
        """
        logger.info("Selecting Staphylococcus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Staphylococcus phage%'
           OR g.genome_name LIKE '%Staphylococcus virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        staphylococcus_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Staphylococcus phages: {len(staphylococcus_phages)}")

        for sp in staphylococcus_phages:
            logger.info(f"    - {sp['genome_name'][:60]}")

        return staphylococcus_phages

    def get_lactobacillus_phages(self, n_target: int = 1) -> List[Dict]:
        """
        Get Lactobacillus bacteriophages.

        Lactobacillus is commensal in urinary tract.

        Genome type: dsDNA
        """
        logger.info("Selecting Lactobacillus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Lactobacillus phage%'
           OR g.genome_name LIKE '%Lactobacillus virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        lactobacillus_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Lactobacillus phages: {len(lactobacillus_phages)}")

        for lp in lactobacillus_phages:
            logger.info(f"    - {lp['genome_name'][:60]}")

        return lactobacillus_phages

    def assign_abundances(self, genomes: List[Dict], categories: Dict[str, int]) -> List[Dict]:
        """
        Assign relative abundances based on urinary virome composition.

        Distribution (literature-based):
        - Papillomaviruses (HPV): 45% (DOMINANT, 95% prevalence)
        - Bacteriophages: 30% (most identifiable viruses)
        - BK/JC polyomaviruses: 15% (critical transplant biomarkers)
        - Anelloviruses: 5% (53% in transplant patients)
        - Others: 5% (adenoviruses, herpesviruses)
        """
        logger.info("\nAssigning relative abundances...")

        # Target abundances
        hpv_abundance = 0.45
        phage_abundance = 0.30
        polyoma_abundance = 0.15
        anello_abundance = 0.05
        other_abundance = 0.05

        for genome in genomes:
            family = genome.get('family') or 'Unknown'
            genome_name = genome.get('genome_name', '')

            # Categorize
            if family == 'Papillomaviridae':
                # HPV: DOMINANT (45%)
                n_hpv = categories['papillomaviruses']
                base_abundance = hpv_abundance / n_hpv if n_hpv > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.2) if n_hpv > 0 else 0

            elif 'phage' in genome_name.lower() or ('virus' in genome_name.lower() and any(bact in genome_name for bact in ['Escherichia', 'Enterobacteria', 'Enterococcus', 'Staphylococcus', 'Lactobacillus'])):
                # Bacteriophages: High (30%)
                n_phages = categories['bacteriophages']
                base_abundance = phage_abundance / n_phages if n_phages > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_phages > 0 else 0

            elif 'BK' in genome_name or 'JC' in genome_name:
                # BK/JC polyomaviruses: Moderate (15%)
                # Boost BK (most important for transplant)
                n_polyoma = categories['polyomaviruses']
                base_abundance = polyoma_abundance / n_polyoma if n_polyoma > 0 else 0

                if 'BK' in genome_name:
                    abundance = base_abundance * np.random.lognormal(0.3, 0.3)  # Higher BK
                else:
                    abundance = base_abundance * np.random.lognormal(0, 0.3)

            elif family == 'Anelloviridae' or 'Torque teno' in genome_name:
                # Anelloviruses: Low (5%)
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello if n_anello > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_anello > 0 else 0

            else:
                # Others (adenoviruses, herpesviruses): Low (5%)
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
        for genome in genomes[:15]:  # Top 15
            logger.info(f"  {genome['abundance_rank']:2d}. {genome['genome_name'][:50]:50s} {genome['relative_abundance']:7.1%}")

        if len(genomes) > 15:
            logger.info(f"  ... and {len(genomes) - 15} more genomes")

        return genomes

    def create_collection(self):
        """Create Collection 28: Urinary Virome."""
        logger.info("=" * 70)
        logger.info("Creating Collection 28: Urinary Virome")
        logger.info("=" * 70)

        # Collect genomes by category
        papillomaviruses = self.get_papillomaviruses(n_target=5)
        bk_virus = self.get_bk_virus(n_target=1)
        jc_virus = self.get_jc_virus(n_target=1)
        anelloviruses = self.get_anelloviruses(n_target=3)
        adenoviruses = self.get_adenoviruses(n_target=2)
        herpesviruses = self.get_herpesviruses(n_target=2)
        escherichia_phages = self.get_escherichia_phages(n_target=2)
        enterococcus_phages = self.get_enterococcus_phages(n_target=2)
        staphylococcus_phages = self.get_staphylococcus_phages(n_target=1)
        lactobacillus_phages = self.get_lactobacillus_phages(n_target=1)

        # Combine all genomes
        all_genomes = (
            papillomaviruses +
            bk_virus +
            jc_virus +
            anelloviruses +
            adenoviruses +
            herpesviruses +
            escherichia_phages +
            enterococcus_phages +
            staphylococcus_phages +
            lactobacillus_phages
        )

        logger.info(f"\nTotal genomes collected: {len(all_genomes)}")

        if len(all_genomes) == 0:
            logger.error("No genomes found! Check database content.")
            return

        # Track category sizes for abundance assignment
        categories = {
            'papillomaviruses': len(papillomaviruses),
            'polyomaviruses': len(bk_virus) + len(jc_virus),
            'anelloviruses': len(anelloviruses),
            'bacteriophages': len(escherichia_phages) + len(enterococcus_phages) + len(staphylococcus_phages) + len(lactobacillus_phages),
            'other': len(adenoviruses) + len(herpesviruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Delete existing Collection 28 if present
        cursor = self.conn.cursor()
        cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 28")
        cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 28")
        self.conn.commit()
        logger.info("\nDeleted existing Collection 28 (if present)")

        # Insert collection
        collection_meta = {
            'collection_id': 28,
            'collection_name': 'Urinary Virome (Healthy)',
            'description': f'Healthy human urinary tract virome with {len(all_genomes)} genomes. Includes papillomaviruses (HPV - DOMINANT 95% prevalence), BK polyomavirus (CRITICAL transplant biomarker 60-100% in complications), JC polyomavirus, anelloviruses (TTV - 53% in transplant patients), adenoviruses (hemorrhagic cystitis), herpesviruses (CMV, EBV), and bacteriophages targeting urinary bacteria (E. coli, Enterococcus, Staphylococcus, Lactobacillus). Host: Homo sapiens, Body site: Urinary tract/Urine. Applications: kidney transplant monitoring (BK/JC viremia), UTI diagnosis, hemorrhagic cystitis detection, immunocompromised patient monitoring.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Santiago-Rodriguez et al. 2015 (Frontiers Microbiol - 20 subjects), Frontiers Medicine 2018 (urinary virome perturbations in kidney transplant), Scientific Reports 2016 (diverse virome in transplant patients). HPV dominant (45% - 95% prevalence), bacteriophages high (30% - most identifiable), BK/JC polyomaviruses moderate (15% - 60-100% in transplant complications), anelloviruses low (5% - 53% in transplant patients), others low (5% - adenoviruses, herpesviruses). Modern ICTV taxonomy (Orthoherpesviridae).',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-10',
            'literature_references': 'Santiago-Rodriguez et al. 2015 Front Microbiol (human urine virome in UTI); Frontiers Medicine 2018 (urinary virome in kidney transplant); Scientific Reports 2016 (diverse virome in transplant patients); ICTV modern taxonomy',
            'version': 1
        }

        # Insert into database
        self._insert_collection(collection_meta, all_genomes)

        logger.info("\n" + "=" * 70)
        logger.info("Collection 28: Urinary Virome - COMPLETE")
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
    curator = UrinaryViromeCurator()

    try:
        curator.create_collection()
    finally:
        curator.close()

    logger.info("\nDone!")


if __name__ == '__main__':
    main()
