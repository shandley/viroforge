#!/usr/bin/env python3
"""
Curate Collection 26: Ocular Surface Virome (Healthy)

Literature basis:
- Doan et al. 2016 (IOVS) - Paucibacterial microbiome and resident DNA virome of healthy conjunctiva
- Clinical ophthalmology literature on viral keratitis

Composition (healthy conjunctiva):
- Torque teno virus (TTV): 65% prevalence (DOMINANT)
- Merkel cell polyomavirus: Present in some individuals
- HPV: Low prevalence, not widespread
- Viral reads: 5% average, 1% median (vs 91% bacterial)
- Paucibacterial environment

Clinically important (infectious keratitis):
- HSV-1: Leading cause of infectious blindness in developed world
- Adenovirus: Most common viral ocular infection worldwide
- VZV, CMV, EBV: Herpesvirus keratitis

Bacteriophages (inferred from bacterial composition):
- Targets: Staphylococcus, Propionibacterium, Corynebacterium (dominant bacteria)

Target: 10-20 genomes
Applications: Ophthalmology, infectious keratitis diagnosis, ocular health
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


class OcularViromeCurator:
    """Curator for ocular surface virome collection."""

    def __init__(self):
        """Initialize connection to viral genomes database."""
        db_path = Path(__file__).parent.parent / "viroforge" / "data" / "viral_genomes.db"
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        logger.info(f"Connected to database: {db_path}")

    def get_anelloviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get anelloviruses (Torque teno viruses).

        DOMINANT virus on ocular surface (65% prevalence).
        TTV found in two-thirds of healthy conjunctiva samples.

        Family: Anelloviridae
        Genome type: Circular ssDNA
        """
        logger.info("Selecting anelloviruses (TTV - dominant at 65% prevalence)...")

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
        logger.info(f"  Anelloviridae (TTV): {len(anelloviruses)}")

        for av in anelloviruses:
            logger.info(f"    - {av['genome_name'][:60]}")

        return anelloviruses

    def get_adenoviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get human adenoviruses.

        Most common viral ocular infection worldwide.
        Causes epidemic keratoconjunctivitis.

        Family: Adenoviridae
        Genome type: dsDNA linear
        """
        logger.info("Selecting adenoviruses (most common viral ocular infection)...")

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

    def get_herpesviruses(self, n_target: int = 4) -> List[Dict]:
        """
        Get human herpesviruses causing ocular infections.

        HSV-1 (HHV-1): Leading cause of infectious blindness in developed world
        VZV (HHV-3): Herpes zoster ophthalmicus
        CMV (HHV-5): CMV retinitis (immunocompromised)
        EBV (HHV-4): Ocular manifestations

        Family: Orthoherpesviridae (modern ICTV name)
        Genome type: dsDNA linear
        """
        logger.info("Selecting herpesviruses (HSV-1, VZV, CMV, EBV for keratitis)...")

        # Specifically select clinically relevant herpesviruses
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Human herpesvirus 1%'        -- HSV-1 (most important!)
           OR g.genome_name LIKE '%Human alphaherpesvirus 1%'    -- HSV-1 alternate name
           OR g.genome_name LIKE '%Human herpesvirus 3%'         -- VZV
           OR g.genome_name LIKE '%Human alphaherpesvirus 3%'    -- VZV alternate name
           OR g.genome_name LIKE '%Human herpesvirus 4%'         -- EBV
           OR g.genome_name LIKE '%Human gammaherpesvirus 4%'    -- EBV alternate name
           OR g.genome_name LIKE '%Human herpesvirus 5%'         -- CMV
           OR g.genome_name LIKE '%Human betaherpesvirus 5%')    -- CMV alternate name
        ORDER BY RANDOM()
        LIMIT ?
        """

        herpesviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Orthoherpesviridae (Human - HSV-1, VZV, EBV, CMV): {len(herpesviruses)}")

        for hv in herpesviruses:
            logger.info(f"    - {hv['genome_name'][:60]}")

        return herpesviruses

    def get_polyomaviruses(self, n_target: int = 1) -> List[Dict]:
        """
        Get Merkel cell polyomavirus.

        Present in some healthy conjunctiva samples (Doan 2016).

        Family: Polyomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting polyomaviruses (Merkel cell polyomavirus)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Polyomaviridae'
           AND g.genome_name LIKE '%Merkel cell polyomavirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        polyomaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Polyomaviridae (Merkel cell): {len(polyomaviruses)}")

        for pv in polyomaviruses:
            logger.info(f"    - {pv['genome_name'][:60]}")

        return polyomaviruses

    def get_papillomaviruses(self, n_target: int = 1) -> List[Dict]:
        """
        Get human papillomaviruses.

        Low prevalence, not widespread in healthy conjunctiva (Doan 2016).

        Family: Papillomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting papillomaviruses (HPV - low prevalence)...")

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

    def get_staphylococcus_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Staphylococcus bacteriophages.

        Staphylococcus is dominant bacterium on ocular surface.
        Bacteriophages inferred from bacterial composition.

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

        staph_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Staphylococcus phages: {len(staph_phages)}")

        for sp in staph_phages:
            logger.info(f"    - {sp['genome_name'][:60]}")

        return staph_phages

    def get_propionibacterium_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Propionibacterium/Cutibacterium bacteriophages.

        Propionibacterium (now Cutibacterium) is common on ocular surface.
        Note: Propionibacterium acnes reclassified to Cutibacterium acnes.

        Genome type: dsDNA
        """
        logger.info("Selecting Propionibacterium/Cutibacterium bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Propionibacterium phage%'
           OR g.genome_name LIKE '%Propionibacterium virus%'
           OR g.genome_name LIKE '%Cutibacterium phage%'
           OR g.genome_name LIKE '%Cutibacterium virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        prop_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Propionibacterium/Cutibacterium phages: {len(prop_phages)}")

        for pp in prop_phages:
            logger.info(f"    - {pp['genome_name'][:60]}")

        return prop_phages

    def get_corynebacterium_phages(self, n_target: int = 1) -> List[Dict]:
        """
        Get Corynebacterium bacteriophages.

        Corynebacterium is common on ocular surface (80% of samples).

        Genome type: dsDNA
        """
        logger.info("Selecting Corynebacterium bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Corynebacterium phage%'
           OR g.genome_name LIKE '%Corynebacterium virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        coryne_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Corynebacterium phages: {len(coryne_phages)}")

        for cp in coryne_phages:
            logger.info(f"    - {cp['genome_name'][:60]}")

        return coryne_phages

    def assign_abundances(self, genomes: List[Dict], categories: Dict[str, int]) -> List[Dict]:
        """
        Assign relative abundances based on ocular surface virome composition.

        Distribution (literature-based):
        - Anelloviruses (TTV): 60% (DOMINANT, 65% prevalence)
        - Bacteriophages: 25% (inferred from bacterial dominance 91%)
        - Adenoviruses: 8% (most common infection)
        - Herpesviruses: 5% (HSV-1 keratitis)
        - Others: 2% (polyomaviruses, HPV - low prevalence)

        Note: Paucibacterial environment, low overall viral abundance
        """
        logger.info("\nAssigning relative abundances...")

        # Target abundances
        anello_abundance = 0.60
        phage_abundance = 0.25
        adeno_abundance = 0.08
        herpes_abundance = 0.05
        other_abundance = 0.02

        for genome in genomes:
            family = genome.get('family') or 'Unknown'
            genome_name = genome.get('genome_name', '')

            # Categorize by family/name
            if family == 'Anelloviridae' or 'Torque teno' in genome_name:
                # TTV: DOMINANT (60%)
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello if n_anello > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.2) if n_anello > 0 else 0

            elif 'phage' in genome_name.lower() or 'virus' in genome_name.lower():
                # Bacteriophages: Moderate (25%)
                n_phages = categories['bacteriophages']
                base_abundance = phage_abundance / n_phages if n_phages > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_phages > 0 else 0

            elif family == 'Adenoviridae':
                # Adenoviruses: Moderate-low (8%)
                n_adeno = categories['adenoviruses']
                base_abundance = adeno_abundance / n_adeno if n_adeno > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_adeno > 0 else 0

            elif family == 'Orthoherpesviridae' or 'herpesvirus' in genome_name.lower():
                # Herpesviruses: Low (5%)
                # Boost HSV-1 (most important for keratitis)
                n_herpes = categories['herpesviruses']
                base_abundance = herpes_abundance / n_herpes if n_herpes > 0 else 0

                if 'herpesvirus 1' in genome_name.lower() or 'HSV-1' in genome_name:
                    abundance = base_abundance * np.random.lognormal(0.3, 0.3)  # Higher HSV-1
                else:
                    abundance = base_abundance * np.random.lognormal(0, 0.3)

            else:
                # Others (polyomaviruses, HPV): Very low (2%)
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
        """Create Collection 26: Ocular Surface Virome."""
        logger.info("=" * 70)
        logger.info("Creating Collection 26: Ocular Surface Virome")
        logger.info("=" * 70)

        # Collect genomes by category
        anelloviruses = self.get_anelloviruses(n_target=3)
        adenoviruses = self.get_adenoviruses(n_target=3)
        herpesviruses = self.get_herpesviruses(n_target=4)
        polyomaviruses = self.get_polyomaviruses(n_target=1)
        papillomaviruses = self.get_papillomaviruses(n_target=1)
        staph_phages = self.get_staphylococcus_phages(n_target=2)
        prop_phages = self.get_propionibacterium_phages(n_target=2)
        coryne_phages = self.get_corynebacterium_phages(n_target=1)

        # Combine all genomes
        all_genomes = (
            anelloviruses +
            adenoviruses +
            herpesviruses +
            polyomaviruses +
            papillomaviruses +
            staph_phages +
            prop_phages +
            coryne_phages
        )

        logger.info(f"\nTotal genomes collected: {len(all_genomes)}")

        if len(all_genomes) == 0:
            logger.error("No genomes found! Check database content.")
            return

        # Track category sizes for abundance assignment
        categories = {
            'anelloviruses': len(anelloviruses),
            'adenoviruses': len(adenoviruses),
            'herpesviruses': len(herpesviruses),
            'bacteriophages': len(staph_phages) + len(prop_phages) + len(coryne_phages),
            'other': len(polyomaviruses) + len(papillomaviruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Delete existing Collection 26 if present
        cursor = self.conn.cursor()
        cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 26")
        cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 26")
        self.conn.commit()
        logger.info("\nDeleted existing Collection 26 (if present)")

        # Insert collection
        collection_meta = {
            'collection_id': 26,
            'collection_name': 'Ocular Surface Virome (Healthy)',
            'description': f'Healthy human ocular surface/conjunctiva virome with {len(all_genomes)} genomes. Paucibacterial environment with resident DNA virome. Includes torque teno virus (TTV - DOMINANT at 65% prevalence), adenoviruses (most common viral ocular infection), herpesviruses (HSV-1 leading cause of infectious blindness, VZV, CMV, EBV), Merkel cell polyomavirus, HPV, and bacteriophages targeting Staphylococcus, Propionibacterium/Cutibacterium, Corynebacterium. Host: Homo sapiens, Body site: Conjunctiva/Ocular surface. Applications: ophthalmology, infectious keratitis diagnosis, ocular health monitoring.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Doan et al. 2016 (IOVS - 107 healthy volunteers, TTV 65% prevalence), clinical ophthalmology literature on viral keratitis. TTV dominant (60% - most abundant), bacteriophages moderate (25% - inferred from 91% bacterial reads), adenoviruses moderate-low (8% - most common infection), herpesviruses low (5% - HSV-1 keratitis), others very low (2% - polyomaviruses, HPV). Paucibacterial environment: 91% bacterial, 5% viral, minimal fungal. Modern ICTV taxonomy (Orthoherpesviridae).',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-09',
            'literature_references': 'Doan et al. 2016 IOVS (paucibacterial microbiome and resident DNA virome); Clinical ophthalmology literature on viral keratitis (HSV-1, adenovirus); ICTV modern taxonomy',
            'version': 1
        }

        # Insert into database
        self._insert_collection(collection_meta, all_genomes)

        logger.info("\n" + "=" * 70)
        logger.info("Collection 26: Ocular Surface Virome - COMPLETE")
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
    curator = OcularViromeCurator()

    try:
        curator.create_collection()
    finally:
        curator.close()

    logger.info("\nDone!")


if __name__ == '__main__':
    main()
