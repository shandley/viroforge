#!/usr/bin/env python3
"""
Curate Collection 19: HIV+ Gut Virome

Comparison to healthy gut (Collection 9):
- Dramatically reduced diversity (40-60 genomes vs 134)
- Severely reduced crAssphage diversity
- Increased eukaryotic viruses (Anelloviridae, Adenoviridae, Herpesviridae)
- Altered phage composition (fewer gut bacteria phages)
- Highly dysbiotic abundance patterns

Literature basis:
- Handley et al. 2012 (Cell): HIV+ patients show dramatically reduced gut virome diversity
- Monaco et al. 2016 (Cell Host Microbe): Altered virome-immune interactions in HIV
- Vujkovic-Cvijin et al. 2013 (Sci Transl Med): Gut dysbiosis in HIV infection

Target size: 40-60 genomes (55-70% reduction from healthy)

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


class HIVGutCurator:
    """Curate HIV+ gut virome collection from database."""

    def __init__(self, db_path: str = 'viroforge/data/viral_genomes.db'):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.random_seed = 42
        np.random.seed(self.random_seed)

    def get_healthy_gut_genomes(self) -> List[Dict]:
        """Get genomes from healthy gut collection (ID 9) for reference."""
        query = """
        SELECT g.genome_id, g.genome_name, t.family, t.genus, cg.relative_abundance
        FROM collection_genomes cg
        JOIN genomes g ON cg.genome_id = g.genome_id
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE cg.collection_id = 9
        ORDER BY cg.relative_abundance DESC
        """
        return [dict(row) for row in self.conn.execute(query)]

    def get_severely_reduced_crassphage(self, n_target: int = 8) -> List[Dict]:
        """
        Get severely reduced crAssphage diversity (HIV+ characteristic).

        Healthy has ~36 crAssphage genomes
        IBD has ~20 (44% reduction)
        HIV+ has ~8 (78% reduction) - dramatic collapse

        Handley et al. 2012 showed severe reduction in gut bacteriophage diversity.
        """
        logger.info("Selecting severely reduced crAssphage diversity...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Intestiviridae', 'Suoliviridae', 'Steigviridae', 'Crevaviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """

        crassphage = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  crAssphage-like: {len(crassphage)} (vs 36 in healthy, 20 in IBD)")
        return crassphage

    def get_increased_eukaryotic_viruses(self, n_target: int = 20) -> List[Dict]:
        """
        Get increased eukaryotic viruses (HIV+ characteristic).

        HIV+ patients show expansion of eukaryotic viruses due to:
        - Immune dysfunction (CD4+ T-cell depletion)
        - Mucosal barrier damage
        - Viral reactivation (especially herpesviruses, anelloviruses)

        Key families:
        - Anelloviridae (torque teno virus - TTV): Highly expanded in HIV+
        - Adenoviridae: Increased due to immune deficiency
        - Herpesviridae: Reactivation of latent viruses (CMV, EBV, HHV-6)
        - Parvoviridae: Opportunistic infections
        """
        logger.info("Selecting increased eukaryotic viruses (HIV+ characteristic)...")

        # Anelloviridae (TTV) - 40% of eukaryotic (highly expanded)
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Anelloviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        anello = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.4),))]
        logger.info(f"  Anelloviridae (TTV): {len(anello)}")

        # Herpesviridae (reactivation) - 30% of eukaryotic
        # Updated to Orthoherpesviridae after taxonomy fix
        # CRITICAL: Select HUMAN herpesviruses for HIV+ collection
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Herpesviridae', 'Orthoherpesviridae', 'Alloherpesviridae')
          AND (g.genome_name LIKE 'Human herpesvirus%'
           OR g.genome_name LIKE 'Human betaherpesvirus%'
           OR g.genome_name LIKE 'Human gammaherpesvirus%'
           OR g.genome_name LIKE 'Human alphaherpesvirus%')
        ORDER BY RANDOM()
        LIMIT ?
        """
        herpes = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.3),))]
        logger.info(f"  Human Orthoherpesviridae (CMV, EBV, HHV): {len(herpes)}")

        # Adenoviridae - 20% of eukaryotic
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Adenoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        adeno = [dict(row) for row in self.conn.execute(query, (int(n_target * 0.2),))]
        logger.info(f"  Adenoviridae: {len(adeno)}")

        # Parvoviridae - 10% of eukaryotic
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family = 'Parvoviridae'
        ORDER BY RANDOM()
        LIMIT ?
        """
        parvo = [dict(row) for row in self.conn.execute(query, (max(int(n_target * 0.1), 2),))]
        logger.info(f"  Parvoviridae: {len(parvo)}")

        eukaryotic = anello + herpes + adeno + parvo
        logger.info(f"  Total eukaryotic viruses: {len(eukaryotic)}")
        return eukaryotic

    def get_reduced_gut_phages(self, n_target: int = 20) -> List[Dict]:
        """
        Get reduced gut bacteriophages (reflecting gut dysbiosis).

        Reduced diversity of gut bacterial phages due to:
        - Altered gut bacterial composition
        - Loss of commensal bacteria
        - Overgrowth of pathobionts
        """
        logger.info("Selecting reduced gut bacteriophages...")

        # Mix of remaining gut phage families
        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Microviridae', 'Inoviridae', 'Ackermannviridae', 'Drexlerviridae')
           OR (t.genus LIKE '%phage%' AND g.genome_name LIKE '%Lactobacillus%')
           OR (t.genus LIKE '%phage%' AND g.genome_name LIKE '%Enterococcus%')
        ORDER BY RANDOM()
        LIMIT ?
        """

        phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Reduced gut phages: {len(phages)}")
        return phages

    def get_pathogen_associated_viruses(self, n_target: int = 7) -> List[Dict]:
        """
        Get viruses associated with opportunistic infections.

        HIV+ patients are susceptible to various viral pathogens.
        """
        logger.info("Selecting pathogen-associated viruses...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE t.family IN ('Papillomaviridae', 'Polyomaviridae')
        ORDER BY RANDOM()
        LIMIT ?
        """

        pathogens = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pathogen-associated: {len(pathogens)}")
        return pathogens

    def assign_severely_dysbiotic_abundances(self, genomes: List[Dict]) -> List[Dict]:
        """
        Assign abundances reflecting severely dysbiotic state.

        HIV+ gut virome characteristics:
        - Extremely uneven distribution (very high dominance)
        - Few viral blooms dominating the community
        - Lowest diversity of all gut conditions

        Uses most skewed log-normal distribution:
        - Healthy: balanced/tiered
        - IBD: skewed (μ=-2.5, σ=2.5)
        - HIV+: extremely skewed (μ=-3.0, σ=3.0)
        """
        n = len(genomes)

        # Most skewed distribution
        mu = -3.0  # Lowest = most skewed
        sigma = 3.0  # Highest = most variable

        raw_abundances = np.random.lognormal(mu, sigma, n)

        # Normalize to sum to 1.0
        normalized = raw_abundances / raw_abundances.sum()

        # Add to genomes
        for genome, abundance in zip(genomes, normalized):
            genome['relative_abundance'] = float(abundance)

        return genomes

    def create_collection(self) -> List[Dict]:
        """Create complete HIV+ gut virome collection."""
        logger.info("=" * 80)
        logger.info("CURATING HIV+ GUT VIROME COLLECTION")
        logger.info("=" * 80)

        # First, show healthy gut for comparison
        healthy = self.get_healthy_gut_genomes()
        logger.info(f"\nHealthy gut (Collection 9): {len(healthy)} genomes")

        # HIV+ characteristics: dramatically reduced diversity
        crassphage = self.get_severely_reduced_crassphage(n_target=8)
        eukaryotic = self.get_increased_eukaryotic_viruses(n_target=20)
        phages = self.get_reduced_gut_phages(n_target=20)
        pathogens = self.get_pathogen_associated_viruses(n_target=7)

        # Combine all
        collection = crassphage + eukaryotic + phages + pathogens

        # Remove any duplicates (shouldn't be any, but safety check)
        seen_ids = set()
        unique_collection = []
        for genome in collection:
            if genome['genome_id'] not in seen_ids:
                unique_collection.append(genome)
                seen_ids.add(genome['genome_id'])

        collection = unique_collection

        # Assign severely dysbiotic abundances
        logger.info("\nAssigning severely dysbiotic abundance pattern...")
        collection = self.assign_severely_dysbiotic_abundances(collection)

        # Sort by abundance (descending)
        collection = sorted(collection, key=lambda x: x['relative_abundance'], reverse=True)

        # Add abundance ranks
        for i, genome in enumerate(collection, 1):
            genome['abundance_rank'] = i

        logger.info(f"\nTotal genomes in HIV+ collection: {len(collection)}")
        logger.info(f"Diversity reduction: {(1 - len(collection)/len(healthy)) * 100:.1f}% fewer genomes than healthy")
        logger.info(f"Total abundance: {sum(g['relative_abundance'] for g in collection):.6f}")

        # Show top 10
        logger.info("\nTop 10 most abundant genomes:")
        for i, genome in enumerate(collection[:10], 1):
            family = genome.get('family', 'Unknown')
            logger.info(f"  {i:2d}. {genome['genome_name'][:50]:50s} {genome['relative_abundance']:.6f} ({family})")

        # Compare diversity to healthy and IBD
        logger.info("\n" + "=" * 80)
        logger.info("COMPARISON TO HEALTHY GUT AND IBD")
        logger.info("=" * 80)
        logger.info(f"Healthy:  {len(healthy)} genomes (100%)")
        logger.info(f"IBD:      ~90 genomes (~67% of healthy)")
        logger.info(f"HIV+:     {len(collection)} genomes ({len(collection)/len(healthy)*100:.1f}% of healthy)")
        logger.info(f"Reduction: {len(healthy) - len(collection)} genomes lost from healthy state")

        return collection

    def insert_collection(self, collection: List[Dict]):
        """Insert collection into database."""
        logger.info("\n" + "=" * 80)
        logger.info("INSERTING COLLECTION INTO DATABASE")
        logger.info("=" * 80)

        cursor = self.conn.cursor()

        collection_meta = {
            'collection_id': 19,
            'collection_name': 'HIV+ Gut Virome',
            'description': (
                'Gut virome from HIV+ patients showing dramatically reduced viral diversity. '
                'Characterized by severe reduction in bacterial phage diversity, '
                'expansion of eukaryotic viruses (Anelloviridae, Herpesviridae), '
                'and highly dysbiotic abundance patterns. '
                'Compare to Collection 9 (healthy gut) and Collection 18 (IBD) to study '
                'progressive dysbiosis and immune dysfunction effects. '
                'Based on Handley et al. 2012, Monaco et al. 2016, and Vujkovic-Cvijin et al. 2013.'
            ),
            'n_genomes': len(collection),
            'selection_criteria': (
                'Dramatically reduced diversity vs healthy (40-60 genomes vs 134). '
                'Severely reduced crAssphage (8 vs 36). Increased eukaryotic viruses '
                '(Anelloviridae, Herpesviridae, Adenoviridae). Reduced gut phages. '
                'Extremely dysbiotic abundance distribution.'
            ),
            'curated_by': 'ViroForge Phase 7',
            'curation_date': '2025-11-09',
            'literature_references': (
                'Handley et al. 2012 (Cell), Monaco et al. 2016 (Cell Host Microbe), '
                'Vujkovic-Cvijin et al. 2013 (Sci Transl Med), '
                'Nganou-Makamdop et al. 2018 (Cell Host Microbe)'
            ),
            'version': 1
        }

        # Check if collection exists
        cursor.execute("SELECT collection_id FROM body_site_collections WHERE collection_id = 19")
        exists = cursor.fetchone()

        if exists:
            logger.info("Collection 19 already exists - DELETING and recreating...")
            cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 19")
            cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 19")

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
                19,
                genome['genome_id'],
                genome['relative_abundance'],
                1.0,  # All genomes present
                genome['abundance_rank']
            ))

        self.conn.commit()
        logger.info(f"✓ Inserted {len(collection)} genome associations")

        logger.info("\n✓ Collection 19 successfully created in database!")

    def close(self):
        """Close database connection."""
        self.conn.close()


def main():
    """Main curation workflow."""
    curator = HIVGutCurator()

    try:
        # Create collection
        collection = curator.create_collection()

        # Insert into database
        curator.insert_collection(collection)

        logger.info("\n" + "=" * 80)
        logger.info("✓ HIV+ GUT VIROME COLLECTION CURATION COMPLETE!")
        logger.info("=" * 80)
        logger.info("\nNext steps:")
        logger.info("  1. Test generation: python scripts/generate_fastq_dataset.py --collection-id 19 --dry-run")
        logger.info("  2. Compare to healthy and IBD: Generate all three collections and compare")
        logger.info("  3. Update documentation")

    finally:
        curator.close()


if __name__ == '__main__':
    main()
