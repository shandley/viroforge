#!/usr/bin/env python3
"""
Curate Collection 27: Lower Respiratory (Lung) Virome (Healthy)

Literature basis:
- Kitsios et al. 2018 - Respiratory microbiome profiling for pneumonia diagnosis
- Frontiers Immunology 2022 - Role of lung virome during respiratory infections
- Viral metagenomics studies on lung transplant recipients

Composition (healthy lower respiratory tract):
- Anelloviruses: MOST ABUNDANT (70% in most tissues, dominant in lungs)
- Bacteriophages: Prevalent in respiratory secretions
  - Target lung bacteria: Propionibacterium, Streptococcus, Haemophilus, Pseudomonas
- Papillomaviruses: Present in healthy lungs
- Herpesviruses: Variable detection (CMV, EBV important in transplant)
- Bunyaviruses: Present in healthy respiratory tract

Clinical relevance (pneumonia, COPD, transplant):
- RSV, Rhinoviruses, Metapneumovirus, Bocavirus
- Influenza, Coronavirus
- Adenovirus

Target: 25-40 genomes
Applications: Pneumonia diagnosis, COPD monitoring, lung transplant surveillance
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


class LungViromeCurator:
    """Curator for lower respiratory/lung virome collection."""

    def __init__(self):
        """Initialize connection to viral genomes database."""
        db_path = Path(__file__).parent.parent / "viroforge" / "data" / "viral_genomes.db"
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        logger.info(f"Connected to database: {db_path}")

    def get_anelloviruses(self, n_target: int = 5) -> List[Dict]:
        """
        Get anelloviruses (Torque teno viruses).

        MOST ABUNDANT in healthy lungs (70% in most tissues).
        Dominant component of lung virome along with bacteriophages.

        Family: Anelloviridae
        Genome type: Circular ssDNA
        """
        logger.info("Selecting anelloviruses (TTV - DOMINANT 70% in tissues)...")

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

    def get_respiratory_syncytial_virus(self, n_target: int = 2) -> List[Dict]:
        """
        Get respiratory syncytial virus (RSV).

        Leading cause of lower respiratory tract infections in infants.
        Important for pneumonia, bronchiolitis.

        NOTE: Generic "Respiratory syncytial virus" in database is the human virus

        Family: Pneumoviridae
        Genome type: ssRNA(-) linear
        """
        logger.info("Selecting respiratory syncytial virus (RSV)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE 'Respiratory syncytial virus%'
           OR g.genome_name LIKE 'Human respiratory syncytial virus%')
           AND g.genome_name NOT LIKE '%Bovine%'
           AND g.genome_name NOT LIKE '%Potato%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        rsv = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  RSV (Human): {len(rsv)}")

        for r in rsv:
            logger.info(f"    - {r['genome_name'][:60]}")

        return rsv

    def get_rhinoviruses(self, n_target: int = 3) -> List[Dict]:
        """
        Get human rhinoviruses.

        Common cause of upper and lower respiratory infections.
        Important for exacerbations in COPD and asthma.

        Family: Picornaviridae
        Genome type: ssRNA(+) linear
        """
        logger.info("Selecting rhinoviruses (common respiratory infections)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Rhinovirus%'
           AND g.genome_name LIKE '%Human%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        rhinoviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Rhinoviruses: {len(rhinoviruses)}")

        for rv in rhinoviruses:
            logger.info(f"    - {rv['genome_name'][:60]}")

        return rhinoviruses

    def get_influenza_viruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get influenza viruses.

        Major cause of pneumonia and severe respiratory infections.
        Exclude parainfluenza (different family).

        Family: Orthomyxoviridae
        Genome type: ssRNA(-) segmented
        """
        logger.info("Selecting influenza viruses (major pneumonia pathogen)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Influenza A virus%'
           AND g.genome_name NOT LIKE '%parainfluenza%'
           AND g.genome_name NOT LIKE '%Parainfluenza%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        influenza = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Influenza A viruses: {len(influenza)}")

        for iv in influenza:
            logger.info(f"    - {iv['genome_name'][:60]}")

        return influenza

    def get_coronaviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human coronaviruses.

        Include SARS-CoV-2, endemic coronaviruses (OC43, 229E, NL63, HKU1).

        Family: Coronaviridae
        Genome type: ssRNA(+) linear
        """
        logger.info("Selecting coronaviruses (SARS-CoV-2, endemic CoVs)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%SARS-CoV-2%'
           OR g.genome_name LIKE '%Human coronavirus%'
           OR (g.genome_name LIKE '%Coronavirus%' AND g.genome_name LIKE '%Human%'))
        ORDER BY RANDOM()
        LIMIT ?
        """

        coronaviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Coronaviruses: {len(coronaviruses)}")

        for cv in coronaviruses:
            logger.info(f"    - {cv['genome_name'][:60]}")

        return coronaviruses

    def get_metapneumovirus(self, n_target: int = 1) -> List[Dict]:
        """
        Get human metapneumovirus.

        Important respiratory pathogen, especially in immunocompromised.

        Family: Pneumoviridae
        Genome type: ssRNA(-) linear
        """
        logger.info("Selecting metapneumovirus...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Human metapneumovirus%'
           OR g.genome_name LIKE '%Metapneumovirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        metapneumovirus = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Metapneumovirus: {len(metapneumovirus)}")

        for mpv in metapneumovirus:
            logger.info(f"    - {mpv['genome_name'][:60]}")

        return metapneumovirus

    def get_bocavirus(self, n_target: int = 1) -> List[Dict]:
        """
        Get human bocavirus.

        Associated with respiratory infections and exacerbations.

        Family: Parvoviridae
        Genome type: ssDNA linear
        """
        logger.info("Selecting bocavirus (Human)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE 'Human bocavirus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        bocavirus = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Bocavirus (Human): {len(bocavirus)}")

        for bv in bocavirus:
            logger.info(f"    - {bv['genome_name'][:60]}")

        return bocavirus

    def get_adenoviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human adenoviruses.

        Cause respiratory infections, particularly in immunocompromised.

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

    def get_herpesviruses(self, n_target: int = 4) -> List[Dict]:
        """
        Get human herpesviruses.

        CMV, EBV particularly important in lung transplant recipients.
        HSV, VZV can cause pneumonia in immunocompromised.

        Family: Orthoherpesviridae (modern ICTV name)
        Genome type: dsDNA linear
        """
        logger.info("Selecting herpesviruses (CMV, EBV for transplant)...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE (g.genome_name LIKE '%Human herpesvirus 1%'        -- HSV-1
           OR g.genome_name LIKE '%Human herpesvirus 3%'         -- VZV
           OR g.genome_name LIKE '%Human herpesvirus 4%'         -- EBV
           OR g.genome_name LIKE '%Human gammaherpesvirus 4%'    -- EBV alternate
           OR g.genome_name LIKE '%Human herpesvirus 5%'         -- CMV
           OR g.genome_name LIKE '%Human betaherpesvirus 5%')    -- CMV alternate
        ORDER BY RANDOM()
        LIMIT ?
        """

        herpesviruses = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Orthoherpesviridae (CMV, EBV, HSV, VZV): {len(herpesviruses)}")

        for hv in herpesviruses:
            logger.info(f"    - {hv['genome_name'][:60]}")

        return herpesviruses

    def get_papillomaviruses(self, n_target: int = 2) -> List[Dict]:
        """
        Get human papillomaviruses.

        Present in healthy respiratory tract (44% of human viruses).

        Family: Papillomaviridae
        Genome type: dsDNA circular
        """
        logger.info("Selecting papillomaviruses (present in healthy lungs)...")

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

    def get_propionibacterium_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Propionibacterium bacteriophages.

        Elevated in children with multiple acute respiratory infections.
        Propionibacterium is common in respiratory tract.

        Genome type: dsDNA
        """
        logger.info("Selecting Propionibacterium bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Propionibacterium phage%'
           OR g.genome_name LIKE '%Cutibacterium phage%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        prop_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Propionibacterium/Cutibacterium phages: {len(prop_phages)}")

        for pp in prop_phages:
            logger.info(f"    - {pp['genome_name'][:60]}")

        return prop_phages

    def get_streptococcus_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Streptococcus bacteriophages.

        Streptococcus pneumoniae is major bacterial pneumonia pathogen.

        Genome type: dsDNA
        """
        logger.info("Selecting Streptococcus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Streptococcus phage%'
           OR g.genome_name LIKE '%Streptococcus virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        strep_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Streptococcus phages: {len(strep_phages)}")

        for sp in strep_phages:
            logger.info(f"    - {sp['genome_name'][:60]}")

        return strep_phages

    def get_haemophilus_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Haemophilus bacteriophages.

        Haemophilus influenzae causes respiratory infections.

        Genome type: dsDNA
        """
        logger.info("Selecting Haemophilus bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Haemophilus phage%'
           OR g.genome_name LIKE '%Haemophilus virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        haemo_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Haemophilus phages: {len(haemo_phages)}")

        for hp in haemo_phages:
            logger.info(f"    - {hp['genome_name'][:60]}")

        return haemo_phages

    def get_pseudomonas_phages(self, n_target: int = 2) -> List[Dict]:
        """
        Get Pseudomonas bacteriophages.

        Pseudomonas aeruginosa important in hospital-acquired pneumonia, CF.

        Genome type: dsDNA
        """
        logger.info("Selecting Pseudomonas bacteriophages...")

        query = """
        SELECT DISTINCT g.genome_id, g.genome_name, t.family, t.genus, t.species, g.length, g.gc_content
        FROM genomes g
        LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        WHERE g.genome_name LIKE '%Pseudomonas phage%'
           OR g.genome_name LIKE '%Pseudomonas virus%'
        ORDER BY RANDOM()
        LIMIT ?
        """

        pseudo_phages = [dict(row) for row in self.conn.execute(query, (n_target,))]
        logger.info(f"  Pseudomonas phages: {len(pseudo_phages)}")

        for pp in pseudo_phages:
            logger.info(f"    - {pp['genome_name'][:60]}")

        return pseudo_phages

    def assign_abundances(self, genomes: List[Dict], categories: Dict[str, int]) -> List[Dict]:
        """
        Assign relative abundances based on lung virome composition.

        Distribution (literature-based):
        - Anelloviruses: 40% (most abundant, 70% in tissues)
        - Bacteriophages: 35% (prevalent in respiratory secretions)
        - Respiratory viruses: 15% (RSV, rhinovirus, influenza, coronavirus, etc.)
        - Herpesviruses: 5% (variable, important in transplant)
        - Others: 5% (papillomaviruses, bunyaviruses)
        """
        logger.info("\nAssigning relative abundances...")

        # Target abundances
        anello_abundance = 0.40
        phage_abundance = 0.35
        respiratory_abundance = 0.15
        herpes_abundance = 0.05
        other_abundance = 0.05

        for genome in genomes:
            family = genome.get('family') or 'Unknown'
            genome_name = genome.get('genome_name', '')

            # Categorize
            if family == 'Anelloviridae' or 'Torque teno' in genome_name:
                # Anelloviruses: DOMINANT (40%)
                n_anello = categories['anelloviruses']
                base_abundance = anello_abundance / n_anello if n_anello > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.2) if n_anello > 0 else 0

            elif 'phage' in genome_name.lower() or ('virus' in genome_name.lower() and any(bact in genome_name for bact in ['Propionibacterium', 'Cutibacterium', 'Streptococcus', 'Haemophilus', 'Pseudomonas'])):
                # Bacteriophages: High (35%)
                n_phages = categories['bacteriophages']
                base_abundance = phage_abundance / n_phages if n_phages > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_phages > 0 else 0

            elif family in ['Pneumoviridae', 'Picornaviridae', 'Orthomyxoviridae', 'Coronaviridae', 'Parvoviridae', 'Adenoviridae'] and 'herpes' not in genome_name.lower():
                # Respiratory viruses: Moderate (15%)
                n_resp = categories['respiratory_viruses']
                base_abundance = respiratory_abundance / n_resp if n_resp > 0 else 0
                abundance = base_abundance * np.random.lognormal(0, 0.3) if n_resp > 0 else 0

            elif family == 'Orthoherpesviridae' or 'herpesvirus' in genome_name.lower():
                # Herpesviruses: Low (5%)
                # Boost CMV and EBV (important in transplant)
                n_herpes = categories['herpesviruses']
                base_abundance = herpes_abundance / n_herpes if n_herpes > 0 else 0

                if 'herpesvirus 5' in genome_name.lower() or 'CMV' in genome_name:  # CMV
                    abundance = base_abundance * np.random.lognormal(0.3, 0.3)
                elif 'herpesvirus 4' in genome_name.lower() or 'gammaherpesvirus 4' in genome_name.lower():  # EBV
                    abundance = base_abundance * np.random.lognormal(0.2, 0.3)
                else:
                    abundance = base_abundance * np.random.lognormal(0, 0.3)

            else:
                # Others (papillomaviruses, etc.): Very low (5%)
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
        """Create Collection 27: Lower Respiratory (Lung) Virome."""
        logger.info("=" * 70)
        logger.info("Creating Collection 27: Lower Respiratory (Lung) Virome")
        logger.info("=" * 70)

        # Collect genomes by category
        anelloviruses = self.get_anelloviruses(n_target=5)
        rsv = self.get_respiratory_syncytial_virus(n_target=2)
        rhinoviruses = self.get_rhinoviruses(n_target=3)
        influenza = self.get_influenza_viruses(n_target=2)
        coronaviruses = self.get_coronaviruses(n_target=2)
        metapneumovirus = self.get_metapneumovirus(n_target=1)
        bocavirus = self.get_bocavirus(n_target=1)
        adenoviruses = self.get_adenoviruses(n_target=2)
        herpesviruses = self.get_herpesviruses(n_target=4)
        papillomaviruses = self.get_papillomaviruses(n_target=2)
        prop_phages = self.get_propionibacterium_phages(n_target=2)
        strep_phages = self.get_streptococcus_phages(n_target=2)
        haemo_phages = self.get_haemophilus_phages(n_target=2)
        pseudo_phages = self.get_pseudomonas_phages(n_target=2)

        # Combine all genomes
        all_genomes = (
            anelloviruses +
            rsv +
            rhinoviruses +
            influenza +
            coronaviruses +
            metapneumovirus +
            bocavirus +
            adenoviruses +
            herpesviruses +
            papillomaviruses +
            prop_phages +
            strep_phages +
            haemo_phages +
            pseudo_phages
        )

        logger.info(f"\nTotal genomes collected: {len(all_genomes)}")

        if len(all_genomes) == 0:
            logger.error("No genomes found! Check database content.")
            return

        # Track category sizes for abundance assignment
        categories = {
            'anelloviruses': len(anelloviruses),
            'respiratory_viruses': len(rsv) + len(rhinoviruses) + len(influenza) + len(coronaviruses) + len(metapneumovirus) + len(bocavirus) + len(adenoviruses),
            'herpesviruses': len(herpesviruses),
            'bacteriophages': len(prop_phages) + len(strep_phages) + len(haemo_phages) + len(pseudo_phages),
            'other': len(papillomaviruses)
        }

        # Assign abundances
        all_genomes = self.assign_abundances(all_genomes, categories)

        # Delete existing Collection 27 if present
        cursor = self.conn.cursor()
        cursor.execute("DELETE FROM collection_genomes WHERE collection_id = 27")
        cursor.execute("DELETE FROM body_site_collections WHERE collection_id = 27")
        self.conn.commit()
        logger.info("\nDeleted existing Collection 27 (if present)")

        # Insert collection
        collection_meta = {
            'collection_id': 27,
            'collection_name': 'Lower Respiratory (Lung) Virome (Healthy)',
            'description': f'Healthy human lower respiratory tract/lung virome with {len(all_genomes)} genomes. Includes anelloviruses (TTV - DOMINANT 70% in tissues), respiratory viruses (RSV, rhinovirus, influenza, coronavirus, metapneumovirus, bocavirus, adenovirus), herpesviruses (CMV, EBV, HSV, VZV - important in transplant), papillomaviruses, and bacteriophages targeting lung bacteria (Propionibacterium, Streptococcus, Haemophilus, Pseudomonas). Host: Homo sapiens, Body site: Lower respiratory tract/Lungs. Applications: pneumonia diagnosis, COPD monitoring, lung transplant surveillance, respiratory infection research.',
            'n_genomes': len(all_genomes),
            'selection_criteria': 'Literature-validated composition from Kitsios et al. 2018 (respiratory microbiome profiling), Frontiers Immunology 2022 (lung virome role), viral metagenomics studies on lung transplant. Anelloviruses dominant (40% - 70% in tissues), bacteriophages high (35% - Propionibacterium elevated in infections), respiratory viruses moderate (15% - RSV, rhinovirus, influenza, coronavirus, metapneumovirus, bocavirus, adenovirus), herpesviruses low (5% - CMV/EBV important in transplant), others low (5% - papillomaviruses). Modern ICTV taxonomy (Orthoherpesviridae, Pneumoviridae).',
            'curated_by': 'ViroForge Development Team',
            'curation_date': '2025-11-10',
            'literature_references': 'Kitsios et al. 2018 Front Microbiol (respiratory microbiome profiling); Frontiers Immunology 2022 (lung virome role); Viral metagenomics in lung transplant recipients; ICTV modern taxonomy',
            'version': 1
        }

        # Insert into database
        self._insert_collection(collection_meta, all_genomes)

        logger.info("\n" + "=" * 70)
        logger.info("Collection 27: Lower Respiratory (Lung) Virome - COMPLETE")
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
    curator = LungViromeCurator()

    try:
        curator.create_collection()
    finally:
        curator.close()

    logger.info("\nDone!")


if __name__ == '__main__':
    main()
