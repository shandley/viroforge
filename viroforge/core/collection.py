"""
ViroForge collection loading utilities.

Provides the CollectionLoader class for loading viral genome collections
from the ViroForge SQLite database, and the label_fastq_headers helper
for annotating FASTQ reads with ground-truth source labels.
"""

import logging
import random
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def label_fastq_headers(
    fastq_path: Path,
    genome_source_map: Dict[str, str],
    output_path: Optional[Path] = None,
) -> Path:
    """Add source type labels to FASTQ read headers.

    ISS encodes the source genome in read headers as @{genome_id}_{index}_{pair}/{read}.
    This function appends source=viral/host_dna/rrna/phix/reagent to each header
    so downstream tools can compute exact classification metrics.

    Args:
        fastq_path: Path to input FASTQ file.
        genome_source_map: Dict mapping genome_id prefixes to source types.
        output_path: Output path (default: overwrite in place).

    Returns:
        Path to labeled FASTQ file.
    """
    if output_path is None:
        output_path = fastq_path

    labeled_lines = []
    with open(fastq_path) as f:
        for line_num, line in enumerate(f):
            if line_num % 4 == 0 and line.startswith("@"):
                # Parse genome ID from ISS header: @{genome_id}_{index}_{pair}/{read}
                header = line.rstrip()
                read_id = header[1:]  # strip @

                # Find matching genome ID (ISS uses genome_id as prefix)
                source = "unknown"
                for genome_id, src_type in genome_source_map.items():
                    if read_id.startswith(genome_id):
                        source = src_type
                        break

                labeled_lines.append(f"{header} source={source}\n")
            else:
                labeled_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(labeled_lines)

    return output_path


class CollectionLoader:
    """Load genomes from curated body site collections."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)

        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

    def list_collections(self) -> List[Dict]:
        """List available collections."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        cursor = conn.execute("""
            SELECT collection_id, collection_name, n_genomes, description
            FROM body_site_collections
            ORDER BY collection_name
        """)

        collections = [dict(row) for row in cursor.fetchall()]
        conn.close()

        return collections

    def load_collection(self, collection_id: int) -> Tuple[Dict, List[Dict]]:
        """
        Load collection metadata and genomes.

        Returns:
            Tuple of (collection_metadata, genomes_with_abundances)
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Get collection metadata
        cursor = conn.execute("""
            SELECT *
            FROM body_site_collections
            WHERE collection_id = ?
        """, (collection_id,))

        collection = cursor.fetchone()
        if not collection:
            conn.close()
            raise ValueError(f"Collection {collection_id} not found")

        collection_meta = dict(collection)

        # Get genomes with abundances
        cursor = conn.execute("""
            SELECT
                cg.genome_id,
                cg.relative_abundance,
                cg.abundance_rank,
                g.genome_name,
                g.length,
                g.gc_content,
                g.genome_type,
                g.sequence,
                t.family,
                t.genus,
                t.species
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE cg.collection_id = ?
            ORDER BY cg.relative_abundance DESC
        """, (collection_id,))

        genomes = [dict(row) for row in cursor.fetchall()]
        conn.close()

        logger.info(f"Loaded collection '{collection_meta['collection_name']}' "
                   f"with {len(genomes)} genomes")

        return collection_meta, genomes

    def load_dark_matter_genomes(
        self,
        n_genomes: int,
        exclude_ids: set,
        random_seed: int = 42
    ) -> List[Dict]:
        """
        Load unclassified viral genomes (dark matter) from the database.

        Selects genomes with family='Unknown' (real viral sequences that lack
        taxonomic classification, mirroring the large unclassifiable fraction of
        real virome datasets), excluding genomes already in the collection and
        genomes that are not true dark matter:
        - Known human viruses that belong in the classified pool (HHV, HIV, HPV,
          norovirus, influenza, hepatitis, etc.)
        - Animal viruses (bovine, porcine, avian, murine, bat, etc.)
        - Insect viruses (baculoviruses, nucleopolyhedroviruses)

        Plant viruses are intentionally kept: dietary plant viruses (PMMoV, ToMV,
        TMV) are a real and abundant component of human gut viromes.

        Args:
            n_genomes: Number of dark matter genomes to select
            exclude_ids: Set of genome_ids already in the collection
            random_seed: Seed for reproducible selection

        Returns:
            List of genome dictionaries (same format as load_collection)
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Genomes in the Unknown family that are not true dark matter. Plant
        # viruses are deliberately not excluded (dietary viruses are real).
        host_exclusions = """
            AND g.genome_name NOT LIKE 'Human herpesvirus%'
            AND g.genome_name NOT LIKE 'Human alpha%herpes%'
            AND g.genome_name NOT LIKE 'Human beta%herpes%'
            AND g.genome_name NOT LIKE 'Human gamma%herpes%'
            AND g.genome_name NOT LIKE 'Human immunodeficiency%'
            AND g.genome_name NOT LIKE 'Human papillomavirus%'
            AND g.genome_name NOT LIKE 'Human bocavirus%'
            AND g.genome_name NOT LIKE 'Human enterovirus%'
            AND g.genome_name NOT LIKE 'Human rhinovirus%'
            AND g.genome_name NOT LIKE 'Human adenovirus%'
            AND g.genome_name NOT LIKE 'Human astrovirus%'
            AND g.genome_name NOT LIKE 'Human rotavirus%'
            AND g.genome_name NOT LIKE 'Hepatitis%'
            AND g.genome_name NOT LIKE 'Norovirus%'
            AND g.genome_name NOT LIKE 'Influenza%'
            AND g.genome_name NOT LIKE '%Gallid%'
            AND g.genome_name NOT LIKE '%Bovine%'
            AND g.genome_name NOT LIKE '%Porcine%'
            AND g.genome_name NOT LIKE '%Canine%'
            AND g.genome_name NOT LIKE '%Feline%'
            AND g.genome_name NOT LIKE '%Murine%'
            AND g.genome_name NOT LIKE '%Equine%'
            AND g.genome_name NOT LIKE '%Avian%'
            AND g.genome_name NOT LIKE '%Simian%'
            AND g.genome_name NOT LIKE '%Bat %'
            AND g.genome_name NOT LIKE '%Alcelaphine%'
            AND g.genome_name NOT LIKE '%Helicoverpa%'
            AND g.genome_name NOT LIKE '%Acyrthosiphon%'
            AND g.genome_name NOT LIKE '%nucleopolyhedrovirus%'
            AND g.genome_name NOT LIKE '%granulovirus%'
            AND g.genome_name NOT LIKE '%baculovirus%'
            AND g.genome_name NOT LIKE '%Drosophila%'
        """

        if exclude_ids:
            placeholders = ','.join('?' for _ in exclude_ids)
            exclude_clause = f"AND g.genome_id NOT IN ({placeholders})"
            params = list(exclude_ids)
        else:
            exclude_clause = ""
            params = []

        # Fetch eligible genome IDs in a deterministic order, then sample with a
        # seeded RNG (SQLite has no seedable RANDOM()). Only the selected rows'
        # sequences are fetched, to avoid loading thousands of sequences.
        id_rows = conn.execute(f"""
            SELECT g.genome_id
            FROM genomes g
            JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE t.family = 'Unknown'
            AND g.sequence IS NOT NULL
            AND length(g.sequence) > 500
            {host_exclusions}
            {exclude_clause}
            ORDER BY g.genome_id
        """, params).fetchall()
        eligible_ids = [r['genome_id'] for r in id_rows]

        if not eligible_ids:
            logger.warning("No dark matter genomes available in database")
            conn.close()
            return []

        n_to_select = min(n_genomes, len(eligible_ids))
        selected_ids = random.Random(random_seed).sample(eligible_ids, n_to_select)

        sel_placeholders = ','.join('?' for _ in selected_ids)
        cursor = conn.execute(f"""
            SELECT g.genome_id, g.genome_name, g.length, g.gc_content,
                   g.genome_type, g.sequence, t.family, t.genus, t.species
            FROM genomes g
            JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE g.genome_id IN ({sel_placeholders})
        """, selected_ids)

        dark_matter = []
        for row in cursor.fetchall():
            genome = dict(row)
            genome['relative_abundance'] = 0.0  # Will be set by caller
            genome['abundance_rank'] = 0
            genome['is_dark_matter'] = True
            dark_matter.append(genome)

        conn.close()

        logger.info(f"Loaded {len(dark_matter)} dark matter genomes "
                   f"(from {len(eligible_ids)} eligible; excluded animal/insect/"
                   f"known-human, kept plant/dietary)")

        return dark_matter
