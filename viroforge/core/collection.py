"""
ViroForge collection loading utilities.

Provides the CollectionLoader class for loading viral genome collections
from the ViroForge SQLite database, the CustomCollectionBuilder for
creating collections from user-provided FASTA files, and the
label_fastq_headers helper for annotating FASTQ reads with ground-truth
source labels.
"""

import csv
import logging
import random
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

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


class CustomCollectionBuilder:
    """Build a ViroForge-compatible collection from user-provided data.

    Accepts a FASTA file of genomes and optionally a TSV of abundances,
    producing the same (collection_meta, genomes) tuple that
    CollectionLoader.load_collection() returns. This lets users run the
    full ViroForge pipeline (VLP enrichment, contamination, artifacts)
    on their own genomes.

    Usage:
        builder = CustomCollectionBuilder("Penguin Gut Virome")
        meta, genomes = builder.from_fasta("my_genomes.fasta")

        # Or with explicit abundances
        meta, genomes = builder.from_fasta(
            "my_genomes.fasta",
            abundances_tsv="abundances.tsv"
        )
    """

    def __init__(
        self,
        name: str,
        description: str = "",
        host: str = "unknown",
        body_site: str = "unknown",
    ):
        self.name = name
        self.description = description or f"Custom collection: {name}"
        self.host = host
        self.body_site = body_site

    def from_fasta(
        self,
        fasta_path: Path,
        abundances_tsv: Optional[Path] = None,
        abundance_distribution: str = "lognormal",
        random_seed: int = 42,
    ) -> Tuple[Dict, List[Dict]]:
        """Load genomes from a FASTA file and build a collection.

        Args:
            fasta_path: Path to FASTA file containing genome sequences.
            abundances_tsv: Optional TSV with columns: genome_id, abundance.
                If not provided, abundances are generated automatically
                using the specified distribution.
            abundance_distribution: How to assign abundances when no TSV
                is provided. Options: "lognormal" (realistic, default),
                "uniform" (equal abundances).
            random_seed: Seed for reproducible abundance generation.

        Returns:
            Tuple of (collection_meta, genomes) matching the format
            returned by CollectionLoader.load_collection().
        """
        from Bio import SeqIO

        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        # Load genomes from FASTA
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in {fasta_path}")

        logger.info(f"Loaded {len(records)} genomes from {fasta_path.name}")

        # Load or generate abundances
        if abundances_tsv is not None:
            abundance_map = self._load_abundances_tsv(abundances_tsv)
        else:
            abundance_map = self._generate_abundances(
                [r.id for r in records],
                distribution=abundance_distribution,
                random_seed=random_seed,
            )

        # Build genome list
        genomes = []
        for rank, record in enumerate(records, 1):
            seq_str = str(record.seq).upper()
            gc_count = seq_str.count("G") + seq_str.count("C")
            gc_content = gc_count / len(seq_str) if len(seq_str) > 0 else 0.0

            # Infer genome type from sequence content
            genome_type = self._infer_genome_type(seq_str)

            # Parse family from FASTA description if present
            family, genus, species = self._parse_taxonomy_from_description(
                record.description
            )

            abundance = abundance_map.get(record.id, 0.0)

            genomes.append({
                "genome_id": record.id,
                "genome_name": record.description.split()[0]
                    if len(record.description.split()) > 1
                    else record.id,
                "length": len(record.seq),
                "gc_content": round(gc_content, 4),
                "genome_type": genome_type,
                "sequence": seq_str,
                "relative_abundance": abundance,
                "abundance_rank": rank,
                "family": family,
                "genus": genus,
                "species": species,
            })

        # Sort by abundance (descending) and reassign ranks
        genomes.sort(key=lambda g: g["relative_abundance"], reverse=True)
        for i, g in enumerate(genomes, 1):
            g["abundance_rank"] = i

        # Build collection metadata
        # Map body_site to host_organism for contamination defaults
        host_map = {
            "human": "human",
            "mouse": "mouse",
            "rat": "rat",
            "penguin": "none",
            "chicken": "none",
            "fish": "none",
        }
        host_organism = host_map.get(self.host.lower(), "none")

        collection_meta = {
            "collection_id": 0,  # custom, not in database
            "collection_name": self.name,
            "description": self.description,
            "n_genomes": len(genomes),
            "host_organism": host_organism,
            "body_site": self.body_site,
            "is_custom": True,
        }

        logger.info(
            f"Custom collection '{self.name}': {len(genomes)} genomes, "
            f"host={self.host}, body_site={self.body_site}"
        )

        return collection_meta, genomes

    def _load_abundances_tsv(self, tsv_path: Path) -> Dict[str, float]:
        """Load genome abundances from a TSV file.

        Expected format: genome_id<tab>abundance (with optional header).
        Abundances are normalized to sum to 1.0.
        """
        tsv_path = Path(tsv_path)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Abundances file not found: {tsv_path}")

        abundance_map = {}
        with open(tsv_path) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                # Skip header
                if row[0].lower() in ("genome_id", "id", "name"):
                    continue
                if len(row) >= 2:
                    try:
                        abundance_map[row[0]] = float(row[1])
                    except ValueError:
                        continue

        # Normalize
        total = sum(abundance_map.values())
        if total > 0:
            for k in abundance_map:
                abundance_map[k] /= total

        logger.info(f"Loaded abundances for {len(abundance_map)} genomes from {tsv_path.name}")
        return abundance_map

    def _generate_abundances(
        self,
        genome_ids: List[str],
        distribution: str = "lognormal",
        random_seed: int = 42,
    ) -> Dict[str, float]:
        """Generate realistic abundance distributions.

        Args:
            genome_ids: List of genome identifiers.
            distribution: "lognormal" for realistic community structure
                (few dominant, long tail of rare), or "uniform" for
                equal abundances.
            random_seed: Random seed for reproducibility.

        Returns:
            Dict mapping genome_id to normalized abundance.
        """
        rng = np.random.default_rng(random_seed)
        n = len(genome_ids)

        if distribution == "uniform":
            abundances = np.ones(n) / n
        elif distribution == "lognormal":
            # Log-normal distribution mimics real virome communities:
            # a few dominant genomes and a long tail of rare ones
            raw = rng.lognormal(mean=0.0, sigma=2.0, size=n)
            abundances = raw / raw.sum()
        else:
            raise ValueError(
                f"Unknown distribution: {distribution}. "
                f"Choose from: lognormal, uniform"
            )

        return dict(zip(genome_ids, abundances.tolist()))

    def _infer_genome_type(self, sequence: str) -> str:
        """Infer genome type (dsDNA, ssDNA, etc.) from sequence content.

        This is a rough heuristic. Most viral genomes in FASTA format
        are dsDNA. Without additional metadata we default to dsDNA.
        """
        # If very short (<10kb) and no ambiguous bases, could be ssDNA
        if len(sequence) < 10000:
            return "ssDNA"
        return "dsDNA"

    def _parse_taxonomy_from_description(
        self, description: str
    ) -> Tuple[str, str, str]:
        """Try to extract family/genus/species from FASTA description.

        Looks for patterns like:
            >genome_id family:Siphoviridae genus:Lambda species:Lambda phage
            >genome_id Siphoviridae
            >genome_id some virus name

        Returns:
            (family, genus, species) — "Unknown" for missing fields.
        """
        family = "Unknown"
        genus = "Unknown"
        species = "Unknown"

        parts = description.split()
        for part in parts[1:]:  # skip the ID
            lower = part.lower()
            if lower.startswith("family:"):
                family = part.split(":", 1)[1]
            elif lower.startswith("genus:"):
                genus = part.split(":", 1)[1]
            elif lower.startswith("species:"):
                species = part.split(":", 1)[1]

        # If no explicit tags, check if any word ends in -viridae (family)
        if family == "Unknown":
            for part in parts[1:]:
                if part.endswith("viridae") or part.endswith("viridae,"):
                    family = part.rstrip(",")
                    break

        return family, genus, species

    def register_in_database(
        self,
        db_path: Path,
        collection_meta: Dict,
        genomes: List[Dict],
    ) -> int:
        """Register a custom collection in the ViroForge database.

        Inserts genomes, taxonomy, and collection entries so the
        collection can be used with --collection-id like any built-in
        collection.

        Args:
            db_path: Path to ViroForge SQLite database.
            collection_meta: Collection metadata dict.
            genomes: List of genome dicts from from_fasta().

        Returns:
            The new collection_id assigned by the database.
        """
        from datetime import datetime

        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA foreign_keys = ON")

        try:
            # Check if collection name already exists
            existing = conn.execute(
                "SELECT collection_id FROM body_site_collections WHERE collection_name = ?",
                (collection_meta["collection_name"],)
            ).fetchone()

            if existing:
                raise ValueError(
                    f"Collection '{collection_meta['collection_name']}' already exists "
                    f"(ID {existing[0]}). Choose a different name."
                )

            # Insert genomes (skip if genome_id already exists in the database)
            now = datetime.now().isoformat()
            for genome in genomes:
                existing_genome = conn.execute(
                    "SELECT genome_id FROM genomes WHERE genome_id = ?",
                    (genome["genome_id"],)
                ).fetchone()

                if not existing_genome:
                    conn.execute("""
                        INSERT INTO genomes (
                            genome_id, genome_name, sequence, length,
                            gc_content, genome_type, source_database,
                            date_added
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        genome["genome_id"],
                        genome["genome_name"],
                        genome["sequence"],
                        genome["length"],
                        genome["gc_content"],
                        genome["genome_type"],
                        "custom",
                        now,
                    ))

                    # Insert taxonomy
                    conn.execute("""
                        INSERT OR IGNORE INTO taxonomy (
                            genome_id, family, genus, species
                        ) VALUES (?, ?, ?, ?)
                    """, (
                        genome["genome_id"],
                        genome.get("family", "Unknown"),
                        genome.get("genus", "Unknown"),
                        genome.get("species", "Unknown"),
                    ))

            # Create the collection with contamination defaults
            host_organism = collection_meta.get("host_organism", "unknown")

            # Body-site-aware contamination defaults
            # These approximate real sample characteristics
            body_site_defaults = {
                "gut": {"host": 5.0, "rrna": 3.0, "bacterial": 70.0, "fungal": 1.0},
                "oral": {"host": 10.0, "rrna": 4.0, "bacterial": 60.0, "fungal": 0.5},
                "skin": {"host": 15.0, "rrna": 3.0, "bacterial": 50.0, "fungal": 2.0},
                "respiratory": {"host": 20.0, "rrna": 4.0, "bacterial": 40.0, "fungal": 0.5},
                "vaginal": {"host": 20.0, "rrna": 3.0, "bacterial": 50.0, "fungal": 1.0},
                "blood": {"host": 40.0, "rrna": 0.5, "bacterial": 5.0, "fungal": 0.5},
                "urinary": {"host": 10.0, "rrna": 3.0, "bacterial": 30.0, "fungal": 1.0},
                "marine": {"host": 0.05, "rrna": 5.0, "bacterial": 60.0, "fungal": 0.5},
                "soil": {"host": 0.05, "rrna": 4.0, "bacterial": 65.0, "fungal": 3.0},
                "freshwater": {"host": 0.05, "rrna": 5.0, "bacterial": 55.0, "fungal": 1.0},
                "wastewater": {"host": 1.0, "rrna": 8.0, "bacterial": 70.0, "fungal": 1.0},
            }
            body_site = collection_meta.get("body_site", "unknown").lower()
            defaults = body_site_defaults.get(body_site, {
                "host": 5.0, "rrna": 3.0, "bacterial": 50.0, "fungal": 1.0,
            })

            cursor = conn.execute("""
                INSERT INTO body_site_collections (
                    collection_name, description, n_genomes,
                    curated_by, curation_date, host_organism,
                    default_host_pct, default_rrna_pct,
                    default_reagent_pct, default_phix_pct,
                    default_bacterial_pct, default_fungal_pct
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                collection_meta["collection_name"],
                collection_meta.get("description", ""),
                len(genomes),
                "custom (user-provided FASTA)",
                now,
                host_organism,
                defaults["host"],
                defaults["rrna"],
                0.5,  # reagent
                0.1,  # phix
                defaults["bacterial"],
                defaults["fungal"],
            ))

            collection_id = cursor.lastrowid

            # Insert collection-genome associations with abundances
            for genome in genomes:
                conn.execute("""
                    INSERT INTO collection_genomes (
                        collection_id, genome_id,
                        relative_abundance, abundance_rank
                    ) VALUES (?, ?, ?, ?)
                """, (
                    collection_id,
                    genome["genome_id"],
                    genome["relative_abundance"],
                    genome["abundance_rank"],
                ))

            conn.commit()
            logger.info(
                f"Registered collection '{collection_meta['collection_name']}' "
                f"as ID {collection_id} with {len(genomes)} genomes"
            )
            return collection_id

        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
