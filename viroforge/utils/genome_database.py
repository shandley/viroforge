"""
Viral genome database management utilities.

This module handles downloading, caching, and validating viral genomes
from NCBI RefSeq using Biopython Entrez.

Features:
- Download genomes from NCBI on first use
- Cache downloaded genomes locally
- Validate genome quality (completeness, taxonomy, sequence quality)
- Track metadata for all genomes
- Support for minimal and full genome sets

Usage:
    from viroforge.utils import get_genome_database

    # Download and get genome database (caches automatically)
    genomes = get_genome_database(dataset='minimal')  # or 'full'

Author: ViroForge Development Team
Date: 2025-01-31
"""

import logging
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
from datetime import datetime

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..data.curated_genomes import MINIMAL_TEST_SET, FULL_PRODUCTION_SET
from ..core.community import ViralGenome

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# IMPORTANT: Set your email for NCBI Entrez
# Users should set this in their environment or code
Entrez.email = "viroforge@example.com"  # Default, should be overridden

# Get data directory path
DATA_DIR = Path(__file__).parent.parent / 'data'
GENOMES_DIR = DATA_DIR / 'genomes'
METADATA_DIR = DATA_DIR / 'metadata'

# Ensure directories exist
GENOMES_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)


class GenomeDownloadError(Exception):
    """Raised when genome download fails."""
    pass


class GenomeValidationError(Exception):
    """Raised when genome validation fails."""
    pass


def set_entrez_email(email: str):
    """
    Set email address for NCBI Entrez requests (required by NCBI).

    Args:
        email: Your email address

    Example:
        >>> set_entrez_email('your.email@institution.edu')
    """
    Entrez.email = email
    logger.info(f"Set Entrez email to: {email}")


def download_genome_from_ncbi(
    accession: str,
    retries: int = 3,
    delay: float = 0.5
) -> SeqRecord:
    """
    Download a single genome from NCBI using Biopython Entrez.

    Args:
        accession: NCBI accession number (e.g., 'NC_001422.1')
        retries: Number of retry attempts if download fails
        delay: Delay between retries (seconds)

    Returns:
        Bio.SeqRecord object

    Raises:
        GenomeDownloadError: If download fails after all retries

    Example:
        >>> record = download_genome_from_ncbi('NC_001422.1')
        >>> print(record.id, len(record.seq))
    """
    for attempt in range(retries):
        try:
            logger.info(f"Downloading {accession} (attempt {attempt + 1}/{retries})")

            # Fetch genome from NCBI nucleotide database
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="gb",  # GenBank format (includes annotations)
                retmode="text"
            )

            # Parse GenBank record
            record = SeqIO.read(handle, "genbank")
            handle.close()

            logger.info(f"✓ Downloaded {accession}: {len(record.seq)} bp")
            return record

        except Exception as e:
            logger.warning(f"✗ Download attempt {attempt + 1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise GenomeDownloadError(
                    f"Failed to download {accession} after {retries} attempts: {e}"
                )


def validate_genome_record(
    record: SeqRecord,
    min_length: int = 1000,
    max_length: int = 1_000_000,
    max_n_percent: float = 5.0,
    min_gc: float = 20.0,
    max_gc: float = 80.0
) -> Tuple[bool, List[str]]:
    """
    Validate genome quality.

    Quality checks:
    - Valid DNA characters only (ATCGN)
    - Length in biological range
    - Not excessive N bases
    - GC content in biological range
    - Has taxonomy annotation

    Args:
        record: BioPython SeqRecord
        min_length: Minimum genome length (bp)
        max_length: Maximum genome length (bp)
        max_n_percent: Maximum percent N bases allowed
        min_gc: Minimum GC content (%)
        max_gc: Maximum GC content (%)

    Returns:
        Tuple of (is_valid: bool, errors: List[str])

    Example:
        >>> is_valid, errors = validate_genome_record(record)
        >>> if not is_valid:
        ...     print("Validation errors:", errors)
    """
    errors = []

    # Check sequence length
    seq_len = len(record.seq)
    if seq_len < min_length:
        errors.append(f"Sequence too short: {seq_len} bp < {min_length} bp")
    if seq_len > max_length:
        errors.append(f"Sequence too long: {seq_len} bp > {max_length} bp")

    # Check for valid DNA characters
    seq_str = str(record.seq).upper()
    invalid_chars = set(seq_str) - set('ATCGN')
    if invalid_chars:
        errors.append(f"Invalid DNA characters found: {invalid_chars}")

    # Check N content
    n_count = seq_str.count('N')
    n_percent = (n_count / seq_len) * 100 if seq_len > 0 else 0
    if n_percent > max_n_percent:
        errors.append(f"Excessive N bases: {n_percent:.1f}% > {max_n_percent}%")

    # Check GC content
    g_count = seq_str.count('G')
    c_count = seq_str.count('C')
    gc_content = ((g_count + c_count) / seq_len) * 100 if seq_len > 0 else 0
    if gc_content < min_gc or gc_content > max_gc:
        errors.append(
            f"GC content out of range: {gc_content:.1f}% "
            f"(expected {min_gc}-{max_gc}%)"
        )

    # Check for taxonomy annotation
    if not record.annotations.get('organism'):
        errors.append("Missing organism annotation")

    is_valid = len(errors) == 0
    return is_valid, errors


def extract_genome_metadata(record: SeqRecord) -> Dict:
    """
    Extract metadata from GenBank record.

    Extracts:
    - Accession
    - Organism/taxonomy
    - Description
    - Length
    - GC content
    - Annotations

    Args:
        record: BioPython SeqRecord from GenBank

    Returns:
        Dictionary of metadata

    Example:
        >>> metadata = extract_genome_metadata(record)
        >>> print(metadata['organism'], metadata['length'])
    """
    seq_str = str(record.seq).upper()
    seq_len = len(seq_str)

    # Calculate GC content
    g_count = seq_str.count('G')
    c_count = seq_str.count('C')
    gc_content = ((g_count + c_count) / seq_len) * 100 if seq_len > 0 else 0.0

    # Extract taxonomy
    organism = record.annotations.get('organism', 'Unknown')
    taxonomy = record.annotations.get('taxonomy', [])

    # Try to extract family from taxonomy
    family = 'Unknown'
    for taxon in taxonomy:
        if taxon.endswith('viridae') or taxon.endswith('virales'):
            family = taxon
            break

    metadata = {
        'accession': record.id,
        'organism': organism,
        'description': record.description,
        'length': seq_len,
        'gc_content': round(gc_content, 2),
        'family': family,
        'taxonomy': ';'.join(taxonomy) if taxonomy else 'Unknown',
        'source': record.annotations.get('source', 'Unknown'),
        'download_date': datetime.now().isoformat(),
        'ncbi_taxid': record.features[0].qualifiers.get('db_xref', [''])[0] if record.features else '',
    }

    return metadata


def download_genome_set(
    genome_set: Dict[str, List[str]],
    output_fasta: Path,
    output_metadata: Path,
    force_download: bool = False
) -> Tuple[int, int]:
    """
    Download a complete set of genomes from NCBI.

    Args:
        genome_set: Dictionary of family -> list of accessions
        output_fasta: Path to output FASTA file
        output_metadata: Path to output metadata TSV file
        force_download: Re-download even if files exist

    Returns:
        Tuple of (successful_downloads, failed_downloads)

    Example:
        >>> from viroforge.data.curated_genomes import MINIMAL_TEST_SET
        >>> success, failed = download_genome_set(
        ...     MINIMAL_TEST_SET,
        ...     Path('genomes.fasta'),
        ...     Path('metadata.tsv')
        ... )
    """
    # Check if already downloaded
    if output_fasta.exists() and output_metadata.exists() and not force_download:
        logger.info(f"✓ Genome set already downloaded: {output_fasta}")
        logger.info("  Use force_download=True to re-download")
        return 0, 0

    logger.info("="*60)
    logger.info("ViroForge Genome Database Download")
    logger.info("="*60)

    # Flatten genome set
    all_accessions = []
    accession_to_family = {}
    for family, accessions in genome_set.items():
        all_accessions.extend(accessions)
        for acc in accessions:
            accession_to_family[acc] = family

    total = len(all_accessions)
    logger.info(f"Total genomes to download: {total}")
    logger.info(f"Families: {len(genome_set)}")
    logger.info("")

    # Download genomes
    downloaded_records = []
    metadata_list = []
    success_count = 0
    failed_count = 0

    for i, accession in enumerate(all_accessions, 1):
        family = accession_to_family[accession]
        logger.info(f"[{i}/{total}] {family}: {accession}")

        try:
            # Download
            record = download_genome_from_ncbi(accession, retries=3, delay=0.5)

            # Validate
            is_valid, errors = validate_genome_record(record)
            if not is_valid:
                logger.warning(f"  ⚠ Validation warnings: {'; '.join(errors)}")
                # Continue anyway - these are warnings, not hard failures

            # Extract metadata
            metadata = extract_genome_metadata(record)
            metadata['curated_family'] = family  # Add our curated family annotation

            # Store
            downloaded_records.append(record)
            metadata_list.append(metadata)
            success_count += 1

            logger.info(f"  ✓ {metadata['organism']} ({metadata['length']} bp, "
                       f"{metadata['gc_content']}% GC)")

            # Be nice to NCBI - rate limit
            time.sleep(0.35)  # ~3 requests/second (NCBI allows 3/sec without API key)

        except Exception as e:
            logger.error(f"  ✗ Failed to download {accession}: {e}")
            failed_count += 1
            continue

    # Write FASTA file
    if downloaded_records:
        logger.info("")
        logger.info(f"Writing FASTA to: {output_fasta}")
        SeqIO.write(downloaded_records, output_fasta, "fasta")
        logger.info(f"✓ Wrote {len(downloaded_records)} genomes")

    # Write metadata file
    if metadata_list:
        logger.info(f"Writing metadata to: {output_metadata}")
        df = pd.DataFrame(metadata_list)
        df.to_csv(output_metadata, sep='\t', index=False)
        logger.info(f"✓ Wrote metadata for {len(metadata_list)} genomes")

    # Summary
    logger.info("")
    logger.info("="*60)
    logger.info("Download Summary")
    logger.info("="*60)
    logger.info(f"Successful: {success_count}/{total}")
    logger.info(f"Failed: {failed_count}/{total}")
    logger.info(f"Success rate: {(success_count/total)*100:.1f}%")
    logger.info("="*60)

    return success_count, failed_count


def get_genome_database(
    dataset: str = 'minimal',
    force_download: bool = False,
    email: Optional[str] = None
) -> List[ViralGenome]:
    """
    Get viral genome database (downloads on first use, then caches).

    This is the main entry point for users to get genomes.

    Args:
        dataset: Which dataset to use ('minimal' or 'full')
        force_download: Re-download even if cached
        email: Email for NCBI Entrez (required)

    Returns:
        List of ViralGenome objects

    Raises:
        ValueError: If dataset is invalid
        RuntimeError: If download fails

    Example:
        >>> from viroforge.utils import get_genome_database
        >>> genomes = get_genome_database('minimal', email='you@example.com')
        >>> print(f"Loaded {len(genomes)} genomes")
    """
    # Set email if provided
    if email:
        set_entrez_email(email)
    elif Entrez.email == "viroforge@example.com":
        logger.warning("⚠ Using default Entrez email. Set your email with:")
        logger.warning("  from viroforge.utils import set_entrez_email")
        logger.warning("  set_entrez_email('your.email@institution.edu')")

    # Select genome set
    if dataset == 'minimal':
        genome_set = MINIMAL_TEST_SET
        fasta_name = 'minimal_test_set.fasta'
        metadata_name = 'minimal_test_set_metadata.tsv'
    elif dataset == 'full':
        genome_set = FULL_PRODUCTION_SET
        fasta_name = 'full_production_set.fasta'
        metadata_name = 'full_production_set_metadata.tsv'
    else:
        raise ValueError(f"Invalid dataset: {dataset}. Choose 'minimal' or 'full'")

    # File paths
    fasta_path = GENOMES_DIR / fasta_name
    metadata_path = METADATA_DIR / metadata_name

    # Download if not cached
    if not fasta_path.exists() or not metadata_path.exists() or force_download:
        logger.info(f"Downloading {dataset} genome set...")
        success, failed = download_genome_set(
            genome_set,
            fasta_path,
            metadata_path,
            force_download=force_download
        )

        if failed > 0:
            logger.warning(f"⚠ {failed} genomes failed to download")

        if success == 0:
            raise RuntimeError(f"Failed to download any genomes for {dataset} set")

    # Load genomes from cached files
    logger.info(f"Loading genomes from cache: {fasta_path}")

    # Read FASTA
    records = list(SeqIO.parse(fasta_path, "fasta"))

    # Read metadata
    metadata_df = pd.read_csv(metadata_path, sep='\t')
    metadata_dict = metadata_df.set_index('accession').to_dict('index')

    # Convert to ViralGenome objects
    viral_genomes = []
    for record in records:
        accession = record.id
        meta = metadata_dict.get(accession, {})

        genome = ViralGenome(
            genome_id=accession,
            sequence=record.seq,
            taxonomy=meta.get('taxonomy', 'Unknown'),
            family=meta.get('curated_family', meta.get('family', 'Unknown')),
            genus='Unknown',  # Can extract from full taxonomy if needed
            species=meta.get('organism', 'Unknown'),
            description=meta.get('description', '')
        )

        viral_genomes.append(genome)

    logger.info(f"✓ Loaded {len(viral_genomes)} viral genomes")

    return viral_genomes


def get_database_info(dataset: str = 'minimal') -> Dict:
    """
    Get information about a genome database without downloading.

    Args:
        dataset: Which dataset ('minimal' or 'full')

    Returns:
        Dictionary with database information

    Example:
        >>> info = get_database_info('full')
        >>> print(f"Database has {info['total_genomes']} genomes")
    """
    if dataset == 'minimal':
        genome_set = MINIMAL_TEST_SET
    elif dataset == 'full':
        genome_set = FULL_PRODUCTION_SET
    else:
        raise ValueError(f"Invalid dataset: {dataset}")

    total = sum(len(accs) for accs in genome_set.values())
    family_counts = {fam: len(accs) for fam, accs in genome_set.items()}

    return {
        'dataset': dataset,
        'total_genomes': total,
        'families': len(genome_set),
        'family_counts': family_counts,
        'accessions': genome_set
    }


if __name__ == '__main__':
    # Example usage
    print("ViroForge Genome Database Utility")
    print("=" * 60)

    # Show database info
    print("\nMinimal Test Set:")
    info = get_database_info('minimal')
    print(f"  Total genomes: {info['total_genomes']}")
    print(f"  Families: {info['families']}")

    print("\nFull Production Set:")
    info = get_database_info('full')
    print(f"  Total genomes: {info['total_genomes']}")
    print(f"  Families: {info['families']}")
    for family, count in info['family_counts'].items():
        print(f"    {family}: {count} genomes")

    print("\n" + "=" * 60)
    print("To download genomes, use:")
    print("  from viroforge.utils import get_genome_database")
    print("  genomes = get_genome_database('minimal', email='your@email.com')")
