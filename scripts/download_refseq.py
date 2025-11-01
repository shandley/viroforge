#!/usr/bin/env python3
"""
Download RefSeq viral genomes for ViroForge database.

This script downloads the RefSeq viral genome catalog and FASTA files
for complete viral genomes.

Usage:
    python scripts/download_refseq.py --output data/refseq --limit 100
    python scripts/download_refseq.py --output data/refseq --all

Author: ViroForge Development Team
Date: October 31, 2025
"""

import argparse
import urllib.request
import urllib.error
import gzip
import shutil
import time
import sys
from pathlib import Path
from typing import List, Dict, Tuple
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# RefSeq FTP URLs
REFSEQ_VIRAL_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/"
ASSEMBLY_SUMMARY_URL = f"{REFSEQ_VIRAL_BASE}assembly_summary.txt"


class RefSeqDownloader:
    """Download and manage RefSeq viral genomes."""

    def __init__(self, output_dir: str):
        """
        Initialize downloader.

        Parameters
        ----------
        output_dir : str
            Directory to save downloaded files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.assembly_summary_path = self.output_dir / "assembly_summary.txt"
        self.genomes_dir = self.output_dir / "genomes"
        self.genomes_dir.mkdir(exist_ok=True)

        self.metadata_dir = self.output_dir / "metadata"
        self.metadata_dir.mkdir(exist_ok=True)

    def download_assembly_summary(self, force: bool = False) -> Path:
        """
        Download RefSeq viral assembly summary.

        Parameters
        ----------
        force : bool
            Force re-download even if file exists

        Returns
        -------
        Path
            Path to downloaded assembly summary
        """
        if self.assembly_summary_path.exists() and not force:
            logger.info(f"Assembly summary already exists: {self.assembly_summary_path}")
            return self.assembly_summary_path

        logger.info(f"Downloading assembly summary from: {ASSEMBLY_SUMMARY_URL}")

        try:
            with urllib.request.urlopen(ASSEMBLY_SUMMARY_URL) as response:
                with open(self.assembly_summary_path, 'wb') as f:
                    shutil.copyfileobj(response, f)

            logger.info(f"✓ Downloaded assembly summary: {self.assembly_summary_path}")
            return self.assembly_summary_path

        except urllib.error.URLError as e:
            logger.error(f"Failed to download assembly summary: {e}")
            raise

    def parse_assembly_summary(
        self,
        complete_only: bool = True,
        reference_only: bool = False
    ) -> List[Dict[str, str]]:
        """
        Parse assembly summary file.

        Parameters
        ----------
        complete_only : bool
            Only include complete genomes
        reference_only : bool
            Only include reference/representative genomes

        Returns
        -------
        List[Dict[str, str]]
            List of genome metadata dictionaries
        """
        if not self.assembly_summary_path.exists():
            raise FileNotFoundError(
                f"Assembly summary not found: {self.assembly_summary_path}"
            )

        logger.info("Parsing assembly summary...")

        genomes = []
        with open(self.assembly_summary_path, 'r') as f:
            # Skip comment lines
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 20:
                    continue

                # Parse fields
                genome_info = {
                    'assembly_accession': fields[0],
                    'bioproject': fields[1],
                    'biosample': fields[2],
                    'wgs_master': fields[3],
                    'refseq_category': fields[4],
                    'taxid': fields[5],
                    'species_taxid': fields[6],
                    'organism_name': fields[7],
                    'infraspecific_name': fields[8],
                    'isolate': fields[9],
                    'version_status': fields[10],
                    'assembly_level': fields[11],
                    'release_type': fields[12],
                    'genome_rep': fields[13],
                    'seq_rel_date': fields[14],
                    'asm_name': fields[15],
                    'submitter': fields[16],
                    'gbrs_paired_asm': fields[17],
                    'paired_asm_comp': fields[18],
                    'ftp_path': fields[19] if len(fields) > 19 else '',
                }

                # Apply filters
                if complete_only and genome_info['assembly_level'] != 'Complete Genome':
                    continue

                if reference_only:
                    if genome_info['refseq_category'] not in ['reference genome', 'representative genome']:
                        continue

                # Only include if FTP path exists
                if not genome_info['ftp_path']:
                    continue

                genomes.append(genome_info)

        logger.info(f"✓ Parsed {len(genomes)} genomes matching criteria")
        logger.info(f"  - Complete only: {complete_only}")
        logger.info(f"  - Reference only: {reference_only}")

        return genomes

    def download_genome(
        self,
        genome_info: Dict[str, str],
        retry: int = 3,
        delay: float = 1.0
    ) -> Tuple[bool, Path]:
        """
        Download a single genome FASTA file.

        Parameters
        ----------
        genome_info : Dict[str, str]
            Genome metadata from assembly summary
        retry : int
            Number of retry attempts
        delay : float
            Delay between retries (seconds)

        Returns
        -------
        Tuple[bool, Path]
            (success, path to downloaded file)
        """
        accession = genome_info['assembly_accession']
        ftp_path = genome_info['ftp_path']

        if not ftp_path:
            logger.warning(f"No FTP path for {accession}")
            return False, None

        # Convert FTP to HTTPS
        https_path = ftp_path.replace('ftp://', 'https://')

        # Construct genomic FASTA URL
        filename = Path(ftp_path).name
        fasta_url = f"{https_path}/{filename}_genomic.fna.gz"

        # Output path
        output_path = self.genomes_dir / f"{accession}_genomic.fna.gz"

        # Skip if already downloaded
        if output_path.exists():
            logger.debug(f"Already downloaded: {accession}")
            return True, output_path

        # Attempt download with retries
        for attempt in range(retry):
            try:
                logger.info(f"Downloading {accession} (attempt {attempt + 1}/{retry})...")
                with urllib.request.urlopen(fasta_url, timeout=30) as response:
                    with open(output_path, 'wb') as f:
                        shutil.copyfileobj(response, f)

                logger.info(f"✓ Downloaded: {accession}")
                return True, output_path

            except (urllib.error.URLError, TimeoutError) as e:
                logger.warning(f"Download failed for {accession}: {e}")
                if attempt < retry - 1:
                    time.sleep(delay)
                else:
                    logger.error(f"✗ Failed after {retry} attempts: {accession}")
                    return False, None

    def download_genomes_batch(
        self,
        genomes: List[Dict[str, str]],
        limit: int = None,
        delay: float = 0.1
    ) -> Dict[str, int]:
        """
        Download multiple genomes.

        Parameters
        ----------
        genomes : List[Dict[str, str]]
            List of genome metadata
        limit : int
            Maximum number of genomes to download (None = all)
        delay : float
            Delay between downloads to be polite to NCBI

        Returns
        -------
        Dict[str, int]
            Statistics: total, success, failed, skipped
        """
        if limit:
            genomes = genomes[:limit]

        stats = {
            'total': len(genomes),
            'success': 0,
            'failed': 0,
            'skipped': 0
        }

        logger.info(f"Starting download of {stats['total']} genomes...")

        for i, genome_info in enumerate(genomes, 1):
            accession = genome_info['assembly_accession']

            # Check if already exists
            output_path = self.genomes_dir / f"{accession}_genomic.fna.gz"
            if output_path.exists():
                stats['skipped'] += 1
                if i % 100 == 0:
                    logger.info(f"Progress: {i}/{stats['total']} ({i/stats['total']*100:.1f}%)")
                continue

            # Download
            success, path = self.download_genome(genome_info)

            if success:
                stats['success'] += 1
            else:
                stats['failed'] += 1

            # Progress update
            if i % 10 == 0:
                logger.info(
                    f"Progress: {i}/{stats['total']} ({i/stats['total']*100:.1f}%) | "
                    f"Success: {stats['success']} | Failed: {stats['failed']} | "
                    f"Skipped: {stats['skipped']}"
                )

            # Be polite to NCBI servers
            time.sleep(delay)

        logger.info("=" * 60)
        logger.info("Download complete!")
        logger.info(f"  Total: {stats['total']}")
        logger.info(f"  Success: {stats['success']}")
        logger.info(f"  Failed: {stats['failed']}")
        logger.info(f"  Skipped: {stats['skipped']}")
        logger.info("=" * 60)

        return stats

    def save_genome_list(self, genomes: List[Dict[str, str]], filename: str = "genome_list.tsv"):
        """
        Save genome metadata to TSV file.

        Parameters
        ----------
        genomes : List[Dict[str, str]]
            List of genome metadata
        filename : str
            Output filename
        """
        output_path = self.metadata_dir / filename

        with open(output_path, 'w') as f:
            # Header
            if genomes:
                f.write('\t'.join(genomes[0].keys()) + '\n')

                # Data
                for genome in genomes:
                    f.write('\t'.join(str(v) for v in genome.values()) + '\n')

        logger.info(f"✓ Saved genome list: {output_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Download RefSeq viral genomes for ViroForge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download assembly summary only
  python scripts/download_refseq.py --output data/refseq --summary-only

  # Download first 100 genomes (testing)
  python scripts/download_refseq.py --output data/refseq --limit 100

  # Download all complete genomes
  python scripts/download_refseq.py --output data/refseq --all

  # Download only reference/representative genomes
  python scripts/download_refseq.py --output data/refseq --reference-only --all
        """
    )

    parser.add_argument(
        '--output',
        type=str,
        default='data/refseq',
        help='Output directory (default: data/refseq)'
    )

    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='Maximum number of genomes to download (default: None = all)'
    )

    parser.add_argument(
        '--all',
        action='store_true',
        help='Download all genomes (same as --limit None)'
    )

    parser.add_argument(
        '--summary-only',
        action='store_true',
        help='Only download assembly summary, do not download genomes'
    )

    parser.add_argument(
        '--reference-only',
        action='store_true',
        help='Only include reference/representative genomes'
    )

    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download even if files exist'
    )

    parser.add_argument(
        '--delay',
        type=float,
        default=0.1,
        help='Delay between downloads in seconds (default: 0.1)'
    )

    args = parser.parse_args()

    # Initialize downloader
    downloader = RefSeqDownloader(args.output)

    # Download assembly summary
    logger.info("Step 1: Downloading assembly summary...")
    downloader.download_assembly_summary(force=args.force)

    # Parse assembly summary
    logger.info("Step 2: Parsing assembly summary...")
    genomes = downloader.parse_assembly_summary(
        complete_only=True,
        reference_only=args.reference_only
    )

    # Save genome list
    logger.info("Step 3: Saving genome metadata...")
    downloader.save_genome_list(genomes)

    # Download genomes if requested
    if not args.summary_only:
        logger.info("Step 4: Downloading genome sequences...")

        limit = None if args.all else args.limit
        if limit is None and not args.all:
            limit = 100  # Default to 100 for safety
            logger.warning(f"No limit specified, defaulting to {limit} genomes")
            logger.warning("Use --all to download all genomes")

        stats = downloader.download_genomes_batch(
            genomes,
            limit=limit,
            delay=args.delay
        )

        logger.info("✓ Download complete!")
        logger.info(f"  Genomes downloaded: {stats['success']}")
        logger.info(f"  Output directory: {downloader.genomes_dir}")

    else:
        logger.info("✓ Summary download complete (--summary-only specified)")
        logger.info(f"  {len(genomes)} genomes available for download")
        logger.info(f"  Run without --summary-only to download genomes")


if __name__ == "__main__":
    main()
