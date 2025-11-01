#!/usr/bin/env python3
"""
Parse RefSeq viral genomes and prepare for database insertion.

This script reads downloaded RefSeq genome FASTA files, extracts sequences,
calculates genome statistics (length, GC content), and prepares data for
insertion into the ViroForge genome database.

Usage:
    python scripts/parse_genomes.py --input data/refseq --output data/parsed
    python scripts/parse_genomes.py --input data/refseq --output data/parsed --limit 100

Author: ViroForge Development Team
Date: November 1, 2025
"""

import argparse
import gzip
import json
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
from datetime import datetime
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GenomeParser:
    """Parse viral genome FASTA files and extract metadata."""

    def __init__(self, input_dir: str, output_dir: str):
        """
        Initialize genome parser.

        Parameters
        ----------
        input_dir : str
            Directory containing downloaded RefSeq data
        output_dir : str
            Directory to save parsed genome data
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.genomes_dir = self.input_dir / "genomes"
        self.metadata_dir = self.input_dir / "metadata"
        self.metadata_file = self.metadata_dir / "genome_list.tsv"

        # Output files
        self.parsed_dir = self.output_dir / "genomes"
        self.parsed_dir.mkdir(exist_ok=True)

        self.stats_file = self.output_dir / "genome_stats.tsv"
        self.failed_file = self.output_dir / "failed_genomes.txt"

    def parse_fasta_header(self, header: str) -> Dict[str, str]:
        """
        Parse FASTA header to extract genome information.

        Parameters
        ----------
        header : str
            FASTA header line (without '>')

        Returns
        -------
        Dict[str, str]
            Parsed header information

        Examples
        --------
        >>> parser.parse_fasta_header("NC_003663.2 Cowpox virus, complete genome")
        {'accession': 'NC_003663.2', 'description': 'Cowpox virus, complete genome'}
        """
        parts = header.split(None, 1)
        if len(parts) == 2:
            return {
                'accession': parts[0],
                'description': parts[1]
            }
        else:
            return {
                'accession': parts[0] if parts else '',
                'description': ''
            }

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content of a sequence.

        Parameters
        ----------
        sequence : str
            DNA/RNA sequence

        Returns
        -------
        float
            GC content (0.0 to 1.0)
        """
        sequence_upper = sequence.upper()
        g_count = sequence_upper.count('G')
        c_count = sequence_upper.count('C')
        total_count = len(sequence_upper)

        if total_count == 0:
            return 0.0

        return (g_count + c_count) / total_count

    def detect_genome_type(self, header_desc: str, sequence: str) -> str:
        """
        Detect genome type (DNA/RNA, ss/ds) from header and sequence.

        Parameters
        ----------
        header_desc : str
            FASTA header description
        sequence : str
            Genome sequence

        Returns
        -------
        str
            Genome type: 'dsDNA', 'ssDNA', 'dsRNA', 'ssRNA', or 'unknown'

        Notes
        -----
        Detection strategy:
        1. Check header for explicit mentions of RNA/DNA and ss/ds
        2. Check for 'U' in sequence (RNA marker)
        3. Default to dsDNA (most common for RefSeq viral)
        """
        header_lower = header_desc.lower()

        # Check for RNA
        if 'rna' in header_lower or 'U' in sequence.upper():
            # Check for single/double stranded
            if 'ssrna' in header_lower or 'single-stranded rna' in header_lower:
                return 'ssRNA'
            elif 'dsrna' in header_lower or 'double-stranded rna' in header_lower:
                return 'dsRNA'
            else:
                return 'ssRNA'  # Default to ssRNA for RNA viruses

        # Check for DNA
        elif 'dna' in header_lower or 'U' not in sequence.upper():
            # Check for single/double stranded
            if 'ssdna' in header_lower or 'single-stranded dna' in header_lower:
                return 'ssDNA'
            elif 'dsdna' in header_lower or 'double-stranded dna' in header_lower:
                return 'dsDNA'
            else:
                return 'dsDNA'  # Default to dsDNA (most common)

        return 'dsDNA'  # Default

    def parse_genome_file(self, genome_file: Path, genome_id: str) -> Optional[Dict[str, any]]:
        """
        Parse a single genome FASTA file.

        Parameters
        ----------
        genome_file : Path
            Path to gzipped FASTA file
        genome_id : str
            Genome RefSeq accession

        Returns
        -------
        Optional[Dict[str, any]]
            Parsed genome data, or None if parsing failed
        """
        try:
            with gzip.open(genome_file, 'rt') as f:
                lines = f.readlines()

            if not lines:
                logger.warning(f"Empty file: {genome_id}")
                return None

            # Parse header
            header_line = lines[0].strip()
            if not header_line.startswith('>'):
                logger.warning(f"Invalid FASTA format (no header): {genome_id}")
                return None

            header_info = self.parse_fasta_header(header_line[1:])

            # Parse sequence
            sequence_lines = [line.strip() for line in lines[1:] if not line.startswith('>')]
            sequence = ''.join(sequence_lines)

            if not sequence:
                logger.warning(f"Empty sequence: {genome_id}")
                return None

            # Calculate statistics
            length = len(sequence)
            gc_content = self.calculate_gc_content(sequence)
            genome_type = self.detect_genome_type(header_info['description'], sequence)

            # Check for multi-segment genomes
            n_segments = len([line for line in lines if line.startswith('>')])

            genome_data = {
                'genome_id': genome_id,
                'genome_name': header_info['description'],
                'sequence': sequence,
                'length': length,
                'gc_content': gc_content,
                'genome_type': genome_type,
                'genome_structure': 'linear',  # RefSeq viruses are typically linear
                'n_segments': n_segments,
                'assembly_level': 'Complete Genome',
                'source_database': 'RefSeq',
                'genbank_accession': header_info['accession'],
                'date_added': datetime.now().isoformat(),
                'version': 1
            }

            return genome_data

        except Exception as e:
            logger.error(f"Error parsing {genome_id}: {e}")
            return None

    def load_metadata(self) -> Dict[str, Dict[str, str]]:
        """
        Load genome metadata from TSV file.

        Returns
        -------
        Dict[str, Dict[str, str]]
            Metadata keyed by assembly accession
        """
        metadata = {}

        if not self.metadata_file.exists():
            logger.warning(f"Metadata file not found: {self.metadata_file}")
            return metadata

        with open(self.metadata_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                accession = row['assembly_accession']
                metadata[accession] = row

        logger.info(f"Loaded metadata for {len(metadata)} genomes")
        return metadata

    def parse_all_genomes(self, limit: Optional[int] = None) -> Tuple[List[Dict], List[str]]:
        """
        Parse all downloaded genome files.

        Parameters
        ----------
        limit : Optional[int]
            Maximum number of genomes to parse (None = all)

        Returns
        -------
        Tuple[List[Dict], List[str]]
            (successfully parsed genomes, list of failed genome IDs)
        """
        logger.info("Starting genome parsing...")

        # Get list of genome files
        genome_files = sorted(self.genomes_dir.glob("*_genomic.fna.gz"))

        if limit:
            genome_files = genome_files[:limit]

        logger.info(f"Found {len(genome_files)} genome files to parse")

        # Load metadata
        metadata = self.load_metadata()

        # Parse genomes
        parsed_genomes = []
        failed_genomes = []

        for i, genome_file in enumerate(genome_files, 1):
            # Extract genome ID from filename
            genome_id = genome_file.name.replace('_genomic.fna.gz', '')

            # Parse genome
            genome_data = self.parse_genome_file(genome_file, genome_id)

            if genome_data:
                # Add metadata if available
                if genome_id in metadata:
                    meta = metadata[genome_id]
                    genome_data['organism_name'] = meta.get('organism_name', '')
                    genome_data['taxid'] = meta.get('taxid', '')
                    genome_data['species_taxid'] = meta.get('species_taxid', '')
                    genome_data['refseq_category'] = meta.get('refseq_category', '')

                parsed_genomes.append(genome_data)
            else:
                failed_genomes.append(genome_id)

            # Progress update
            if i % 100 == 0:
                logger.info(
                    f"Progress: {i}/{len(genome_files)} ({i/len(genome_files)*100:.1f}%) | "
                    f"Success: {len(parsed_genomes)} | Failed: {len(failed_genomes)}"
                )

        logger.info("=" * 60)
        logger.info("Parsing complete!")
        logger.info(f"  Total files: {len(genome_files)}")
        logger.info(f"  Successfully parsed: {len(parsed_genomes)}")
        logger.info(f"  Failed: {len(failed_genomes)}")
        logger.info("=" * 60)

        return parsed_genomes, failed_genomes

    def save_parsed_genomes(self, genomes: List[Dict]) -> None:
        """
        Save parsed genome data to JSON file.

        Parameters
        ----------
        genomes : List[Dict]
            List of parsed genome data dictionaries
        """
        output_file = self.output_dir / "parsed_genomes.json"

        # Don't save full sequences in JSON (too large)
        # Save only metadata
        metadata_only = []
        for genome in genomes:
            meta = genome.copy()
            # Save sequence to separate file
            seq_file = self.parsed_dir / f"{genome['genome_id']}.fasta"
            with open(seq_file, 'w') as f:
                f.write(f">{genome['genome_id']} {genome['genome_name']}\n")
                # Write sequence in 80-character lines
                sequence = genome['sequence']
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

            # Remove sequence from metadata JSON
            meta.pop('sequence', None)
            metadata_only.append(meta)

        with open(output_file, 'w') as f:
            json.dump(metadata_only, f, indent=2)

        logger.info(f"✓ Saved parsed genome metadata: {output_file}")
        logger.info(f"✓ Saved {len(genomes)} genome sequences to: {self.parsed_dir}")

    def save_statistics(self, genomes: List[Dict]) -> None:
        """
        Save genome statistics to TSV file.

        Parameters
        ----------
        genomes : List[Dict]
            List of parsed genome data dictionaries
        """
        with open(self.stats_file, 'w', newline='') as f:
            fieldnames = [
                'genome_id', 'genome_name', 'organism_name', 'length',
                'gc_content', 'genome_type', 'n_segments', 'taxid'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')

            writer.writeheader()
            for genome in genomes:
                row = {field: genome.get(field, '') for field in fieldnames}
                writer.writerow(row)

        logger.info(f"✓ Saved genome statistics: {self.stats_file}")

    def save_failed_list(self, failed_genomes: List[str]) -> None:
        """
        Save list of failed genome IDs.

        Parameters
        ----------
        failed_genomes : List[str]
            List of genome IDs that failed parsing
        """
        if not failed_genomes:
            return

        with open(self.failed_file, 'w') as f:
            for genome_id in failed_genomes:
                f.write(f"{genome_id}\n")

        logger.info(f"✓ Saved failed genome list: {self.failed_file}")

    def get_parsing_summary(self, genomes: List[Dict]) -> Dict[str, any]:
        """
        Generate parsing summary statistics.

        Parameters
        ----------
        genomes : List[Dict]
            List of parsed genome data dictionaries

        Returns
        -------
        Dict[str, any]
            Summary statistics
        """
        if not genomes:
            return {}

        lengths = [g['length'] for g in genomes]
        gc_contents = [g['gc_content'] for g in genomes]

        # Count genome types
        genome_types = {}
        for genome in genomes:
            gtype = genome['genome_type']
            genome_types[gtype] = genome_types.get(gtype, 0) + 1

        summary = {
            'total_genomes': len(genomes),
            'total_bases': sum(lengths),
            'mean_length': sum(lengths) / len(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_gc_content': sum(gc_contents) / len(gc_contents),
            'genome_type_distribution': genome_types
        }

        return summary


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Parse RefSeq viral genomes for ViroForge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse all downloaded genomes
  python scripts/parse_genomes.py --input data/refseq --output data/parsed

  # Parse first 100 genomes
  python scripts/parse_genomes.py --input data/refseq --output data/parsed --limit 100

  # Parse with verbose output
  python scripts/parse_genomes.py --input data/refseq --output data/parsed --verbose
        """
    )

    parser.add_argument(
        '--input',
        type=str,
        default='data/refseq',
        help='Input directory containing RefSeq data (default: data/refseq)'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='data/parsed',
        help='Output directory for parsed genomes (default: data/parsed)'
    )

    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='Maximum number of genomes to parse (default: None = all)'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Initialize parser
    genome_parser = GenomeParser(args.input, args.output)

    # Parse genomes
    logger.info("Step 1: Parsing genome files...")
    parsed_genomes, failed_genomes = genome_parser.parse_all_genomes(limit=args.limit)

    if not parsed_genomes:
        logger.error("No genomes were successfully parsed!")
        return

    # Save results
    logger.info("Step 2: Saving parsed genomes...")
    genome_parser.save_parsed_genomes(parsed_genomes)

    logger.info("Step 3: Saving statistics...")
    genome_parser.save_statistics(parsed_genomes)

    if failed_genomes:
        logger.info("Step 4: Saving failed genome list...")
        genome_parser.save_failed_list(failed_genomes)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("PARSING SUMMARY")
    logger.info("=" * 60)

    summary = genome_parser.get_parsing_summary(parsed_genomes)
    logger.info(f"Total genomes parsed: {summary['total_genomes']}")
    logger.info(f"Total bases: {summary['total_bases']:,}")
    logger.info(f"Mean genome length: {summary['mean_length']:.0f} bp")
    logger.info(f"Length range: {summary['min_length']:,} - {summary['max_length']:,} bp")
    logger.info(f"Mean GC content: {summary['mean_gc_content']:.3f}")
    logger.info("\nGenome type distribution:")
    for gtype, count in sorted(summary['genome_type_distribution'].items()):
        logger.info(f"  {gtype}: {count}")

    if failed_genomes:
        logger.info(f"\nFailed genomes: {len(failed_genomes)}")

    logger.info("=" * 60)
    logger.info("✓ Parsing complete!")
    logger.info(f"  Output directory: {genome_parser.output_dir}")


if __name__ == "__main__":
    main()
