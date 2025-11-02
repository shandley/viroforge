#!/usr/bin/env python3
"""
Validate ViroForge-generated FASTQ datasets

Comprehensive validation suite for testing generated datasets independent of any
specific analysis pipeline. Useful for quality control and verification.

Usage:
    python scripts/validate_fastq_dataset.py --dataset data/benchmark_datasets/collection_16_cov10x_novaseq_vlp

    # Verbose output
    python scripts/validate_fastq_dataset.py --dataset my_dataset --verbose

    # Save report
    python scripts/validate_fastq_dataset.py --dataset my_dataset --report validation_report.txt
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional
import gzip

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class DatasetValidator:
    """Validate ViroForge-generated datasets"""

    def __init__(self, dataset_dir: Path, verbose: bool = False):
        self.dataset_dir = Path(dataset_dir)
        self.verbose = verbose
        self.errors = []
        self.warnings = []
        self.stats = {}

        # Expected directories
        self.fasta_dir = self.dataset_dir / 'fasta'
        self.fastq_dir = self.dataset_dir / 'fastq'
        self.metadata_dir = self.dataset_dir / 'metadata'

    def log_error(self, message: str):
        """Log an error"""
        self.errors.append(message)
        logger.error(message)

    def log_warning(self, message: str):
        """Log a warning"""
        self.warnings.append(message)
        logger.warning(message)

    def log_info(self, message: str):
        """Log info if verbose"""
        if self.verbose:
            logger.info(message)

    def validate_all(self) -> bool:
        """Run all validation tests"""
        logger.info(f"Validating dataset: {self.dataset_dir}")
        logger.info("=" * 80)

        # Test 1: Directory structure
        logger.info("Test 1: Validating directory structure...")
        if not self.validate_directory_structure():
            return False

        # Test 2: File presence
        logger.info("Test 2: Checking required files...")
        files = self.check_required_files()
        if not files:
            return False

        # Test 3: Metadata integrity
        logger.info("Test 3: Validating metadata integrity...")
        metadata = self.validate_metadata()
        if not metadata:
            return False

        # Test 4: FASTA validation
        logger.info("Test 4: Validating FASTA reference...")
        fasta_genomes = self.validate_fasta(files['fasta'])
        if not fasta_genomes:
            return False

        # Test 5: FASTQ format validation
        logger.info("Test 5: Validating FASTQ format...")
        fastq_stats = self.validate_fastq_format(files['r1'], files['r2'])
        if not fastq_stats:
            return False

        # Test 6: Abundance validation
        logger.info("Test 6: Validating abundances...")
        if not self.validate_abundances(metadata, fasta_genomes):
            return False

        # Test 7: Coverage/read count validation
        logger.info("Test 7: Validating coverage and read counts...")
        if not self.validate_coverage(metadata, fasta_genomes, fastq_stats):
            return False

        # Test 8: Quality score validation
        logger.info("Test 8: Validating quality scores...")
        if not self.validate_quality_scores(files['r1'], files['r2']):
            return False

        # Test 9: Compositional checks
        logger.info("Test 9: Checking composition...")
        if not self.validate_composition(metadata):
            return False

        logger.info("=" * 80)
        return True

    def validate_directory_structure(self) -> bool:
        """Validate expected directory structure"""
        if not self.dataset_dir.exists():
            self.log_error(f"Dataset directory not found: {self.dataset_dir}")
            return False

        missing = []
        for dir_path in [self.fasta_dir, self.fastq_dir, self.metadata_dir]:
            if not dir_path.exists():
                missing.append(dir_path.name)

        if missing:
            self.log_error(f"Missing directories: {', '.join(missing)}")
            return False

        self.log_info("✓ Directory structure valid")
        return True

    def check_required_files(self) -> Optional[Dict[str, Path]]:
        """Check for required files"""
        files = {}

        # Find FASTA file
        fasta_files = list(self.fasta_dir.glob('*.fasta')) + list(self.fasta_dir.glob('*.fa'))
        if not fasta_files:
            self.log_error("No FASTA file found in fasta/")
            return None
        if len(fasta_files) > 1:
            self.log_warning(f"Multiple FASTA files found, using: {fasta_files[0].name}")
        files['fasta'] = fasta_files[0]

        # Find FASTQ files
        r1_files = list(self.fastq_dir.glob('*_R1.fastq')) + list(self.fastq_dir.glob('*_R1.fq'))
        r2_files = list(self.fastq_dir.glob('*_R2.fastq')) + list(self.fastq_dir.glob('*_R2.fq'))

        if not r1_files:
            self.log_error("No R1 FASTQ file found in fastq/")
            return None
        if not r2_files:
            self.log_error("No R2 FASTQ file found in fastq/")
            return None

        files['r1'] = r1_files[0]
        files['r2'] = r2_files[0]

        # Find metadata files
        json_files = list(self.metadata_dir.glob('*_metadata.json'))
        tsv_files = list(self.metadata_dir.glob('*_composition.tsv'))
        abund_files = list(self.metadata_dir.glob('*_abundances.txt'))

        if not json_files:
            self.log_error("No metadata JSON file found")
            return None
        if not tsv_files:
            self.log_warning("No composition TSV file found")
        if not abund_files:
            self.log_warning("No abundances file found")

        files['metadata'] = json_files[0]
        files['composition'] = tsv_files[0] if tsv_files else None
        files['abundances'] = abund_files[0] if abund_files else None

        self.log_info(f"✓ Found all required files")
        return files

    def validate_metadata(self) -> Optional[Dict]:
        """Validate metadata JSON"""
        json_files = list(self.metadata_dir.glob('*_metadata.json'))
        if not json_files:
            self.log_error("No metadata JSON found")
            return None

        try:
            with open(json_files[0]) as f:
                metadata = json.load(f)
        except json.JSONDecodeError as e:
            self.log_error(f"Invalid JSON: {e}")
            return None

        # Check required fields
        required_fields = ['genomes']  # Minimal required field
        missing = [f for f in required_fields if f not in metadata]
        if missing:
            self.log_error(f"Missing metadata fields: {', '.join(missing)}")
            return None

        # Validate genomes list
        if not isinstance(metadata['genomes'], list):
            self.log_error("'genomes' field is not a list")
            return None

        if len(metadata['genomes']) == 0:
            self.log_error("No genomes in metadata")
            return None

        # Check genome fields (handle both 'abundance' and 'relative_abundance')
        required_genome_fields = ['genome_id', 'length']
        for i, genome in enumerate(metadata['genomes']):
            missing = [f for f in required_genome_fields if f not in genome]
            if missing:
                self.log_error(f"Genome {i} missing fields: {', '.join(missing)}")
                return None

            # Check for abundance field (either name works)
            if 'abundance' not in genome and 'relative_abundance' not in genome:
                self.log_error(f"Genome {i} missing abundance field")
                return None

        self.log_info(f"✓ Metadata valid ({len(metadata['genomes'])} genomes)")
        return metadata

    def validate_fasta(self, fasta_path: Path) -> Optional[Dict[str, int]]:
        """Validate FASTA and return genome lengths"""
        genomes = {}
        current_id = None
        current_seq = []

        try:
            with open(fasta_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous
                        if current_id:
                            genomes[current_id] = len(''.join(current_seq))
                        # Start new
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)

                # Save last
                if current_id:
                    genomes[current_id] = len(''.join(current_seq))

        except Exception as e:
            self.log_error(f"Error reading FASTA: {e}")
            return None

        if not genomes:
            self.log_error("No sequences found in FASTA")
            return None

        # Check for empty sequences
        empty = [gid for gid, length in genomes.items() if length == 0]
        if empty:
            self.log_error(f"Empty sequences found: {', '.join(empty[:5])}")
            return None

        self.log_info(f"✓ FASTA valid ({len(genomes)} sequences)")
        return genomes

    def validate_fastq_format(self, r1_path: Path, r2_path: Path) -> Optional[Dict]:
        """Validate FASTQ format and return statistics"""
        stats = {
            'r1_reads': 0,
            'r2_reads': 0,
            'read_lengths': [],
            'quality_lengths': []
        }

        # Validate R1
        if not self._validate_single_fastq(r1_path, stats, 'r1'):
            return None

        # Validate R2
        if not self._validate_single_fastq(r2_path, stats, 'r2'):
            return None

        # Check read counts match
        if stats['r1_reads'] != stats['r2_reads']:
            self.log_error(f"Read count mismatch: R1={stats['r1_reads']}, R2={stats['r2_reads']}")
            return None

        # Calculate statistics
        if stats['read_lengths']:
            stats['mean_read_length'] = sum(stats['read_lengths']) / len(stats['read_lengths'])
            stats['min_read_length'] = min(stats['read_lengths'])
            stats['max_read_length'] = max(stats['read_lengths'])

        self.log_info(f"✓ FASTQ format valid ({stats['r1_reads']:,} read pairs)")
        if self.verbose and stats['read_lengths']:
            self.log_info(f"  Read length: {stats['mean_read_length']:.1f} bp (range: {stats['min_read_length']}-{stats['max_read_length']})")

        return stats

    def _validate_single_fastq(self, fastq_path: Path, stats: Dict, read_type: str) -> bool:
        """Validate a single FASTQ file"""
        read_count = 0
        line_num = 0

        try:
            open_func = gzip.open if str(fastq_path).endswith('.gz') else open
            mode = 'rt' if str(fastq_path).endswith('.gz') else 'r'

            with open_func(fastq_path, mode) as f:
                while True:
                    # Read 4-line record
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline().strip()
                    plus = f.readline()
                    qual = f.readline().strip()

                    line_num += 4
                    read_count += 1

                    # Validate header
                    if not header.startswith('@'):
                        self.log_error(f"{fastq_path.name} line {line_num-3}: Invalid header")
                        return False

                    # Validate plus line
                    if not plus.startswith('+'):
                        self.log_error(f"{fastq_path.name} line {line_num-1}: Invalid plus line")
                        return False

                    # Validate lengths match
                    if len(seq) != len(qual):
                        self.log_error(f"{fastq_path.name} line {line_num-3}: Length mismatch (seq={len(seq)}, qual={len(qual)})")
                        return False

                    # Check for empty
                    if len(seq) == 0:
                        self.log_error(f"{fastq_path.name} line {line_num-3}: Empty sequence")
                        return False

                    # Store lengths (sample first 1000)
                    if read_count <= 1000:
                        stats['read_lengths'].append(len(seq))
                        stats['quality_lengths'].append(len(qual))

        except Exception as e:
            self.log_error(f"Error reading {fastq_path.name}: {e}")
            return False

        stats[f'{read_type}_reads'] = read_count
        return True

    def validate_abundances(self, metadata: Dict, fasta_genomes: Dict[str, int]) -> bool:
        """Validate abundance values"""
        # Check all genomes have abundances
        abundances = []
        genome_ids_in_metadata = set()

        for genome in metadata['genomes']:
            genome_id = genome['genome_id']
            # Handle both 'abundance' and 'relative_abundance'
            abundance = genome.get('abundance') or genome.get('relative_abundance')
            genome_ids_in_metadata.add(genome_id)

            # Check abundance range
            if abundance < 0:
                self.log_error(f"Negative abundance for {genome_id}: {abundance}")
                return False
            if abundance > 1:
                self.log_error(f"Abundance > 1 for {genome_id}: {abundance}")
                return False

            abundances.append(abundance)

        # Check abundances sum to ~1.0
        total = sum(abundances)
        if not (0.99 <= total <= 1.01):
            self.log_warning(f"Abundances sum to {total:.6f} (expected ~1.0)")

        # Check metadata genomes match FASTA
        fasta_ids = set(fasta_genomes.keys())
        if genome_ids_in_metadata != fasta_ids:
            missing_in_metadata = fasta_ids - genome_ids_in_metadata
            missing_in_fasta = genome_ids_in_metadata - fasta_ids

            if missing_in_metadata:
                self.log_error(f"Genomes in FASTA but not metadata: {len(missing_in_metadata)}")
                if self.verbose:
                    for gid in list(missing_in_metadata)[:5]:
                        self.log_info(f"  {gid}")
                return False

            if missing_in_fasta:
                self.log_error(f"Genomes in metadata but not FASTA: {len(missing_in_fasta)}")
                if self.verbose:
                    for gid in list(missing_in_fasta)[:5]:
                        self.log_info(f"  {gid}")
                return False

        self.log_info(f"✓ Abundances valid (sum={total:.6f})")
        return True

    def validate_coverage(self, metadata: Dict, fasta_genomes: Dict[str, int],
                         fastq_stats: Dict) -> bool:
        """Validate coverage calculations"""
        # Get generation params (handle both 'generation_params' and 'configuration')
        params = metadata.get('generation_params') or metadata.get('configuration', {})
        target_coverage = params.get('coverage', None)
        read_length = params.get('read_length', None)

        if target_coverage is None:
            self.log_warning("No target coverage in metadata, skipping coverage validation")
            return True

        if read_length is None:
            self.log_warning("No read length in metadata, skipping coverage validation")
            return True

        # Calculate total genome length
        total_length = sum(fasta_genomes.values())

        # Calculate expected reads
        # n_reads = (total_length * coverage) / (2 * read_length)
        expected_reads = (total_length * target_coverage) / (2 * read_length)

        # Get actual reads
        actual_reads = fastq_stats['r1_reads']

        # Allow 5% tolerance
        ratio = actual_reads / expected_reads
        if not (0.95 <= ratio <= 1.05):
            self.log_warning(f"Read count deviation: expected={expected_reads:,.0f}, actual={actual_reads:,} (ratio={ratio:.3f})")
        else:
            self.log_info(f"✓ Read count matches coverage (expected={expected_reads:,.0f}, actual={actual_reads:,})")

        # Calculate actual coverage
        total_bases = actual_reads * 2 * read_length  # Both R1 and R2
        actual_coverage = total_bases / total_length

        self.log_info(f"✓ Actual coverage: {actual_coverage:.2f}x (target: {target_coverage}x)")

        return True

    def validate_quality_scores(self, r1_path: Path, r2_path: Path) -> bool:
        """Validate quality score distributions"""
        qual_scores = []
        reads_sampled = 0
        max_sample = 10000

        try:
            open_func = gzip.open if str(r1_path).endswith('.gz') else open
            mode = 'rt' if str(r1_path).endswith('.gz') else 'r'

            with open_func(r1_path, mode) as f:
                while reads_sampled < max_sample:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline()
                    plus = f.readline()
                    qual = f.readline().strip()

                    # Convert quality string to scores
                    for q in qual:
                        qual_scores.append(ord(q))

                    reads_sampled += 1

        except Exception as e:
            self.log_error(f"Error reading quality scores: {e}")
            return False

        if not qual_scores:
            self.log_error("No quality scores found")
            return False

        # Check quality score range (Phred33: 33-126)
        min_qual = min(qual_scores)
        max_qual = max(qual_scores)

        if min_qual < 33:
            self.log_error(f"Quality scores below Phred33 minimum (33): {min_qual}")
            return False

        if max_qual > 126:
            self.log_error(f"Quality scores above Phred33 maximum (126): {max_qual}")
            return False

        mean_qual = sum(qual_scores) / len(qual_scores)

        self.log_info(f"✓ Quality scores valid (range: {min_qual}-{max_qual}, mean: {mean_qual:.1f})")

        return True

    def validate_composition(self, metadata: Dict) -> bool:
        """Validate taxonomic composition"""
        # Calculate Shannon diversity (handle both abundance field names)
        abundances = []
        for g in metadata['genomes']:
            abund = g.get('abundance') or g.get('relative_abundance', 0)
            if abund > 0:
                abundances.append(abund)

        if abundances:
            import math
            shannon = -sum(a * math.log(a) for a in abundances if a > 0)

            # Expected ranges based on collection type
            collection_info = metadata.get('collection_info') or metadata.get('collection', {})
            collection_name = collection_info.get('name', '')

            self.log_info(f"✓ Shannon diversity: {shannon:.3f}")

            if 'skin' in collection_name.lower() and shannon > 3.0:
                self.log_warning(f"High diversity for skin collection: {shannon:.3f}")
            elif 'marine' in collection_name.lower() and shannon < 3.0:
                self.log_warning(f"Low diversity for marine collection: {shannon:.3f}")

        # Count families/genera if taxonomy available
        families = Counter()
        genera = Counter()

        for genome in metadata['genomes']:
            tax = genome.get('taxonomy', {})
            if tax.get('family'):
                families[tax['family']] += 1
            if tax.get('genus'):
                genera[tax['genus']] += 1

        if families:
            self.log_info(f"  Families: {len(families)}")
            if self.verbose:
                for family, count in families.most_common(5):
                    self.log_info(f"    {family}: {count}")

        if genera:
            self.log_info(f"  Genera: {len(genera)}")

        return True

    def generate_report(self) -> str:
        """Generate validation report"""
        lines = []
        lines.append("=" * 80)
        lines.append("VIROFORGE DATASET VALIDATION REPORT")
        lines.append("=" * 80)
        lines.append(f"Dataset: {self.dataset_dir}")
        lines.append("")

        if not self.errors and not self.warnings:
            lines.append("STATUS: ✓ PASSED - All validation tests passed")
        elif self.errors:
            lines.append(f"STATUS: ✗ FAILED - {len(self.errors)} error(s) found")
        else:
            lines.append(f"STATUS: ⚠ WARNING - {len(self.warnings)} warning(s) found")

        lines.append("")

        if self.errors:
            lines.append("ERRORS:")
            for error in self.errors:
                lines.append(f"  ✗ {error}")
            lines.append("")

        if self.warnings:
            lines.append("WARNINGS:")
            for warning in self.warnings:
                lines.append(f"  ⚠ {warning}")
            lines.append("")

        lines.append("=" * 80)

        return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Validate ViroForge-generated FASTQ datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate a dataset
  python scripts/validate_fastq_dataset.py --dataset data/benchmark_datasets/collection_16_cov10x_novaseq_vlp

  # Verbose output
  python scripts/validate_fastq_dataset.py --dataset my_dataset --verbose

  # Save report to file
  python scripts/validate_fastq_dataset.py --dataset my_dataset --report validation_report.txt
        """
    )

    parser.add_argument(
        '--dataset',
        type=str,
        required=True,
        help='Path to dataset directory'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose output'
    )

    parser.add_argument(
        '--report',
        type=str,
        help='Save validation report to file'
    )

    args = parser.parse_args()

    # Validate dataset
    validator = DatasetValidator(args.dataset, verbose=args.verbose)
    success = validator.validate_all()

    # Generate report
    report = validator.generate_report()
    print("\n" + report)

    # Save report if requested
    if args.report:
        with open(args.report, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to: {args.report}")

    # Exit with appropriate code
    sys.exit(0 if success and not validator.errors else 1)


if __name__ == '__main__':
    main()
