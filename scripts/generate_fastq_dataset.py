#!/usr/bin/env python3
"""
ViroForge FASTQ Dataset Generator

Generate realistic synthetic virome FASTQ files from curated body site collections.
Integrates with existing ViroForge infrastructure for VLP enrichment, amplification
bias, and platform artifacts.

Usage:
    # Generate gut virome with standard VLP enrichment
    python generate_fastq_dataset.py \\
        --collection gut \\
        --output data/fastq/gut_standard \\
        --coverage 10 \\
        --vlp standard \\
        --amplification rdab \\
        --platform novaseq

    # Generate marine virome without VLP (bulk metag comparison)
    python generate_fastq_dataset.py \\
        --collection marine \\
        --output data/fastq/marine_bulk \\
        --coverage 10 \\
        --no-vlp \\
        --platform miseq

Author: ViroForge Development Team
Date: 2025-11-01
"""

import argparse
import sys
import logging
import sqlite3
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
from datetime import datetime
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: Biopython not installed. Install with: pip install biopython")
    sys.exit(1)

# Import ViroForge modules
try:
    from viroforge.enrichment.vlp import VLPEnrichment, VLPProtocol
    from viroforge.core.contamination import (
        create_contamination_profile,
        create_rna_contamination_profile,
        ContaminationProfile
    )
    from viroforge.amplification import (
        rdab_40_cycles,
        rdab_30_cycles,
        mda_standard,
        mda_overnight,
        linker_standard,
        no_amplification
    )
    from viroforge.workflows.rna_virome import (
        RNAViromeWorkflow,
        ReverseTranscription,
        RiboDepletion,
        RNADegradation,
        RNAVirusType,
        PrimerType,
        RiboDepleteMethod,
        infer_virus_type_from_taxonomy
    )
    from viroforge.core.collection import CollectionLoader, label_fastq_headers
    from viroforge.generator import FASTQGenerator, run_generation
except ImportError as e:
    print(f"Error importing ViroForge modules: {e}")
    print("Make sure you're running from the ViroForge root directory")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge FASTQ Dataset Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate gut virome with VLP enrichment
  python generate_fastq_dataset.py \\
      --collection-id 9 \\
      --output data/fastq/gut_standard \\
      --coverage 10 \\
      --platform novaseq

  # List available collections
  python generate_fastq_dataset.py --list-collections

  # Generate without VLP (bulk comparison)
  python generate_fastq_dataset.py \\
      --collection-id 13 \\
      --output data/fastq/marine_bulk \\
      --no-vlp \\
      --coverage 10
        """
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database'
    )

    parser.add_argument(
        '--list-collections',
        action='store_true',
        help='List available collections and exit'
    )

    parser.add_argument(
        '--collection-id',
        type=int,
        help='Collection ID to generate FASTQs from'
    )

    parser.add_argument(
        '--output',
        help='Output directory for generated files'
    )

    parser.add_argument(
        '--coverage',
        type=float,
        default=10.0,
        help='Mean coverage depth (default: 10x)'
    )

    parser.add_argument(
        '--n-reads',
        type=int,
        help='Number of reads (overrides --coverage)'
    )

    parser.add_argument(
        '--read-length',
        type=int,
        default=150,
        help='Read length in bp (default: 150)'
    )

    parser.add_argument(
        '--insert-size',
        type=int,
        default=350,
        help='Insert size for paired-end (default: 350)'
    )

    parser.add_argument(
        '--platform',
        choices=['novaseq', 'miseq', 'hiseq', 'pacbio-hifi', 'nanopore'],
        default='novaseq',
        help='Sequencing platform (default: novaseq). ' +
             'Short-read: novaseq, miseq, hiseq (NOTE: these run InSilicoSeq in ' +
             'basic mode and currently produce identical 125 bp reads regardless ' +
             'of choice or --read-length; they are interchangeable). ' +
             'Long-read: pacbio-hifi, nanopore'
    )

    # Long-read specific options
    lr_group = parser.add_argument_group('long-read options (PacBio HiFi, Nanopore)')
    lr_group.add_argument(
        '--depth',
        type=float,
        default=10.0,
        help='Sequencing depth for long reads (default: 10x). ' +
             'Note: --coverage is for short reads, --depth is for long reads'
    )
    lr_group.add_argument(
        '--pacbio-passes',
        type=int,
        default=10,
        help='Number of CCS passes for PacBio HiFi (default: 10)'
    )
    lr_group.add_argument(
        '--pacbio-read-length',
        type=int,
        default=15000,
        help='Mean read length for PacBio HiFi in bp (default: 15000)'
    )
    lr_group.add_argument(
        '--ont-chemistry',
        choices=['R9.4', 'R10.4'],
        default='R10.4',
        help='Nanopore chemistry version (default: R10.4)'
    )
    lr_group.add_argument(
        '--ont-read-length',
        type=int,
        default=20000,
        help='Mean read length for Nanopore in bp (default: 20000)'
    )

    parser.add_argument(
        '--no-vlp',
        action='store_true',
        help='Skip VLP enrichment simulation (generate bulk metagenome)'
    )

    parser.add_argument(
        '--vlp-protocol',
        choices=['tangential_flow', 'tangential_flow_045', 'tangential_flow_01',
                 'syringe', 'syringe_045',
                 'ultracentrifugation', 'norgen'],
        default='tangential_flow',
        help='VLP enrichment protocol (default: tangential_flow). '
             'Filtration pore sizes: tangential_flow=0.2um, '
             'tangential_flow_045=0.45um (giant viruses), '
             'tangential_flow_01=0.1um (tight, phage-focused), '
             'syringe=0.2um, syringe_045=0.45um. '
             'Or use --pore-size for custom values.'
    )

    parser.add_argument(
        '--pore-size',
        type=float,
        default=None,
        help='Custom filtration pore size in micrometers (overrides protocol default). '
             'Common values: 0.1, 0.2, 0.22, 0.45. '
             'Only applies to filtration protocols (tangential_flow, syringe).'
    )

    parser.add_argument(
        '--contamination-level',
        choices=['clean', 'realistic', 'heavy'],
        default='realistic',
        help='Contamination level (default: realistic)'
    )

    parser.add_argument(
        '--molecule-type',
        choices=['dna', 'rna'],
        default='dna',
        help='Molecule type: dna (default) or rna. ' +
             'RNA viromes include reverse transcription, rRNA depletion (Ribo-Zero), ' +
             'and RNA degradation modeling.'
    )

    parser.add_argument(
        '--rna-primer',
        choices=['random_hexamer', 'random_octamer', 'oligo_dt', 'specific'],
        default='random_hexamer',
        help='RT primer type for RNA viromes (default: random_hexamer)'
    )

    parser.add_argument(
        '--rna-depletion',
        choices=['ribo_zero', 'ribominus', 'none'],
        default='ribo_zero',
        help='rRNA depletion method for RNA viromes (default: ribo_zero)'
    )

    parser.add_argument(
        '--amplification',
        choices=['none', 'rdab', 'rdab-30', 'mda', 'mda-long', 'linker'],
        default='linker',
        help='Library preparation amplification method (default: linker). ' +
             'rdab=RdAB 40 cycles (standard), rdab-30=RdAB 30 cycles (moderate), ' +
             'mda=MDA 4h (low biomass), mda-long=MDA 16h (overnight), ' +
             'linker=Linker-based (minimal bias), none=No amplification'
    )

    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be generated without running ISS'
    )

    parser.add_argument(
        '--enable-benchmarking',
        action='store_true',
        default=True,
        help='Enable benchmarking metadata export (contamination manifest, expected coverage). ' +
             'Default: True. Use --no-enable-benchmarking to disable.'
    )

    parser.add_argument(
        '--no-enable-benchmarking',
        action='store_false',
        dest='enable_benchmarking',
        help='Disable benchmarking metadata export (smaller metadata files)'
    )

    # Contamination realism options
    contam_group = parser.add_argument_group('contamination realism')
    contam_group.add_argument(
        '--adapter-rate',
        type=float,
        default=0.03,
        help='Fraction of reads with adapter read-through (0.0-1.0, default: 0.03). '
             'Set to 0.0 to disable.'
    )
    contam_group.add_argument(
        '--adapter-type',
        choices=['truseq', 'nextera'],
        default='truseq',
        help='Illumina adapter type (default: truseq)'
    )
    contam_group.add_argument(
        '--mean-insert-size',
        type=int,
        default=None,
        help='Mean insert size for insert-size-driven adapter contamination (bp). '
             'When set, adapter rate is determined naturally by the insert size '
             'distribution vs read length. Overrides --adapter-rate. '
             'Example: --mean-insert-size 200 with 150bp reads gives ~2%% adapter rate.'
    )
    contam_group.add_argument(
        '--insert-size-sd',
        type=int,
        default=50,
        help='Standard deviation of insert size distribution (default: 50bp). '
             'Used with --mean-insert-size.'
    )
    contam_group.add_argument(
        '--chimera-rate',
        type=float,
        default=0.0,
        help='Fraction of reads with internal adapter chimeras (0.0-1.0, '
             'default: 0.0). Models chimeric ligation events where adapter '
             'sequence appears inside a read, not just at the 3\' end.'
    )
    contam_group.add_argument(
        '--host-genome',
        type=str,
        default=None,
        help='Path to full host genome FASTA for host DNA contamination '
             '(overrides bundled fragments)'
    )
    contam_group.add_argument(
        '--rrna-database',
        type=str,
        default=None,
        help='Path to rRNA database FASTA (overrides bundled rRNA references)'
    )
    contam_group.add_argument(
        '--no-real-contaminants',
        action='store_true',
        help='Force synthetic contamination sequences (old behavior). '
             'By default, ViroForge uses bundled real reference sequences.'
    )
    contam_group.add_argument(
        '--low-complexity-rate',
        type=float,
        default=0.005,
        help='Fraction of reads replaced with low-complexity artifacts (0.0-1.0, '
             'default: 0.005). Set to 0.0 to disable. '
             'Models homopolymer runs, dinucleotide repeats, simple repeats, '
             'and low-entropy sequences from adapter dimers and PCR failures.'
    )
    contam_group.add_argument(
        '--entropy-range',
        type=str,
        default=None,
        help='Generate low-complexity reads with controlled Shannon entropy '
             'instead of predefined artifact types. Format: MIN-MAX (e.g., '
             '0.3-0.7). Entropy is in bits (0.0=homopolymer, 2.0=random). '
             'Useful for testing complexity filter threshold sensitivity.'
    )
    contam_group.add_argument(
        '--duplicate-rate',
        type=float,
        default=0.10,
        help='Fraction of reads that become PCR duplicate templates (0.0-1.0, '
             'default: 0.10). Set to 0.0 to disable. '
             'Each template gets 1-5 copies (geometric distribution).'
    )
    contam_group.add_argument(
        '--duplicate-max-copies',
        type=int,
        default=5,
        help='Maximum copies per PCR duplicate template (default: 5).'
    )
    contam_group.add_argument(
        '--duplicate-error-rate',
        type=float,
        default=0.001,
        help='Per-base substitution rate in PCR copies (default: 0.001). '
             'Models PCR polymerase errors. Set to 0 for exact duplicates.'
    )
    contam_group.add_argument(
        '--mda-chimera-rate',
        type=float,
        default=0.0,
        help='Fraction of reads replaced with MDA chimeric artifacts '
             '(0.0-0.3, default: 0.0, or 0.15 when --amplification is mda/mda-long '
             'and this is left at 0). Models phi29 branch-migration chimeras where '
             'a read joins sequence from two different genomic regions.'
    )

    # ERV injection options
    erv_group = parser.add_argument_group('retroviral read injection')
    erv_group.add_argument(
        '--erv-endogenous-rate',
        type=float,
        default=0.0,
        help='Abundance fraction for endogenous retroviral reads (0.0-1.0, '
             'default: 0.0). Models polymorphic ERV insertions from HERV '
             'consensus sequences. Requires --herv-fasta or VIROFORGE_HERV_DB.'
    )
    erv_group.add_argument(
        '--erv-exogenous-rate',
        type=float,
        default=0.0,
        help='Abundance fraction for exogenous retroviral reads (0.0-1.0, '
             'default: 0.0). Models active retroviral infections from '
             'Retroviridae in the ViroForge database (44 genomes).'
    )
    erv_group.add_argument(
        '--herv-fasta',
        type=str,
        default=None,
        help='Path to HERV consensus FASTA for endogenous ERV injection. '
             'Can also be set via VIROFORGE_HERV_DB environment variable.'
    )
    erv_group.add_argument(
        '--erv-exogenous-viruses',
        type=str,
        nargs='+',
        default=None,
        help='Specific exogenous retroviruses to include (name patterns). '
             'Default: random sample from all non-endogenous Retroviridae. '
             'Example: --erv-exogenous-viruses "Human immunodeficiency" "Primate T-lymphotropic"'
    )

    # Dark matter (unclassified viral sequences)
    dm_group = parser.add_argument_group('viral dark matter')
    dm_group.add_argument(
        '--dark-matter-fraction',
        type=float,
        default=0.30,
        help='Fraction of reads from unclassified viral genomes (0.0-1.0, '
             'default: 0.30). Adds real but taxonomically unclassified '
             'sequences from RefSeq to simulate the 60-90%% of reads in '
             'real viromes that do not match known references. '
             'Set to 0.0 to disable.'
    )
    dm_group.add_argument(
        '--dark-matter-count',
        type=int,
        default=None,
        help='Number of dark matter genomes to include (default: auto, '
             'scales with collection size). Overrides automatic calculation.'
    )

    args = parser.parse_args()
    return run_generation(args)


if __name__ == '__main__':
    sys.exit(main() or 0)
