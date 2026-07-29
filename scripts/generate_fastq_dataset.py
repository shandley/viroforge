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
    from viroforge.generator import FASTQGenerator, build_parser, run_generation
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
    args = build_parser().parse_args()
    return run_generation(args)


if __name__ == '__main__':
    sys.exit(main() or 0)
