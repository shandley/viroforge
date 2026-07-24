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
    from viroforge.generator import FASTQGenerator
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

    # Load collections
    loader = CollectionLoader(args.database)

    # List collections if requested
    if args.list_collections:
        collections = loader.list_collections()
        print("\nAvailable Collections:")
        print("=" * 80)
        for coll in collections:
            print(f"ID: {coll['collection_id']}")
            print(f"  Name: {coll['collection_name']}")
            print(f"  Genomes: {coll['n_genomes']}")
            print(f"  Description: {coll.get('description', 'N/A')}")
            print()
        return

    # Validate required arguments
    if not args.collection_id:
        parser.error("--collection-id required (use --list-collections to see options)")

    if not args.output:
        parser.error("--output required")

    # Load collection
    collection_meta, genomes = loader.load_collection(args.collection_id)

    # Add viral dark matter (unclassified genomes)
    dark_matter_stats = {}
    if args.dark_matter_fraction > 0:
        dm_fraction = min(args.dark_matter_fraction, 0.95)  # Cap at 95%

        # Determine how many dark matter genomes to add
        n_collection = len(genomes)
        if args.dark_matter_count is not None:
            n_dark_matter = args.dark_matter_count
        else:
            # Scale with collection size: roughly 0.5x the collection size
            # e.g., 90 genomes → ~45 dark matter genomes
            n_dark_matter = max(10, int(n_collection * 0.5))

        # Load dark matter genomes
        existing_ids = {g['genome_id'] for g in genomes}
        dark_matter_genomes = loader.load_dark_matter_genomes(
            n_genomes=n_dark_matter,
            exclude_ids=existing_ids,
            random_seed=args.seed
        )

        if dark_matter_genomes:
            # Rescale existing abundances: known viral fraction = (1 - dm_fraction)
            known_fraction = 1.0 - dm_fraction
            for g in genomes:
                g['relative_abundance'] *= known_fraction

            # Distribute dark matter abundance using log-normal distribution
            rng = np.random.default_rng(args.seed + 1000)
            dm_raw = rng.lognormal(mean=0.0, sigma=1.5, size=len(dark_matter_genomes))
            dm_raw /= dm_raw.sum()
            dm_abundances = dm_raw * dm_fraction

            for i, dm_genome in enumerate(dark_matter_genomes):
                dm_genome['relative_abundance'] = float(dm_abundances[i])
                dm_genome['abundance_rank'] = n_collection + i + 1
                genomes.append(dm_genome)

            dark_matter_stats = {
                'dark_matter_fraction': dm_fraction,
                'dark_matter_genomes': len(dark_matter_genomes),
                'known_viral_genomes': n_collection,
                'total_genomes': len(genomes),
            }

            logger.info(f"Added {len(dark_matter_genomes)} dark matter genomes "
                       f"({dm_fraction*100:.0f}% of total abundance)")
        else:
            logger.warning("No dark matter genomes available — proceeding without")
    else:
        logger.info("Dark matter disabled (--dark-matter-fraction 0.0)")

    # Create generator (sanitize collection name for file system)
    import re
    collection_name = collection_meta['collection_name']
    # Remove/replace invalid filename characters
    collection_name = re.sub(r'[^\w\s-]', '', collection_name)  # Remove special chars
    collection_name = collection_name.replace(' ', '_').lower()  # Replace spaces, lowercase

    generator = FASTQGenerator(
        output_dir=args.output,
        collection_name=collection_name,
        molecule_type=args.molecule_type,
        random_seed=args.seed
    )

    # Prepare RNA workflow parameters (if RNA virome)
    rna_workflow_params = None
    if args.molecule_type == 'rna':
        # Map string args to enum values
        primer_map = {
            'random_hexamer': PrimerType.RANDOM_HEXAMER,
            'random_octamer': PrimerType.RANDOM_OCTAMER,
            'oligo_dt': PrimerType.OLIGO_DT,
            'specific': PrimerType.SPECIFIC
        }
        ribo_map = {
            'ribo_zero': RiboDepleteMethod.RIBO_ZERO,
            'ribominus': RiboDepleteMethod.RIBOMINUS,
            'none': RiboDepleteMethod.NONE
        }

        rna_workflow_params = {
            'primer_type': primer_map[args.rna_primer],
            'ribo_method': ribo_map[args.rna_depletion],
            'degradation_rate': 0.20  # Default 20% degradation
        }

    # Determine read type based on platform
    is_long_read = args.platform in ['pacbio-hifi', 'nanopore']
    read_type = "long" if is_long_read else "short"

    logger.info(f"Platform: {args.platform} (read_type={read_type})")

    # Prepare genomes with VLP enrichment
    vlp_protocol = None if args.no_vlp else args.vlp_protocol

    # Handle custom pore size override — map to matching protocol variant
    if args.pore_size is not None and vlp_protocol is not None:
        base = vlp_protocol.split('_0')[0]  # Strip existing pore size suffix
        if base in ('tangential_flow', 'syringe'):
            pore_map = {0.1: '_01', 0.2: '', 0.22: '', 0.45: '_045'}
            suffix = pore_map.get(args.pore_size)
            if suffix is not None:
                vlp_protocol = base + suffix
                logger.info(f"Pore size {args.pore_size} μm → protocol: {vlp_protocol}")
            else:
                # For non-standard pore sizes, dynamically create a config
                # and register it in the protocol map
                from viroforge.enrichment.vlp import VLPProtocol as _VLP
                custom_config = _VLP.with_custom_pore_size(base, args.pore_size)
                custom_key = f"{base}_custom"
                # Monkey-patch: we'll add it to the map inside _apply_vlp_enrichment
                # by setting vlp_protocol to the custom key and storing the config
                vlp_protocol = custom_key
                generator._custom_vlp_configs = getattr(generator, '_custom_vlp_configs', {})
                generator._custom_vlp_configs[custom_key] = custom_config
                logger.info(f"Custom pore size: {args.pore_size} μm for {base}")
        else:
            logger.warning(f"--pore-size ignored: only applies to tangential_flow or syringe, "
                         f"not {vlp_protocol}")

    use_real = not args.no_real_contaminants

    # Build ERV kwargs if any ERV rates are specified
    erv_kwargs = {}
    if args.erv_endogenous_rate > 0:
        erv_kwargs['erv_endogenous_pct'] = args.erv_endogenous_rate * 100
        if args.herv_fasta:
            erv_kwargs['herv_fasta_path'] = Path(args.herv_fasta)
    if args.erv_exogenous_rate > 0:
        erv_kwargs['erv_exogenous_pct'] = args.erv_exogenous_rate * 100
        erv_kwargs['db_path'] = Path(args.database)
        if args.erv_exogenous_viruses:
            erv_kwargs['erv_exogenous_viruses'] = args.erv_exogenous_viruses

    # Use the collection's sample-type contamination baseline when it has one
    # (blood host-heavy, marine host-free, ...); contamination-level then scales it.
    # Falls back to the global preset when the collection predates the defaults.
    collection_defaults = collection_meta if collection_meta.get('default_host_pct') is not None else None

    sequences, abundances, enrichment_stats, contamination_profile = generator.prepare_genomes(
        genomes,
        vlp_protocol=vlp_protocol,
        contamination_level=args.contamination_level,
        rna_workflow_params=rna_workflow_params,
        read_type=read_type,
        use_real_references=use_real,
        erv_kwargs=erv_kwargs or None,
        collection_defaults=collection_defaults,
    )

    # Apply amplification bias
    abundances, amplification_stats = generator.apply_amplification(
        sequences,
        abundances,
        args.amplification
    )

    # Preserve the configured dark-matter fraction through VLP enrichment and
    # amplification bias. Both reweight genomes by their properties, which pulls
    # the realized dark-matter share off the target (observed 0.11-0.48 for a
    # 0.30 target). Rescale the dark-matter and known-viral blocks to the
    # configured ratio while preserving total viral mass and each block's
    # internal (enriched, amplified) structure, then record the realized value.
    if dark_matter_stats:
        dm_ids = {g['genome_id'] for g in genomes if g.get('is_dark_matter')}
        viral_ids = {g['genome_id'] for g in genomes}
        abund = np.asarray(abundances, dtype=float)
        dark_idx = [i for i, s in enumerate(sequences) if s.id in dm_ids]
        known_idx = [i for i, s in enumerate(sequences)
                     if s.id in viral_ids and s.id not in dm_ids]
        d_sum = float(abund[dark_idx].sum()) if dark_idx else 0.0
        k_sum = float(abund[known_idx].sum()) if known_idx else 0.0
        target = dark_matter_stats['dark_matter_fraction']
        if d_sum > 0 and k_sum > 0:
            viral_mass = d_sum + k_sum
            abund[dark_idx] *= (target * viral_mass) / d_sum
            abund[known_idx] *= ((1.0 - target) * viral_mass) / k_sum
            abundances = abund.tolist()
            d2 = float(abund[dark_idx].sum())
            k2 = float(abund[known_idx].sum())
            dark_matter_stats['realized_fraction_of_viral'] = d2 / (d2 + k2)
        else:
            dark_matter_stats['realized_fraction_of_viral'] = (
                d_sum / (d_sum + k_sum) if (d_sum + k_sum) > 0 else None
            )

    # Write FASTA
    fasta_path = generator.write_fasta(sequences, abundances)

    # Configuration for metadata
    config = {
        'coverage': args.coverage,
        'n_reads': args.n_reads,
        'read_length': args.read_length,
        'insert_size': args.insert_size,
        'platform': args.platform,
        'molecule_type': args.molecule_type,
        'vlp_protocol': vlp_protocol if vlp_protocol else 'none',
        'contamination_level': args.contamination_level,
        'amplification': args.amplification,
        'dark_matter': dark_matter_stats if dark_matter_stats else None,
    }

    # Add RNA-specific parameters if applicable
    if args.molecule_type == 'rna':
        config['rna_primer'] = args.rna_primer
        config['rna_depletion'] = args.rna_depletion

    # Export metadata with enrichment and amplification stats
    generator.export_metadata(
        collection_meta,
        sequences,
        abundances,
        config,
        enrichment_stats,
        amplification_stats,
        contamination_profile,
        args.enable_benchmarking,
        db_path=args.database
    )

    if args.dry_run:
        logger.info("Dry run complete - FASTQ generation skipped")
        return

    # Route to appropriate sequencing simulator based on platform
    if is_long_read:
        # ===============================================================
        # LONG-READ SEQUENCING (PacBio HiFi, Nanopore)
        # ===============================================================
        logger.info(f"Generating long-read FASTQ files ({args.platform})...")

        # Import long-read config classes
        from viroforge.simulators.longread import (
            PacBioHiFiConfig,
            NanoporeConfig
        )
        import tempfile

        output_prefix = str(generator.fastq_dir / collection_name)

        # Calculate per-genome depths based on abundance
        # PBSIM3 only accepts a single --depth value, so we must run it
        # per-genome to preserve community composition
        per_genome_depths = {
            seq.id: args.depth * abundance
            for seq, abundance in zip(sequences, abundances)
        }
        min_depth_threshold = 0.01  # Skip genomes below this depth

        # Determine platform and run simulation
        if args.platform == 'pacbio-hifi':
            platform_config = PacBioHiFiConfig(
                passes=args.pacbio_passes,
                read_length_mean=args.pacbio_read_length
            )
            logger.info(f"  PacBio HiFi config: {platform_config.passes} passes, "
                       f"{platform_config.read_length_mean}bp mean read length")
            logger.info(f"  Depth: {args.depth}x")

            # Step 1: Generate CLR with PBSIM3 (per-genome for correct abundances)
            logger.info("  Step 1/2: Generating CLR reads with PBSIM3...")

            # Resolve full path to QSHMM model file
            import shutil as _shutil
            qshmm_name = platform_config.accuracy_model
            if not qshmm_name.endswith('.model'):
                qshmm_name += '.model'
            qshmm_path = qshmm_name
            pbsim_path = _shutil.which('pbsim')
            if pbsim_path:
                conda_data = Path(pbsim_path).parent.parent / 'data' / qshmm_name
                if conda_data.exists():
                    qshmm_path = str(conda_data)

            # Run PBSIM3 per-genome with abundance-weighted depth
            skipped = 0
            genome_bam_files = []
            for i, (seq, abundance) in enumerate(zip(sequences, abundances)):
                genome_depth = per_genome_depths[seq.id]
                if genome_depth < min_depth_threshold:
                    skipped += 1
                    continue

                # Write single-genome FASTA
                genome_fasta = Path(f"{output_prefix}_genome_{i:04d}.fasta")
                SeqIO.write([seq], str(genome_fasta), "fasta")
                genome_prefix = f"{output_prefix}_genome_{i:04d}_clr"

                pbsim_cmd = [
                    'pbsim',
                    '--strategy', 'wgs',
                    '--method', 'qshmm',
                    '--qshmm', qshmm_path,
                    '--depth', str(genome_depth),
                    '--genome', str(genome_fasta),
                    '--pass-num', str(platform_config.passes),
                    '--length-mean', str(platform_config.read_length_mean),
                    '--length-sd', str(platform_config.read_length_sd),
                    '--accuracy-mean', str(1.0 - platform_config.clr_error_rate),
                    '--prefix', genome_prefix
                ]

                if args.seed:
                    pbsim_cmd.extend(['--seed', str(args.seed + i)])

                try:
                    result = subprocess.run(pbsim_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    if result.returncode != 0 and result.returncode != -13:
                        logger.warning(f"  PBSIM3 failed for genome {seq.id} (exit {result.returncode})")
                        genome_fasta.unlink(missing_ok=True)
                        continue
                except FileNotFoundError:
                    logger.error("PBSIM3 (pbsim) not found in PATH")
                    logger.error("Install with: conda install -c bioconda pbsim3")
                    sys.exit(1)

                # Collect per-genome BAM output (PBSIM3 outputs {prefix}_0001.bam)
                import glob as _glob_clr
                genome_bams = sorted(_glob_clr.glob(f"{genome_prefix}_*.bam"))
                genome_bam_files.extend(genome_bams)

                # Also check for SAM output (older PBSIM3 versions)
                genome_sam = Path(f"{genome_prefix}.sam")
                if genome_sam.exists() and not genome_bams:
                    # Convert SAM to BAM
                    genome_bam_path = Path(f"{genome_prefix}.bam")
                    sam_conv = subprocess.run(
                        ['samtools', 'view', '-b', '-o', str(genome_bam_path), str(genome_sam)],
                        capture_output=True, text=True
                    )
                    if sam_conv.returncode == 0:
                        genome_bam_files.append(str(genome_bam_path))
                    genome_sam.unlink()

                # Clean up per-genome temp files
                genome_fasta.unlink(missing_ok=True)
                for f in _glob_clr.glob(f"{genome_prefix}*.ref"):
                    Path(f).unlink(missing_ok=True)
                for f in _glob_clr.glob(f"{genome_prefix}*.maf"):
                    Path(f).unlink(missing_ok=True)
                for f in _glob_clr.glob(f"{genome_prefix}*.maf.gz"):
                    Path(f).unlink(missing_ok=True)

                if (i + 1) % 50 == 0:
                    logger.info(f"    Processed {i + 1}/{len(sequences)} genomes")

            if skipped > 0:
                logger.info(f"  Skipped {skipped} genomes with depth < {min_depth_threshold}x")
            logger.info(f"  PBSIM3 CLR generation complete ({len(genome_bam_files)} BAM files)")

            # Step 2: Merge BAMs and run ccs
            logger.info("  Step 2/2: Generating HiFi consensus with ccs...")
            clr_bam = Path(output_prefix + '_clr.bam')
            hifi_fastq = Path(output_prefix + '_hifi.fastq.gz')

            if not genome_bam_files:
                logger.error("No CLR BAM files generated by PBSIM3")
                sys.exit(1)

            # Merge per-genome BAMs into single BAM
            try:
                if len(genome_bam_files) == 1:
                    import shutil
                    shutil.move(str(genome_bam_files[0]), str(clr_bam))
                else:
                    merge_cmd = ['samtools', 'merge', '-f', str(clr_bam)] + [str(b) for b in genome_bam_files]
                    subprocess.run(merge_cmd, capture_output=True, text=True, check=True)
                logger.info(f"  Merged {len(genome_bam_files)} BAM files")
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools merge failed: {e.stderr}")
                raise
            except FileNotFoundError:
                logger.error("samtools not found in PATH")
                logger.error("Install with: conda install -c bioconda samtools")
                sys.exit(1)

            # Clean up per-genome BAMs
            for bam_file in genome_bam_files:
                Path(bam_file).unlink(missing_ok=True)

            # Run ccs
            ccs_cmd = [
                'ccs',
                str(clr_bam),
                str(hifi_fastq),
                '--min-passes', str(platform_config.min_passes),
                '--min-rq', '0.99',
                '--log-level', 'INFO'
            ]

            try:
                subprocess.run(ccs_cmd, capture_output=True, text=True, check=True)
                logger.info("  PacBio ccs complete")
            except subprocess.CalledProcessError as e:
                logger.error(f"PacBio ccs failed: {e.stderr}")
                raise
            except FileNotFoundError:
                logger.error("PacBio ccs not found in PATH")
                logger.error("Install with: conda install -c bioconda pbccs")
                import platform as _platform
                if _platform.machine() != 'x86_64':
                    logger.error(f"Note: pbccs is only available for Linux x86-64 (your system: {_platform.machine()})")
                    logger.error("PacBio HiFi generation must be run on a Linux x86-64 machine or cluster")
                sys.exit(1)

            reads_path = hifi_fastq

            # Clean up intermediate files
            clr_bam.unlink(missing_ok=True)

        else:  # nanopore
            platform_config = NanoporeConfig(
                chemistry=args.ont_chemistry,
                read_length_mean=args.ont_read_length
            )
            logger.info(f"  Nanopore config: {platform_config.chemistry} chemistry, "
                       f"{platform_config.read_length_mean}bp mean read length")
            logger.info(f"  Depth: {args.depth}x")

            # Generate Nanopore reads with PBSIM3 (per-genome for correct abundances)
            logger.info("  Generating Nanopore reads with PBSIM3...")
            nanopore_fastq = Path(output_prefix + '.fastq')

            # Map chemistry version to PBSIM3 error model file
            errhmm_map = {
                'R10.4': 'ERRHMM-ONT-HQ.model',  # R10.4 = high quality
                'R9.4': 'ERRHMM-ONT.model',       # R9.4 = standard
            }
            errhmm_name = errhmm_map.get(platform_config.chemistry, 'ERRHMM-ONT-HQ.model')

            # Find the model file: check conda data dir, then fall back to just the name
            import shutil
            pbsim_path = shutil.which('pbsim')
            errhmm_path = errhmm_name
            if pbsim_path:
                conda_data = Path(pbsim_path).parent.parent / 'data' / errhmm_name
                if conda_data.exists():
                    errhmm_path = str(conda_data)

            # Run PBSIM3 per-genome with abundance-weighted depth
            skipped = 0
            total_reads = 0
            with open(nanopore_fastq, 'w') as fq_out:
                for i, (seq, abundance) in enumerate(zip(sequences, abundances)):
                    genome_depth = per_genome_depths[seq.id]
                    if genome_depth < min_depth_threshold:
                        skipped += 1
                        continue

                    # Write single-genome FASTA
                    genome_fasta = Path(f"{output_prefix}_genome_{i:04d}.fasta")
                    SeqIO.write([seq], str(genome_fasta), "fasta")
                    genome_prefix = f"{output_prefix}_genome_{i:04d}"

                    pbsim_cmd = [
                        'pbsim',
                        '--strategy', 'wgs',
                        '--method', 'errhmm',
                        '--errhmm', errhmm_path,
                        '--depth', str(genome_depth),
                        '--genome', str(genome_fasta),
                        '--length-mean', str(platform_config.read_length_mean),
                        '--length-sd', str(platform_config.read_length_sd),
                        '--accuracy-mean', str(1.0 - platform_config.error_rate),
                        '--hp-del-bias', str(platform_config.hp_del_bias),
                        '--prefix', genome_prefix
                    ]

                    if args.seed:
                        pbsim_cmd.extend(['--seed', str(args.seed + i)])

                    try:
                        result = subprocess.run(pbsim_cmd, capture_output=True, text=True)
                        if result.returncode != 0 and result.returncode != -13:
                            logger.warning(f"  PBSIM3 failed for genome {seq.id} (exit {result.returncode})")
                            continue
                    except FileNotFoundError:
                        logger.error("PBSIM3 (pbsim) not found in PATH")
                        logger.error("Install with: conda install -c bioconda pbsim3")
                        sys.exit(1)

                    # Collect FASTQ output for this genome
                    # PBSIM3 outputs .fq.gz (gzipped) files
                    import glob as _glob_nano
                    import gzip
                    nano_fastqs = _glob_nano.glob(f"{genome_prefix}*.fq.gz")
                    if not nano_fastqs:
                        logger.debug(f"  No .fq.gz for genome {i} ({seq.id}), depth={genome_depth:.3f}, rc={result.returncode}")
                    if not nano_fastqs:
                        # Fall back to uncompressed .fastq
                        nano_fastqs = _glob_nano.glob(f"{genome_prefix}*.fastq")
                    for fq in sorted(nano_fastqs):
                        opener = gzip.open if fq.endswith('.gz') else open
                        with opener(fq, 'rt') as inf:
                            for line in inf:
                                fq_out.write(line)
                                if line.startswith('@'):
                                    total_reads += 1
                        Path(fq).unlink()

                    # Clean up per-genome temp files
                    genome_fasta.unlink(missing_ok=True)
                    for f in _glob_nano.glob(f"{genome_prefix}*.ref"):
                        Path(f).unlink(missing_ok=True)
                    for f in _glob_nano.glob(f"{genome_prefix}*.maf*"):
                        Path(f).unlink(missing_ok=True)

                    if (i + 1) % 50 == 0:
                        logger.info(f"    Processed {i + 1}/{len(sequences)} genomes ({total_reads} reads)")

            if skipped > 0:
                logger.info(f"  Skipped {skipped} genomes with depth < {min_depth_threshold}x")
            logger.info(f"  PBSIM3 Nanopore generation complete ({total_reads} reads)")

            reads_path = nanopore_fastq

        # Create ground truth file
        ground_truth_path = generator.fastq_dir / f"{collection_name}_ground_truth.tsv"
        ground_truth_data = []
        for i, (seq, abundance) in enumerate(zip(sequences, abundances)):
            seq_type = 'viral' if i < enrichment_stats['n_viral_genomes'] else 'contaminant'
            ground_truth_data.append({
                'genome_id': seq.id,
                'genome_type': seq_type,
                'length': len(seq.seq),
                'relative_abundance': abundance,
                'platform': args.platform,
                'read_type': 'long'
            })

        import pandas as pd
        pd.DataFrame(ground_truth_data).to_csv(ground_truth_path, sep='\t', index=False)

        logger.info(f"Long-read generation complete!")
        logger.info(f"  Reads: {reads_path}")
        logger.info(f"  Ground truth: {ground_truth_path}")

        # Update output message for long reads
        r1_path = reads_path
        r2_path = None  # No R2 for long reads

    else:
        # ===============================================================
        # SHORT-READ SEQUENCING (Illumina: NovaSeq, MiSeq, HiSeq)
        # ===============================================================
        logger.info(f"Generating short-read FASTQ files ({args.platform})...")

        r1_path, r2_path = generator.generate_fastq_with_iss(
            fasta_path=fasta_path,
            sequences=sequences,
            abundances=abundances,
            coverage=args.coverage,
            read_length=args.read_length,
            insert_size=args.insert_size,
            platform=args.platform,
            n_reads=args.n_reads
        )

        # Post-process: add source type labels to read headers
        if contamination_profile is not None:
            genome_source_map = {}
            # Viral genomes (check if dark matter)
            n_viral = enrichment_stats.get('n_viral_genomes', 0)
            # Build set of dark matter genome IDs for labeling
            dark_matter_ids = {g['genome_id'] for g in genomes if g.get('is_dark_matter')}
            for seq in sequences[:n_viral]:
                if seq.id in dark_matter_ids:
                    genome_source_map[seq.id] = "dark_matter"
                else:
                    genome_source_map[seq.id] = "viral"
            # Contaminant genomes
            for contaminant in contamination_profile.contaminants:
                genome_source_map[contaminant.genome_id] = contaminant.contaminant_type.value

            logger.info("Labeling read headers with source types...")
            label_fastq_headers(r1_path, genome_source_map)
            label_fastq_headers(r2_path, genome_source_map)

        # Post-process: inject low-complexity artifact reads
        if args.low_complexity_rate > 0:
            from viroforge.simulators.low_complexity import add_low_complexity_reads

            # Parse entropy range if provided
            entropy_range = None
            if args.entropy_range:
                try:
                    lo, hi = args.entropy_range.split("-")
                    entropy_range = (float(lo), float(hi))
                except ValueError:
                    logger.error(f"Invalid --entropy-range format: {args.entropy_range} (expected MIN-MAX, e.g. 0.3-0.7)")
                    raise

            lc_manifest = generator.metadata_dir / f"{collection_name}_low_complexity_manifest.tsv"
            lc_stats = add_low_complexity_reads(
                r1_path=r1_path,
                r2_path=r2_path,
                rate=args.low_complexity_rate,
                entropy_range=entropy_range,
                random_seed=args.seed,
                in_place=True,
                manifest_path=lc_manifest,
            )
            logger.info(
                f"Low-complexity artifacts: {lc_stats['reads_modified']}/{lc_stats['reads_total']} "
                f"reads ({args.low_complexity_rate:.1%})"
            )

            # Append to metadata JSON
            import glob as _glob
            for mf in _glob.glob(str(generator.metadata_dir / "*_metadata.json")):
                with open(mf) as f:
                    metadata = json.load(f)
                metadata["low_complexity_stats"] = {
                    "rate": args.low_complexity_rate,
                    "reads_modified": lc_stats["reads_modified"],
                    "artifact_counts": lc_stats["artifact_counts"],
                    "manifest_file": str(lc_manifest),
                }
                with open(mf, "w") as f:
                    json.dump(metadata, f, indent=2)

        # Post-process: add adapter contamination
        # Prefer insert-size-driven mode when --mean-insert-size is set
        use_insert_size_mode = getattr(args, 'mean_insert_size', None) is not None
        use_adapter_rate_mode = args.adapter_rate > 0 and not use_insert_size_mode

        if use_insert_size_mode:
            from viroforge.simulators.adapters import add_insert_size_adapters

            manifest = generator.metadata_dir / f"{collection_name}_adapter_manifest.tsv"
            adapter_stats = add_insert_size_adapters(
                r1_path=r1_path,
                r2_path=r2_path,
                adapter_type=args.adapter_type,
                mean_insert=args.mean_insert_size,
                std_insert=getattr(args, 'insert_size_sd', 50),
                chimera_rate=getattr(args, 'chimera_rate', 0.0),
                random_seed=args.seed,
                in_place=True,
                manifest_path=manifest,
            )
            logger.info(
                f"Insert-size adapter contamination: "
                f"{adapter_stats['reads_readthrough']} read-through "
                f"({adapter_stats['emergent_adapter_rate']:.1%}), "
                f"{adapter_stats['reads_chimera']} chimeras"
            )

            # Append adapter stats to metadata JSON
            import glob
            for mf in glob.glob(str(generator.metadata_dir / "*_metadata.json")):
                with open(mf) as f:
                    metadata = json.load(f)
                metadata["adapter_stats"] = {
                    "mode": "insert_size_driven",
                    "adapter_type": args.adapter_type,
                    "mean_insert_size": args.mean_insert_size,
                    "insert_size_sd": getattr(args, 'insert_size_sd', 50),
                    "chimera_rate": getattr(args, 'chimera_rate', 0.0),
                    "reads_total": adapter_stats["reads_total"],
                    "reads_readthrough": adapter_stats["reads_readthrough"],
                    "reads_chimera": adapter_stats["reads_chimera"],
                    "emergent_adapter_rate": adapter_stats["emergent_adapter_rate"],
                    "mean_adapter_length": adapter_stats["mean_adapter_length"],
                    "manifest_file": str(manifest),
                }
                with open(mf, "w") as f:
                    json.dump(metadata, f, indent=2)

        elif use_adapter_rate_mode:
            from viroforge.simulators.adapters import add_adapter_readthrough

            manifest = generator.metadata_dir / f"{collection_name}_adapter_manifest.tsv"
            adapter_stats = add_adapter_readthrough(
                r1_path=r1_path,
                r2_path=r2_path,
                output_r1=r1_path,
                output_r2=r2_path,
                adapter_rate=args.adapter_rate,
                adapter_type=args.adapter_type,
                random_seed=args.seed,
                in_place=True,
                manifest_path=manifest,
            )
            logger.info(
                f"Adapter contamination: {adapter_stats['reads_modified']}/{adapter_stats['reads_total']} "
                f"reads ({args.adapter_rate:.1%}), mean adapter length: "
                f"{adapter_stats['mean_adapter_length']:.1f} bp"
            )

            # Append adapter stats to metadata JSON
            import glob
            for mf in glob.glob(str(generator.metadata_dir / "*_metadata.json")):
                with open(mf) as f:
                    metadata = json.load(f)
                metadata["adapter_stats"] = {
                    "mode": "fixed_rate",
                    "adapter_rate": args.adapter_rate,
                    "adapter_type": args.adapter_type,
                    "reads_total": adapter_stats["reads_total"],
                    "reads_modified": adapter_stats["reads_modified"],
                    "mean_adapter_length": adapter_stats["mean_adapter_length"],
                    "manifest_file": str(manifest),
                }
                with open(mf, "w") as f:
                    json.dump(metadata, f, indent=2)

        # Post-process: inject PCR duplicates
        if args.duplicate_rate > 0:
            from viroforge.simulators.duplicates import add_pcr_duplicates

            dup_manifest = generator.metadata_dir / f"{collection_name}_duplicate_manifest.tsv"
            # MDA (phi29) produces a different duplicate profile than PCR
            # (power-law copies, stronger GC bias); map the amplification method.
            dup_amp_method = "mda" if args.amplification in ("mda", "mda-long") else "pcr"
            dup_stats = add_pcr_duplicates(
                r1_path=r1_path,
                r2_path=r2_path,
                duplicate_rate=args.duplicate_rate,
                max_copies=args.duplicate_max_copies,
                error_rate=args.duplicate_error_rate,
                amplification_method=dup_amp_method,
                random_seed=args.seed,
                in_place=True,
                manifest_path=dup_manifest,
            )
            logger.info(
                f"PCR duplicates: {dup_stats['copies_generated']} copies from "
                f"{dup_stats['templates']} templates "
                f"({dup_stats['duplicate_fraction']:.1%} duplicate fraction)"
            )

            # Append to metadata JSON
            import glob as _glob2
            for mf in _glob2.glob(str(generator.metadata_dir / "*_metadata.json")):
                with open(mf) as f:
                    metadata = json.load(f)
                # Record the values actually used, not the raw CLI args: MDA
                # auto-upgrades max_copies/error_rate/distribution inside
                # add_pcr_duplicates, so args would misreport them for MDA runs.
                metadata["duplicate_stats"] = {
                    "duplicate_rate": args.duplicate_rate,
                    "amplification_method": dup_stats["amplification_method"],
                    "max_copies": dup_stats["max_copies_used"],
                    "error_rate": dup_stats["error_rate_used"],
                    "copy_distribution": dup_stats["copy_distribution"],
                    "reads_original": dup_stats["reads_original"],
                    "reads_output": dup_stats["reads_output"],
                    "templates": dup_stats["templates"],
                    "copies_generated": dup_stats["copies_generated"],
                    "duplicate_fraction": dup_stats["duplicate_fraction"],
                    "copy_count_distribution": dup_stats["copy_count_distribution"],
                    "manifest_file": str(dup_manifest),
                }
                with open(mf, "w") as f:
                    json.dump(metadata, f, indent=2)

        # Post-process: inject MDA chimeras (only for MDA amplification)
        chimera_rate = args.mda_chimera_rate
        if chimera_rate == 0.0 and args.amplification in ('mda', 'mda-long'):
            chimera_rate = 0.15  # Default 15% for MDA if not set
        if chimera_rate > 0 and args.amplification in ('mda', 'mda-long'):
            from viroforge.simulators.duplicates import add_mda_chimeras

            chimera_manifest = generator.metadata_dir / f"{collection_name}_chimera_manifest.tsv"
            chimera_stats = add_mda_chimeras(
                r1_path=r1_path,
                r2_path=r2_path,
                chimera_rate=chimera_rate,
                random_seed=args.seed,
                in_place=True,
                manifest_path=chimera_manifest,
            )
            logger.info(
                f"MDA chimeras: {chimera_stats['chimeras_created']} chimeric reads "
                f"({chimera_stats['chimera_fraction']:.1%}), "
                f"{chimera_stats['inter_genomic']} inter-genomic, "
                f"{chimera_stats['intra_genomic']} intra-genomic"
            )

            import glob as _glob3
            for mf in _glob3.glob(str(generator.metadata_dir / "*_metadata.json")):
                with open(mf) as f:
                    metadata = json.load(f)
                metadata["chimera_stats"] = {
                    "chimera_rate": chimera_rate,
                    "chimeras_created": chimera_stats["chimeras_created"],
                    "chimera_fraction": chimera_stats["chimera_fraction"],
                    "inter_genomic": chimera_stats["inter_genomic"],
                    "intra_genomic": chimera_stats["intra_genomic"],
                    "manifest_file": str(chimera_manifest),
                }
                with open(mf, "w") as f:
                    json.dump(metadata, f, indent=2)

    # Format output message
    if is_long_read:
        output_msg = f"""
✓ Long-read FASTQ generation complete!
  Collection: {collection_meta['collection_name']}
  Viral genomes: {enrichment_stats['n_viral_genomes']}
  Contaminants: {enrichment_stats['n_contaminants']}
  Total sequences: {enrichment_stats['n_viral_genomes'] + enrichment_stats['n_contaminants']}
  Depth: {args.depth}x
  Platform: {args.platform}
  VLP protocol: {vlp_protocol if vlp_protocol else 'none (bulk metagenome)'}
  Contamination level: {args.contamination_level}
  Viral fraction: {enrichment_stats['viral_fraction']*100:.2f}%
  Contamination: {enrichment_stats['contamination_fraction']*100:.2f}%

  Output files:
    - Reads: {r1_path}
    - Ground truth: {ground_truth_path}
    - FASTA: {fasta_path}
    - Metadata: {generator.metadata_dir}
"""
    else:
        output_msg = f"""
✓ Short-read FASTQ generation complete!
  Collection: {collection_meta['collection_name']}
  Viral genomes: {enrichment_stats['n_viral_genomes']}
  Contaminants: {enrichment_stats['n_contaminants']}
  Total sequences: {enrichment_stats['n_viral_genomes'] + enrichment_stats['n_contaminants']}
  Coverage: {args.coverage}x
  Platform: {args.platform}
  VLP protocol: {vlp_protocol if vlp_protocol else 'none (bulk metagenome)'}
  Contamination level: {args.contamination_level}
  Viral fraction: {enrichment_stats['viral_fraction']*100:.2f}%
  Contamination: {enrichment_stats['contamination_fraction']*100:.2f}%

  Output files:
    - R1: {r1_path}
    - R2: {r2_path}
    - FASTA: {fasta_path}
    - Metadata: {generator.metadata_dir}
"""
    logger.info(output_msg)


if __name__ == '__main__':
    main()
