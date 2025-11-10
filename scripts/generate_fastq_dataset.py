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


class FASTQGenerator:
    """Generate FASTQ files from collection genomes."""

    def __init__(
        self,
        output_dir: Path,
        collection_name: str,
        molecule_type: str = 'dna',
        random_seed: int = 42
    ):
        self.output_dir = Path(output_dir)
        self.collection_name = collection_name
        self.molecule_type = molecule_type.lower()
        self.random_seed = random_seed

        # Validate molecule type
        if self.molecule_type not in ['dna', 'rna']:
            raise ValueError(f"molecule_type must be 'dna' or 'rna', got: {molecule_type}")

        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fasta_dir = self.output_dir / 'fasta'
        self.fastq_dir = self.output_dir / 'fastq'
        self.metadata_dir = self.output_dir / 'metadata'

        for dir in [self.fasta_dir, self.fastq_dir, self.metadata_dir]:
            dir.mkdir(exist_ok=True)

        np.random.seed(random_seed)

        logger.info(f"Initialized FASTQGenerator for {molecule_type.upper()} virome")

    def prepare_genomes(
        self,
        genomes: List[Dict],
        vlp_protocol: Optional[str] = 'tangential_flow',
        contamination_level: str = 'realistic',
        rna_workflow_params: Optional[Dict] = None
    ) -> Tuple[List[SeqRecord], List[float], Dict]:
        """
        Prepare genome sequences with VLP enrichment and contamination.

        For RNA viromes, also applies:
        - Reverse transcription (RNA → cDNA)
        - rRNA depletion (Ribo-Zero/RiboMinus)
        - RNA degradation modeling

        Args:
            genomes: List of genome dictionaries
            vlp_protocol: VLP protocol name (tangential_flow, syringe, ultracentrifugation,
                         norgen, none) or None for no VLP
            contamination_level: Contamination level (clean, realistic, heavy)
            rna_workflow_params: Optional parameters for RNA workflow (primer_type, ribo_method, etc.)

        Returns:
            Tuple of (sequence_records, abundances, enrichment_stats)
        """
        logger.info(f"Preparing {len(genomes)} viral genomes...")

        # Convert genomes to viral sequences
        viral_sequences = []
        viral_abundances = []

        for genome in genomes:
            record = SeqRecord(
                Seq(genome['sequence']),
                id=genome['genome_id'],
                description=f"{genome['genome_name']} | {genome.get('family', 'Unknown')}"
            )
            viral_sequences.append(record)
            viral_abundances.append(genome['relative_abundance'])

        viral_abundances = np.array(viral_abundances)

        # ===================================================================
        # RNA VIROME WORKFLOW (if molecule_type='rna')
        # ===================================================================
        rna_workflow_stats = None
        if self.molecule_type == 'rna':
            logger.info("=" * 80)
            logger.info("APPLYING RNA VIROME WORKFLOW")
            logger.info("=" * 80)

            # Set default RNA workflow parameters
            if rna_workflow_params is None:
                rna_workflow_params = {}

            primer_type = rna_workflow_params.get('primer_type', PrimerType.RANDOM_HEXAMER)
            ribo_method = rna_workflow_params.get('ribo_method', RiboDepleteMethod.RIBO_ZERO)
            degradation_rate = rna_workflow_params.get('degradation_rate', 0.20)

            # Create virus type mapping for RNA workflow
            virus_types = {}
            for genome in genomes:
                taxonomy = {
                    'family': genome.get('family', 'Unknown'),
                    'genus': genome.get('genus', 'Unknown'),
                    'species': genome.get('species', 'Unknown'),
                    'genome_name': genome.get('genome_name', '')
                }
                virus_types[genome['genome_id']] = infer_virus_type_from_taxonomy(taxonomy)

            # Create RNA virome workflow
            rna_workflow = RNAViromeWorkflow(
                reverse_transcription=ReverseTranscription(
                    primer_type=primer_type,
                    random_seed=self.random_seed
                ),
                ribo_depletion=RiboDepletion(
                    method=ribo_method,
                    random_seed=self.random_seed
                ),
                rna_degradation=RNADegradation(
                    degradation_rate=degradation_rate,
                    random_seed=self.random_seed
                ),
                random_seed=self.random_seed
            )

            # Apply RNA workflow (degradation, RT, rRNA depletion)
            viral_sequences, rna_workflow_stats = rna_workflow.apply(
                sequences=viral_sequences,
                virus_types=virus_types,
                rrna_abundance_before=0.90  # 90% rRNA before Ribo-Zero
            )

            logger.info("RNA workflow complete")
            logger.info("=" * 80)

        # If no VLP protocol specified, return viral genomes only
        if vlp_protocol is None or vlp_protocol == 'none':
            logger.info("No VLP enrichment applied (bulk metagenome mode)")

            # Create contamination profile (RNA or DNA specific)
            if self.molecule_type == 'rna':
                contam_profile = create_rna_contamination_profile(
                    contamination_level,
                    ribo_depletion=True,  # Assume Ribo-Zero applied
                    microbiome_rich=True,
                    random_seed=self.random_seed
                )
            else:
                contam_profile = create_contamination_profile(
                    contamination_level,
                    random_seed=self.random_seed
                )

            # Combine viral + contamination
            sequences, abundances = self._combine_viral_and_contamination(
                viral_sequences,
                viral_abundances,
                contam_profile
            )

            stats = {
                'molecule_type': self.molecule_type,
                'vlp_protocol': 'none',
                'contamination_level': contamination_level,
                'viral_fraction': float(viral_abundances.sum()),
                'contamination_fraction': float(contam_profile.get_total_abundance()),
                'n_viral_genomes': len(viral_sequences),
                'n_contaminants': len(contam_profile),
                'rna_workflow': rna_workflow_stats if rna_workflow_stats else None
            }

            return sequences, abundances, stats

        # Apply VLP enrichment workflow
        logger.info(f"Applying VLP enrichment: {vlp_protocol}")
        sequences, abundances, stats = self._apply_vlp_enrichment(
            viral_sequences,
            viral_abundances,
            genomes,
            vlp_protocol,
            contamination_level
        )

        return sequences, abundances, stats

    def _apply_vlp_enrichment(
        self,
        viral_sequences: List[SeqRecord],
        viral_abundances: np.ndarray,
        genomes: List[Dict],
        vlp_protocol: str,
        contamination_level: str
    ) -> Tuple[List[SeqRecord], List[float], Dict]:
        """
        Apply VLP enrichment protocol with size-based enrichment and contamination reduction.

        Args:
            viral_sequences: Viral genome SeqRecords
            viral_abundances: Viral relative abundances
            genomes: Genome metadata dictionaries
            vlp_protocol: VLP protocol name
            contamination_level: Contamination level

        Returns:
            Tuple of (combined_sequences, combined_abundances, enrichment_stats)
        """
        # Map protocol name to VLPProtocol
        protocol_map = {
            'tangential_flow': VLPProtocol.tangential_flow_standard(),
            'syringe': VLPProtocol.syringe_filter_standard(),
            'ultracentrifugation': VLPProtocol.ultracentrifugation(),
            'norgen': VLPProtocol.norgen_kit()
        }

        if vlp_protocol not in protocol_map:
            raise ValueError(
                f"Unknown VLP protocol: {vlp_protocol}. "
                f"Choose from: {', '.join(protocol_map.keys())}"
            )

        protocol_config = protocol_map[vlp_protocol]

        # Initialize VLP enrichment
        vlp = VLPEnrichment(
            protocol=protocol_config,
            random_seed=self.random_seed
        )

        # Apply size-based enrichment to viral genomes
        logger.info("Applying size-based viral enrichment...")
        enriched_viral_abundances, viral_stats = vlp.apply_enrichment(
            genomes=genomes,
            abundances=viral_abundances
        )

        # Create contamination profile (RNA or DNA specific)
        logger.info(f"Creating contamination profile: {contamination_level} ({self.molecule_type.upper()})")
        if self.molecule_type == 'rna':
            contam_profile = create_rna_contamination_profile(
                contamination_level,
                ribo_depletion=True,  # Assume Ribo-Zero applied
                microbiome_rich=True,
                random_seed=self.random_seed
            )
        else:
            contam_profile = create_contamination_profile(
                contamination_level,
                random_seed=self.random_seed
            )

        # Apply contamination reduction
        logger.info("Applying VLP contamination reduction...")
        reduced_contam_profile, contam_stats = vlp.apply_contamination_reduction(
            contam_profile
        )

        # Combine viral genomes + reduced contamination
        sequences, abundances = self._combine_viral_and_contamination(
            viral_sequences,
            enriched_viral_abundances,
            reduced_contam_profile
        )

        # Calculate FINAL viral and contamination fractions after normalization
        n_viral = len(viral_sequences)
        final_viral_fraction = sum(abundances[:n_viral])
        final_contam_fraction = sum(abundances[n_viral:])

        # Compile statistics
        stats = {
            'vlp_protocol': vlp_protocol,
            'contamination_level': contamination_level,
            'viral_enrichment': viral_stats,
            'contamination_reduction': contam_stats,
            'viral_fraction': float(final_viral_fraction),
            'contamination_fraction': float(final_contam_fraction),
            'n_viral_genomes': len(viral_sequences),
            'n_contaminants': len(reduced_contam_profile)
        }

        logger.info(f"VLP enrichment complete:")
        logger.info(f"  Final viral fraction: {stats['viral_fraction']*100:.2f}%")
        logger.info(f"  Final contamination: {stats['contamination_fraction']*100:.2f}%")
        logger.info(f"  Total genomes: {len(sequences)}")

        return sequences, abundances, stats

    def _combine_viral_and_contamination(
        self,
        viral_sequences: List[SeqRecord],
        viral_abundances: np.ndarray,
        contam_profile: ContaminationProfile
    ) -> Tuple[List[SeqRecord], List[float]]:
        """
        Combine viral genomes and contamination into single set of sequences.

        Args:
            viral_sequences: Viral genome SeqRecords
            viral_abundances: Viral relative abundances (sum may be < 1.0)
            contam_profile: Contamination profile with contaminants

        Returns:
            Tuple of (combined_sequences, combined_abundances)
        """
        sequences = list(viral_sequences)
        abundances = list(viral_abundances)

        # Add contaminants
        for contaminant in contam_profile.contaminants:
            record = SeqRecord(
                contaminant.sequence,
                id=contaminant.genome_id,
                description=f"{contaminant.organism} | {contaminant.contaminant_type.value}"
            )
            sequences.append(record)
            abundances.append(contaminant.abundance)

        # Normalize to sum to 1.0
        abundances = np.array(abundances)
        total_abundance = abundances.sum()

        # Validate total abundance
        if total_abundance == 0 or not np.isfinite(total_abundance):
            raise ValueError(
                "Total abundance is zero or invalid. This indicates a problem with "
                "VLP enrichment or contamination profile generation."
            )

        abundances = abundances / total_abundance

        # Validate final abundances
        if not np.allclose(abundances.sum(), 1.0, atol=1e-6):
            logger.warning(
                f"Abundances do not sum to 1.0 (sum={abundances.sum():.8f}). "
                "Renormalizing..."
            )
            abundances = abundances / abundances.sum()

        if np.any(~np.isfinite(abundances)):
            raise ValueError("Abundances contain NaN or infinite values")

        logger.info(f"Combined {len(viral_sequences)} viral + {len(contam_profile)} contaminants")
        logger.info(f"Final abundance range: {abundances.min():.6f} - {abundances.max():.6f}")

        return sequences, abundances.tolist()

    def apply_amplification(
        self,
        sequences: List[SeqRecord],
        abundances: List[float],
        amplification_method: str
    ) -> Tuple[List[float], Dict]:
        """
        Apply library preparation amplification bias to abundances.

        Args:
            sequences: List of SeqRecord objects
            abundances: List of relative abundances
            amplification_method: Method name (none, rdab, rdab-30, mda, mda-long, linker)

        Returns:
            Tuple of (modified_abundances, amplification_stats)
        """
        if amplification_method == 'none':
            logger.info("No amplification bias applied")
            return abundances, {'method': 'none', 'bias_applied': False}

        logger.info(f"Applying amplification bias: {amplification_method}")

        # Map method name to amplification function
        method_map = {
            'rdab': rdab_40_cycles,
            'rdab-30': rdab_30_cycles,
            'mda': mda_standard,
            'mda-long': mda_overnight,
            'linker': linker_standard
        }

        if amplification_method not in method_map:
            raise ValueError(f"Unknown amplification method: {amplification_method}")

        # Get amplification method instance
        amp_method = method_map[amplification_method]()

        # Convert sequences and abundances to format amplification code expects
        # The amplification methods work on genome objects with length, gc_content, abundance
        # We need to create temporary genome-like objects

        # Calculate GC content for each sequence
        abundances_array = np.array(abundances)
        original_abundances = abundances_array.copy()

        # Create genome-like dicts for amplification calculation
        temp_genomes = []
        for seq, ab in zip(sequences, abundances_array):
            gc_content = (str(seq.seq).count('G') + str(seq.seq).count('C')) / len(seq.seq)
            temp_genomes.append({
                'genome_id': seq.id,
                'length': len(seq.seq),
                'gc_content': gc_content,
                'abundance': ab
            })

        # Apply bias based on method type
        if amplification_method in ['rdab', 'rdab-30']:
            # RdAB: length and GC bias
            from viroforge.amplification import RdABAmplification
            biased_abundances = self._apply_rdab_bias(temp_genomes, amp_method)
        elif amplification_method in ['mda', 'mda-long']:
            # MDA: extreme GC bias and stochasticity
            from viroforge.amplification import MDAAmplification
            biased_abundances = self._apply_mda_bias(temp_genomes, amp_method)
        elif amplification_method == 'linker':
            # Linker: minimal GC bias only
            from viroforge.amplification import LinkerAmplification
            biased_abundances = self._apply_linker_bias(temp_genomes, amp_method)
        else:
            biased_abundances = abundances_array

        # Renormalize to sum to 1.0
        biased_abundances = biased_abundances / biased_abundances.sum()

        # Calculate statistics
        fold_changes = biased_abundances / (original_abundances + 1e-10)
        stats = {
            'method': amplification_method,
            'bias_applied': True,
            'mean_fold_change': float(np.mean(fold_changes)),
            'max_fold_change': float(np.max(fold_changes)),
            'min_fold_change': float(np.min(fold_changes)),
            'abundance_correlation': float(np.corrcoef(original_abundances, biased_abundances)[0, 1])
        }

        logger.info(f"Amplification bias applied: mean fold change = {stats['mean_fold_change']:.2f}x")

        return biased_abundances.tolist(), stats

    def _apply_rdab_bias(self, genomes, amp_method):
        """Apply RdAB amplification bias (length + GC)."""
        abundances = np.array([g['abundance'] for g in genomes])

        for i, genome in enumerate(genomes):
            # Use the amplification method's own efficiency calculations
            length_eff = amp_method.calculate_length_efficiency(genome['length'])
            gc_eff = amp_method.calculate_gc_efficiency(genome['gc_content'])

            # Combined efficiency per cycle, raised to number of cycles
            cycle_efficiency = length_eff * gc_eff
            amplification_factor = cycle_efficiency ** amp_method.cycles

            abundances[i] *= amplification_factor

        return abundances

    def _apply_mda_bias(self, genomes, amp_method):
        """Apply MDA amplification bias (extreme GC + stochasticity)."""
        abundances = np.array([g['abundance'] for g in genomes])

        for i, genome in enumerate(genomes):
            # Use the amplification method's GC efficiency calculation
            gc_eff = amp_method.calculate_gc_efficiency(genome['gc_content'])

            # Add stochastic variation using the method's own function
            final_eff = amp_method.add_stochastic_variation(gc_eff)

            # Amplification scales with time
            time_scaling = amp_method.amplification_time / 4.0  # Normalize to 4h standard
            amplification_factor = final_eff ** time_scaling

            abundances[i] *= amplification_factor

        return abundances

    def _apply_linker_bias(self, genomes, amp_method):
        """Apply linker amplification bias (minimal GC only)."""
        abundances = np.array([g['abundance'] for g in genomes])

        for i, genome in enumerate(genomes):
            # Use the amplification method's GC efficiency calculation
            gc_eff = amp_method.calculate_gc_efficiency(genome['gc_content'])

            # Amplification factor scales with cycles (no length bias for linker)
            amplification_factor = gc_eff ** amp_method.cycles

            abundances[i] *= amplification_factor

        return abundances

    def write_fasta(
        self,
        sequences: List[SeqRecord],
        abundances: List[float]
    ) -> Path:
        """Write genomes to FASTA file with abundance annotations."""
        fasta_path = self.fasta_dir / f"{self.collection_name}.fasta"

        with open(fasta_path, 'w') as f:
            for record, abundance in zip(sequences, abundances):
                # Add abundance to description
                record.description = f"{record.description} | abundance={abundance:.8f}"
                SeqIO.write(record, f, 'fasta')

        logger.info(f"Wrote {len(sequences)} sequences to: {fasta_path}")
        return fasta_path

    def generate_fastq_with_iss(
        self,
        fasta_path: Path,
        sequences: List[SeqRecord],
        abundances: List[float],
        coverage: float = 10.0,
        read_length: int = 150,
        insert_size: int = 350,
        platform: str = 'novaseq',
        n_reads: Optional[int] = None
    ) -> Tuple[Path, Path]:
        """
        Generate FASTQ files using InSilicoSeq.

        Args:
            fasta_path: Path to input FASTA
            sequences: List of SeqRecord objects (for getting IDs)
            abundances: Relative abundances (must sum to 1.0)
            coverage: Mean coverage depth
            read_length: Read length (bp)
            insert_size: Insert size for paired-end (bp)
            platform: Sequencing platform (novaseq, miseq, hiseq)
            n_reads: Number of reads (overrides coverage if set)

        Returns:
            Tuple of (R1_path, R2_path)
        """
        import subprocess

        # Validate coverage
        if coverage <= 0:
            raise ValueError(f"Coverage must be positive, got: {coverage}")
        if coverage > 100:
            logger.warning(
                f"Very high coverage requested ({coverage}x). "
                f"This may take hours and produce large files (>10 GB)."
            )
            if coverage > 500:
                raise ValueError(
                    f"Coverage {coverage}x exceeds maximum allowed (500x). "
                    f"If you really need this, use --n-reads directly."
                )

        # Calculate total genome length for coverage calculation
        if n_reads is None:
            total_length = sum(len(record.seq) for record in sequences)
            # Calculate reads needed: (total_length * coverage) / (2 * read_length)
            # Factor of 2 because paired-end reads
            n_reads = int((total_length * coverage) / (2 * read_length))

            # Validate calculated reads
            if n_reads == 0:
                raise ValueError(
                    f"Calculated zero reads for {coverage}x coverage "
                    f"(total genome length: {total_length:,} bp). "
                    "Check your coverage and genome sizes."
                )

            if n_reads > 1_000_000_000:  # 1 billion reads = ~150 GB for 2x150bp
                raise ValueError(
                    f"Calculated read count ({n_reads:,}) exceeds reasonable limit "
                    f"(1 billion reads). Reduce coverage or use smaller collection."
                )

            logger.info(f"Calculated {n_reads:,} reads needed for {coverage}x coverage")
            logger.info(f"  Total genome length: {total_length:,} bp")
            logger.info(f"  Estimated output size: ~{(n_reads * read_length * 4 / 1e9):.1f} GB")
        else:
            # Validate user-provided n_reads
            if n_reads <= 0:
                raise ValueError(f"n_reads must be positive, got: {n_reads}")
            if n_reads > 1_000_000_000:
                logger.warning(
                    f"Very high read count requested ({n_reads:,}). "
                    f"Estimated output size: ~{(n_reads * read_length * 4 / 1e9):.1f} GB"
                )

        # Create abundance file for ISS using actual genome IDs
        abundance_file = self.metadata_dir / f"{self.collection_name}_abundances.txt"
        with open(abundance_file, 'w') as f:
            for record, abundance in zip(sequences, abundances):
                f.write(f"{record.id}\t{abundance:.8f}\n")

        # Prepare ISS command
        output_prefix = self.fastq_dir / self.collection_name

        # Map our platform names to ISS error models
        iss_models = {
            'novaseq': 'novaseq',
            'miseq': 'miseq',
            'hiseq': 'hiseq'
        }
        error_model = iss_models.get(platform, 'novaseq')

        cmd = [
            'iss', 'generate',
            '--genomes', str(fasta_path),
            '--abundance_file', str(abundance_file),
            '--model', error_model,
            '--output', str(output_prefix),
            '--n_reads', str(n_reads),
            '--mode', 'basic',  # Use basic mode for faster generation
            '--seed', str(self.random_seed)
        ]

        logger.info(f"Running InSilicoSeq...")
        logger.info(f"Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("InSilicoSeq completed successfully")

            # ISS creates files with _R1.fastq and _R2.fastq suffixes
            r1_path = Path(str(output_prefix) + '_R1.fastq')
            r2_path = Path(str(output_prefix) + '_R2.fastq')

            if not r1_path.exists() or not r2_path.exists():
                raise FileNotFoundError(f"ISS output files not found: {output_prefix}")

            # Validate file sizes and format
            for fastq_path in [r1_path, r2_path]:
                file_size = fastq_path.stat().st_size
                if file_size == 0:
                    raise ValueError(f"ISS output file is empty: {fastq_path}")

                # Quick read count check (first few lines)
                with open(fastq_path) as f:
                    try:
                        # Check FASTQ format (first record)
                        header = f.readline()
                        if not header.startswith('@'):
                            raise ValueError(
                                f"Invalid FASTQ format in {fastq_path}: "
                                f"Expected '@' header, got: {header[:50]}"
                            )
                    except Exception as e:
                        raise ValueError(
                            f"Failed to validate FASTQ format in {fastq_path}: {e}"
                        )

            logger.info(f"R1 file size: {r1_path.stat().st_size:,} bytes")
            logger.info(f"R2 file size: {r2_path.stat().st_size:,} bytes")

            return r1_path, r2_path

        except subprocess.CalledProcessError as e:
            logger.error(f"InSilicoSeq failed: {e.stderr}")
            raise
        except FileNotFoundError:
            logger.error("InSilicoSeq (iss) not found in PATH")
            logger.error("Install with: conda install -c bioconda insilicoseq")
            raise

    def export_metadata(
        self,
        collection_meta: Dict,
        sequences: List[SeqRecord],
        abundances: List[float],
        config: Dict,
        enrichment_stats: Optional[Dict] = None,
        amplification_stats: Optional[Dict] = None
    ):
        """Export complete ground truth metadata including all sequences (viral + contaminants)."""

        # Validate that sequences and abundances match
        if len(sequences) != len(abundances):
            raise ValueError(
                f"Sequence count ({len(sequences)}) doesn't match "
                f"abundance count ({len(abundances)})"
            )

        # Get counts from enrichment stats
        n_viral = enrichment_stats.get('n_viral_genomes', len(sequences)) if enrichment_stats else len(sequences)
        n_contaminants = enrichment_stats.get('n_contaminants', 0) if enrichment_stats else 0

        metadata = {
            'generation_info': {
                'timestamp': datetime.now().isoformat(),
                'viroforge_version': '0.4.0',
                'random_seed': self.random_seed
            },
            'collection': {
                'id': collection_meta['collection_id'],
                'name': collection_meta['collection_name'],
                'description': collection_meta.get('description', ''),
                'n_viral_genomes': n_viral,
                'n_contaminants': n_contaminants,
                'total_sequences': len(sequences)
            },
            'configuration': config,
            'enrichment_stats': enrichment_stats,
            'amplification_stats': amplification_stats,
            'sequences': []
        }

        # Add ALL sequence details (viral + contaminants)
        for i, (seq, abundance) in enumerate(zip(sequences, abundances)):
            seq_type = 'viral' if i < n_viral else 'contaminant'

            # Parse description for metadata
            # Format: "name | family | ..." or "name" for contaminants
            desc_parts = seq.description.split('|')
            seq_name = desc_parts[0].strip() if desc_parts else seq.id

            seq_info = {
                'genome_id': seq.id,
                'genome_name': seq_name,
                'sequence_type': seq_type,
                'length': len(seq.seq),
                'relative_abundance': abundance
            }

            # Add viral-specific taxonomy if available
            if seq_type == 'viral' and len(desc_parts) >= 2:
                # Try to parse taxonomy from description
                for part in desc_parts[1:]:
                    part = part.strip()
                    if 'family:' in part.lower():
                        seq_info['family'] = part.split(':')[-1].strip()
                    elif 'genus:' in part.lower():
                        seq_info['genus'] = part.split(':')[-1].strip()
                    elif 'species:' in part.lower():
                        seq_info['species'] = part.split(':')[-1].strip()
                    elif not part.startswith('abundance='):
                        # Assume it's family if no prefix
                        seq_info['family'] = part

            metadata['sequences'].append(seq_info)

        # Write metadata JSON
        metadata_path = self.metadata_dir / f"{self.collection_name}_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        logger.info(f"Exported metadata to: {metadata_path}")
        logger.info(f"  {n_viral} viral genomes + {n_contaminants} contaminants = {len(sequences)} total sequences")

        # Also write TSV for easy viewing
        import pandas as pd
        df = pd.DataFrame(metadata['sequences'])
        tsv_path = self.metadata_dir / f"{self.collection_name}_composition.tsv"
        df.to_csv(tsv_path, sep='\t', index=False)
        logger.info(f"Exported composition table to: {tsv_path}")


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
        choices=['novaseq', 'miseq', 'hiseq'],
        default='novaseq',
        help='Sequencing platform (default: novaseq)'
    )

    parser.add_argument(
        '--no-vlp',
        action='store_true',
        help='Skip VLP enrichment simulation (generate bulk metagenome)'
    )

    parser.add_argument(
        '--vlp-protocol',
        choices=['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen'],
        default='tangential_flow',
        help='VLP enrichment protocol (default: tangential_flow)'
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
        default='none',
        help='Library preparation amplification method (default: none). ' +
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

    # Prepare genomes with VLP enrichment
    vlp_protocol = None if args.no_vlp else args.vlp_protocol
    sequences, abundances, enrichment_stats = generator.prepare_genomes(
        genomes,
        vlp_protocol=vlp_protocol,
        contamination_level=args.contamination_level,
        rna_workflow_params=rna_workflow_params
    )

    # Apply amplification bias
    abundances, amplification_stats = generator.apply_amplification(
        sequences,
        abundances,
        args.amplification
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
        'amplification': args.amplification
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
        amplification_stats
    )

    if args.dry_run:
        logger.info("Dry run complete - FASTQ generation skipped")
        return

    # Generate FASTQs
    logger.info("Generating FASTQ files...")
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

    # Format output message
    output_msg = f"""
✓ FASTQ generation complete!
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
