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
        rna_workflow_params: Optional[Dict] = None,
        read_type: str = "short",
        use_real_references: bool = True,
    ) -> Tuple[List[SeqRecord], List[float], Dict, Optional[ContaminationProfile]]:
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
            read_type: Read type ("short" or "long") for VLP size bias modeling

        Returns:
            Tuple of (sequence_records, abundances, enrichment_stats, contamination_profile)
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

            # Resync abundances with modified sequence list.
            # RNA workflow can fragment (add) and drop (remove) sequences.
            # Map each surviving sequence back to its parent genome's abundance.
            original_abundance_map = {}
            for genome in genomes:
                original_abundance_map[genome['genome_id']] = genome['relative_abundance']

            new_abundances = []
            for seq_record in viral_sequences:
                # Fragments have IDs like "GCF_xxx_frag0"; strip suffix to find parent
                parent_id = seq_record.id.split('_frag')[0]
                parent_abundance = original_abundance_map.get(parent_id, 0.0)
                new_abundances.append(parent_abundance)

            # Renormalize
            viral_abundances = np.array(new_abundances)
            total = viral_abundances.sum()
            if total > 0:
                viral_abundances = viral_abundances / total

            logger.info(
                f"RNA workflow complete: {len(viral_sequences)} sequences "
                f"(was {len(genomes)})"
            )
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
                    random_seed=self.random_seed,
                    use_real_references=use_real_references,
                )
            else:
                contam_profile = create_contamination_profile(
                    contamination_level,
                    random_seed=self.random_seed,
                    use_real_references=use_real_references,
                )

            # Combine viral + contamination
            sequences, abundances = self._combine_viral_and_contamination(
                viral_sequences,
                viral_abundances,
                contam_profile
            )

            # Compute fractions from combined (normalized) abundances
            total_abundance = sum(abundances)
            viral_total = sum(abundances[:len(viral_sequences)])
            contam_total = sum(abundances[len(viral_sequences):])

            stats = {
                'molecule_type': self.molecule_type,
                'vlp_protocol': 'none',
                'contamination_level': contamination_level,
                'viral_fraction': float(viral_total / total_abundance) if total_abundance > 0 else 0.0,
                'contamination_fraction': float(contam_total / total_abundance) if total_abundance > 0 else 0.0,
                'n_viral_genomes': len(viral_sequences),
                'n_contaminants': len(contam_profile),
                'rna_workflow': rna_workflow_stats if rna_workflow_stats else None
            }

            return sequences, abundances, stats, contam_profile

        # Apply VLP enrichment workflow
        logger.info(f"Applying VLP enrichment: {vlp_protocol}")
        sequences, abundances, stats, contam_profile = self._apply_vlp_enrichment(
            viral_sequences,
            viral_abundances,
            genomes,
            vlp_protocol,
            contamination_level,
            read_type,
            use_real_references=use_real_references,
        )

        return sequences, abundances, stats, contam_profile

    def _apply_vlp_enrichment(
        self,
        viral_sequences: List[SeqRecord],
        viral_abundances: np.ndarray,
        genomes: List[Dict],
        vlp_protocol: str,
        contamination_level: str,
        read_type: str = "short",
        use_real_references: bool = True,
    ) -> Tuple[List[SeqRecord], List[float], Dict, ContaminationProfile]:
        """
        Apply VLP enrichment protocol with size-based enrichment and contamination reduction.

        Args:
            viral_sequences: Viral genome SeqRecords
            viral_abundances: Viral relative abundances
            genomes: Genome metadata dictionaries
            vlp_protocol: VLP protocol name
            contamination_level: Contamination level
            read_type: Read type ("short" or "long") for size bias modeling

        Returns:
            Tuple of (combined_sequences, combined_abundances, enrichment_stats, contamination_profile)
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
        logger.info(f"Applying size-based viral enrichment (read_type={read_type})...")
        enriched_viral_abundances, viral_stats = vlp.apply_enrichment(
            genomes=genomes,
            abundances=viral_abundances,
            read_type=read_type
        )

        # Create contamination profile (RNA or DNA specific)
        logger.info(f"Creating contamination profile: {contamination_level} ({self.molecule_type.upper()})")
        if self.molecule_type == 'rna':
            contam_profile = create_rna_contamination_profile(
                contamination_level,
                ribo_depletion=True,  # Assume Ribo-Zero applied
                microbiome_rich=True,
                random_seed=self.random_seed,
                use_real_references=use_real_references,
            )
        else:
            contam_profile = create_contamination_profile(
                contamination_level,
                random_seed=self.random_seed,
                use_real_references=use_real_references,
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

        return sequences, abundances, stats, reduced_contam_profile

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

    def _categorize_coverage(self, coverage: float) -> str:
        """
        Categorize coverage into quality bins for benchmarking.

        Coverage categories match assembly quality standards:
        - complete: ≥20x (expect ≥95% genome recovery)
        - high_quality: ≥10x (expect ≥75% recovery)
        - partial: ≥5x (expect ≥50% recovery)
        - fragmented: ≥1x (expect <50% recovery)
        - missing: <1x (unlikely to recover)

        Args:
            coverage: Expected coverage depth

        Returns:
            Coverage category string
        """
        if coverage >= 20:
            return 'complete'
        elif coverage >= 10:
            return 'high_quality'
        elif coverage >= 5:
            return 'partial'
        elif coverage >= 1:
            return 'fragmented'
        else:
            return 'missing'

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
            sys.exit(1)

    def export_metadata(
        self,
        collection_meta: Dict,
        sequences: List[SeqRecord],
        abundances: List[float],
        config: Dict,
        enrichment_stats: Optional[Dict] = None,
        amplification_stats: Optional[Dict] = None,
        contamination_profile: Optional[ContaminationProfile] = None,
        enable_benchmarking: bool = True
    ):
        """
        Export complete ground truth metadata including all sequences (viral + contaminants).

        Args:
            collection_meta: Collection metadata dictionary
            sequences: List of all sequences (viral + contaminants)
            abundances: Relative abundances for all sequences
            config: Configuration parameters
            enrichment_stats: VLP enrichment statistics
            amplification_stats: Amplification bias statistics
            contamination_profile: ContaminationProfile object (for benchmarking metadata)
            enable_benchmarking: Enable benchmarking metadata export (default: True)
        """

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
            'metadata_version': '1.1' if enable_benchmarking else '1.0',
            'generation_info': {
                'timestamp': datetime.now().isoformat(),
                'viroforge_version': '0.11.0',
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

        # ===================================================================
        # BENCHMARKING METADATA (v1.1)
        # ===================================================================
        if enable_benchmarking:
            logger.info("Generating benchmarking metadata (contamination manifest + expected coverage)...")

            benchmarking = {}

            # 1. Contamination Manifest
            if contamination_profile is not None:
                contamination_manifest = []
                for contaminant in contamination_profile.contaminants:
                    contam_dict = contaminant.to_dict()
                    # Remove sequence (too large for JSON)
                    contam_dict.pop('sequence', None)
                    contamination_manifest.append(contam_dict)

                benchmarking['contamination_manifest'] = {
                    'profile_name': contamination_profile.name,
                    'total_contamination_pct': float(contamination_profile.get_total_abundance() * 100),
                    'n_contaminants': len(contamination_profile.contaminants),
                    'contaminants': contamination_manifest
                }
                logger.info(f"  Contamination manifest: {len(contamination_manifest)} contaminants")
            else:
                benchmarking['contamination_manifest'] = None
                logger.info("  Contamination manifest: None (no contamination profile)")

            # 2. Expected Coverage per Genome
            # Calculate expected coverage based on relative abundance and sequencing parameters
            expected_coverage_list = []
            total_genome_length = sum(len(seq.seq) for seq in sequences)

            # Extract coverage/depth from config
            coverage_param = config.get('coverage', config.get('depth', 10.0))
            n_reads = config.get('n_reads')
            read_length = config.get('read_length', 150)
            platform = config.get('platform', 'novaseq')
            is_long_read = platform in ['pacbio-hifi', 'nanopore']

            # Calculate total sequencing output (in bp)
            if n_reads is not None:
                # User specified read count
                if is_long_read:
                    # Long-read: single-end, variable length
                    # Approximate: use 80% of mean read length (accounting for length variation)
                    effective_read_length = read_length * 0.8
                    total_bp = n_reads * effective_read_length
                else:
                    # Short-read: paired-end
                    total_bp = n_reads * read_length * 2
            else:
                # Calculate from coverage parameter
                total_bp = total_genome_length * coverage_param

            for i, (seq, abundance) in enumerate(zip(sequences, abundances)):
                genome_length = len(seq.seq)
                seq_type = 'viral' if i < n_viral else 'contaminant'

                # Calculate expected reads for this genome
                genome_bp = total_bp * abundance

                # Expected coverage = (genome_bp) / genome_length
                expected_coverage = genome_bp / genome_length if genome_length > 0 else 0.0

                # Expected completeness (what fraction of genome we expect to recover)
                # Using Lander-Waterman formula: C = 1 - e^(-coverage)
                # For viruses, this is reasonable approximation
                expected_completeness = 1.0 - np.exp(-expected_coverage)

                expected_coverage_list.append({
                    'genome_id': seq.id,
                    'sequence_type': seq_type,
                    'length': genome_length,
                    'relative_abundance': float(abundance),
                    'expected_coverage': float(expected_coverage),
                    'expected_completeness': float(expected_completeness),
                    'coverage_category': self._categorize_coverage(expected_coverage)
                })

            benchmarking['expected_coverage'] = {
                'total_genome_length': total_genome_length,
                'total_sequencing_bp': int(total_bp),
                'coverage_parameter': float(coverage_param),
                'read_length': read_length,
                'platform': platform,
                'per_genome': expected_coverage_list
            }
            logger.info(f"  Expected coverage: calculated for {len(expected_coverage_list)} genomes")

            # 3. Notes on additional benchmarking metadata (for future phases)
            benchmarking['notes'] = {
                'read_manifest': 'Not implemented in Phase 13A (optional for advanced benchmarking)',
                'gene_annotations': 'Not implemented in Phase 13A (planned for Phase 13D)'
            }

            metadata['benchmarking'] = benchmarking
            logger.info("Benchmarking metadata generation complete")

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
        choices=['novaseq', 'miseq', 'hiseq', 'pacbio-hifi', 'nanopore'],
        default='novaseq',
        help='Sequencing platform (default: novaseq). ' +
             'Short-read: novaseq, miseq, hiseq. ' +
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
        default=0.0,
        help='Fraction of reads with adapter read-through (0.0-1.0, default: 0.0). '
             'Set to 0.05 for typical adapter contamination.'
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
        default=0.0,
        help='Fraction of reads replaced with low-complexity artifacts (0.0-1.0, '
             'default: 0.0). Set to 0.01 for typical artifact rate. '
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
        default=0.0,
        help='Fraction of reads that become PCR duplicate templates (0.0-1.0, '
             'default: 0.0). Set to 0.10-0.30 for typical PCR duplicate rates. '
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

    # Determine read type based on platform
    is_long_read = args.platform in ['pacbio-hifi', 'nanopore']
    read_type = "long" if is_long_read else "short"

    logger.info(f"Platform: {args.platform} (read_type={read_type})")

    # Prepare genomes with VLP enrichment
    vlp_protocol = None if args.no_vlp else args.vlp_protocol
    use_real = not args.no_real_contaminants
    sequences, abundances, enrichment_stats, contamination_profile = generator.prepare_genomes(
        genomes,
        vlp_protocol=vlp_protocol,
        contamination_level=args.contamination_level,
        rna_workflow_params=rna_workflow_params,
        read_type=read_type,
        use_real_references=use_real,
    )

    # Add ERV sequences (before ISS, so they get realistic error profiles)
    if args.erv_endogenous_rate > 0 or args.erv_exogenous_rate > 0:
        from viroforge.core.contamination import (
            add_erv_endogenous,
            add_erv_exogenous,
            ContaminationProfile as _CP,
        )

        # Create a temporary profile for ERV contaminants
        erv_profile = _CP(name="erv_injection")

        if args.erv_endogenous_rate > 0:
            herv_path = Path(args.herv_fasta) if args.herv_fasta else None
            add_erv_endogenous(
                erv_profile,
                abundance_pct=args.erv_endogenous_rate * 100,
                herv_fasta_path=herv_path,
                random_seed=args.seed,
            )

        if args.erv_exogenous_rate > 0:
            add_erv_exogenous(
                erv_profile,
                abundance_pct=args.erv_exogenous_rate * 100,
                db_path=Path(args.database),
                virus_names=args.erv_exogenous_viruses,
                random_seed=args.seed,
            )

        if erv_profile.contaminants:
            # Add ERV sequences to the combined sequence/abundance lists
            for contaminant in erv_profile.contaminants:
                erv_record = SeqRecord(
                    Seq(str(contaminant.sequence)),
                    id=contaminant.genome_id,
                    description=contaminant.description,
                )
                sequences.append(erv_record)
                abundances.append(contaminant.abundance)

                # Also add to the main contamination profile for metadata
                if contamination_profile is not None:
                    contamination_profile.add_contaminant(contaminant)

            # Renormalize abundances
            total = sum(abundances)
            abundances = [a / total for a in abundances]

            logger.info(
                f"Added {len(erv_profile.contaminants)} ERV sequences "
                f"({sum(1 for c in erv_profile.contaminants if c.contaminant_type.value == 'erv_endogenous')} endogenous, "
                f"{sum(1 for c in erv_profile.contaminants if c.contaminant_type.value == 'erv_exogenous')} exogenous)"
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
        amplification_stats,
        contamination_profile,
        args.enable_benchmarking
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

        # Calculate average depth for PBSIM3
        total_length = sum(len(seq.seq) for seq in sequences)
        weighted_depth = sum(args.depth * abundance for abundance in abundances)

        output_prefix = str(generator.fastq_dir / collection_name)

        # Determine platform and run simulation
        if args.platform == 'pacbio-hifi':
            platform_config = PacBioHiFiConfig(
                passes=args.pacbio_passes,
                read_length_mean=args.pacbio_read_length
            )
            logger.info(f"  PacBio HiFi config: {platform_config.passes} passes, "
                       f"{platform_config.read_length_mean}bp mean read length")
            logger.info(f"  Depth: {args.depth}x")

            # Step 1: Generate CLR with PBSIM3
            logger.info("  Step 1/2: Generating CLR reads with PBSIM3...")
            clr_sam = Path(output_prefix + '_clr.sam')

            pbsim_cmd = [
                'pbsim',
                '--strategy', 'wgs',
                '--method', 'qshmm',
                '--qshmm', platform_config.accuracy_model,
                '--depth', str(weighted_depth),
                '--genome', str(fasta_path),
                '--pass-num', str(platform_config.passes),
                '--length-mean', str(platform_config.read_length_mean),
                '--length-sd', str(platform_config.read_length_sd),
                '--accuracy-mean', str(1.0 - platform_config.clr_error_rate),
                '--prefix', output_prefix + '_clr'
            ]

            if args.seed:
                pbsim_cmd.extend(['--seed', str(args.seed)])

            try:
                result = subprocess.run(pbsim_cmd, capture_output=True, text=True, check=True)
                logger.info("  PBSIM3 CLR generation complete")
            except subprocess.CalledProcessError as e:
                logger.error(f"PBSIM3 failed: {e.stderr}")
                raise
            except FileNotFoundError:
                logger.error("PBSIM3 (pbsim) not found in PATH")
                logger.error("Install with: conda install -c bioconda pbsim3")
                sys.exit(1)

            # Step 2: Convert SAM to BAM and run ccs
            logger.info("  Step 2/2: Generating HiFi consensus with ccs...")
            clr_bam = Path(output_prefix + '_clr.bam')
            hifi_fastq = Path(output_prefix + '_hifi.fastq.gz')

            # SAM → BAM
            sam_to_bam_cmd = ['samtools', 'view', '-b', '-o', str(clr_bam), str(clr_sam)]
            try:
                subprocess.run(sam_to_bam_cmd, capture_output=True, text=True, check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools failed: {e.stderr}")
                raise
            except FileNotFoundError:
                logger.error("samtools not found in PATH")
                logger.error("Install with: conda install -c bioconda samtools")
                sys.exit(1)

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
                sys.exit(1)

            reads_path = hifi_fastq

            # Clean up intermediate files
            clr_sam.unlink(missing_ok=True)
            clr_bam.unlink(missing_ok=True)

        else:  # nanopore
            platform_config = NanoporeConfig(
                chemistry=args.ont_chemistry,
                read_length_mean=args.ont_read_length
            )
            logger.info(f"  Nanopore config: {platform_config.chemistry} chemistry, "
                       f"{platform_config.read_length_mean}bp mean read length")
            logger.info(f"  Depth: {args.depth}x")

            # Generate Nanopore reads with PBSIM3
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

            pbsim_cmd = [
                'pbsim',
                '--strategy', 'wgs',
                '--method', 'errhmm',
                '--errhmm', errhmm_path,
                '--depth', str(weighted_depth),
                '--genome', str(fasta_path),
                '--length-mean', str(platform_config.read_length_mean),
                '--length-sd', str(platform_config.read_length_sd),
                '--accuracy-mean', str(1.0 - platform_config.error_rate),
                '--hp-del-bias', str(platform_config.hp_del_bias),
                '--prefix', output_prefix
            ]

            if args.seed:
                pbsim_cmd.extend(['--seed', str(args.seed)])

            try:
                result = subprocess.run(pbsim_cmd, capture_output=True, text=True, check=True)
                logger.info("  PBSIM3 Nanopore generation complete")
            except subprocess.CalledProcessError as e:
                logger.error(f"PBSIM3 failed: {e.stderr}")
                raise
            except FileNotFoundError:
                logger.error("PBSIM3 (pbsim) not found in PATH")
                logger.error("Install with: conda install -c bioconda pbsim3")
                sys.exit(1)

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
            # Viral genomes
            n_viral = enrichment_stats.get('n_viral_genomes', 0)
            for seq in sequences[:n_viral]:
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
            dup_stats = add_pcr_duplicates(
                r1_path=r1_path,
                r2_path=r2_path,
                duplicate_rate=args.duplicate_rate,
                max_copies=args.duplicate_max_copies,
                error_rate=args.duplicate_error_rate,
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
                metadata["duplicate_stats"] = {
                    "duplicate_rate": args.duplicate_rate,
                    "max_copies": args.duplicate_max_copies,
                    "error_rate": args.duplicate_error_rate,
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
