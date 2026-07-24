"""
ViroForge FASTQ Generator Engine.

Core ``FASTQGenerator`` class responsible for preparing viral genome sequences,
applying VLP enrichment, contamination, amplification bias, and generating
FASTQ files via InSilicoSeq.

Extracted from ``scripts/generate_fastq_dataset.py`` so that it can be imported
as a library component (e.g. ``from viroforge.generator import FASTQGenerator``).
"""

import json
import logging
import sqlite3
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viroforge.enrichment.vlp import VLPEnrichment, VLPProtocol
from viroforge.core.contamination import (
    create_contamination_profile,
    create_rna_contamination_profile,
    ContaminationProfile,
)
from viroforge.amplification import (
    rdab_40_cycles,
    rdab_30_cycles,
    mda_standard,
    mda_overnight,
    linker_standard,
    no_amplification,
)
from viroforge.workflows.rna_virome import (
    RNAViromeWorkflow,
    ReverseTranscription,
    RiboDepletion,
    RNADegradation,
    RNAVirusType,
    PrimerType,
    RiboDepleteMethod,
    infer_virus_type_from_taxonomy,
)
from viroforge.core.collection import CollectionLoader, label_fastq_headers

logger = logging.getLogger(__name__)


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
        erv_kwargs: Optional[Dict] = None,
        collection_defaults: Optional[Dict] = None,
    ) -> Tuple[List[SeqRecord], List[float], Dict, Optional[ContaminationProfile]]:
        """
        Prepare genome sequences with VLP enrichment and contamination.

        For RNA viromes, also applies:
        - Reverse transcription (RNA -> cDNA)
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
            original_genome_count = len(genomes)
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

            # Rebuild genome metadata list to mirror post-RNA-workflow sequences.
            # VLP apply_enrichment requires a list parallel to viral_sequences
            # (uses 'length' and 'genome_type' per entry). RNA degradation can
            # fragment a genome into multiple shorter sequences, changing the count.
            genome_meta_map = {g['genome_id']: g for g in genomes}
            genomes = []
            for seq_record in viral_sequences:
                parent_id = seq_record.id.split('_frag')[0]
                parent_meta = genome_meta_map.get(parent_id, {})
                genomes.append({
                    **parent_meta,
                    'genome_id': seq_record.id,
                    'length': len(seq_record.seq),
                })

            logger.info(
                f"RNA workflow complete: {len(viral_sequences)} sequences "
                f"(was {len(genome_meta_map)})"
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
                    **(erv_kwargs or {}),
                )
            else:
                contam_profile = create_contamination_profile(
                    contamination_level,
                    random_seed=self.random_seed,
                    use_real_references=use_real_references,
                    collection_defaults=collection_defaults,
                    **(erv_kwargs or {}),
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
            erv_kwargs=erv_kwargs,
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
        erv_kwargs: Optional[Dict] = None,
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
            'tangential_flow_045': VLPProtocol.tangential_flow_045(),
            'tangential_flow_01': VLPProtocol.tangential_flow_01(),
            'syringe': VLPProtocol.syringe_filter_standard(),
            'syringe_045': VLPProtocol.syringe_filter_045(),
            'ultracentrifugation': VLPProtocol.ultracentrifugation(),
            'norgen': VLPProtocol.norgen_kit()
        }

        # Check for dynamically registered custom pore size configs
        custom_configs = getattr(self, '_custom_vlp_configs', {})
        if vlp_protocol in custom_configs:
            protocol_config = custom_configs[vlp_protocol]
        elif vlp_protocol not in protocol_map:
            raise ValueError(
                f"Unknown VLP protocol: {vlp_protocol}. "
                f"Choose from: {', '.join(protocol_map.keys())}"
            )
        else:
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
                **(erv_kwargs or {}),
            )
        else:
            contam_profile = create_contamination_profile(
                contamination_level,
                random_seed=self.random_seed,
                use_real_references=use_real_references,
                collection_defaults=collection_defaults,
                **(erv_kwargs or {}),
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

    def _export_taxonomy(self, db_path: str, genome_ids: List[str]) -> Dict[str, Dict]:
        """Look up ICTV lineage + NCBI taxid for the given viral genomes.

        Returns {genome_id: {ncbi_taxid, realm, kingdom, phylum, class, order,
        family, subfamily, genus, species, is_known}}. is_known is False when the
        family is Unknown (no ICTV match), which the taxonomy benchmark uses to
        stratify classifiable versus dark/novel content.
        """
        if not genome_ids:
            return {}
        out: Dict[str, Dict] = {}
        conn = sqlite3.connect(db_path)
        try:
            cur = conn.cursor()
            qmarks = ",".join("?" * len(genome_ids))
            cur.execute(
                f"""
                SELECT genome_id, ncbi_taxid, realm, kingdom, phylum, class,
                       order_name, family, subfamily, genus, species
                FROM taxonomy WHERE genome_id IN ({qmarks})
                """,
                genome_ids,
            )
            for row in cur.fetchall():
                gid, taxid, realm, kingdom, phylum, cls, order_name, family, subfamily, genus, species = row
                out[gid] = {
                    "ncbi_taxid": int(taxid) if taxid not in (None, "") else None,
                    "realm": realm, "kingdom": kingdom, "phylum": phylum,
                    "class": cls, "order": order_name, "family": family,
                    "subfamily": subfamily, "genus": genus, "species": species,
                    "is_known": family not in (None, "", "Unknown"),
                }
        finally:
            conn.close()
        return out

    def _categorize_coverage(self, coverage: float) -> str:
        """
        Categorize coverage into quality bins for benchmarking.

        Coverage categories match assembly quality standards:
        - complete: >=20x (expect >=95% genome recovery)
        - high_quality: >=10x (expect >=75% recovery)
        - partial: >=5x (expect >=50% recovery)
        - fragmented: >=1x (expect <50% recovery)
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
        # Enforce minimum abundance floor so every genome gets at least a few reads.
        # Without this, VLP enrichment + amplification bias can reduce rare genomes
        # to zero abundance, causing ISS to skip them entirely.
        min_abundance = 1e-6  # Floor: ~1 read per million
        abundances_arr = np.array(abundances, dtype=float)
        below_floor = abundances_arr < min_abundance
        if np.any(below_floor & (abundances_arr > 0)):
            n_boosted = int(np.sum(below_floor & (abundances_arr > 0)))
            abundances_arr[below_floor & (abundances_arr > 0)] = min_abundance
            # Renormalize to maintain total = 1.0
            abundances_arr /= abundances_arr.sum()
            logger.info(f"Applied minimum abundance floor ({min_abundance}) to {n_boosted} "
                       f"rare genomes to ensure read generation")
            abundances = list(abundances_arr)

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
        enable_benchmarking: bool = True,
        db_path: Optional[str] = None
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
                'viroforge_version': '0.13.0',
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

            # 2b. Per-genome taxonomy for taxonomy benchmarking (Module 4).
            # Export ncbi_taxid + full ICTV lineage so the taxonomy benchmark is
            # self-contained (no 500 MB database needed at benchmark time).
            if db_path:
                viral_ids = [seq.id for i, seq in enumerate(sequences) if i < n_viral]
                benchmarking['taxonomy'] = self._export_taxonomy(db_path, viral_ids)
                logger.info(f"  Taxonomy: exported for {len(benchmarking['taxonomy'])} viral genomes")

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
