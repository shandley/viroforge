#!/usr/bin/env python3
"""
Curate reference sequences for realistic contamination modeling.

Downloads and bundles small reference FASTA files that ViroForge uses
to generate biologically realistic contaminant reads (detectable by
real QC tools like SortMeRNA, BBDuk, fastp, Kraken2).

Reference sources:
- PhiX174: NCBI RefSeq NC_001422.1 (5,386 bp)
- rRNA: NCBI RefSeq RNA database (representative 16S/18S/23S/28S/5S/5.8S)
- Host fragments: NCBI T2T-CHM13v2.0 (50 x 10kb fragments from all chromosomes)
- Adapters: Illumina published sequences (TruSeq, Nextera)

All sequences are public domain (NCBI) or freely redistributable.

Usage:
    python scripts/curate_reference_sequences.py [--output-dir viroforge/data/references]
    python scripts/curate_reference_sequences.py --skip-host  # skip large host download
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).parent.parent / "viroforge" / "data" / "references"


# ---------------------------------------------------------------------------
# Illumina adapter sequences (published, widely redistributed)
# ---------------------------------------------------------------------------

ILLUMINA_ADAPTERS = {
    # TruSeq adapters
    "TruSeq_Universal_Adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "TruSeq_Adapter_Read1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "TruSeq_Adapter_Read2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "TruSeq_Index_Adapter_i7": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "TruSeq_Index_Adapter_i5": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    # Nextera adapters
    "Nextera_Transposase_Adapter1": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "Nextera_Transposase_Adapter2": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "Nextera_Read1_Adapter": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "Nextera_Read2_Adapter": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    # Common short adapter contaminants (what appears in read-through)
    "TruSeq_R1_Readthrough": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "TruSeq_R2_Readthrough": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
}


# ---------------------------------------------------------------------------
# rRNA accessions from NCBI RefSeq RNA
# ---------------------------------------------------------------------------

RRNA_ACCESSIONS = {
    # E. coli (primary model organism for rRNA contamination)
    "NR_024570.1": "Escherichia_coli_16S_rRNA",
    "NR_103073.1": "Escherichia_coli_23S_rRNA",
    "NR_132729.1": "Escherichia_coli_5S_rRNA",
    # Human rRNA (host contamination)
    "NR_003286.4": "Homo_sapiens_18S_rRNA",
    "NR_003287.4": "Homo_sapiens_28S_rRNA",
    "NR_003285.3": "Homo_sapiens_5.8S_rRNA",
    "NR_023379.1": "Homo_sapiens_5S_rRNA",
    # Common gut bacteria 16S
    "NR_074548.2": "Bacteroides_fragilis_16S_rRNA",
    "NR_028961.1": "Faecalibacterium_prausnitzii_16S_rRNA",
    "NR_042311.1": "Bifidobacterium_longum_16S_rRNA",
    "NR_112558.1": "Lactobacillus_acidophilus_16S_rRNA",
    "NR_113347.1": "Clostridium_perfringens_16S_rRNA",
    "NR_074596.1": "Prevotella_copri_16S_rRNA",
    # Common environmental bacteria 16S
    "NR_113599.1": "Pseudomonas_aeruginosa_16S_rRNA",
    "NR_114042.1": "Staphylococcus_aureus_16S_rRNA",
    "NR_074910.1": "Bacillus_subtilis_16S_rRNA",
    # Archaea 16S
    "NR_042783.1": "Methanobrevibacter_smithii_16S_rRNA",
    # Mouse rRNA (common host)
    "NR_003278.3": "Mus_musculus_18S_rRNA",
    "NR_003279.1": "Mus_musculus_28S_rRNA",
    # Additional common 16S for diversity
    "NR_112116.2": "Enterococcus_faecalis_16S_rRNA",
    "NR_074540.1": "Streptococcus_pneumoniae_16S_rRNA",
    "NR_025082.1": "Salmonella_enterica_16S_rRNA",
    "NR_102783.2": "Akkermansia_muciniphila_16S_rRNA",
    "NR_113250.1": "Ruminococcus_bromii_16S_rRNA",
}

# PhiX174 accession
PHIX_ACCESSION = "NC_001422.1"

# Human genome assembly for host fragments
# T2T-CHM13v2.0 chromosome accessions
HUMAN_CHROMOSOMES = {
    "chr1": "NC_060925.1",
    "chr2": "NC_060926.1",
    "chr3": "NC_060927.1",
    "chr4": "NC_060928.1",
    "chr5": "NC_060929.1",
    "chr6": "NC_060930.1",
    "chr7": "NC_060931.1",
    "chr8": "NC_060932.1",
    "chr9": "NC_060933.1",
    "chr10": "NC_060934.1",
    "chr11": "NC_060935.1",
    "chr12": "NC_060936.1",
    "chr13": "NC_060937.1",
    "chr14": "NC_060938.1",
    "chr15": "NC_060939.1",
    "chr16": "NC_060940.1",
    "chr17": "NC_060941.1",
    "chr18": "NC_060942.1",
    "chr19": "NC_060943.1",
    "chr20": "NC_060944.1",
    "chr21": "NC_060945.1",
    "chr22": "NC_060946.1",
    "chrX": "NC_060947.1",
}

# Positions chosen to sample diverse genomic regions (gene-rich, gene-poor,
# repetitive). Each entry is (chromosome, start_position, length).
# These are deterministic so the bundled file is reproducible.
HOST_FRAGMENT_REGIONS = [
    ("chr1", 1_000_000, 10_000),
    ("chr1", 50_000_000, 10_000),
    ("chr1", 150_000_000, 10_000),
    ("chr2", 10_000_000, 10_000),
    ("chr2", 100_000_000, 10_000),
    ("chr3", 20_000_000, 10_000),
    ("chr3", 120_000_000, 10_000),
    ("chr4", 30_000_000, 10_000),
    ("chr4", 90_000_000, 10_000),
    ("chr5", 15_000_000, 10_000),
    ("chr5", 80_000_000, 10_000),
    ("chr6", 25_000_000, 10_000),
    ("chr6", 100_000_000, 10_000),
    ("chr7", 5_000_000, 10_000),
    ("chr7", 75_000_000, 10_000),
    ("chr8", 10_000_000, 10_000),
    ("chr8", 70_000_000, 10_000),
    ("chr9", 20_000_000, 10_000),
    ("chr9", 65_000_000, 10_000),
    ("chr10", 15_000_000, 10_000),
    ("chr10", 80_000_000, 10_000),
    ("chr11", 5_000_000, 10_000),
    ("chr11", 60_000_000, 10_000),
    ("chr12", 10_000_000, 10_000),
    ("chr12", 90_000_000, 10_000),
    ("chr13", 20_000_000, 10_000),
    ("chr13", 50_000_000, 10_000),
    ("chr14", 25_000_000, 10_000),
    ("chr14", 60_000_000, 10_000),
    ("chr15", 20_000_000, 10_000),
    ("chr15", 55_000_000, 10_000),
    ("chr16", 10_000_000, 10_000),
    ("chr16", 50_000_000, 10_000),
    ("chr17", 5_000_000, 10_000),
    ("chr17", 45_000_000, 10_000),
    ("chr18", 15_000_000, 10_000),
    ("chr18", 40_000_000, 10_000),
    ("chr19", 5_000_000, 10_000),
    ("chr19", 30_000_000, 10_000),
    ("chr20", 10_000_000, 10_000),
    ("chr20", 35_000_000, 10_000),
    ("chr21", 10_000_000, 10_000),
    ("chr21", 25_000_000, 10_000),
    ("chr22", 15_000_000, 10_000),
    ("chr22", 30_000_000, 10_000),
    ("chrX", 10_000_000, 10_000),
    ("chrX", 75_000_000, 10_000),
    ("chrX", 120_000_000, 10_000),
]


def fetch_sequences_from_ncbi(accessions: dict[str, str], output_path: Path) -> None:
    """Fetch sequences from NCBI Entrez and write to FASTA."""
    from Bio import Entrez, SeqIO

    Entrez.email = "viroforge@example.com"

    records = []
    for accession, label in accessions.items():
        logger.info(f"Fetching {accession} ({label})...")
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="fasta", retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            record.id = label
            record.description = f"{label} [{accession}]"
            records.append(record)
        except Exception as e:
            logger.warning(f"Failed to fetch {accession}: {e}")

    SeqIO.write(records, output_path, "fasta")
    logger.info(f"Wrote {len(records)} sequences to {output_path}")


def fetch_phix(output_path: Path) -> None:
    """Fetch PhiX174 genome from NCBI."""
    fetch_sequences_from_ncbi({PHIX_ACCESSION: "PhiX174"}, output_path)


def fetch_rrna(output_path: Path) -> None:
    """Fetch representative rRNA sequences from NCBI."""
    fetch_sequences_from_ncbi(RRNA_ACCESSIONS, output_path)


def write_adapters(output_path: Path) -> None:
    """Write Illumina adapter sequences to FASTA."""
    with open(output_path, "w") as f:
        for name, seq in ILLUMINA_ADAPTERS.items():
            f.write(f">{name}\n{seq}\n")
    logger.info(f"Wrote {len(ILLUMINA_ADAPTERS)} adapter sequences to {output_path}")


def fetch_host_fragments(output_path: Path) -> None:
    """Fetch human genome fragments from NCBI T2T-CHM13v2.0.

    Downloads only the specific regions defined in HOST_FRAGMENT_REGIONS,
    not entire chromosomes.
    """
    from Bio import Entrez, SeqIO

    Entrez.email = "viroforge@example.com"

    records = []
    # Group regions by chromosome to minimize fetches
    from collections import defaultdict

    chrom_regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, length in HOST_FRAGMENT_REGIONS:
        chrom_regions[chrom].append((start, length))

    for chrom, regions in sorted(chrom_regions.items()):
        accession = HUMAN_CHROMOSOMES.get(chrom)
        if not accession:
            logger.warning(f"No accession for {chrom}, skipping")
            continue

        for start, length in regions:
            end = start + length
            label = f"human_{chrom}_{start}_{end}"
            logger.info(f"Fetching {label} from {accession}...")
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="fasta",
                    retmode="text",
                    seq_start=start + 1,  # Entrez uses 1-based coords
                    seq_stop=end,
                )
                record = SeqIO.read(handle, "fasta")
                handle.close()
                record.id = label
                record.description = f"{label} [T2T-CHM13v2.0 {accession}:{start+1}-{end}]"
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to fetch {label}: {e}")

    SeqIO.write(records, output_path, "fasta")
    logger.info(f"Wrote {len(records)} host fragments to {output_path}")


def verify_output(output_path: Path, expected_min: int) -> bool:
    """Verify a FASTA file has at least expected_min sequences."""
    from Bio import SeqIO

    count = sum(1 for _ in SeqIO.parse(output_path, "fasta"))
    total_bp = sum(len(r.seq) for r in SeqIO.parse(output_path, "fasta"))
    logger.info(f"  {output_path.name}: {count} sequences, {total_bp:,} bp")
    if count < expected_min:
        logger.error(
            f"  Expected at least {expected_min} sequences, got {count}"
        )
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Curate reference sequences for ViroForge contamination modeling"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Output directory for reference FASTA files",
    )
    parser.add_argument(
        "--skip-host",
        action="store_true",
        help="Skip host fragment download (requires many NCBI fetches)",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Adapters (no download needed)
    logger.info("=== Writing adapter sequences ===")
    write_adapters(args.output_dir / "adapters.fasta")

    # PhiX174
    logger.info("=== Fetching PhiX174 ===")
    fetch_phix(args.output_dir / "phix174.fasta")

    # rRNA representatives
    logger.info("=== Fetching rRNA representatives ===")
    fetch_rrna(args.output_dir / "rrna_representatives.fasta")

    # Host fragments
    if not args.skip_host:
        logger.info("=== Fetching human genome fragments (T2T-CHM13v2.0) ===")
        logger.info("This downloads 47 x 10kb regions from NCBI (may take a few minutes)")
        fetch_host_fragments(args.output_dir / "host_fragments.fasta")
    else:
        logger.info("=== Skipping host fragments (--skip-host) ===")

    # Verification
    logger.info("=== Verification ===")
    ok = True
    ok &= verify_output(args.output_dir / "adapters.fasta", 10)
    ok &= verify_output(args.output_dir / "phix174.fasta", 1)
    ok &= verify_output(args.output_dir / "rrna_representatives.fasta", 15)
    if not args.skip_host:
        ok &= verify_output(args.output_dir / "host_fragments.fasta", 40)

    if ok:
        logger.info("All reference files curated successfully.")
    else:
        logger.error("Some files failed verification.")
        sys.exit(1)


if __name__ == "__main__":
    main()
