#!/usr/bin/env python3
"""
Build bacterial reference fragments for ViroForge.

Downloads representative genomes for the top 50 human gut bacteria,
extracts random 10kb fragments, and saves them as a FASTA file
with abundance metadata.

Species list and abundances based on:
- Qin et al. 2010 (MetaHIT) - Nature 464:59-65
- Lloyd-Price et al. 2019 (HMP2/IBDMDB) - Nature 569:655-662
- Almeida et al. 2021 (UHGG v2) - Nature Biotechnol 39:105-114
"""

import os
import sys
import json
import random
import subprocess
from pathlib import Path

# Top 50 human gut bacterial species with:
# - RefSeq accession (representative genome)
# - Typical relative abundance in healthy adult gut (Western diet)
# - Prevalence (fraction of samples where detected)
GUT_BACTERIA = [
    # Bacteroidetes (dominant in Western diet)
    {"species": "Bacteroides vulgatus", "accession": "GCF_000012825.1", "abundance": 0.065, "prevalence": 0.90},
    {"species": "Bacteroides uniformis", "accession": "GCF_000154205.1", "abundance": 0.055, "prevalence": 0.85},
    {"species": "Bacteroides thetaiotaomicron", "accession": "GCF_000011065.1", "abundance": 0.045, "prevalence": 0.80},
    {"species": "Bacteroides fragilis", "accession": "GCF_000009925.1", "abundance": 0.035, "prevalence": 0.75},
    {"species": "Bacteroides ovatus", "accession": "GCF_000154125.1", "abundance": 0.030, "prevalence": 0.70},
    {"species": "Bacteroides caccae", "accession": "GCF_000154845.1", "abundance": 0.020, "prevalence": 0.60},
    {"species": "Parabacteroides distasonis", "accession": "GCF_000012845.1", "abundance": 0.025, "prevalence": 0.70},
    {"species": "Parabacteroides merdae", "accession": "GCF_000154465.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Prevotella copri", "accession": "GCF_000157935.1", "abundance": 0.040, "prevalence": 0.35},
    {"species": "Alistipes putredinis", "accession": "GCF_000154525.1", "abundance": 0.020, "prevalence": 0.65},

    # Firmicutes - Lachnospiraceae
    {"species": "Roseburia intestinalis", "accession": "GCF_000169035.1", "abundance": 0.035, "prevalence": 0.75},
    {"species": "Eubacterium rectale", "accession": "GCF_000020605.1", "abundance": 0.040, "prevalence": 0.80},
    {"species": "Coprococcus comes", "accession": "GCF_000155815.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Dorea formicigenerans", "accession": "GCF_000156015.1", "abundance": 0.012, "prevalence": 0.60},
    {"species": "Dorea longicatena", "accession": "GCF_000155875.1", "abundance": 0.012, "prevalence": 0.55},
    {"species": "Blautia obeum", "accession": "GCF_000154065.1", "abundance": 0.018, "prevalence": 0.65},

    # Firmicutes - Ruminococcaceae
    {"species": "Faecalibacterium prausnitzii", "accession": "GCF_000162015.1", "abundance": 0.055, "prevalence": 0.90},
    {"species": "Ruminococcus bromii", "accession": "GCF_000156375.1", "abundance": 0.025, "prevalence": 0.70},
    {"species": "Subdoligranulum variabile", "accession": "GCF_000177015.2", "abundance": 0.015, "prevalence": 0.55},

    # Firmicutes - Clostridiaceae / other
    {"species": "Clostridium bolteae", "accession": "GCF_000155435.1", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Clostridium symbiosum", "accession": "GCF_000155205.1", "abundance": 0.008, "prevalence": 0.40},
    {"species": "Dialister invisus", "accession": "GCF_000157955.2", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Phascolarctobacterium succinatutens", "accession": "GCF_000424985.1", "abundance": 0.012, "prevalence": 0.50},

    # Firmicutes - Erysipelotrichaceae
    {"species": "Holdemanella biformis", "accession": "GCF_000156495.1", "abundance": 0.008, "prevalence": 0.35},

    # Firmicutes - Lactobacillaceae / Streptococcaceae
    {"species": "Streptococcus thermophilus", "accession": "GCF_000011825.1", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Lactobacillus ruminis", "accession": "GCF_000224985.2", "abundance": 0.003, "prevalence": 0.25},

    # Actinobacteria
    {"species": "Bifidobacterium longum", "accession": "GCF_000007525.1", "abundance": 0.030, "prevalence": 0.75},
    {"species": "Bifidobacterium adolescentis", "accession": "GCF_000010425.1", "abundance": 0.020, "prevalence": 0.60},
    {"species": "Bifidobacterium bifidum", "accession": "GCF_000167455.2", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Collinsella aerofaciens", "accession": "GCF_000169035.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Eggerthella lenta", "accession": "GCF_000024265.1", "abundance": 0.008, "prevalence": 0.40},

    # Verrucomicrobia
    {"species": "Akkermansia muciniphila", "accession": "GCF_000020225.1", "abundance": 0.025, "prevalence": 0.65},

    # Proteobacteria (low abundance in healthy gut)
    {"species": "Escherichia coli", "accession": "GCF_000005845.2", "abundance": 0.010, "prevalence": 0.70},
    {"species": "Klebsiella pneumoniae", "accession": "GCF_000240185.1", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Sutterella wadsworthensis", "accession": "GCF_000154245.1", "abundance": 0.008, "prevalence": 0.40},
    {"species": "Bilophila wadsworthia", "accession": "GCF_000185705.2", "abundance": 0.006, "prevalence": 0.35},
    {"species": "Desulfovibrio piger", "accession": "GCF_000156375.1", "abundance": 0.004, "prevalence": 0.25},
    {"species": "Haemophilus parainfluenzae", "accession": "GCF_000191405.1", "abundance": 0.003, "prevalence": 0.20},

    # Fusobacteria
    {"species": "Fusobacterium nucleatum", "accession": "GCF_000007325.1", "abundance": 0.003, "prevalence": 0.20},

    # Archaea (methanogenic)
    {"species": "Methanobrevibacter smithii", "accession": "GCF_000016525.1", "abundance": 0.015, "prevalence": 0.50},

    # Additional common species
    {"species": "Odoribacter splanchnicus", "accession": "GCF_000190535.1", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Barnesiella intestinihominis", "accession": "GCF_000311925.1", "abundance": 0.008, "prevalence": 0.35},
    {"species": "Megamonas hypermegale", "accession": "GCF_000428645.1", "abundance": 0.006, "prevalence": 0.25},
    {"species": "Mitsuokella multacida", "accession": "GCF_000178935.2", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Veillonella parvula", "accession": "GCF_000024945.1", "abundance": 0.005, "prevalence": 0.35},
    {"species": "Megasphaera elsdenii", "accession": "GCF_000284795.1", "abundance": 0.004, "prevalence": 0.20},
    {"species": "Citrobacter freundii", "accession": "GCF_000648515.1", "abundance": 0.003, "prevalence": 0.15},
    {"species": "Enterococcus faecalis", "accession": "GCF_000007785.1", "abundance": 0.003, "prevalence": 0.25},
    {"species": "Ruminococcus gnavus", "accession": "GCF_000156515.1", "abundance": 0.018, "prevalence": 0.65},
    {"species": "Coprococcus eutactus", "accession": "GCF_000155435.1", "abundance": 0.012, "prevalence": 0.50},
    {"species": "Anaerostipes hadrus", "accession": "GCF_000332875.2", "abundance": 0.014, "prevalence": 0.55},
]

FRAGMENT_LENGTH = 10000
FRAGMENTS_PER_SPECIES = 50
RANDOM_SEED = 42


def download_genome(accession: str, output_dir: Path) -> Path:
    """Download genome from NCBI using datasets CLI or efetch."""
    fasta_path = output_dir / f"{accession}.fasta"
    if fasta_path.exists():
        return fasta_path

    # Try NCBI datasets CLI first
    try:
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(output_dir / f"{accession}.zip")
        ]
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)

        # Extract
        import zipfile
        zip_path = output_dir / f"{accession}.zip"
        with zipfile.ZipFile(zip_path) as zf:
            for name in zf.namelist():
                if name.endswith(".fna"):
                    with zf.open(name) as src, open(fasta_path, 'wb') as dst:
                        dst.write(src.read())
                    break
        zip_path.unlink()
        return fasta_path
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Fallback: efetch via Entrez
    try:
        from Bio import Entrez, SeqIO
        Entrez.email = "viroforge@example.com"

        # Get assembly info to find nucleotide accession
        handle = Entrez.esearch(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            assembly_id = record["IdList"][0]
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()

            # Get FTP path
            ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0].get('FtpPath_RefSeq', '')
            if ftp_path:
                import urllib.request
                fname = ftp_path.split('/')[-1]
                url = f"{ftp_path}/{fname}_genomic.fna.gz"
                gz_path = output_dir / f"{accession}.fna.gz"
                urllib.request.urlretrieve(url, gz_path)

                import gzip
                with gzip.open(gz_path, 'rt') as gz, open(fasta_path, 'w') as out:
                    out.write(gz.read())
                gz_path.unlink()
                return fasta_path
    except Exception as e:
        print(f"  Warning: Could not download {accession}: {e}")

    return None


def extract_fragments(fasta_path: Path, n_fragments: int, fragment_length: int,
                      species_name: str, rng: random.Random) -> list:
    """Extract random fragments from a genome FASTA."""
    from Bio import SeqIO

    # Read all sequences
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) >= fragment_length:
            sequences.append(record)

    if not sequences:
        print(f"  Warning: No sequences >= {fragment_length}bp in {fasta_path}")
        return []

    fragments = []
    species_tag = species_name.replace(" ", "_")

    for i in range(n_fragments):
        # Pick a random sequence weighted by length
        weights = [len(s.seq) - fragment_length for s in sequences]
        total_w = sum(weights)
        if total_w <= 0:
            continue

        seq_record = rng.choices(sequences, weights=weights, k=1)[0]
        max_start = len(seq_record.seq) - fragment_length
        start = rng.randint(0, max_start)
        end = start + fragment_length

        frag_seq = str(seq_record.seq[start:end])
        frag_id = f"bact_{species_tag}_{i:04d}"
        frag_desc = f"[{species_name}] {seq_record.id}:{start+1}-{end}"

        fragments.append((frag_id, frag_desc, frag_seq))

    return fragments


def main():
    output_dir = Path("viroforge/data/references/bacterial_fragments")
    genome_cache = output_dir / "genome_cache"
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_cache.mkdir(parents=True, exist_ok=True)

    rng = random.Random(RANDOM_SEED)

    print(f"Building bacterial reference fragments for human gut microbiome")
    print(f"  Species: {len(GUT_BACTERIA)}")
    print(f"  Fragments per species: {FRAGMENTS_PER_SPECIES}")
    print(f"  Fragment length: {FRAGMENT_LENGTH} bp")
    print()

    all_fragments = []
    metadata = []
    failed = []

    for i, bact in enumerate(GUT_BACTERIA):
        species = bact["species"]
        accession = bact["accession"]
        print(f"  [{i+1}/{len(GUT_BACTERIA)}] {species} ({accession})...", end=" ", flush=True)

        fasta_path = download_genome(accession, genome_cache)
        if fasta_path is None or not fasta_path.exists():
            print("FAILED (download)")
            failed.append(species)
            continue

        fragments = extract_fragments(fasta_path, FRAGMENTS_PER_SPECIES,
                                      FRAGMENT_LENGTH, species, rng)
        if not fragments:
            print(f"FAILED (no fragments)")
            failed.append(species)
            continue

        all_fragments.extend(fragments)
        metadata.append({
            "species": species,
            "accession": accession,
            "abundance": bact["abundance"],
            "prevalence": bact["prevalence"],
            "n_fragments": len(fragments),
            "fragment_length": FRAGMENT_LENGTH,
        })
        print(f"OK ({len(fragments)} fragments)")

    # Write combined FASTA
    fasta_out = output_dir / "gut_human_fragments.fasta"
    with open(fasta_out, 'w') as f:
        for frag_id, frag_desc, frag_seq in all_fragments:
            f.write(f">{frag_id} {frag_desc}\n")
            # Write sequence in 80-char lines
            for j in range(0, len(frag_seq), 80):
                f.write(frag_seq[j:j+80] + "\n")

    # Write metadata
    meta_out = output_dir / "gut_human_metadata.json"
    with open(meta_out, 'w') as f:
        json.dump({
            "body_site": "gut",
            "host_organism": "human",
            "description": "Top 50 human gut bacteria (healthy Western diet)",
            "source": "MetaHIT/HMP2/UHGG",
            "n_species": len(metadata),
            "n_fragments_total": len(all_fragments),
            "fragment_length": FRAGMENT_LENGTH,
            "species": metadata,
        }, f, indent=2)

    # Write abundance table (TSV for easy use)
    abund_out = output_dir / "gut_human_abundances.tsv"
    with open(abund_out, 'w') as f:
        f.write("species\taccession\tabundance\tprevalence\n")
        for m in metadata:
            f.write(f"{m['species']}\t{m['accession']}\t{m['abundance']}\t{m['prevalence']}\n")

    print(f"\n  Output:")
    print(f"    FASTA: {fasta_out} ({len(all_fragments)} fragments)")
    print(f"    Metadata: {meta_out}")
    print(f"    Abundances: {abund_out}")
    print(f"    Successfully processed: {len(metadata)}/{len(GUT_BACTERIA)} species")
    if failed:
        print(f"    Failed: {', '.join(failed)}")

    fasta_size = os.path.getsize(fasta_out)
    print(f"    FASTA size: {fasta_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
