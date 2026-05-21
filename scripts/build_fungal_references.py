#!/usr/bin/env python3
"""
Build fungal reference fragments for ViroForge.

Downloads representative genomes for common human gut/body site fungi,
extracts random 50kb fragments, and saves them as a FASTA file.

Species based on:
- Nash et al. 2017 (HMP mycobiome) - mSphere 2:e00351-17
- Sokol et al. 2017 (IBD mycobiome) - Gut 66:1039-1048
- Auchtung et al. 2018 - mSystems 3:e00092-18
"""

import os
import sys
import json
import random
from pathlib import Path

# Import shared download/extract functions
sys.path.insert(0, str(Path(__file__).parent))
from build_bacterial_references import download_genome, extract_fragments

# Common human-associated fungi
# Gut mycobiome is dominated by Candida + Saccharomyces
GUT_FUNGI = [
    # Candida (dominant in gut, especially IBD)
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.30, "prevalence": 0.70},
    {"species": "Candida tropicalis", "accession": "GCF_000006335.3", "abundance": 0.10, "prevalence": 0.35},
    {"species": "Candida parapsilosis", "accession": "GCF_000182765.1", "abundance": 0.05, "prevalence": 0.25},
    {"species": "Candida glabrata", "accession": "GCF_000002545.3", "abundance": 0.05, "prevalence": 0.20},

    # Saccharomyces (dietary, very common)
    {"species": "Saccharomyces cerevisiae", "accession": "GCF_000146045.2", "abundance": 0.20, "prevalence": 0.65},

    # Malassezia (skin commensal, detected in gut)
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.08, "prevalence": 0.30},

    # Cladosporium (environmental, common in gut surveys)
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.05, "prevalence": 0.25},

    # Aspergillus (environmental/food-borne)
    {"species": "Aspergillus niger", "accession": "GCF_000002855.4", "abundance": 0.04, "prevalence": 0.20},

    # Penicillium (environmental/food-borne)
    {"species": "Penicillium rubens", "accession": "GCF_000226395.1", "abundance": 0.03, "prevalence": 0.15},

    # Debaryomyces (food-associated, cheese/meat)
    {"species": "Debaryomyces hansenii", "accession": "GCF_000006445.2", "abundance": 0.05, "prevalence": 0.20},

    # Cryptococcus (environmental, occasional gut detection)
    {"species": "Cryptococcus neoformans", "accession": "GCF_000149245.1", "abundance": 0.02, "prevalence": 0.10},

    # Fusarium (environmental)
    {"species": "Fusarium oxysporum", "accession": "GCF_013085055.1", "abundance": 0.03, "prevalence": 0.15},
]

FRAGMENT_LENGTH = 50000
FRAGMENTS_PER_SPECIES = 20  # fewer than bacteria since fungi are a small fraction
RANDOM_SEED = 42


def main():
    output_dir = Path("viroforge/data/references/fungal_fragments")
    genome_cache = output_dir / "genome_cache"
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_cache.mkdir(parents=True, exist_ok=True)

    rng = random.Random(RANDOM_SEED)

    print(f"Building fungal reference fragments for human gut mycobiome")
    print(f"  Species: {len(GUT_FUNGI)}")
    print(f"  Fragments per species: {FRAGMENTS_PER_SPECIES}")
    print(f"  Fragment length: {FRAGMENT_LENGTH} bp")
    print()

    all_fragments = []
    metadata = []
    failed = []

    for i, fungus in enumerate(GUT_FUNGI):
        species = fungus["species"]
        accession = fungus["accession"]
        print(f"  [{i+1}/{len(GUT_FUNGI)}] {species} ({accession})...", end=" ", flush=True)

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
            "abundance": fungus["abundance"],
            "prevalence": fungus["prevalence"],
            "n_fragments": len(fragments),
            "fragment_length": FRAGMENT_LENGTH,
        })
        print(f"OK ({len(fragments)} fragments)")

    # Write combined FASTA
    fasta_out = output_dir / "gut_human_fragments.fasta"
    with open(fasta_out, 'w') as f:
        for frag_id, frag_desc, frag_seq in all_fragments:
            f.write(f">{frag_id} {frag_desc}\n")
            for j in range(0, len(frag_seq), 80):
                f.write(frag_seq[j:j+80] + "\n")

    # Write metadata
    meta_out = output_dir / "gut_human_metadata.json"
    with open(meta_out, 'w') as f:
        json.dump({
            "body_site": "gut",
            "host_organism": "human",
            "description": "Common human gut fungi (mycobiome)",
            "source": "Nash et al. 2017, Sokol et al. 2017",
            "n_species": len(metadata),
            "n_fragments_total": len(all_fragments),
            "fragment_length": FRAGMENT_LENGTH,
            "species": metadata,
        }, f, indent=2)

    # Write abundance table
    abund_out = output_dir / "gut_human_abundances.tsv"
    with open(abund_out, 'w') as f:
        f.write("species\taccession\tabundance\tprevalence\n")
        for m in metadata:
            f.write(f"{m['species']}\t{m['accession']}\t{m['abundance']}\t{m['prevalence']}\n")

    print(f"\n  Output:")
    print(f"    FASTA: {fasta_out} ({len(all_fragments)} fragments)")
    print(f"    Metadata: {meta_out}")
    print(f"    Abundances: {abund_out}")
    print(f"    Successfully processed: {len(metadata)}/{len(GUT_FUNGI)} species")
    if failed:
        print(f"    Failed: {', '.join(failed)}")

    fasta_size = os.path.getsize(fasta_out)
    print(f"    FASTA size: {fasta_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
