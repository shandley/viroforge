#!/usr/bin/env python3
"""
Build fungal reference fragments for ViroForge.

Downloads representative genomes for common human body site fungi,
extracts random 50kb fragments, and saves them as FASTA files.

Usage:
    python scripts/build_fungal_references.py                # Build all body sites
    python scripts/build_fungal_references.py --site gut_human
    python scripts/build_fungal_references.py --site vaginal_human
    python scripts/build_fungal_references.py --list-sites

Literature references:
- Gut: Nash et al. 2017 (mSphere 2:e00351); Sokol et al. 2017 (Gut 66:1039)
- Vaginal: Bradford & Ravel 2017 (Curr Opin Microbiol 40:58); Drell et al. 2013
- Skin: Findley et al. 2013 (Nature 498:367); Byrd et al. 2018
- Oral: Ghannoum et al. 2010 (PLoS Pathog 6:e1000713)
- Respiratory: Nguyen et al. 2015 (PLoS Pathog 11:e1004625)
"""

import argparse
import os
import sys
import json
import random
from pathlib import Path

# Import shared download/extract functions
sys.path.insert(0, str(Path(__file__).parent))
from build_bacterial_references import download_genome, extract_fragments

# ============================================================
# Body site fungal species lists
# ============================================================

# Gut mycobiome — Candida + Saccharomyces dominated
GUT_FUNGI = [
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.30, "prevalence": 0.70},
    {"species": "Candida tropicalis", "accession": "GCF_000006335.3", "abundance": 0.10, "prevalence": 0.35},
    {"species": "Candida parapsilosis", "accession": "GCF_000182765.1", "abundance": 0.05, "prevalence": 0.25},
    {"species": "Candida glabrata", "accession": "GCF_000002545.3", "abundance": 0.05, "prevalence": 0.20},
    {"species": "Saccharomyces cerevisiae", "accession": "GCF_000146045.2", "abundance": 0.20, "prevalence": 0.65},
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.08, "prevalence": 0.30},
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.05, "prevalence": 0.25},
    {"species": "Aspergillus niger", "accession": "GCF_000002855.4", "abundance": 0.04, "prevalence": 0.20},
    {"species": "Penicillium rubens", "accession": "GCF_000226395.1", "abundance": 0.03, "prevalence": 0.15},
    {"species": "Debaryomyces hansenii", "accession": "GCF_000006445.2", "abundance": 0.05, "prevalence": 0.20},
    {"species": "Cryptococcus neoformans", "accession": "GCF_000149245.1", "abundance": 0.02, "prevalence": 0.10},
    {"species": "Fusarium oxysporum", "accession": "GCF_013085055.1", "abundance": 0.03, "prevalence": 0.15},
]

# Vaginal mycobiome — Candida dominated (especially C. albicans)
# Bradford & Ravel 2017; Drell et al. 2013
VAGINAL_FUNGI = [
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.55, "prevalence": 0.60},
    {"species": "Candida glabrata", "accession": "GCF_000002545.3", "abundance": 0.15, "prevalence": 0.20},
    {"species": "Candida parapsilosis", "accession": "GCF_000182765.1", "abundance": 0.08, "prevalence": 0.15},
    {"species": "Candida tropicalis", "accession": "GCF_000006335.3", "abundance": 0.05, "prevalence": 0.10},
    {"species": "Candida krusei", "accession": "GCF_000956635.1", "abundance": 0.04, "prevalence": 0.08},
    {"species": "Saccharomyces cerevisiae", "accession": "GCF_000146045.2", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Rhodotorula mucilaginosa", "accession": "GCF_011057735.1", "abundance": 0.02, "prevalence": 0.08},
]

# Skin mycobiome — Malassezia dominated
# Findley et al. 2013 (Nature 498:367)
SKIN_FUNGI = [
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.40, "prevalence": 0.90},
    {"species": "Malassezia globosa", "accession": "GCF_000181695.1", "abundance": 0.25, "prevalence": 0.80},
    {"species": "Malassezia sympodialis", "accession": "GCF_001600215.1", "abundance": 0.10, "prevalence": 0.50},
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.05, "prevalence": 0.20},
    {"species": "Cryptococcus neoformans", "accession": "GCF_000149245.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Aspergillus niger", "accession": "GCF_000002855.4", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.04, "prevalence": 0.15},
    {"species": "Rhodotorula mucilaginosa", "accession": "GCF_011057735.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Penicillium rubens", "accession": "GCF_000226395.1", "abundance": 0.02, "prevalence": 0.08},
    {"species": "Epicoccum nigrum", "accession": "GCF_014839945.1", "abundance": 0.02, "prevalence": 0.08},
]

# Oral mycobiome — Candida + diverse environmental
# Ghannoum et al. 2010 (PLoS Pathog 6:e1000713)
ORAL_FUNGI = [
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.35, "prevalence": 0.75},
    {"species": "Candida parapsilosis", "accession": "GCF_000182765.1", "abundance": 0.08, "prevalence": 0.25},
    {"species": "Candida glabrata", "accession": "GCF_000002545.3", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.10, "prevalence": 0.40},
    {"species": "Aspergillus niger", "accession": "GCF_000002855.4", "abundance": 0.06, "prevalence": 0.20},
    {"species": "Fusarium oxysporum", "accession": "GCF_013085055.1", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Penicillium rubens", "accession": "GCF_000226395.1", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Saccharomyces cerevisiae", "accession": "GCF_000146045.2", "abundance": 0.08, "prevalence": 0.30},
    {"species": "Cryptococcus neoformans", "accession": "GCF_000149245.1", "abundance": 0.03, "prevalence": 0.08},
    {"species": "Rhodotorula mucilaginosa", "accession": "GCF_011057735.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.04, "prevalence": 0.12},
    {"species": "Debaryomyces hansenii", "accession": "GCF_000006445.2", "abundance": 0.03, "prevalence": 0.10},
]

# Respiratory/Lung mycobiome — Aspergillus + Candida + environmental
# Nguyen et al. 2015 (PLoS Pathog 11:e1004625)
RESPIRATORY_FUNGI = [
    {"species": "Aspergillus fumigatus", "accession": "GCF_000002655.1", "abundance": 0.20, "prevalence": 0.35},
    {"species": "Aspergillus niger", "accession": "GCF_000002855.4", "abundance": 0.08, "prevalence": 0.20},
    {"species": "Candida albicans", "accession": "GCF_000182965.3", "abundance": 0.20, "prevalence": 0.50},
    {"species": "Candida glabrata", "accession": "GCF_000002545.3", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Cladosporium herbarum", "accession": "GCF_013340325.1", "abundance": 0.10, "prevalence": 0.30},
    {"species": "Penicillium rubens", "accession": "GCF_000226395.1", "abundance": 0.08, "prevalence": 0.20},
    {"species": "Saccharomyces cerevisiae", "accession": "GCF_000146045.2", "abundance": 0.05, "prevalence": 0.15},
    {"species": "Pneumocystis jirovecii", "accession": "GCF_001477535.1", "abundance": 0.08, "prevalence": 0.15},
    {"species": "Malassezia restricta", "accession": "GCF_003290485.1", "abundance": 0.04, "prevalence": 0.10},
    {"species": "Cryptococcus neoformans", "accession": "GCF_000149245.1", "abundance": 0.04, "prevalence": 0.08},
    {"species": "Fusarium oxysporum", "accession": "GCF_013085055.1", "abundance": 0.03, "prevalence": 0.10},
]

# ============================================================
# Body site registry
# ============================================================
BODY_SITE_FUNGI = {
    "gut_human": {
        "species_list": GUT_FUNGI,
        "description": "Common human gut fungi (mycobiome)",
        "source": "Nash2017/Sokol2017",
    },
    "vaginal_human": {
        "species_list": VAGINAL_FUNGI,
        "description": "Vaginal mycobiome (Candida-dominated)",
        "source": "Bradford_Ravel2017/Drell2013",
    },
    "skin_human": {
        "species_list": SKIN_FUNGI,
        "description": "Skin mycobiome (Malassezia-dominated, sebaceous sites)",
        "source": "Findley2013/Byrd2018",
    },
    "oral_human": {
        "species_list": ORAL_FUNGI,
        "description": "Oral mycobiome (Candida + environmental fungi)",
        "source": "Ghannoum2010",
    },
    "respiratory_human": {
        "species_list": RESPIRATORY_FUNGI,
        "description": "Respiratory/lung mycobiome (Aspergillus + Candida)",
        "source": "Nguyen2015",
    },
    # Low-biomass sites share respiratory/gut fungi at very low levels
    "blood_human": {
        "species_list": GUT_FUNGI,  # Translocated from gut
        "description": "Blood mycobiome (translocated, very low biomass)",
        "source": "Proxy: gut fungi",
    },
    "ocular_human": {
        "species_list": SKIN_FUNGI,  # Similar to periocular skin
        "description": "Ocular mycobiome (periocular skin-like)",
        "source": "Proxy: skin fungi",
    },
    "urinary_human": {
        "species_list": VAGINAL_FUNGI,  # Anatomically adjacent
        "description": "Urinary mycobiome (Candida-dominated)",
        "source": "Proxy: vaginal fungi",
    },
    "lung_human": {
        "species_list": RESPIRATORY_FUNGI,
        "description": "Lower respiratory mycobiome",
        "source": "Nguyen2015",
    },
}

FRAGMENT_LENGTH = 50000
FRAGMENTS_PER_SPECIES = 20  # fewer than bacteria since fungi are a small fraction
RANDOM_SEED = 42


def build_site(body_site: str, output_base: Path, genome_cache: Path):
    """Build fungal reference fragments for a single body site."""
    if body_site not in BODY_SITE_FUNGI:
        print(f"ERROR: Unknown body site '{body_site}'")
        print(f"Available: {', '.join(BODY_SITE_FUNGI.keys())}")
        return False

    site_info = BODY_SITE_FUNGI[body_site]
    species_list = site_info["species_list"]

    output_base.mkdir(parents=True, exist_ok=True)
    genome_cache.mkdir(parents=True, exist_ok=True)

    rng = random.Random(RANDOM_SEED)

    print(f"\nBuilding fungal reference fragments: {body_site}")
    print(f"  Description: {site_info['description']}")
    print(f"  Species: {len(species_list)}")
    print(f"  Fragments per species: {FRAGMENTS_PER_SPECIES}")
    print(f"  Fragment length: {FRAGMENT_LENGTH} bp")
    print()

    all_fragments = []
    metadata = []
    failed = []

    for i, fungus in enumerate(species_list):
        species = fungus["species"]
        accession = fungus["accession"]
        print(f"  [{i+1}/{len(species_list)}] {species} ({accession})...", end=" ", flush=True)

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

    if not all_fragments:
        print(f"  ERROR: No fragments generated for {body_site}")
        return False

    # Write combined FASTA
    fasta_out = output_base / f"{body_site}_fragments.fasta"
    with open(fasta_out, 'w') as f:
        for frag_id, frag_desc, frag_seq in all_fragments:
            f.write(f">{frag_id} {frag_desc}\n")
            for j in range(0, len(frag_seq), 80):
                f.write(frag_seq[j:j+80] + "\n")

    # Write metadata
    meta_out = output_base / f"{body_site}_metadata.json"
    with open(meta_out, 'w') as f:
        json.dump({
            "body_site": body_site,
            "host_organism": "human",
            "description": site_info["description"],
            "source": site_info["source"],
            "n_species": len(metadata),
            "n_fragments_total": len(all_fragments),
            "fragment_length": FRAGMENT_LENGTH,
            "species": metadata,
        }, f, indent=2)

    # Write abundance table
    abund_out = output_base / f"{body_site}_abundances.tsv"
    with open(abund_out, 'w') as f:
        f.write("species\taccession\tabundance\tprevalence\n")
        for m in metadata:
            f.write(f"{m['species']}\t{m['accession']}\t{m['abundance']}\t{m['prevalence']}\n")

    print(f"\n  Output:")
    print(f"    FASTA: {fasta_out} ({len(all_fragments)} fragments)")
    print(f"    Metadata: {meta_out}")
    print(f"    Abundances: {abund_out}")
    print(f"    Successfully processed: {len(metadata)}/{len(species_list)} species")
    if failed:
        print(f"    Failed: {', '.join(failed)}")

    fasta_size = os.path.getsize(fasta_out)
    print(f"    FASTA size: {fasta_size / 1024 / 1024:.1f} MB")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Build fungal reference fragments for ViroForge body sites"
    )
    parser.add_argument(
        "--site", type=str, default=None,
        help="Body site to build (e.g., gut_human, vaginal_human). "
             "Default: build all sites."
    )
    parser.add_argument(
        "--list-sites", action="store_true",
        help="List available body sites and exit"
    )
    args = parser.parse_args()

    if args.list_sites:
        print("Available body sites:")
        for site, info in BODY_SITE_FUNGI.items():
            n_species = len(info["species_list"])
            print(f"  {site:25s} {n_species:3d} species  {info['description']}")
        return

    output_base = Path("viroforge/data/references/fungal_fragments")
    genome_cache = output_base / "genome_cache"

    if args.site:
        build_site(args.site, output_base, genome_cache)
    else:
        print(f"Building fungal references for {len(BODY_SITE_FUNGI)} body sites")
        results = {}
        for site in BODY_SITE_FUNGI:
            ok = build_site(site, output_base, genome_cache)
            results[site] = "OK" if ok else "FAILED"

        print("\n" + "=" * 60)
        print("Summary:")
        for site, status in results.items():
            print(f"  {site:25s} {status}")


if __name__ == "__main__":
    main()
