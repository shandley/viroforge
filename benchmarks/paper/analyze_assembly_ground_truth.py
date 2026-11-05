#!/usr/bin/env python3
"""Analyze assembly chimeras and missing genomes with ground truth."""
import json
import sqlite3
from collections import Counter
from pathlib import Path

BENCH = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/assembly/all_assembly_benchmarks.json")
DB = Path("/scratch/sahlab/leran/viroforge/viroforge/data/viral_genomes.db")

with open(BENCH) as f:
    data = json.load(f)

conn = sqlite3.connect(str(DB))

def get_genome_info(genome_id):
    row = conn.execute(
        "SELECT g.genome_name, g.length, g.gc_content, t.family, t.genus "
        "FROM genomes g JOIN taxonomy t ON g.genome_id = t.genome_id "
        "WHERE g.genome_id = ?", (genome_id,)
    ).fetchone()
    if row:
        return {"name": row[0][:60], "length": row[1], "gc": row[2], "family": row[3], "genus": row[4]}
    return None

def get_abundance(genome_id, collection_id=1):
    row = conn.execute(
        "SELECT relative_abundance FROM collection_genomes "
        "WHERE genome_id = ? AND collection_id = ?", (genome_id, collection_id)
    ).fetchone()
    return row[0] if row else None

print("="*80)
print("CHIMERA ANALYSIS")
print("="*80)

for assembler in ["megahit", "metaspades"]:
    d = data[assembler]
    chimeras = d.get("chimeras", {})
    examples = chimeras.get("examples", [])

    print(f"\n--- {assembler.upper()} chimeras: {chimeras.get('n', 0)} ({chimeras.get('rate', 0):.4f}) ---")

    # Classify chimeras: viral-viral, viral-bacterial, bacterial-bacterial
    viral_viral = []
    viral_bacterial = []
    bacterial_bacterial = []
    other = []

    for ch in examples:
        segments = ch.get("segments", [])
        types = []
        for seg in segments:
            genome = seg["genome"]
            if genome.startswith("bacterial_bg_"):
                types.append("bacterial")
            elif genome.startswith("GCF_") or genome.startswith("GCA_"):
                types.append("viral")
            else:
                types.append("other")

        types_set = set(types)
        if types_set == {"viral"}:
            viral_viral.append(ch)
        elif types_set == {"bacterial"}:
            bacterial_bacterial.append(ch)
        elif "viral" in types_set and "bacterial" in types_set:
            viral_bacterial.append(ch)
        else:
            other.append(ch)

    print(f"  Viral-viral: {len(viral_viral)}")
    print(f"  Bacterial-bacterial: {len(bacterial_bacterial)}")
    print(f"  Viral-bacterial: {len(viral_bacterial)}")
    print(f"  Other: {len(other)}")

    if viral_viral:
        print(f"\n  Viral-viral chimeras:")
        for ch in viral_viral[:5]:
            contig = ch["contig"]
            length = ch["length"]
            genomes = [s["genome"] for s in ch["segments"]]
            genome_infos = []
            for g in genomes:
                info = get_genome_info(g)
                if info:
                    genome_infos.append(f"{g} ({info['name'][:30]}, {info['family']})")
                else:
                    genome_infos.append(g)
            print(f"    {contig} ({length} bp): {' + '.join(genome_infos)}")

    if viral_bacterial:
        print(f"\n  Viral-bacterial chimeras:")
        for ch in viral_bacterial[:5]:
            contig = ch["contig"]
            length = ch["length"]
            genomes = [s["genome"] for s in ch["segments"]]
            genome_infos = []
            for g in genomes:
                if g.startswith("bacterial_bg_"):
                    genome_infos.append(f"{g} (bacterial)")
                else:
                    info = get_genome_info(g)
                    if info:
                        genome_infos.append(f"{g} ({info['name'][:30]}, {info['family']})")
                    else:
                        genome_infos.append(g)
            print(f"    {contig} ({length} bp): {' + '.join(genome_infos)}")

    if bacterial_bacterial:
        print(f"\n  Bacterial-bacterial chimeras: {len(bacterial_bacterial)} (all between bacterial background fragments)")

print("\n" + "="*80)
print("MISSING GENOME ANALYSIS")
print("="*80)

for assembler in ["megahit", "metaspades"]:
    d = data[assembler]
    per_genome = d["per_genome"]

    missing = []
    for gid, info in per_genome.items():
        if info["category"] == "missing":
            gi = get_genome_info(gid)
            ab = get_abundance(gid)
            missing.append({"id": gid, "gi": gi, "abundance": ab or 0, **info})

    print(f"\n--- {assembler.upper()}: {len(missing)} missing genomes ---")

    # By family
    fam_counts = Counter()
    for m in missing:
        if m["gi"]:
            fam_counts[m["gi"]["family"]] += 1
    print(f"By family: {dict(fam_counts.most_common(10))}")

    # By size
    sizes = {"<5kb": 0, "5-10kb": 0, "10-50kb": 0, "50-100kb": 0, ">100kb": 0}
    for m in missing:
        if m["gi"]:
            l = m["gi"]["length"]
            if l < 5000: sizes["<5kb"] += 1
            elif l < 10000: sizes["5-10kb"] += 1
            elif l < 50000: sizes["10-50kb"] += 1
            elif l < 100000: sizes["50-100kb"] += 1
            else: sizes[">100kb"] += 1
    print(f"By size: {sizes}")

    # By GC content
    gc_bins = {"<35%": 0, "35-45%": 0, "45-55%": 0, "55-65%": 0, ">65%": 0}
    for m in missing:
        if m["gi"] and m["gi"]["gc"] is not None:
            gc = m["gi"]["gc"]
            if gc < 0.35: gc_bins["<35%"] += 1
            elif gc < 0.45: gc_bins["35-45%"] += 1
            elif gc < 0.55: gc_bins["45-55%"] += 1
            elif gc < 0.65: gc_bins["55-65%"] += 1
            else: gc_bins[">65%"] += 1
    print(f"By GC: {gc_bins}")

    # List all with details
    missing.sort(key=lambda x: x["abundance"], reverse=True)
    print(f"\nAll missing (by abundance):")
    for m in missing:
        gi = m["gi"]
        if gi:
            print(f"  {m['id']}: {gi['name']} | {gi['length']:,} bp | GC={gi['gc']:.1%} | {gi['family']} | abund={m['abundance']:.6f}")
        else:
            print(f"  {m['id']}: NOT IN DB | abund={m['abundance']:.6f}")

conn.close()
print("\nDone.")
