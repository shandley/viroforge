#!/usr/bin/env python3
"""Analyze taxonomy benchmark by source label - v3 fixed Kaiju parser."""
import json
from collections import Counter, defaultdict
from pathlib import Path

TAXDIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/taxonomy")
TMPDIR = Path("/scratch/sahlab/leran/tmp_analysis")

# Load ground truth from preprocessed headers
print("Loading preprocessed headers...")
read_source = {}
read_genome = {}
with open(TMPDIR / "preproc_headers.txt") as f:
    for line in f:
        parts = line.strip().split()
        name = parts[0].lstrip("@")
        source = None
        for p in parts:
            if p.startswith("source="):
                source = p.split("=", 1)[1]
                break
        read_source[name] = source
        if name.startswith("GCF_") or name.startswith("GCA_"):
            segs = name.split("_")
            genome_id = f"{segs[0]}_{segs[1]}_{segs[2]}"
        elif "_bg_" in name:
            genome_id = "bacterial_background"
        elif name.startswith("rrna_"):
            genome_id = "rrna_reference"
        elif name.startswith("host_"):
            if "herv" in name.lower():
                genome_id = "host_dna_herv"
            else:
                genome_id = "host_dna"
        elif name.startswith("phix_") or name.startswith("NC_001422"):
            genome_id = "phix174"
        elif name.startswith("reagent_"):
            genome_id = "reagent"
        elif name.startswith("fungal_"):
            genome_id = "fungal_background"
        elif name.startswith("archaeal_"):
            genome_id = "archaeal_background"
        elif name.startswith("dietary_"):
            genome_id = "dietary_background"
        else:
            genome_id = "unknown"
        read_genome[name] = genome_id
print(f"  Loaded {len(read_source)} reads")

VIRAL_SOURCES = {"viral", "dark_matter", "viral_prophage", "viral_provirus",
                 "viral_endogenous", "viral_satellite", "viral_viroid"}

def parse_classifier(path, fmt):
    """Parse Kraken2 or Kaiju output."""
    results = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            status = parts[0]  # C or U
            read_name = parts[1]
            taxid_str = parts[2] if len(parts) > 2 else "0"
            try:
                taxid = int(taxid_str)
            except ValueError:
                taxid = 0
            results[read_name] = (status == "C", taxid)
    return results

def analyze_classifier(clf_name, results):
    print(f"\n{'='*80}")
    print(f"  {clf_name}")
    print(f"{'='*80}")

    # Match reads between ground truth and classifier output
    matched = 0
    # Per source label: classified/unclassified counts
    clf_by_source = defaultdict(lambda: {"classified": 0, "unclassified": 0})
    # Non-viral reads that got classified
    nonviral_classified = Counter()
    # Viral reads that got unclassified
    viral_unclassified_by_source = Counter()
    viral_unclassified_genomes = Counter()
    # Per viral genome: classification rate
    genome_clf = defaultdict(int)
    genome_uncl = defaultdict(int)
    genome_total_reads = defaultdict(int)
    # Taxid diversity per genome
    genome_taxids = defaultdict(lambda: Counter())

    for read_name in read_source:
        if read_name not in results:
            # Try stripping /1 or /2 suffix for PE
            base = read_name.rsplit("/", 1)[0] if "/" in read_name else read_name
            if base not in results:
                continue
            classified, taxid = results[base]
        else:
            classified, taxid = results[read_name]

        matched += 1
        source = read_source[read_name]
        genome_id = read_genome[read_name]
        is_viral = source in VIRAL_SOURCES

        if classified:
            clf_by_source[source]["classified"] += 1
        else:
            clf_by_source[source]["unclassified"] += 1

        if not is_viral and classified:
            nonviral_classified[source] += 1

        if is_viral and not classified:
            viral_unclassified_by_source[source] += 1
            viral_unclassified_genomes[genome_id] += 1

        if genome_id.startswith("GC"):
            genome_total_reads[genome_id] += 1
            if classified:
                genome_clf[genome_id] += 1
                genome_taxids[genome_id][taxid] += 1
            else:
                genome_uncl[genome_id] += 1

    print(f"  Matched {matched:,} / {len(read_source):,} ground-truth reads")

    print(f"\n  Classification rate by source label:")
    for source in sorted(clf_by_source.keys()):
        stats = clf_by_source[source]
        total = stats["classified"] + stats["unclassified"]
        rate = stats["classified"] / total if total > 0 else 0
        print(f"    {source}: {stats['classified']:,}/{total:,} ({rate*100:.1f}%)")

    total_nv = sum(nonviral_classified.values())
    print(f"\n  Non-viral reads classified: {total_nv:,}")
    for src, cnt in nonviral_classified.most_common(10):
        print(f"    {src}: {cnt:,}")

    total_vu = sum(viral_unclassified_by_source.values())
    print(f"\n  Viral reads UNCLASSIFIED: {total_vu:,}")
    for src, cnt in viral_unclassified_by_source.most_common():
        print(f"    {src}: {cnt:,}")

    # Top genomes with unclassified viral reads
    print(f"\n  Top 30 viral genomes with MOST unclassified reads:")
    for genome_id, cnt in viral_unclassified_genomes.most_common(30):
        total = genome_total_reads.get(genome_id, cnt)
        clf = genome_clf.get(genome_id, 0)
        rate = clf / total * 100 if total > 0 else 0
        print(f"    {genome_id}: {cnt:,} unclassified / {total:,} total ({rate:.1f}% classified)")

    # Lowest classification rate genomes (min 10 reads)
    genome_rates = []
    for gid in genome_total_reads:
        total = genome_total_reads[gid]
        clf = genome_clf[gid]
        if total >= 10:
            genome_rates.append((gid, clf/total, total, clf))
    genome_rates.sort(key=lambda x: x[1])

    print(f"\n  Genomes with LOWEST classification rate (min 10 reads, bottom 30):")
    for gid, rate, total, clf in genome_rates[:30]:
        print(f"    {gid}: {clf}/{total} ({rate*100:.1f}%)")

    # Most diverse taxid assignments
    diverse = [(gid, len(tids), tids) for gid, tids in genome_taxids.items() if len(tids) > 1]
    diverse.sort(key=lambda x: x[1], reverse=True)
    print(f"\n  Genomes with MOST DIVERSE taxid assignments (top 15):")
    for gid, n_taxids, tids in diverse[:15]:
        top3 = tids.most_common(3)
        top3_str = ", ".join(f"taxid={t}:{c}" for t, c in top3)
        total = genome_total_reads[gid]
        print(f"    {gid}: {n_taxids} taxids / {total} reads | top: {top3_str}")

# Process each classifier
for clf_name, path, fmt in [
    ("kraken2_standard", TAXDIR / "kraken2_standard.out", "kraken2"),
    ("kraken2_viral", TAXDIR / "kraken2_viral.out", "kraken2"),
    ("kaiju_viral", TAXDIR / "kaiju.out", "kaiju"),
]:
    print(f"\nParsing {clf_name}...")
    results = parse_classifier(path, fmt)
    print(f"  Parsed {len(results):,} classifications")
    analyze_classifier(clf_name, results)

print("\nDone.")
