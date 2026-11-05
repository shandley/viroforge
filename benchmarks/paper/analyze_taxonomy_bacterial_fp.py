#!/usr/bin/env python3
"""Analyze which bacterial background reads get classified by viral databases."""
from collections import Counter, defaultdict
from pathlib import Path

TAXDIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/taxonomy")
TMPDIR = Path("/scratch/sahlab/leran/tmp_analysis")

# Load ground truth
print("Loading preprocessed headers...")
bacterial_reads = {}  # name -> genome_subid (e.g., "bacterial_bg_gut_0075")
with open(TMPDIR / "preproc_headers.txt") as f:
    for line in f:
        parts = line.strip().split()
        name = parts[0].lstrip("@")
        source = None
        for p in parts:
            if p.startswith("source="):
                source = p.split("=", 1)[1]
                break
        if source == "bacterial_background":
            # Extract the bacterial fragment ID: bacterial_bg_<community>_<number>
            name_parts = name.split("_")
            # e.g., bacterial_bg_gut_0075_790_1 -> bacterial_bg_gut_0075
            if len(name_parts) >= 4:
                frag_id = "_".join(name_parts[:4])
            else:
                frag_id = name
            bacterial_reads[name] = frag_id

print(f"  Loaded {len(bacterial_reads)} bacterial background reads")

# Count reads per fragment
frag_counts = Counter(bacterial_reads.values())
print(f"  From {len(frag_counts)} distinct fragments")
print(f"  Top fragments by read count:")
for frag, cnt in frag_counts.most_common(10):
    print(f"    {frag}: {cnt}")

# Parse Kraken2 viral output
print("\nParsing Kraken2 viral output...")
k2v_classified_bact = {}  # name -> taxid for classified bacterial reads
with open(TAXDIR / "kraken2_viral.out") as f:
    for line in f:
        parts = line.strip().split("\t")
        status = parts[0]
        read_name = parts[1]
        taxid = int(parts[2])
        if read_name in bacterial_reads and status == "C":
            k2v_classified_bact[read_name] = taxid

print(f"  {len(k2v_classified_bact)} bacterial reads classified by Kraken2 viral")

# Which fragments contribute classified reads?
frag_classified = Counter()
frag_taxids = defaultdict(lambda: Counter())
for name, taxid in k2v_classified_bact.items():
    frag = bacterial_reads[name]
    frag_classified[frag] += 1
    frag_taxids[frag][taxid] += 1

print(f"\nBacterial fragments with MOST reads classified as viral:")
for frag, cnt in frag_classified.most_common(20):
    total = frag_counts[frag]
    rate = cnt / total * 100
    top_taxid = frag_taxids[frag].most_common(1)[0]
    print(f"  {frag}: {cnt}/{total} ({rate:.1f}%) | top taxid: {top_taxid[0]} ({top_taxid[1]} reads)")

# Check if challenge regions are over-represented
print(f"\nChallenge region analysis:")
challenge_frags = [f for f in frag_classified if "challenge" in f or int(f.split("_")[-1]) >= 120]
generic_frags = [f for f in frag_classified if f not in challenge_frags]
# The 15 challenge regions are fragments 120-134 (bacterial_bg_gut_0120 through 0134)
challenge_ids = set()
generic_ids = set()
for frag in frag_classified:
    parts = frag.split("_")
    try:
        idx = int(parts[-1])
        if idx >= 120:
            challenge_ids.add(frag)
        else:
            generic_ids.add(frag)
    except ValueError:
        generic_ids.add(frag)

challenge_clf = sum(frag_classified[f] for f in challenge_ids)
challenge_total = sum(frag_counts[f] for f in challenge_ids)
generic_clf = sum(frag_classified[f] for f in generic_ids)
generic_total = sum(frag_counts[f] for f in generic_ids)

print(f"  Challenge regions (idx >= 120): {len(challenge_ids)} fragments, {challenge_clf}/{challenge_total} reads classified ({challenge_clf/challenge_total*100:.1f}%)" if challenge_total > 0 else "  No challenge fragments found")
print(f"  Generic fragments (idx < 120): {len(generic_ids)} fragments, {generic_clf}/{generic_total} reads classified ({generic_clf/generic_total*100:.1f}%)" if generic_total > 0 else "  No generic fragments found")

# Now do same for Kaiju
print("\n\nParsing Kaiju output...")
kaiju_classified_bact = {}
with open(TAXDIR / "kaiju.out") as f:
    for line in f:
        parts = line.strip().split("\t")
        status = parts[0]
        read_name = parts[1]
        taxid = int(parts[2]) if len(parts) > 2 and parts[2].strip() else 0
        if read_name in bacterial_reads and status == "C":
            kaiju_classified_bact[read_name] = taxid

print(f"  {len(kaiju_classified_bact)} bacterial reads classified by Kaiju viral")

frag_classified_kaiju = Counter()
for name, taxid in kaiju_classified_bact.items():
    frag = bacterial_reads[name]
    frag_classified_kaiju[frag] += 1

print(f"\nBacterial fragments with MOST reads classified as viral (Kaiju):")
for frag, cnt in frag_classified_kaiju.most_common(20):
    total = frag_counts[frag]
    rate = cnt / total * 100
    print(f"  {frag}: {cnt}/{total} ({rate:.1f}%)")

# Kaiju challenge vs generic
kaiju_challenge_clf = sum(frag_classified_kaiju[f] for f in challenge_ids)
kaiju_generic_clf = sum(frag_classified_kaiju[f] for f in generic_ids)
print(f"\n  Challenge regions: {kaiju_challenge_clf}/{challenge_total} reads classified ({kaiju_challenge_clf/challenge_total*100:.1f}%)" if challenge_total > 0 else "")
print(f"  Generic fragments: {kaiju_generic_clf}/{generic_total} reads classified ({kaiju_generic_clf/generic_total*100:.1f}%)" if generic_total > 0 else "")

print("\nDone.")
