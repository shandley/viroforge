#!/usr/bin/env python3
"""Analyze false positives and false negatives for each QC benchmark step."""
import json, gzip
from collections import Counter
from pathlib import Path

BASE = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/qc")
TMPDIR = Path("/scratch/sahlab/leran/tmp_analysis")

print("Loading preprocessed headers...")
read_info = {}
with open(TMPDIR / "preproc_headers.txt") as f:
    for line in f:
        parts = line.strip().split()
        name = parts[0].lstrip("@")
        source = None
        for p in parts:
            if p.startswith("source="):
                source = p.split("=", 1)[1]
                break
        read_info[name] = source
print(f"  Loaded {len(read_info)} reads")

policies = {
    "step2_host": {"remove": {"host_dna", "host_dna_herv"}, "label": "Host DNA removal"},
    "step3_rrna": {"remove": {"rrna"}, "label": "rRNA removal"},
    "step4_phix": {"remove": {"phix"}, "label": "PhiX removal"},
    "step6_complexity": {"remove": {"artifact_low_complexity"}, "label": "Complexity filtering"},
}

tool_dirs = {
    "step2_host": {
        "bowtie2_viroforge_refs": "bowtie2_viroforge_refs",
        "bowtie2_grch38": "bowtie2_grch38",
    },
    "step3_rrna": {
        "sortmerna_viroforge_refs": "sortmerna_viroforge_refs",
        "bbduk_viroforge_refs": "bbduk_viroforge_refs",
        "bowtie2_silva": "bowtie2_silva",
    },
    "step4_phix": {
        "bowtie2": "bowtie2",
        "bbduk": "bbduk",
    },
    "step6_complexity": {
        "bbduk": "bbduk",
    },
}

def load_clean_names(step, tool_dir):
    clean_path = BASE / step / tool_dir / "clean_R1.fastq.gz"
    if not clean_path.exists():
        clean_path = BASE / step / tool_dir / "clean_R1.fastq"
        if not clean_path.exists():
            return None
    names = set()
    opener = gzip.open if str(clean_path).endswith(".gz") else open
    with opener(clean_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                names.add(line.strip().split()[0].lstrip("@"))
    return names

def extract_genome(name):
    if name.startswith("GCF_") or name.startswith("GCA_"):
        parts = name.split("_")
        return f"{parts[0]}_{parts[1]}_{parts[2]}"
    elif "_bg_" in name:
        prefix = name.split("_bg_")[0]
        return f"{prefix}_background"
    elif name.startswith("rrna_"):
        return "rrna_reference"
    elif name.startswith("host_"):
        if "herv" in name.lower():
            return "host_dna_herv"
        return "host_dna_fragment"
    elif name.startswith("phix_") or name.startswith("NC_001422"):
        return "phix174"
    elif name.startswith("reagent_"):
        return "reagent"
    else:
        return name.split("_")[0] + "_..."

for step_name, policy in policies.items():
    print(f"\n{'='*80}")
    print(f"  {policy['label']} ({step_name})")
    print(f"{'='*80}")
    remove_labels = policy["remove"]

    for tool_name, tool_dir in tool_dirs.get(step_name, {}).items():
        clean_names = load_clean_names(step_name, tool_dir)
        if clean_names is None:
            print(f"\n  [{tool_name}] - cleaned file not found")
            continue

        fp_sources = Counter()
        fp_genomes = Counter()
        fn_sources = Counter()
        fn_genomes = Counter()
        tp = 0
        tn = 0

        for name, source in read_info.items():
            was_removed = name not in clean_names
            should_remove = source in remove_labels

            if was_removed and should_remove:
                tp += 1
            elif was_removed and not should_remove:
                fp_sources[source] += 1
                fp_genomes[extract_genome(name)] += 1
            elif not was_removed and should_remove:
                fn_sources[source] += 1
                fn_genomes[extract_genome(name)] += 1
            else:
                tn += 1

        total_fp = sum(fp_sources.values())
        total_fn = sum(fn_sources.values())

        print(f"\n  [{tool_name}]")
        print(f"  TP={tp}, FP={total_fp}, FN={total_fn}, TN={tn}")

        if fp_sources:
            print(f"\n  FALSE POSITIVES ({total_fp} reads removed that should have been kept):")
            for src, cnt in fp_sources.most_common():
                print(f"    {src}: {cnt:,} ({cnt/total_fp*100:.1f}%)")
            print(f"  Top FP genomes:")
            for genome, cnt in fp_genomes.most_common(10):
                print(f"    {genome}: {cnt:,}")

        if fn_sources:
            print(f"\n  FALSE NEGATIVES ({total_fn} reads kept that should have been removed):")
            for src, cnt in fn_sources.most_common():
                print(f"    {src}: {cnt:,} ({cnt/total_fn*100:.1f}%)")
            print(f"  Top FN genomes:")
            for genome, cnt in fn_genomes.most_common(10):
                print(f"    {genome}: {cnt:,}")

        if not fp_sources and not fn_sources:
            print("  Perfect - no FP or FN.")

print("\nDone.")
