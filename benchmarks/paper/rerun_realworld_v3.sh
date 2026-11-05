#!/bin/bash
#SBATCH --job-name=realworld_v3
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch/sahlab/leran/paper_benchmark/realworld_v3_%j.out
#SBATCH --error=/scratch/sahlab/leran/paper_benchmark/realworld_v3_%j.err

set -euo pipefail

# ── Environment ──────────────────────────────────────────────────────────────
QC_ENV=/scratch/sahlab/leran/miniforge3/envs/qc_bench/bin
export PATH="$QC_ENV:$PATH"
export LD_LIBRARY_PATH=/scratch/sahlab/leran/miniforge3/envs/qc_bench/lib:${LD_LIBRARY_PATH:-}

THREADS=16
cd /scratch/sahlab/leran/viroforge

DATASET=/scratch/sahlab/leran/viroforge/data/gut_10x
QC_DIR=/scratch/sahlab/leran/paper_benchmark/gut_10x/qc
TAX_DIR=/scratch/sahlab/leran/paper_benchmark/gut_10x/taxonomy
PREP_R1=$QC_DIR/preprocessed/preprocessed_R1.fastq.gz
PREP_R2=$QC_DIR/preprocessed/preprocessed_R2.fastq.gz

# All databases on /scratch
INDEXES=/scratch/sahlab/leran/paper_benchmark/indexes
GRCH38_INDEX=$INDEXES/GRCh38_noalt_as
SILVA_INDEX=$INDEXES/SILVA_128_LSUParc_SSUParc_ribosomal_RNA
KAIJU_FMI=$INDEXES/kaiju_db_viruses.fmi
KAIJU_NODES=/scratch/sahlab/leran/paper_benchmark/kaiju_nodes.dmp
KAIJU_NAMES=$INDEXES/kaiju_names.dmp
TAXDUMP=$INDEXES/taxdump

echo "============================================================"
echo "REAL-WORLD DATABASE BENCHMARK v3 — GUT 10x"
echo "============================================================"
echo "Date: $(date)"
echo "Fixes: pipe bowtie2 directly to BAM; handle SILVA dups; fix singleton loss"
echo ""

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: HOST REMOVAL WITH FULL GRCH38
# ═══════════════════════════════════════════════════════════════════════════
echo "=== STEP 1: Host removal with GRCh38 ==="
HOST_DIR=$QC_DIR/step2_host/bowtie2_grch38
mkdir -p $HOST_DIR
rm -f $HOST_DIR/clean_R1.fastq.gz $HOST_DIR/clean_R2.fastq.gz

echo "  Aligning and extracting in one pipe (no intermediate SAM)..."
# Use bowtie2 --un-conc-gz to directly write unmapped read pairs
# This avoids the samtools name-sort singleton-loss problem entirely
bowtie2 -x $GRCH38_INDEX \
    -1 $PREP_R1 -2 $PREP_R2 \
    --very-sensitive \
    --threads $THREADS \
    --un-conc-gz $HOST_DIR/clean_R%.fastq.gz \
    -S /dev/null \
    2>$HOST_DIR/bowtie2.log

# bowtie2 --un-conc-gz writes clean_R1.fastq.gz and clean_R2.fastq.gz
echo "  Bowtie2 alignment stats:"
cat $HOST_DIR/bowtie2.log

NCLEAN=$(zcat $HOST_DIR/clean_R1.fastq.gz | wc -l | awk '{print $1/4}')
echo "  Clean read pairs: $NCLEAN"
echo "  DONE"
echo ""

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: rRNA REMOVAL WITH SILVA (bowtie2, kneaddata-style)
# ═══════════════════════════════════════════════════════════════════════════
echo "=== STEP 2: rRNA removal with bowtie2 + SILVA_128 ==="
RRNA_DIR=$QC_DIR/step3_rrna/bowtie2_silva
mkdir -p $RRNA_DIR
rm -f $RRNA_DIR/clean_R1.fastq.gz $RRNA_DIR/clean_R2.fastq.gz

echo "  Aligning to SILVA_128 and extracting unmapped pairs directly..."
# Same approach: --un-conc-gz avoids the SAM header dup issue entirely
bowtie2 -x $SILVA_INDEX \
    -1 $PREP_R1 -2 $PREP_R2 \
    --very-sensitive \
    --threads $THREADS \
    --large-index \
    --un-conc-gz $RRNA_DIR/clean_R%.fastq.gz \
    -S /dev/null \
    2>$RRNA_DIR/bowtie2.log

echo "  Bowtie2 SILVA alignment stats:"
cat $RRNA_DIR/bowtie2.log

NCLEAN=$(zcat $RRNA_DIR/clean_R1.fastq.gz | wc -l | awk '{print $1/4}')
echo "  Clean read pairs: $NCLEAN"
echo "  DONE"
echo ""

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: KAIJU TAXONOMY CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════
echo "=== STEP 3: Kaiju taxonomy classification ==="

echo "  Running Kaiju with RefSeq_virus_prot DB..."
kaiju \
    -t $KAIJU_NODES \
    -f $KAIJU_FMI \
    -i $PREP_R1 \
    -j $PREP_R2 \
    -z $THREADS \
    -o $TAX_DIR/kaiju.out \
    -v 2>$TAX_DIR/kaiju.log

CLASSIFIED=$(grep -c '^C' $TAX_DIR/kaiju.out || echo 0)
TOTAL=$(wc -l < $TAX_DIR/kaiju.out)
echo "  Kaiju classified: $CLASSIFIED / $TOTAL reads"

kaiju2table \
    -t $KAIJU_NODES \
    -n $KAIJU_NAMES \
    -r species \
    -o $TAX_DIR/kaiju.report \
    $TAX_DIR/kaiju.out 2>$TAX_DIR/kaiju2table.log || echo "  kaiju2table failed (non-critical)"

echo "  DONE"
echo ""

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: SCORE ALL QC RESULTS
# ═══════════════════════════════════════════════════════════════════════════
echo "============================================================"
echo "SCORING ALL RESULTS"
echo "============================================================"

python3 << 'SCORE_EOF'
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/sahlab/leran/viroforge")
from viroforge.benchmarking.parsers import read_labels, read_names
from viroforge.benchmarking.qc import benchmark_qc

RAW_R1 = Path("/scratch/sahlab/leran/viroforge/data/gut_10x/fastq/gut_virome_-_adult_healthy_western_diet_R1.fastq")
QC_DIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/qc")
RESULTS = QC_DIR / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

raw = read_labels(RAW_R1)
print(f"Raw reads loaded: {len(raw):,}")

POLICIES = {
    "step2_host": {
        "viral": "keep", "dark_matter": "keep", "erv_exogenous": "keep",
        "host_dna": "remove", "host_dna_herv": "remove",
        "rrna": "keep", "phix": "keep", "reagent_bacteria": "keep",
        "bacterial_background": "keep", "fungal_background": "keep",
        "archaeal_background": "keep", "artifact_low_complexity": "keep",
        "erv_endogenous": "keep", "dietary_background": "keep",
    },
    "step3_rrna": {
        "viral": "keep", "dark_matter": "keep", "erv_exogenous": "keep",
        "host_dna": "keep", "host_dna_herv": "keep",
        "rrna": "remove",
        "phix": "keep", "reagent_bacteria": "keep",
        "bacterial_background": "keep", "fungal_background": "keep",
        "archaeal_background": "keep", "artifact_low_complexity": "keep",
        "erv_endogenous": "keep", "dietary_background": "keep",
    },
}

TOOLS = {
    "step2_host": {
        "bowtie2_viroforge_refs": QC_DIR / "step2_host/bowtie2/clean_R1.fastq.gz",
        "bowtie2_grch38": QC_DIR / "step2_host/bowtie2_grch38/clean_R1.fastq.gz",
    },
    "step3_rrna": {
        "sortmerna_viroforge_refs": QC_DIR / "step3_rrna/sortmerna/clean_fwd.fq.gz",
        "bbduk_viroforge_refs": QC_DIR / "step3_rrna/bbduk/clean_R1.fastq.gz",
        "bowtie2_silva": QC_DIR / "step3_rrna/bowtie2_silva/clean_R1.fastq.gz",
    },
}

all_results = {}
for step, tools in TOOLS.items():
    policy = POLICIES[step]
    step_results = {}
    for tool_name, clean_path in tools.items():
        if not clean_path.exists():
            print(f"  SKIP: {step}/{tool_name} -- not found at {clean_path}")
            continue
        print(f"  Scoring: {step}/{tool_name}")
        kept = read_names(clean_path)
        print(f"    Kept reads: {len(kept):,}")
        metrics = benchmark_qc(raw, kept, keep_remove=policy)
        out_file = RESULTS / f"{step}_{tool_name}.json"
        with open(out_file, "w") as f:
            json.dump(metrics, f, indent=2, default=str)
        vr = metrics.get("viral_retention")
        cont = metrics.get("contamination", {})
        rc, pr, f1 = cont.get("recall"), cont.get("precision"), cont.get("f1")
        if all(v is not None for v in [vr, rc, pr, f1]):
            print(f"    viral_ret={vr:.4f}  recall={rc:.4f}  precision={pr:.4f}  f1={f1:.4f}")
        step_results[tool_name] = metrics
    all_results[step] = step_results

# Merge existing phix/dedup/complexity
existing_file = QC_DIR / "all_qc_benchmarks.json"
if existing_file.exists():
    with open(existing_file) as f:
        existing = json.load(f)
    for step in ["step4_phix", "step5_dedup", "step6_complexity"]:
        if step in existing:
            all_results[step] = existing[step]

with open(QC_DIR / "all_qc_benchmarks_realworld.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

print("\n" + "=" * 110)
print("QC BENCHMARK RESULTS -- GUT 10x (REAL-WORLD + VIROFORGE DATABASES)")
print("=" * 110)
for step, tools in sorted(all_results.items()):
    if not tools:
        continue
    print(f"\n-- {step} --")
    print(f"  {'Tool':<30s}  {'Viral Ret':>10s}  {'Recall':>10s}  {'Precision':>10s}  {'F1':>10s}")
    print(f"  {'-'*30}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
    for tool, m in tools.items():
        def fmt(v): return f"{v:.4f}" if isinstance(v, float) else "N/A"
        vr = m.get("viral_retention")
        cont = m.get("contamination", {})
        print(f"  {tool:<30s}  {fmt(vr):>10s}  "
              f"{fmt(cont.get('recall')):>10s}  "
              f"{fmt(cont.get('precision')):>10s}  "
              f"{fmt(cont.get('f1')):>10s}")
        if "dedup" in m and m.get("dedup", {}).get("n_duplicates"):
            dd = m["dedup"]
            print(f"    dedup: {dd['n_duplicates']:,} duplicates, "
                  f"removal rate: {fmt(dd.get('dedup_rate'))}")
print("\nQC scoring complete.")
SCORE_EOF

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: SCORE KAIJU TAXONOMY
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Scoring Kaiju taxonomy ==="

python3 << 'TAXSCORE_EOF'
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/sahlab/leran/viroforge")
from viroforge.benchmarking.taxonomy import benchmark_taxonomy

METADATA = Path("/scratch/sahlab/leran/viroforge/data/gut_10x/metadata/gut_virome_-_adult_healthy_western_diet_metadata.json")
RAW_R1 = Path("/scratch/sahlab/leran/viroforge/data/gut_10x/fastq/gut_virome_-_adult_healthy_western_diet_R1.fastq")
TAX_DIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/taxonomy")
RESULTS = TAX_DIR / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
TAXDUMP = Path("/scratch/sahlab/leran/paper_benchmark/indexes/taxdump")

with open(METADATA) as f:
    metadata = json.load(f)

kaiju_out = TAX_DIR / "kaiju.out"
if kaiju_out.exists() and kaiju_out.stat().st_size > 0:
    print("Scoring Kaiju (RefSeq_virus_prot)...")
    result = benchmark_taxonomy(
        classifier_output=str(kaiju_out),
        raw_reads=str(RAW_R1),
        metadata=metadata,
        taxdump_dir=str(TAXDUMP),
        tool="kaiju",
    )
    with open(RESULTS / "kaiju_viral.json", "w") as f:
        json.dump(result, f, indent=2, default=str)

    rl = result.get("read_level", {})
    print(f"  Sensitivity: {rl.get('sensitivity', 'N/A')}")
    print(f"  Species F1:  {rl.get('species_f1', 'N/A')}")
    print(f"  Genus F1:    {rl.get('genus_f1', 'N/A')}")
    print(f"  Family F1:   {rl.get('family_f1', 'N/A')}")
else:
    print("ERROR: kaiju.out not found or empty")
    result = None

# Merge with existing taxonomy results
existing_file = TAX_DIR / "all_taxonomy_benchmarks.json"
all_tax = {}
if existing_file.exists():
    with open(existing_file) as f:
        all_tax = json.load(f)
if result:
    all_tax["kaiju_viral"] = result

with open(TAX_DIR / "all_taxonomy_benchmarks_realworld.json", "w") as f:
    json.dump(all_tax, f, indent=2, default=str)

print("\n" + "=" * 100)
print("TAXONOMY BENCHMARK RESULTS -- GUT 10x (ALL TOOLS)")
print("=" * 100)
print(f"  {'Tool':<25s}  {'Sensitivity':>12s}  {'Species F1':>12s}  {'Genus F1':>12s}  {'Family F1':>12s}")
print(f"  {'-'*25}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")
for name, r in all_tax.items():
    rl = r.get("read_level", {})
    def fmt(v):
        if v is None: return "N/A"
        if isinstance(v, float): return f"{v:.4f}"
        return str(v)
    print(f"  {name:<25s}  {fmt(rl.get('sensitivity')):>12s}  "
          f"{fmt(rl.get('species_f1')):>12s}  "
          f"{fmt(rl.get('genus_f1')):>12s}  "
          f"{fmt(rl.get('family_f1')):>12s}")
print("\nAll scoring complete.")
TAXSCORE_EOF

echo ""
echo "============================================================"
echo "REAL-WORLD BENCHMARK v3 COMPLETE: $(date)"
echo "============================================================"
echo "Output:"
echo "  QC:       $QC_DIR/all_qc_benchmarks_realworld.json"
echo "  Taxonomy: $TAX_DIR/all_taxonomy_benchmarks_realworld.json"
