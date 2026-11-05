#!/bin/bash
#SBATCH --job-name=bench_10x
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/sahlab/leran/paper_benchmark/bench_10x_%j.out
#SBATCH --error=/scratch/sahlab/leran/paper_benchmark/bench_10x_%j.err

set -uo pipefail
# Note: -e removed so individual tool failures don't kill the whole pipeline

# ── Environment ──────────────────────────────────────────────────────────────
QC_ENV=/scratch/sahlab/leran/miniforge3/envs/qc_bench/bin
export PATH="$QC_ENV:$PATH"
export LD_LIBRARY_PATH=/scratch/sahlab/leran/miniforge3/envs/qc_bench/lib:${LD_LIBRARY_PATH:-}
FASTP=/scratch/sahlab/leran/fastp  # static binary

# Shared tools and databases on the cluster
KRAKEN2_BIN=/ref/sahlab/software/kraken2/kraken2
KRAKEN2_STANDARD_DB=/ref/sahlab/data/kraken2_db/kraken2_standard
KRAKEN2_VIRAL_DB=/ref/sahlab/data/viral_analysis_DBs/kraken_viral_Nucl_DB
TAXDUMP_DIR=/ref/sahlab/data/kraken2_db/kraken2_standard/taxonomy
KAIJU_BIN=/ref/sahlab/software/miniforge3/envs/kaiju/bin/kaiju
KAIJU_DB=/ref/gdlab/data/kaiju_db/RefSeq_virus_prot/kaiju_db_viruses.fmi
KAIJU_NODES=/ref/gdlab/data/kaiju_db/RefSeq_virus_prot/nodes.dmp

# Paths
DATASET=/scratch/sahlab/leran/viroforge/data/gut_10x
RAW_R1=$DATASET/fastq/gut_virome_-_adult_healthy_western_diet_R1.fastq
RAW_R2=$DATASET/fastq/gut_virome_-_adult_healthy_western_diet_R2.fastq
METADATA=$DATASET/metadata/gut_virome_-_adult_healthy_western_diet_metadata.json
GENOME_FASTA=$DATASET/fasta/gut_virome_-_adult_healthy_western_diet.fasta

WORKDIR=/scratch/sahlab/leran/paper_benchmark/gut_10x
IDX_DIR=/scratch/sahlab/leran/paper_benchmark/qc/indexes  # reuse existing indices
THREADS=8

cd /scratch/sahlab/leran/viroforge
mkdir -p $WORKDIR/{qc,assembly,taxonomy}/{results}

echo "=== GUT 10x FULL BENCHMARK PIPELINE: $(date) ==="
echo "Dataset: $DATASET"
echo "Reads: $(grep -c '^@' $RAW_R1) read pairs"

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: QC BENCHMARK
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "================================================================"
echo "PART 1: QC BENCHMARK"
echo "================================================================"

QC_DIR=$WORKDIR/qc
PREP_DIR=$QC_DIR/preprocessed
mkdir -p $PREP_DIR

# ── Step 1: fastp preprocessing (adapter + quality only) ─────────────────
echo ""
echo "=== Step 1: fastp preprocessing ==="
PREP_R1=$PREP_DIR/preprocessed_R1.fastq.gz
PREP_R2=$PREP_DIR/preprocessed_R2.fastq.gz

if [ ! -f "$PREP_R1" ]; then
    $FASTP \
        -i $RAW_R1 -I $RAW_R2 \
        -o $PREP_R1 -O $PREP_R2 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        --length_required 50 \
        --dont_eval_duplication \
        -j $PREP_DIR/fastp_preprocess.json \
        -h $PREP_DIR/fastp_preprocess.html \
        -w $THREADS
    echo "fastp preprocessing done: $(zcat $PREP_R1 | grep -c '^@') read pairs"
else
    echo "Preprocessed reads already exist, skipping"
fi

# ── Step 2: Host DNA removal ─────────────────────────────────────────────
echo ""
echo "=== Step 2: Host DNA removal ==="

STEP2=$QC_DIR/step2_host
mkdir -p $STEP2/{bowtie2,bwa_mem2}

# bowtie2
echo "  Testing: bowtie2"
BT2_CLEAN=$STEP2/bowtie2/clean_R1.fastq.gz
if [ ! -f "$BT2_CLEAN" ]; then
    bowtie2 -x ${IDX_DIR}/host_bt2 \
        -1 $PREP_R1 -2 $PREP_R2 \
        --very-sensitive -p $THREADS \
        --un-conc-gz $STEP2/bowtie2/clean_R%.fastq.gz \
        -S /dev/null 2>$STEP2/bowtie2/bowtie2.log
    echo "  bowtie2 host removal done"
else
    echo "  bowtie2 output exists, skipping"
fi

# bwa-mem2 (if available) — wrapped because ViroForge read names have /1 /2
# suffixes that bwa-mem2 treats as mismatched pair names
BWA_CLEAN=$STEP2/bwa_mem2/clean_R1.fastq.gz
if command -v bwa-mem2 &>/dev/null; then
    echo "  Testing: bwa-mem2"
    if [ ! -f "$BWA_CLEAN" ]; then
        BWA_IDX=$QC_DIR/indexes/host_bwa
        mkdir -p $BWA_IDX
        HOST_REF=$(python3 -c "import viroforge; from pathlib import Path; print(Path(viroforge.__file__).parent / 'data/references/host_fragments.fasta')")
        if [ ! -f "$BWA_IDX/host.fasta.bwt.2bit.64" ]; then
            cp $HOST_REF $BWA_IDX/host.fasta
            bwa-mem2 index $BWA_IDX/host.fasta
        fi
        # Use -5SP and pipe through fixmate to handle /1 /2 suffixes
        if bwa-mem2 mem -t $THREADS -5SP $BWA_IDX/host.fasta $PREP_R1 $PREP_R2 2>$STEP2/bwa_mem2/bwa.log \
            | samtools view -b -@ 4 -o $STEP2/bwa_mem2/mapped.bam; then
            samtools sort -n -@ 4 -o $STEP2/bwa_mem2/namesorted.bam $STEP2/bwa_mem2/mapped.bam
            samtools fastq -f 12 -F 256 $STEP2/bwa_mem2/namesorted.bam \
                -1 $BWA_CLEAN -2 $STEP2/bwa_mem2/clean_R2.fastq.gz \
                -0 /dev/null -s /dev/null 2>$STEP2/bwa_mem2/extract.log
            echo "  bwa-mem2 host removal done"
        else
            echo "  bwa-mem2 FAILED (read name mismatch) — skipping"
        fi
    else
        echo "  bwa-mem2 output exists, skipping"
    fi
else
    echo "  bwa-mem2 not available, skipping"
fi

# ── Step 3: rRNA removal ─────────────────────────────────────────────────
echo ""
echo "=== Step 3: rRNA removal ==="

STEP3=$QC_DIR/step3_rrna
mkdir -p $STEP3/{sortmerna,bbduk,ribodetector}

RRNA_REF=$(python3 -c "import viroforge; from pathlib import Path; print(Path(viroforge.__file__).parent / 'data/references/rrna_representatives.fasta')")

# SortMeRNA
echo "  Testing: SortMeRNA"
SMR_CLEAN=$STEP3/sortmerna/clean_fwd.fq.gz
if [ ! -f "$SMR_CLEAN" ]; then
    rm -rf $STEP3/sortmerna/kvdb 2>/dev/null
    sortmerna \
        --ref $RRNA_REF \
        --reads $PREP_R1 --reads $PREP_R2 \
        --aligned $STEP3/sortmerna/rrna --other $STEP3/sortmerna/clean \
        --paired_in --fastx --out2 \
        --threads $THREADS --workdir $STEP3/sortmerna/kvdb
    echo "  SortMeRNA done"
else
    echo "  SortMeRNA output exists, skipping"
fi

# bbduk
echo "  Testing: bbduk (rRNA)"
BBDUK_RRNA_CLEAN=$STEP3/bbduk/clean_R1.fastq.gz
if [ ! -f "$BBDUK_RRNA_CLEAN" ]; then
    bbduk.sh -Xmx8g in1=$PREP_R1 in2=$PREP_R2 \
        out1=$BBDUK_RRNA_CLEAN out2=$STEP3/bbduk/clean_R2.fastq.gz \
        ref=$RRNA_REF k=31 hdist=1 \
        stats=$STEP3/bbduk/bbduk_stats.txt
    echo "  bbduk rRNA done"
else
    echo "  bbduk rRNA output exists, skipping"
fi

# ribodetector
echo "  Testing: ribodetector"
RIBO_CLEAN=$STEP3/ribodetector/clean_R1.fastq.gz
if [ ! -f "$RIBO_CLEAN" ]; then
    ribodetector_cpu -t $THREADS \
        -i $PREP_R1 $PREP_R2 \
        -o $RIBO_CLEAN $STEP3/ribodetector/clean_R2.fastq.gz \
        -r $STEP3/ribodetector/rrna_R1.fastq.gz $STEP3/ribodetector/rrna_R2.fastq.gz \
        --len 150 -e rrna || echo "  ribodetector FAILED"
    echo "  ribodetector done"
else
    echo "  ribodetector output exists, skipping"
fi

# ── Step 4: PhiX removal ─────────────────────────────────────────────────
echo ""
echo "=== Step 4: PhiX removal ==="

STEP4=$QC_DIR/step4_phix
mkdir -p $STEP4/{bowtie2,bbduk}

# bowtie2
echo "  Testing: bowtie2 (PhiX)"
BT2_PHIX_CLEAN=$STEP4/bowtie2/clean_R1.fastq.gz
if [ ! -f "$BT2_PHIX_CLEAN" ]; then
    bowtie2 -x ${IDX_DIR}/phix_bt2 \
        -1 $PREP_R1 -2 $PREP_R2 \
        --very-sensitive -p $THREADS \
        --un-conc-gz $STEP4/bowtie2/clean_R%.fastq.gz \
        -S /dev/null 2>$STEP4/bowtie2/bowtie2.log
    echo "  bowtie2 PhiX done"
else
    echo "  bowtie2 PhiX output exists, skipping"
fi

# bbduk
PHIX_REF=$(python3 -c "import viroforge; from pathlib import Path; print(Path(viroforge.__file__).parent / 'data/references/phix174.fasta')")
echo "  Testing: bbduk (PhiX)"
BBDUK_PHIX_CLEAN=$STEP4/bbduk/clean_R1.fastq.gz
if [ ! -f "$BBDUK_PHIX_CLEAN" ]; then
    bbduk.sh -Xmx8g in1=$PREP_R1 in2=$PREP_R2 \
        out1=$BBDUK_PHIX_CLEAN out2=$STEP4/bbduk/clean_R2.fastq.gz \
        ref=$PHIX_REF k=31 hdist=1 \
        stats=$STEP4/bbduk/bbduk_stats.txt
    echo "  bbduk PhiX done"
else
    echo "  bbduk PhiX output exists, skipping"
fi

# ── Step 5: Deduplication ─────────────────────────────────────────────────
echo ""
echo "=== Step 5: Deduplication ==="

STEP5=$QC_DIR/step5_dedup
mkdir -p $STEP5/{fastp,clumpify}

# fastp --dedup
echo "  Testing: fastp --dedup"
FASTP_DEDUP_CLEAN=$STEP5/fastp/clean_R1.fastq.gz
if [ ! -f "$FASTP_DEDUP_CLEAN" ]; then
    $FASTP \
        -i $PREP_R1 -I $PREP_R2 \
        -o $FASTP_DEDUP_CLEAN -O $STEP5/fastp/clean_R2.fastq.gz \
        --dedup \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        -j $STEP5/fastp/fastp_dedup.json \
        -h $STEP5/fastp/fastp_dedup.html \
        -w $THREADS
    echo "  fastp dedup done"
else
    echo "  fastp dedup output exists, skipping"
fi

# clumpify
echo "  Testing: clumpify"
CLUMP_CLEAN=$STEP5/clumpify/clean_R1.fastq.gz
if [ ! -f "$CLUMP_CLEAN" ]; then
    clumpify.sh -Xmx8g \
        in1=$PREP_R1 in2=$PREP_R2 \
        out1=$CLUMP_CLEAN out2=$STEP5/clumpify/clean_R2.fastq.gz \
        dedupe=t optical=f
    echo "  clumpify done"
else
    echo "  clumpify output exists, skipping"
fi

# ── Step 6: Low-complexity filtering ─────────────────────────────────────
echo ""
echo "=== Step 6: Low-complexity filtering ==="

STEP6=$QC_DIR/step6_complexity
mkdir -p $STEP6/bbduk

# bbduk entropy
echo "  Testing: bbduk (entropy)"
BBDUK_LC_CLEAN=$STEP6/bbduk/clean_R1.fastq.gz
if [ ! -f "$BBDUK_LC_CLEAN" ]; then
    bbduk.sh -Xmx8g in1=$PREP_R1 in2=$PREP_R2 \
        out1=$BBDUK_LC_CLEAN out2=$STEP6/bbduk/clean_R2.fastq.gz \
        entropy=0.5 entropywindow=50 \
        stats=$STEP6/bbduk/bbduk_stats.txt
    echo "  bbduk entropy done"
else
    echo "  bbduk entropy output exists, skipping"
fi

# ── Score all QC steps ───────────────────────────────────────────────────
echo ""
echo "=== Scoring QC steps with viroforge benchmark qc ==="

python3 << 'SCORE_QC_EOF'
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/sahlab/leran/viroforge")
from viroforge.benchmarking.parsers import read_labels, read_names
from viroforge.benchmarking.qc import benchmark_qc

RAW_R1 = Path("/scratch/sahlab/leran/viroforge/data/gut_10x/fastq/gut_virome_-_adult_healthy_western_diet_R1.fastq")
QC_DIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/qc")
RESULTS = QC_DIR / "results"

raw = read_labels(RAW_R1)
print(f"Raw reads loaded: {len(raw):,}")

# Step-specific policies
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
    "step4_phix": {
        "viral": "keep", "dark_matter": "keep", "erv_exogenous": "keep",
        "host_dna": "keep", "host_dna_herv": "keep",
        "rrna": "keep",
        "phix": "remove",
        "reagent_bacteria": "keep",
        "bacterial_background": "keep", "fungal_background": "keep",
        "archaeal_background": "keep", "artifact_low_complexity": "keep",
        "erv_endogenous": "keep", "dietary_background": "keep",
    },
    "step5_dedup": None,  # uses default policy
    "step6_complexity": {
        "viral": "keep", "dark_matter": "keep", "erv_exogenous": "keep",
        "host_dna": "keep", "host_dna_herv": "keep",
        "rrna": "keep", "phix": "keep", "reagent_bacteria": "keep",
        "bacterial_background": "keep", "fungal_background": "keep",
        "archaeal_background": "keep",
        "artifact_low_complexity": "remove",
        "erv_endogenous": "keep", "dietary_background": "keep",
    },
}

# Tool -> cleaned R1 file paths
TOOLS = {
    "step2_host": {
        "bowtie2": QC_DIR / "step2_host/bowtie2/clean_R1.fastq.gz",
        "bwa_mem2": QC_DIR / "step2_host/bwa_mem2/clean_R1.fastq.gz",
    },
    "step3_rrna": {
        "sortmerna": QC_DIR / "step3_rrna/sortmerna/clean_fwd.fq.gz",
        "bbduk": QC_DIR / "step3_rrna/bbduk/clean_R1.fastq.gz",
        "ribodetector": QC_DIR / "step3_rrna/ribodetector/clean_R1.fastq.gz",
    },
    "step4_phix": {
        "bowtie2": QC_DIR / "step4_phix/bowtie2/clean_R1.fastq.gz",
        "bbduk": QC_DIR / "step4_phix/bbduk/clean_R1.fastq.gz",
    },
    "step5_dedup": {
        "fastp_dedup": QC_DIR / "step5_dedup/fastp/clean_R1.fastq.gz",
        "clumpify": QC_DIR / "step5_dedup/clumpify/clean_R1.fastq.gz",
    },
    "step6_complexity": {
        "bbduk": QC_DIR / "step6_complexity/bbduk/clean_R1.fastq.gz",
    },
}

all_results = {}
for step, tools in TOOLS.items():
    policy = POLICIES[step]
    step_results = {}
    for tool_name, clean_path in tools.items():
        if not clean_path.exists():
            print(f"  SKIP: {step}/{tool_name} — output not found")
            continue
        print(f"  Scoring: {step}/{tool_name}")
        kept = read_names(clean_path)
        metrics = benchmark_qc(raw, kept, keep_remove=policy)
        out_file = RESULTS / f"{step}_{tool_name}.json"
        with open(out_file, "w") as f:
            json.dump(metrics, f, indent=2, default=str)
        vr = metrics.get("viral_retention")
        rc = metrics["contamination"].get("recall")
        pr = metrics["contamination"].get("precision")
        f1 = metrics["contamination"].get("f1")
        print(f"    viral_retention={vr:.4f}  recall={rc:.4f}  precision={pr:.4f}  f1={f1:.4f}" if all(v is not None for v in [vr,rc,pr,f1]) else f"    metrics incomplete")
        step_results[tool_name] = metrics
    all_results[step] = step_results

with open(QC_DIR / "all_qc_benchmarks.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

# Print summary table
print("\n" + "=" * 100)
print("QC BENCHMARK RESULTS — GUT 10x")
print("=" * 100)
for step, tools in all_results.items():
    if not tools:
        continue
    print(f"\n── {step} ──")
    print(f"  {'Tool':<18s}  {'Viral Ret':>10s}  {'Recall':>10s}  {'Precision':>10s}  {'F1':>10s}")
    print(f"  {'─'*18}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}")
    for tool, m in tools.items():
        def fmt(v): return f"{v:.4f}" if v is not None else "N/A"
        print(f"  {tool:<18s}  {fmt(m.get('viral_retention')):>10s}  "
              f"{fmt(m['contamination'].get('recall')):>10s}  "
              f"{fmt(m['contamination'].get('precision')):>10s}  "
              f"{fmt(m['contamination'].get('f1')):>10s}")
        if "dedup" in m and m["dedup"]["n_duplicates"]:
            dd = m["dedup"]
            print(f"    └─ dedup: {dd['n_duplicates']:,} duplicates, "
                  f"removal rate: {fmt(dd['dedup_rate'])}")

SCORE_QC_EOF

echo "QC Benchmark scoring complete: $(date)"

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: ASSEMBLY BENCHMARK
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "================================================================"
echo "PART 2: ASSEMBLY BENCHMARK"
echo "================================================================"

ASM_DIR=$WORKDIR/assembly

# MEGAHIT
echo "=== Running MEGAHIT ==="
MEGAHIT_DIR=$ASM_DIR/megahit
if [ ! -f "$MEGAHIT_DIR/final.contigs.fa" ]; then
    rm -rf $MEGAHIT_DIR 2>/dev/null  # megahit needs empty output dir
    megahit -1 $RAW_R1 -2 $RAW_R2 \
        -o $MEGAHIT_DIR \
        -t $THREADS \
        --min-contig-len 500
    echo "MEGAHIT done: $(grep -c '^>' $MEGAHIT_DIR/final.contigs.fa) contigs"
else
    echo "MEGAHIT output exists, skipping"
fi

# metaSPAdes
echo ""
echo "=== Running metaSPAdes ==="
SPADES_DIR=$ASM_DIR/metaspades
if [ ! -f "$SPADES_DIR/contigs.fasta" ]; then
    spades.py --meta \
        -1 $RAW_R1 -2 $RAW_R2 \
        -o $SPADES_DIR \
        -t $THREADS \
        --memory 64 \
        --only-assembler
    echo "metaSPAdes done: $(grep -c '^>' $SPADES_DIR/contigs.fasta) contigs"
else
    echo "metaSPAdes output exists, skipping"
fi

# Score assemblies
echo ""
echo "=== Scoring assemblies ==="

python3 << 'SCORE_ASM_EOF'
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/sahlab/leran/viroforge")

DATASET = Path("/scratch/sahlab/leran/viroforge/data/gut_10x")
ASM_DIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/assembly")
RESULTS = ASM_DIR / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

GENOME_FASTA = DATASET / "fasta/gut_virome_-_adult_healthy_western_diet.fasta"
META_JSON = DATASET / "metadata/gut_virome_-_adult_healthy_western_diet_metadata.json"

with open(META_JSON) as f:
    metadata = json.load(f)

from viroforge.benchmarking.assembly import benchmark_assembly

assemblies = {
    "megahit": ASM_DIR / "megahit/final.contigs.fa",
    "metaspades": ASM_DIR / "metaspades/contigs.fasta",
}

all_results = {}
for name, contigs_path in assemblies.items():
    if not contigs_path.exists():
        print(f"  SKIP: {name} — contigs not found")
        continue
    print(f"  Scoring: {name}")
    metrics = benchmark_assembly(
        contigs_fasta=str(contigs_path),
        reference_fasta=str(GENOME_FASTA),
        metadata=metadata,
    )
    out_file = RESULTS / f"assembly_{name}.json"
    with open(out_file, "w") as f:
        json.dump(metrics, f, indent=2, default=str)

    gs = metrics.get("genome_stats", {})
    print(f"    Contigs: {metrics.get('n_contigs', 'N/A')}")
    print(f"    Recovery rate: {gs.get('recovery_rate', 'N/A')}")
    print(f"    N50: {metrics.get('n50', 'N/A')}")
    print(f"    Chimeric: {metrics.get('n_chimeric', 'N/A')}")
    all_results[name] = metrics

with open(ASM_DIR / "all_assembly_benchmarks.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

print("\nAssembly benchmark scoring complete")
SCORE_ASM_EOF

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: TAXONOMY BENCHMARK
# ══════════════════════════════════════════════════════════════════════════════

echo ""
echo "================================================================"
echo "PART 3: TAXONOMY BENCHMARK"
echo "================================================================"

TAX_DIR=$WORKDIR/taxonomy

# Kraken2 standard DB
echo "=== Running Kraken2 (standard DB) ==="
K2_STD_OUT=$TAX_DIR/kraken2_standard.out
K2_STD_RPT=$TAX_DIR/kraken2_standard.report
if [ ! -f "$K2_STD_OUT" ] || [ ! -s "$K2_STD_OUT" ]; then
    $KRAKEN2_BIN \
        --db $KRAKEN2_STANDARD_DB \
        --paired $RAW_R1 $RAW_R2 \
        --output $K2_STD_OUT \
        --report $K2_STD_RPT \
        --threads $THREADS
    echo "Kraken2 standard done"
else
    echo "Kraken2 standard output exists, skipping"
fi

# Kraken2 viral DB
echo ""
echo "=== Running Kraken2 (viral DB) ==="
K2_VIR_OUT=$TAX_DIR/kraken2_viral.out
K2_VIR_RPT=$TAX_DIR/kraken2_viral.report
if [ ! -f "$K2_VIR_OUT" ] || [ ! -s "$K2_VIR_OUT" ]; then
    $KRAKEN2_BIN \
        --db $KRAKEN2_VIRAL_DB \
        --paired $RAW_R1 $RAW_R2 \
        --output $K2_VIR_OUT \
        --report $K2_VIR_RPT \
        --threads $THREADS
    echo "Kraken2 viral done"
else
    echo "Kraken2 viral output exists, skipping"
fi

# Kaiju
echo ""
echo "=== Running Kaiju ==="
KAIJU_OUT=$TAX_DIR/kaiju.out
if [ ! -f "$KAIJU_OUT" ] || [ ! -s "$KAIJU_OUT" ]; then
    $KAIJU_BIN \
        -t $KAIJU_NODES \
        -f $KAIJU_DB \
        -i $RAW_R1 -j $RAW_R2 \
        -o $KAIJU_OUT \
        -z $THREADS \
        -v
    echo "Kaiju done"
else
    echo "Kaiju output exists, skipping"
fi

# Score taxonomy
echo ""
echo "=== Scoring taxonomy ==="

python3 << 'SCORE_TAX_EOF'
import json
import sys
from pathlib import Path

sys.path.insert(0, "/scratch/sahlab/leran/viroforge")

DATASET = Path("/scratch/sahlab/leran/viroforge/data/gut_10x")
TAX_DIR = Path("/scratch/sahlab/leran/paper_benchmark/gut_10x/taxonomy")
RESULTS = TAX_DIR / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

TAXDUMP_DIR = "/ref/sahlab/data/kraken2_db/kraken2_standard/taxonomy"

META_JSON = DATASET / "metadata/gut_virome_-_adult_healthy_western_diet_metadata.json"
with open(META_JSON) as f:
    metadata = json.load(f)

taxonomy_gt = metadata.get("benchmarking", {}).get("taxonomy", {})
if not taxonomy_gt:
    print("ERROR: No benchmarking.taxonomy in metadata")
    sys.exit(1)

from viroforge.benchmarking.taxonomy import benchmark_taxonomy, detect_format, PARSERS
from viroforge.benchmarking.ncbi_tree import NcbiTree

ncbi_tree = NcbiTree.from_dir(TAXDUMP_DIR)
print(f"NCBI tree loaded: {len(ncbi_tree._parents)} taxa")

classifiers = {
    "kraken2_standard": TAX_DIR / "kraken2_standard.out",
    "kraken2_viral": TAX_DIR / "kraken2_viral.out",
    "kaiju": TAX_DIR / "kaiju.out",
}

all_results = {}
for name, output_path in classifiers.items():
    if not output_path.exists() or output_path.stat().st_size == 0:
        print(f"  SKIP: {name} — output not found")
        continue

    print(f"  Scoring: {name}")
    fmt = detect_format(str(output_path))
    print(f"    Format detected: {fmt}")
    assignments = PARSERS[fmt](str(output_path))
    print(f"    Parsed {len(assignments)} read classifications")

    metrics = benchmark_taxonomy(
        assignments=assignments,
        taxonomy_gt=taxonomy_gt,
        ncbi_tree=ncbi_tree,
    )

    out_file = RESULTS / f"taxonomy_{name}.json"
    with open(out_file, "w") as f:
        json.dump(metrics, f, indent=2, default=str)

    known = metrics.get("known_viruses", {})
    dark = metrics.get("dark_matter", {})
    per_rank = metrics.get("per_rank", {})
    ab = metrics.get("abundance_profile", {})

    print(f"    Known: sensitivity={known.get('sensitivity', 0):.4f}, precision={known.get('precision', 0):.4f}")
    print(f"    Dark matter: unclassified={dark.get('unclassified_rate', 0):.4f}")
    for rank in ["species", "genus", "family"]:
        r = per_rank.get(rank, {})
        print(f"    {rank}: recall={r.get('recall', 0):.4f}, precision={r.get('precision', 0):.4f}, f1={r.get('f1', 0):.4f}")
    print(f"    Abundance: BC={ab.get('bray_curtis', 0):.4f}, Spearman={ab.get('spearman', 0):.4f}")

    all_results[name] = metrics

with open(TAX_DIR / "all_taxonomy_benchmarks.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

print("\nTaxonomy benchmark scoring complete")
SCORE_TAX_EOF

echo ""
echo "================================================================"
echo "ALL BENCHMARKS COMPLETE: $(date)"
echo "================================================================"
echo "Results in: $WORKDIR"
echo "  QC:       $WORKDIR/qc/all_qc_benchmarks.json"
echo "  Assembly: $WORKDIR/assembly/all_assembly_benchmarks.json"
echo "  Taxonomy: $WORKDIR/taxonomy/all_taxonomy_benchmarks.json"
