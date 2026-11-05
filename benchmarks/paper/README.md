# Paper Benchmark: Gut Virome 10x Coverage

This directory contains the scripts used to generate and benchmark the gut virome
dataset reported in the ViroForge paper. All benchmarks were run on the WashU HTCF
cluster.

## Dataset generation

The 10x gut virome dataset was generated with ViroForge v0.22.0 on the
`release/clean-v0.22.0` branch:

```bash
viroforge generate \
    --collection-id 1 \
    --output data/gut_10x \
    --platform novaseq \
    --coverage 10 \
    --vlp-protocol tangential_flow \
    --amplification linker \
    --contamination-level realistic \
    --read-length 150 \
    --insert-size 350 \
    --seed 42
```

**Dataset location on HTCF:** `/scratch/sahlab/leran/viroforge/data/gut_10x/`

**Dataset properties:**
- Collection 1 (Gut Virome - Adult Healthy, Western Diet), 133 known + 66 dark matter = 199 genomes
- NovaSeq 6000, 150 bp paired-end, 350 bp insert size
- Linker amplification, tangential flow filtration (0.2 um), realistic contamination
- ~747,929 read pairs after fastp preprocessing

## Benchmark pipeline scripts

### `run_10x_benchmark.sh`

Main SLURM script that runs the complete benchmark pipeline in three parts:

1. **QC benchmark** (Part 1): fastp preprocessing, then five QC steps each tested
   with multiple tools:
   - Host removal: Bowtie2 (ViroForge bundled refs), BWA-MEM2
   - rRNA removal: SortMeRNA, BBDuk, RiboDetector (all with ViroForge bundled refs)
   - PhiX removal: Bowtie2, BBDuk
   - Deduplication: fastp --dedup, Clumpify
   - Complexity filtering: BBDuk entropy

2. **Assembly benchmark** (Part 2): MEGAHIT and metaSPAdes, scored with
   `viroforge benchmark assembly`

3. **Taxonomy benchmark** (Part 3): Kraken2 (standard DB, viral DB) and Kaiju
   (viral DB), scored with `viroforge benchmark taxonomy`

```bash
sbatch run_10x_benchmark.sh
```

### `rerun_realworld_v3.sh`

Re-runs QC and taxonomy benchmarks with full-size real-world databases:
- Host removal: Bowtie2 with GRCh38 (full human genome)
- rRNA removal: Bowtie2 with SILVA 128
- Taxonomy: Kaiju with RefSeq virus protein DB

These results are reported alongside the ViroForge-bundled-reference results in
the paper (Tables 1a, 1c).

```bash
sbatch rerun_realworld_v3.sh
```

## Ground-truth analysis scripts

These scripts trace benchmark results back to per-read source labels to produce
the findings in Sections 3.5-3.7 of the paper. They run on the cluster against
the benchmark output files.

### `analyze_qc_fp_fn.py`

Traces every QC false positive and false negative back to its source label and
genome of origin. Produces the data for Table 3 and Section 3.5 (per-read source
label analysis of QC errors).

### `analyze_assembly_ground_truth.py`

Classifies chimeric contigs by domain (viral-viral, bacterial-bacterial,
viral-bacterial) and analyzes missing genomes by family, size, and GC content.
Produces the data for Table 4 and Section 3.6.

### `analyze_taxonomy_by_source.py`

Joins taxonomy classifier output with per-read source labels to compute
classification rates by source label. Produces the data for Table 5 and
Section 3.7. Identifies genome-specific false negatives and taxonomic ambiguity.

### `analyze_taxonomy_bacterial_fp.py`

Deep-dives into bacterial background false positives in viral-only databases.
Compares false-positive rates between curated challenge regions (prophages,
CRISPR, ICEs) and generic bacterial fragments.

## Cluster paths

| Resource | Path |
|----------|------|
| Generated dataset | `/scratch/sahlab/leran/viroforge/data/gut_10x/` |
| Benchmark results | `/scratch/sahlab/leran/paper_benchmark/gut_10x/` |
| QC results JSON | `.../gut_10x/qc/all_qc_benchmarks_realworld.json` |
| Assembly results JSON | `.../gut_10x/assembly/all_assembly_benchmarks.json` |
| Taxonomy results JSON | `.../gut_10x/taxonomy/all_taxonomy_benchmarks_realworld.json` |
| Kraken2 standard DB | `/ref/sahlab/data/kraken2_db/kraken2_standard` |
| Kraken2 viral DB | `/ref/sahlab/data/viral_analysis_DBs/kraken_viral_Nucl_DB` |
| Kaiju viral DB | `/ref/gdlab/data/kaiju_db/RefSeq_virus_prot/kaiju_db_viruses.fmi` |
| GRCh38 index | `/scratch/sahlab/leran/paper_benchmark/indexes/GRCh38_noalt_as` |
| SILVA 128 index | `.../indexes/SILVA_128_LSUParc_SSUParc_ribosomal_RNA` |

## Reproducing from scratch

1. Install ViroForge: `pip install -e ".[benchmark]"` from the repo root
2. Set up the database: `viroforge setup-db`
3. Generate the dataset with the command above
4. Adapt cluster paths in the SLURM scripts to your environment
5. Submit: `sbatch run_10x_benchmark.sh`, then `sbatch rerun_realworld_v3.sh`
6. Run analysis scripts against the output
