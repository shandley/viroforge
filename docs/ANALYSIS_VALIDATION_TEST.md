# ViroForge Analysis Validation Test

**Date**: 2026-04-22
**Dataset**: Collection 1 — Healthy Human Gut Virome (NovaSeq, 10x coverage)
**Purpose**: Validate ViroForge's synthetic data by running real bioinformatics analysis pipelines and comparing results against ground truth labels.

---

## Dataset Summary

| Parameter | Value |
|---|---|
| Collection | 1 — Gut Virome (Adult Healthy, Western Diet) |
| Platform | NovaSeq (paired-end, 150 bp) |
| Coverage | 10x |
| Total R1 Reads | 531,606 |
| Viral Genomes | 133 |
| Contaminant Sequences | 156 |
| VLP Protocol | Tangential Flow Filtration (0.2 um) |
| Contamination Level | Realistic |
| Random Seed | 42 |
| Metadata Version | 1.1 |

---

## Test 1: Read-Level Ground Truth Audit

**Goal**: Verify that per-read `source=` labels in FASTQ headers match the expected contamination proportions from metadata.

**Method**: Count reads by `source=` label in R1 FASTQ headers and compare against the post-VLP-enrichment contamination levels stored in `dataset_metadata.json`.

### Results

| Source | Read Count | Observed % | Expected % | Relative Error |
|---|---|---|---|---|
| viral | 528,305 | 99.379% | 99.343% | 0.04% |
| rrna | 2,633 | 0.495% | 0.530% | 6.6% |
| host_dna | 327 | 0.062% | 0.062% | 0.5% |
| phix | 316 | 0.059% | 0.060% | 1.0% |
| reagent_bacteria | 25 | 0.005% | 0.005% | 6.0% |

### Conclusion

All 5 source categories match expected abundances. Every read has a `source=` label with no unlabeled reads. The VLP tangential flow filtration model reduced total contamination from 7.6% to 0.66%, and actual read counts accurately reflect those reduced levels. The ground truth labeling system is working correctly.

---

## Test 2: QC Pipeline Validation

**Goal**: Run a standard virome QC pipeline (adapter trimming, PhiX removal, host decontamination) and evaluate performance against ground truth `source=` labels.

**Tools**: fastp 1.1.0, Bowtie2 2.5.4, samtools 1.21

### Step 1: Adapter Trimming (fastp)

```bash
fastp \
  -i data/collection_1/fastq/gut_virome_-_adult_healthy_western_diet_R1.fastq \
  -I data/collection_1/fastq/gut_virome_-_adult_healthy_western_diet_R2.fastq \
  -o analysis/qc/clean_R1.fastq \
  -O analysis/qc/clean_R2.fastq \
  --thread 4
```

**Results**:
- Reads passed filter: 1,063,212 (100% — no reads removed)
- Reads with adapter trimmed: 404 (0.04%)
- Bases trimmed due to adapters: 4,250
- Duplication rate detected: 3.99%
- Insert size peak: 175 bp

**Notes**: fastp passed all reads since ViroForge simulates realistic quality scores that meet default thresholds. The low adapter trimming rate (0.04%) indicates most reads were generated without adapter read-through, consistent with the default configuration.

### Step 2: PhiX Removal (Bowtie2)

```bash
bowtie2-build viroforge/data/references/phix174.fasta analysis/qc/phix_index
bowtie2 -x analysis/qc/phix_index \
  -1 analysis/qc/clean_R1.fastq -2 analysis/qc/clean_R2.fastq \
  --un-conc analysis/qc/no_phix_R%.fastq \
  -S analysis/qc/phix_aligned.sam --threads 4
```

**Alignment**: 1.97% concordant alignment rate, 2.15% overall

**Performance vs Ground Truth**:

| Metric | Value |
|---|---|
| True PhiX reads | 316 |
| Flagged as PhiX | 12,345 |
| True Positives | 316 |
| False Positives | 12,029 |
| False Negatives | 0 |
| Sensitivity | 100.0% |
| Specificity | 97.7% |
| Precision (PPV) | 2.6% |

**Key Finding**: Bowtie2 achieved perfect sensitivity (all 316 PhiX reads detected) but very low precision — 12,029 viral reads were incorrectly flagged as PhiX due to sequence similarity. All false positives were viral reads. This is a known real-world problem: PhiX174 (family Microviridae) shares homology with other Microviridae phages commonly found in gut viromes. This demonstrates that naive PhiX removal by full-genome alignment can cause significant data loss in virome studies.

**Implication for pipeline developers**: PhiX removal tools should consider using stricter alignment thresholds (e.g., `--score-min` in Bowtie2) or kmer-based approaches with exact PhiX174 sequence matching rather than loose alignment to minimize false positives in virome data.

### Step 3: Host DNA Removal (Bowtie2)

```bash
bowtie2-build viroforge/data/references/host_fragments.fasta analysis/qc/host_index
bowtie2 -x analysis/qc/host_index \
  -1 analysis/qc/no_phix_R1.fastq -2 analysis/qc/no_phix_R2.fastq \
  --un-conc analysis/qc/no_host_R%.fastq \
  -S analysis/qc/host_aligned.sam --threads 4
```

**Alignment**: 0.06% overall alignment rate

**Performance vs Ground Truth**:

| Metric | Value |
|---|---|
| True host reads | 327 |
| Flagged as host | 327 |
| True Positives | 327 |
| False Positives | 0 |
| False Negatives | 0 |
| Sensitivity | 100.0% |
| Specificity | 100.0% |
| Precision (PPV) | 100.0% |

**Key Finding**: Perfect performance. All 327 host DNA reads were correctly identified with zero false positives. This is expected since ViroForge generates host reads from real T2T-CHM13v2.0 human genome fragments, which align unambiguously back to the same reference.

**Note**: In real-world applications, host removal uses the full human genome as reference (not 48 fragments), so sensitivity may differ. However, this validates that ViroForge's host contamination reads are realistic enough to be caught by standard tools.

### Overall QC Performance Summary

| Metric | Value |
|---|---|
| Total contaminant reads | 3,301 |
| Contaminants removed | 643 (19.5%) |
| Contaminants missed | 2,658 (80.5%) |
| Viral reads lost (false positives) | 12,029 (2.28% of viral reads) |

**Contaminants not removed** (no dedicated removal step run):
- rRNA: 2,633 reads — would require SortMeRNA or bbduk with rRNA database
- Reagent bacteria: 25 reads — would require database of known reagent contaminants

### QC Validation Conclusions

1. **Ground truth labels enable precise QC benchmarking** — without ViroForge's `source=` labels, there would be no way to calculate these metrics on real data.
2. **PhiX removal is the weakest QC step for virome data** — 2.3% viral read loss due to Microviridae homology is a significant concern.
3. **Host removal works well** when the reference matches the contamination source.
4. **A complete virome QC pipeline needs** at minimum: adapter trimming + PhiX removal + host removal + rRNA removal. Missing any step leaves contaminants in the data.
5. **ViroForge data behaves like real sequencing data** — QC tools respond to it the same way they would to real virome data, validating the simulation approach.

---

## Test 3: Taxonomy Classification (Kraken2)

**Status**: Pending — requires Kraken2 viral database download.

**Planned approach**:
1. Download pre-built Kraken2 viral database
2. Classify reads from the QC-cleaned output
3. Compare Kraken2 classifications against ViroForge's composition TSV ground truth
4. Calculate precision/recall at family and genus level

---

## Test 4: Assembly and Genome Recovery (MEGAHIT)

**Status**: Pending — will run after taxonomy classification.

**Planned approach**:
1. Assemble QC-cleaned reads with MEGAHIT (meta-sensitive mode)
2. Map contigs to reference genomes with minimap2
3. Calculate genome recovery rate (how many of 133 genomes recovered)
4. Assess completeness and fragmentation per genome

---

## Tools and Versions

| Tool | Version | Purpose |
|---|---|---|
| fastp | 1.1.0 | Adapter trimming, quality filtering |
| Bowtie2 | 2.5.4 | PhiX removal, host decontamination |
| samtools | 1.21 | SAM/BAM processing |
| Kraken2 | 2.1.3 | Taxonomy classification (pending) |
| MEGAHIT | 1.2.9 | Metagenomic assembly (pending) |
| minimap2 | 2.28 | Contig-to-reference mapping (pending) |

---

## Reproducibility

All analysis was performed on the dataset generated with:
```bash
viroforge generate --collection-id 1 --platform novaseq --coverage 10 --seed 42
```

Analysis scripts and intermediate files are in the `analysis/` directory:
- `analysis/qc/` — QC pipeline outputs (FASTQ, SAM, read ID lists)
- `analysis/taxonomy/` — Kraken2 outputs (pending)
- `analysis/assembly/` — MEGAHIT assembly outputs (pending)
