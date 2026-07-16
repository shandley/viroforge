# QC Benchmarking (Phase 13B, Module 1)

Validate a QC pipeline's contamination removal and viral retention against
ViroForge per-read ground truth.

## Idea

ViroForge writes a `source=` label on every read header (viral, dark_matter,
host_dna, rrna, phix, reagent_bacteria, artifact_low_complexity). Run any QC tool
(fastp, bbduk, host removal, complexity/dedup filters) on a ViroForge dataset,
then compare the tool's surviving reads against the raw reads to measure how much
contamination it removed and, critically, how many viral reads it kept.

The per-read labels were independently verified accurate in
`docs/DATA_QUALITY_EVALUATION.md`, so they are trustworthy ground truth.

## Usage

```bash
viroforge benchmark qc \
  --raw-reads   data/gut/fastq/gut_R1.fastq \
  --cleaned-reads results/qc/gut_cleaned_R1.fastq.gz \
  --output  reports/qc_benchmark.json \
  --markdown reports/qc_benchmark.md
```

- `--raw-reads`: the raw ViroForge R1 FASTQ (the ground truth). Plain or `.gz`.
- `--cleaned-reads`: the QC tool's R1 output (the surviving reads). Plain or `.gz`.
- `--output` / `--markdown`: optional JSON and markdown reports. The markdown
  summary is also printed to stdout.
- `--keep-remove SOURCE:keep|remove ...`: override the default keep/remove policy.

Exit code is non-zero when the results are unreliable (see match-rate gate below).

## Metrics

The report leads with the metric that matters most for low-biomass viromes:

- Viral retention: fraction of viral reads the QC tool kept.
- Over-filtering rate: fraction of viral reads wrongly removed.
- Removal rate per contaminant type (host_dna, rrna, phix, reagent_bacteria,
  artifact_low_complexity).
- Aggregate confusion matrix (positive class = should-be-removed) with precision,
  recall, F1. For virome samples with few contaminant reads, aggregate precision
  is noisy; use viral retention and per-type removal as the primary read.

## Ground-truth semantics

- Match-rate gate: read names in the cleaned file are matched back to the raw
  reads. If under 99 percent match, the report is flagged unreliable and the exit
  code is non-zero. This catches tools that rename or reorder reads and wrong
  raw/cleaned file pairings. A trailing `/1` or `/2` mate suffix is tolerated;
  accession versions (`.1`) and duplicate suffixes (`_dupN`) are never stripped.
- PCR duplicates are orthogonal to contamination. A duplicate inherits its
  template's source (a duplicated host read is still host), so contamination
  metrics are computed over non-duplicate reads only, and duplicate removal is
  reported separately as a dedup rate.
- Artifact reads carry two labels (`source=viral source=artifact_low_complexity`);
  the benchmark uses the last one, the read's current nature, which is what QC
  should act on.
- Keep/remove policy is explicit and overridable. Defaults: viral, dark_matter,
  and erv_exogenous are kept; host_dna, rrna, phix, reagent_bacteria,
  artifact_low_complexity, and erv_endogenous are removed.
- Scope (v1): measures read removal, R1 only (each fragment counted once, mates
  share the label). It does not measure trimming quality; adapter read-through
  reads keep `source=viral` and correctly count as keep-and-trim.

## Testing

`tests/test_benchmark_qc.py` builds a cleaned set with fixed, recorded drop counts
per class and asserts the engine recovers the hand-computed rates, so the metrics
are checked against an independent oracle rather than the benchmark's own logic.
