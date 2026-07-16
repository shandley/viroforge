# Taxonomy Benchmarking (Module 4, read-based v1)

Compare a read classifier's per-read taxonomic assignments against ViroForge's
true per-genome taxonomy.

## Idea

ViroForge writes read names that embed the source genome accession
(`GCF_000819615.1_991_0/1`). A classifier (Kraken2, etc.) preserves that read
name and assigns each read an NCBI taxid. So the ground truth for each read is
recoverable from the read name plus the per-genome taxonomy exported in the
dataset metadata: read name -> genome accession -> true NCBI taxid. No raw reads
or database are needed at benchmark time.

The metadata now exports `benchmarking.taxonomy` = {genome_id: {ncbi_taxid,
realm..species, is_known}} for every viral genome (regenerate a dataset with a
current version to get it). `ncbi_taxid` is populated for all genomes, which is
what makes taxid-level comparison possible.

## Usage

```bash
# Kraken2 per-read output
viroforge benchmark taxonomy \
  --pipeline-output results/kraken2.out \
  --ground-truth data/gut/metadata/gut_metadata.json \
  --format kraken2 \
  --output reports/taxonomy.json --markdown reports/taxonomy.md

# Generic TSV (read_id <TAB> taxid)
viroforge benchmark taxonomy --pipeline-output calls.tsv \
  --ground-truth ... --format generic

# With per-rank (genus/family) metrics: supply the NCBI taxdump
viroforge benchmark taxonomy --pipeline-output results/kraken2.out \
  --ground-truth ... --taxdump-dir /path/to/taxdump
```

The taxdump is nodes.dmp + names.dmp from
ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz (Kraken2 users already have it).

## What it measures

Reads are stratified into two groups, which is the crux of scoring viromes fairly:

- Known viruses (ICTV family assigned, likely present in reference databases).
  Reported: sensitivity (correct / all), precision (correct / classified), and
  the unclassified and misclassified counts. This is the classifier's real
  performance on classifiable viral content.
- Dark matter (family Unknown / novel). Reported: the fraction correctly left
  unclassified (expected to be high, since a classifier should not force a call
  on novel content) and how many were classified anyway. A high unclassified rate
  here is correct behavior, not a failure, so it is never folded into a penalty.

Per-rank accuracy (species, genus, family) is reported when an NCBI taxdump is
supplied via `--taxdump-dir`. Both the true taxid and the classifier's assigned
taxid are resolved to their ancestor at each rank and compared, so a genus-level
call is credited as correct at genus and family (but not species). This is
computed over the known-virus stratum only, and a call too shallow to reach a
rank (or an unclassified read) counts against recall at that rank. Recall
therefore rises up the ranks (family easiest, species hardest), which is the
expected shape.

Plus an abundance profile comparison over NCBI taxids (Bray-Curtis dissimilarity,
Pearson, Spearman, mean absolute error) between the classifier's and the true
per-taxid read fractions.

## Scope and semantics

- Taxid-exact (always): a read is exactly correct only if the assigned taxid
  equals the true genome's species-level taxid.
- Per-rank (with `--taxdump-dir`): genus/family precision/recall/F1 that credit
  correct-but-higher-rank calls, resolved through the NCBI taxonomy tree.
  Without the taxdump, only taxid-exact is reported.
- Read-based only. Contig-based taxonomy is a later mode.
- Formats: Kraken2 per-read output and a generic `read_id\ttaxid` TSV. Kraken2
  `--use-names` output (`Name (taxid N)`) is parsed. Centrifuge, MMseqs2, and
  DIAMOND are not yet supported.
- Reads whose name does not map to a viral genome (contaminants) are counted
  separately and not scored. If no reads map to viral genomes the result is
  flagged unreliable.

## Testing

`tests/test_benchmark_taxonomy.py` runs assignments with a fixed, recorded mix of
correct / unclassified / misclassified across the known and dark strata and
asserts the engine recovers the hand-computed per-stratum metrics.
