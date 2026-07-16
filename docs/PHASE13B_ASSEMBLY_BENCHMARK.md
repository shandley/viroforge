# Assembly Benchmarking (Phase 13B, Module 2)

Validate a metagenome assembly against the true ViroForge genomes: how much of
each genome was recovered, at what identity, whether any contigs are chimeric,
assembly contiguity, and how accurately per-genome abundance can be estimated.

## Idea

Assemble a ViroForge dataset with any assembler (SPAdes, MEGAHIT), then align the
resulting contigs to the true genomes and score recovery against the per-genome
expectations ViroForge already records. Viral genomes are small, so complete
recovery is realistic; this module measures how close an assembly gets.

Requires the `benchmark` extra (minimap2 via mappy): `pip install -e ".[benchmark]"`.

## Usage

```bash
viroforge benchmark assembly \
  --contigs results/assembly/contigs.fasta \
  --genomes data/gut/fasta/gut.fasta \
  --ground-truth data/gut/metadata/gut_metadata.json \
  --output reports/assembly.json --markdown reports/assembly.md
```

- `--contigs`: the assembler's output.
- `--genomes`: the ViroForge source genome FASTA for the dataset (`fasta/*.fasta`).
- `--ground-truth`: the dataset metadata JSON (per-genome length, type, abundance,
  and expected completeness). Optional; without it, genome length/type/abundance
  are read from the FASTA headers and the expected-completeness comparison is skipped.

## Metrics

- Genome recovery: per-genome completeness (fraction of the genome covered by
  contigs, from merged alignment intervals), bucketed complete (>=95%) /
  high_quality (>=75%) / partial (>=50%) / fragmented / missing; recovery rate;
  mean completeness; per-genome identity.
- Observed vs expected completeness: compares recovery to what is achievable at
  the dataset's coverage (ViroForge's Lander-Waterman `expected_completeness`), so
  a genome sequenced too shallowly to assemble is not scored as a failure.
- Chimeras: contigs where two or more different genomes each cover a distinct,
  largely non-overlapping segment of the contig.
- Contiguity: N50, L50, longest contig, total bp.
- Abundance accuracy (when per-contig coverage is available): see below.

## Abundance accuracy

Per-contig coverage is read from SPAdes headers (`..._cov_45.2`) or MEGAHIT
headers (`multi=45.2`). Coverage is sequencing depth, which for a genome of
relative abundance a and length L is proportional to a/L (reads track abundance;
depth divides by length). So per-genome abundance is estimated as depth times
length, aggregated over the contigs assigned to that genome, then normalized and
compared to the true relative abundance (Spearman and mean absolute error).

Abundance accuracy is bounded by recovery: a genome that is missing or only
partially assembled cannot be quantified accurately, so this metric reflects both
the quantification and the assembly quality together.

## Scope and thresholds

- v1 covers genome recovery, chimeras, contiguity, and abundance. Coverage
  uniformity / GC-bias plots and HTML visualizations are deferred.
- Alignment uses minimap2 `asm5` (contigs are near-identical to their source
  genome). Alignments below 90% identity are dropped as spurious. A chimera
  requires each genome to cover at least 20% of the contig with less than 50%
  span overlap. These are constants at the top of `viroforge/benchmarking/assembly.py`.

## Testing

`tests/test_benchmark_assembly.py` carves contigs from synthetic genomes with
fixed known completeness (complete / partial / fragmented / missing) plus a
deliberate two-genome chimera, and asserts the engine recovers those exact
figures, so the metrics are checked against an independent oracle.
