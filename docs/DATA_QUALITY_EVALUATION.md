# ViroForge Data Quality Evaluation

**Date**: 2026-07-16
**Version evaluated**: v0.13.0 (canonical seeded database, commit ddd9d38 + working fixes)
**Scope**: Illumina short-read output, DNA and RNA workflows, both read mates.
Long-read (PacBio HiFi, Nanopore) not evaluated: pbccs is Linux x86-64 only and
PBSIM3 was not built in this environment.
**Tool**: `scripts/evaluate_dataset.py` (re-runnable)

## Executive summary

The ground truth that a QC benchmark actually consumes, the per-read `source=`
labels and the artifact tags, is accurate, traceable to the claimed source
genome, and reproducible. On that basis the data is high quality for QC pipeline
benchmarking and you can proceed.

The evaluation also found four fidelity issues that do not affect QC label
correctness but do matter for composition realism and for taxonomy or assembly
benchmarking later. The most important two are that the default 30 percent dark
matter is not actually delivered (realized 11 to 48 percent depending on
collection and seed), and that the three Illumina platforms produce byte for
byte identical reads, so platform selection currently does nothing for short
reads. Both have since been addressed: dark matter now delivers about 0.30, and
the platform behavior is documented. Details, evidence, and resolutions below.

## Method

Nine datasets were generated to span the configuration space, then evaluated by
re-deriving every quantity from the reads and the source FASTA rather than
trusting the metadata:

| dataset | collection | config |
|---|---|---|
| default_s42 / default_s42b | Mouse Gut (8) | v0.13.0 defaults, seed 42 (generated twice) |
| default_s7 | Mouse Gut (8) | defaults, seed 7 |
| mda_s42 | Mouse Gut (8) | `--amplification mda --mda-chimera-rate 0.02` |
| clean_s42 | Mouse Gut (8) | `--contamination-level clean`, all artifacts off |
| heavy_s42 | Mouse Gut (8) | `--contamination-level heavy` |
| miseq_s42 / hiseq_s42 | Mouse Gut (8) | `--platform miseq` / `hiseq` |
| gut_default | Gut Adult Healthy (1) | defaults, 60k reads (larger collection) |
| rna_s42 | Fecal RNA (15) | `--molecule-type rna --rna-depletion ribo_zero` |

Checks: k-mer containment (k=21, canonical, so strand-independent) of each read
against its claimed genome, stratified by source label with a per-label expected
match model; Shannon entropy of low-complexity reads versus normal reads;
duplicate copy-number distributions and template identity; MDA chimera two-genome
content; realized-versus-intended abundance (Spearman); realized-versus-claimed
contamination fraction; read length, quality, GC, N content, and pairing; and
byte-level reproducibility. Pass/fail thresholds were fixed before any numbers
were seen (see the constants at the top of the evaluation script).

## Verdict grid

All ten datasets pass all gated dimensions:

```
dataset        A1   A2   A3   D1   D2   B    E
default_s42    P    P    P    P    P    P    P
default_s42b   P    P    P    P    P    P    P
default_s7     P    P    P    P    P    P    P
mda_s42        P    P    P    P    P    P    P
clean_s42      P    P    P    P    P    P    P
heavy_s42      P    P    P    P    P    P    P
miseq_s42      P    P    P    P    P    P    P
hiseq_s42      P    P    P    P    P    P    P
gut_default    P    P    P    P    P    P    P
rna_s42        P    P    P    P    P    P    P
```

A1 label coverage, A2 label-versus-genome consistency, A3 k-mer traceability,
D1 low-complexity entropy, D2 duplicate identity, B viral abundance rank,
E pairing and format.

## What passed, with evidence

**Read-level labels and traceability (dimension A).** Every read carries a known
`source=` label (100 percent across all datasets). No read is labeled as a
contaminant type whose genome is viral, or vice versa (0 mismatches). k-mer
containment against the claimed genome is 1.00 (median) for every real source
type (viral, dark_matter, host_dna, rrna, phix, reagent_bacteria) and 0.00 for
`artifact_low_complexity` reads, which is exactly right: real reads contain their
genome's sequence, synthetic artifact reads contain nothing. This is the surface
a QC benchmark scores against, and it is correct and precise. It holds for both
read mates (R2 checked on default_s42: 100 percent labeled, median viral
containment 1.00) and for the RNA workflow with rRNA depletion (rna_s42:
viral, dark_matter, and rrna reads all trace at 1.00, artifacts at 0.00).

**Artifacts (dimension D).** Low-complexity reads have clearly lower base entropy
than normal reads (median about 1.4 versus 1.98, and as low as 0.44 when
homopolymers dominate the sample). Duplicate reads are identical to their
templates (486 of 486 sampled for MDA, similar for PCR). MDA chimeras were
verified to contain k-mers from two different genomes. The only duplicates that
do not match their template are ones whose template was itself chimerized after
the copy was made (chimera injection runs after duplicate injection); excluding
those, duplicate identity is exact. Duplicate copy-number distributions match the model (geometric and
capped at 5 for PCR, power-law out to 20 for MDA). Adapter read-through is
present and specific: 235 of 300 adapter-tagged reads contain the TruSeq adapter
(15 bp or more), versus 1 of 500 untagged reads; the roughly 65 not caught have
short read-through below the 15 bp detection floor, which is expected since
read-through length depends on insert size.

**Abundance fidelity (dimension B).** When measured over viral and dark-matter
genomes above the rare-genome floor, realized read counts track intended
abundance with Spearman 0.99 to 0.996. (Including enrichment-suppressed
contaminants and floored-tie genomes drags this down to about 0.85, which is a
measurement artifact, not a generation problem. The evaluation script scores
viral genomes only for this reason.)

**Format and reproducibility (dimensions E and F).** R1 and R2 counts always
match, reads are well formed, quality is Q30.5 mean. Same seed gives byte-for-byte
identical FASTQ (`default_s42` and `default_s42b` have identical md5 sums for R1
and R2), which is foundational for a benchmark generator.

## Findings

### 1. The default 30 percent dark matter is not delivered (high priority)

`--dark-matter-fraction 0.30` is documented and defaulted as "30 percent of the
viral community is dark matter." The realized fraction of viral reads is:

| dataset | collection | realized dark matter | target |
|---|---|---|---|
| default_s42 / _s42b | Mouse Gut | 0.114 | 0.30 |
| default_s7 | Mouse Gut (seed 7) | 0.185 | 0.30 |
| mda_s42 | Mouse Gut | 0.138 | 0.30 |
| gut_default | Gut Adult Healthy | 0.484 | 0.30 |

Root cause, established by controlled runs: the fraction is correct until
amplification bias is applied. With `--no-vlp --amplification none` the dark
matter share of viral mass is exactly 0.300. Adding the default `linker`
amplification drops it to 0.128, and VLP shaves it to 0.116. ISS then realizes
reads faithfully in proportion to that already-diluted abundance (realized reads
0.114 match FASTA mass 0.116; genome lengths are equal so length weighting is not
involved). Because amplification reweights by genome properties and the known
viral genomes differ per collection and per seed, the realized fraction is not
just low, it is uncontrolled: it ranges from 0.11 to 0.48 across the runs above.

The metadata always reports the configured target (`dark_matter_fraction: 0.30`),
so it overstates or understates the realized value with no indication.

Recommendation: apply the dark-matter fraction after amplification bias (or
renormalize the dark-matter and known-viral blocks to preserve the intended
ratio through amplification), and record the realized post-amplification fraction
in the metadata alongside the target. This is a composition-realism and
taxonomy-benchmarking issue, not a QC label issue.

Resolution (applied 2026-07-16): the dark-matter and known-viral blocks are now
renormalized to the configured ratio after VLP enrichment and amplification, with
total viral mass and each block's internal structure preserved. Realized dark
matter is now 0.295 to 0.301 across the collections and seeds above (was 0.11 to
0.48), and the metadata records `realized_fraction_of_viral`.

### 2. Illumina platform selection is a no-op for short reads (high priority)

NovaSeq, MiSeq, and HiSeq runs with the same seed produce byte-for-byte identical
FASTQ (identical md5). Read length is fixed at 125 bp for all three, and mean
quality is Q30.5 for all three. Two consequences:

- The configured `--read-length 150` is ignored; output is 125 bp.
- MiSeq (which should be 300 bp reads) and the platform distinction generally have
  no effect on short-read output.

Cause: ISS is invoked with `--mode basic`, which uses a generic error model and
ignores the platform-specific error model and read length. Long-read platforms
use a different simulator (PBSIM3) and are not affected by this, but were not
evaluated here.

Recommendation: if platform realism matters for benchmarking, drive ISS with the
kde error models (drop `--mode basic`) so NovaSeq, MiSeq, and HiSeq differ in read
length and error profile, and honor `--read-length`. If basic mode is a
deliberate speed or determinism choice, document that the three Illumina
platforms are interchangeable and that reads are 125 bp regardless of settings, so
users do not benchmark "platform effects" that do not exist.

Resolution (applied 2026-07-16): documented as interchangeable. The `--platform`
help text and the README now state that the three Illumina platforms produce
identical 125 bp reads in basic mode. Real per-platform differentiation (kde
error models, honoring `--read-length`) is deferred as a larger change.

### 3. Contamination fraction is approximate, not calibrated (medium priority)

Realized contamination scales correctly with the contamination level
(clean 0.0008, realistic 0.0049, heavy 0.0144) but does not match the metadata's
`contamination_fraction` precisely. For the mouse gut runs realized is about 75
percent of claimed (realistic: 0.0049 versus 0.0065; heavy: 0.0144 versus 0.0195).
For the larger gut collection the direction reverses and realized exceeds claimed
(0.0192 versus 0.0066). So `contamination_fraction` in the metadata is a
reasonable order-of-magnitude figure but should not be treated as the exact
realized value. For QC, the per-read labels are exact regardless, so contamination
removal metrics computed from labels are unaffected; only the reported summary
fraction is loose.

### 4. MDA chimera manifest genome names are inconsistent (low priority)

The chimera manifest sometimes records the donor or target genome with an ISS
read-index suffix (for example `GCF_004146805.1_2170`) rather than the bare
accession. The chimera reads themselves are correct (they contain two genomes'
sequence); only the manifest's genome identifier is occasionally the read name
rather than the genome id, which complicates programmatic use.

## What this means for QC benchmarking

Green light. The two ground-truth surfaces a QC benchmark reads, per-read
`source=` labels and artifact tags, are accurate, traceable, and reproducible.
Contamination-removal precision and recall computed from the labels will be
correct because the labels are correct, independent of the loose summary
fractions in finding 3.

Carry these caveats into later modules: dark-matter fraction (finding 1) and
platform realism (finding 2) affect taxonomy, assembly, and novel-discovery
benchmarking and any claim about composition realism, so they should be fixed or
clearly documented before those modules rely on them.

## Reproducing this evaluation

```bash
# per-dataset evaluation (JSON + human summary)
python scripts/evaluate_dataset.py --dataset <dataset_dir> --json out.json
```

The nine datasets were generated with `validation/eval/generate_eval_datasets.sh`
(collection 8 unless noted, seed 42, 20k reads). Pass/fail thresholds are the
constants at the top of `scripts/evaluate_dataset.py` and were fixed before
results were seen.
