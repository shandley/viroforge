# ViroForge

Synthetic virome dataset generator for benchmarking and validating virome analysis pipelines.

ViroForge creates realistic FASTQ sequencing data from curated viral genome collections with complete ground truth metadata. Use it to test QC tools, assemblers, taxonomic classifiers, and end-to-end virome workflows against datasets where every read's origin is known.

## Features

- 28 curated virome collections (human body sites, disease states, environmental)
- 14,423 RefSeq viral genomes with ICTV taxonomy
- 5 sequencing platforms (NovaSeq, MiSeq, HiSeq, PacBio HiFi, Oxford Nanopore)
- DNA and RNA virome workflows
- VLP enrichment modeling (5 protocols)
- Real reference contamination sequences (rRNA, host DNA, PhiX, adapters)
- Complete ground truth metadata for every dataset

## Installation

```bash
git clone https://github.com/shandley/viroforge.git
cd viroforge
pip install -e .
```

For long-read support (PacBio HiFi, Oxford Nanopore), you also need pbsim3, pbccs, and samtools. See [Long-Read Tutorial](docs/LONGREAD_TUTORIAL.md) for details.

For the web interface:

```bash
pip install -e ".[web]"
```

## Quick start

### Browse available collections

```bash
viroforge browse
```

### Generate a dataset using a preset

```bash
# List presets
viroforge presets list

# Generate
viroforge generate --preset gut-standard

# Override parameters
viroforge generate --preset gut-standard --seed 123 --output my_data
```

### Generate a dataset with full control

```bash
viroforge generate \
    --collection-id 9 \
    --output data/gut \
    --platform novaseq \
    --coverage 30 \
    --vlp-protocol tangential_flow
```

### Generate an RNA virome dataset

```bash
viroforge generate \
    --collection-id 21 \
    --output data/respiratory_rna \
    --molecule-type rna \
    --rna-depletion ribo_zero \
    --platform novaseq \
    --coverage 10
```

### Generate long-read data

```bash
# PacBio HiFi
viroforge generate \
    --collection-id 9 \
    --output data/gut_hifi \
    --platform pacbio-hifi \
    --depth 15

# Oxford Nanopore
viroforge generate \
    --collection-id 7 \
    --output data/soil_nanopore \
    --platform nanopore \
    --depth 20
```

### Generate matched short + long reads for hybrid assembly

```bash
python scripts/generate_hybrid_dataset.py \
    --collection-id 9 \
    --output data/gut_hybrid \
    --short-platform novaseq \
    --long-platform pacbio-hifi \
    --coverage 30 \
    --depth 15 \
    --seed 42
```

### Batch generation

```bash
# Run from a YAML config
viroforge batch examples/batch_configs/technology_comparison.yaml

# Run in parallel
viroforge batch examples/batch_configs/coverage_sweep.yaml --parallel 4
```

### Reporting and comparison

```bash
# View a dataset report
viroforge report data/gut-standard

# Compare multiple datasets
viroforge compare data/gut_novaseq data/gut_hiseq data/gut_pacbio_hifi
```

### Web interface

```bash
viroforge web
# Opens at http://127.0.0.1:5000
```

## Presets

| Preset | Description |
|--------|-------------|
| `gut-standard` | Human gut virome, NovaSeq, 30x, VLP enrichment |
| `gut-bulk` | Gut bulk metagenome, NovaSeq, 50x, no VLP |
| `marine-standard` | Marine virome, MiSeq, 30x |
| `respiratory-rna` | Respiratory RNA virome, NovaSeq, 40x |
| `quick-test-short` | Fast test dataset, 5x coverage |
| `quick-test-long` | Fast long-read test, 5x depth |
| `hybrid-standard` | Hybrid assembly, NovaSeq 30x + HiFi 15x |
| `assembly-high-coverage` | High coverage for assembly, 100x |

## Collections

### Host-associated

| ID | Collection | Genomes | Notes |
|----|-----------|---------|-------|
| 1 | Healthy Human Gut | 134 | Western diet, adult |
| 2 | Healthy Human Skin | 42 | Sebaceous sites |
| 3 | Healthy Human Oral | 67 | Saliva |
| 4 | Healthy Human Urogenital | 31 | Vaginal virome |
| 5 | Healthy Human Respiratory | 58 | Nasopharynx |
| 24 | Vaginal (Healthy) | 26 | Women's health |
| 25 | Blood/Plasma (Healthy) | 21 | Viremia, blood safety |
| 26 | Ocular Surface (Healthy) | 17 | Ophthalmology |
| 27 | Lower Respiratory/Lung (Healthy) | 31 | Pneumonia, transplant |
| 28 | Urinary (Healthy) | 20 | Kidney transplant |

### Disease states

| ID | Collection | Genomes | Notes |
|----|-----------|---------|-------|
| 18 | IBD Gut | 90 | Inflammatory bowel disease |
| 19 | HIV+ Gut | 55 | Includes herpesviruses |
| 20 | CF Respiratory | 81 | Cystic fibrosis lung |

### Environmental

| ID | Collection | Genomes | Notes |
|----|-----------|---------|-------|
| 6 | Marine | 78 | Coastal surface water |
| 7 | Soil | 82 | Agricultural |
| 8 | Freshwater | 71 | Lake surface water |
| 17 | Wastewater | 352 | Epidemiological surveillance |

### RNA viromes

| ID | Collection | Genomes | Notes |
|----|-----------|---------|-------|
| 21 | Human Respiratory RNA | 56 | Influenza, RSV, coronaviruses |
| 22 | Arbovirus Environmental | 39 | Flaviviruses, alphaviruses |
| 23 | Fecal RNA | 58 | Rotavirus, norovirus |

### VLP protocol comparisons (9-15)

Collections 9-15 provide matched compositions under different VLP enrichment conditions (baseline, high/low bacterial lysis, prophage induction, eukaryotic virus shedding).

### Amplification control (16)

Collection 16 is a pre-amplification control for method comparisons (100 genomes).

## VLP enrichment protocols

| Protocol | Method | Contamination reduction | Viral recovery |
|----------|--------|------------------------|----------------|
| `tangential_flow` | 0.2 um TFF | 91% | 85% |
| `ultracentrifugation` | Density gradient | 88% | 90% |
| `norgen` | Column-based | 87% | 70% |
| `syringe` | 0.22 um syringe | 86% | 60% |
| `none` | Bulk metagenome | 0% | 100% |

## RNA virome workflow

ViroForge models complete RNA virome library preparation:

- **Reverse transcription** with virus-type specific efficiency (40-90% depending on ssRNA+, ssRNA-, dsRNA)
- **rRNA depletion** via Ribo-Zero or RiboMinus (90% rRNA reduced to 10%, with 10-20x viral enrichment)
- **RNA degradation** and fragmentation modeling

Use `--molecule-type rna` with `--rna-depletion ribo_zero` to enable the full RNA workflow.

## Contamination modeling

ViroForge generates contamination reads from real reference sequences, making them detectable by standard QC tools (SortMeRNA, fastp, BBDuk, Kraken2).

**Bundled references (used by default):**

| Type | Source | Size |
|------|--------|------|
| rRNA | 23 sequences from NCBI RefSeq (E. coli 16S/23S, human 18S/28S, gut bacteria) | 41 KB |
| Host DNA | 48 fragments from T2T-CHM13v2.0 (all chromosomes, 10 kb each) | 481 KB |
| PhiX174 | NC_001422.1 | 5.4 KB |
| Adapters | TruSeq and Nextera (11 sequences) | <1 KB |

**Adapter read-through** can be added to simulate insert-size-dependent adapter contamination:

```bash
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/gut \
    --adapter-rate 0.05 \
    --adapter-type truseq
```

**Override references** with your own databases:

```bash
--host-genome /path/to/GRCh38.fasta
--rrna-database /path/to/SILVA_138.fasta
--no-real-contaminants   # revert to synthetic sequences
```

## Output structure

```
output/
  fasta/
    collection.fasta          # Reference genomes
  fastq/
    collection_R1.fastq       # Forward reads (Illumina)
    collection_R2.fastq       # Reverse reads (Illumina)
    collection_hifi.fastq.gz  # PacBio HiFi reads
    collection.fastq          # Nanopore reads
  metadata/
    metadata.json             # Ground truth (composition, taxonomy, workflow stats)
    composition.tsv           # Abundance table
    abundances.txt            # InSilicoSeq abundance file
```

## Requirements

- Python 3.9+
- numpy, pandas, biopython, scipy, pyyaml, rich
- InSilicoSeq (for Illumina read simulation)
- PBSIM3, pbccs, samtools (for long-read simulation, optional)
- Flask (for web interface, optional)

## Testing

```bash
pytest tests/ -v
```

## Citation

```bibtex
@software{viroforge2025,
  title = {ViroForge: A Synthetic Virome Data Generator},
  author = {Handley, Scott and contributors},
  year = {2025},
  url = {https://github.com/shandley/viroforge},
  version = {0.11.0}
}
```

## License

MIT. See [LICENSE](LICENSE).

## Links

- [Handley Lab](https://www.handleylab.org), Washington University in St. Louis
- [Hecatomb](https://github.com/shandley/hecatomb) - Viral metagenome assembly pipeline
- [Bug reports](https://github.com/shandley/viroforge/issues)
