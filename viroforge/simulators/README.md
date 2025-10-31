# ViroForge Simulators Module

This module provides sequencing simulators for generating realistic reads from mock virome compositions.

## Illumina Sequencing (`illumina.py`)

Generate realistic Illumina paired-end reads using InSilicoSeq as the underlying simulator.

### Installation

Install InSilicoSeq before using this module:

```bash
# Via conda (recommended)
conda install -c bioconda insilicoseq

# Via pip
pip install InSilicoSeq
```

### Quick Start

```python
from viroforge.simulators import quick_generate

# Generate complete dataset in one line
output = quick_generate(
    body_site='gut',
    contamination_level='realistic',
    n_reads=1_000_000,
    random_seed=42
)

print(f"R1: {output['r1']}")
print(f"R2: {output['r2']}")
print(f"Ground truth: {output['ground_truth']}")
```

### Full Control

```python
from viroforge.utils import create_mock_virome
from viroforge.simulators import generate_reads

# Create custom composition
composition = create_mock_virome(
    name='my_virome',
    body_site='oral',
    contamination_level='clean',
    n_viral_genomes=50,
    viral_fraction=0.95
)

# Generate reads with specific parameters
output = generate_reads(
    composition=composition,
    output_prefix='output/my_dataset',
    n_reads=10_000_000,
    model='NovaSeq',  # Options: NovaSeq, HiSeq, MiSeq
    cpus=8,
    compress=True,  # Gzip compression
    gc_bias=True,   # Enable GC bias simulation
    validate_output=True,  # Validate FASTQ files
    random_seed=42
)
```

### Available Error Models

InSilicoSeq provides realistic error models for different Illumina platforms:

| Model | Read Length | Insert Size | Use Case |
|-------|-------------|-------------|----------|
| **NovaSeq** | 150bp PE | ~350bp | Current standard (2-channel chemistry) |
| **HiSeq** | 125bp PE | ~350bp | Previous generation |
| **MiSeq** | 250bp PE | ~450bp | Longer reads, smaller runs |
| **NextSeq** | 150bp PE | ~350bp | Alternative to HiSeq |

### File Size Estimation

```python
from viroforge.simulators import estimate_file_size

# Estimate before generating
print(estimate_file_size(1_000_000))
# Output: R1: ~95 MB, R2: ~95 MB, Total: ~190 MB

print(estimate_file_size(1_000_000, compress=True))
# Output: R1: ~24 MB, R2: ~24 MB, Total: ~48 MB (compressed)
```

### Output Files

The `generate_reads()` function creates:

1. **`*_R1.fastq`** - Forward reads (or .fastq.gz if compressed)
2. **`*_R2.fastq`** - Reverse reads (or .fastq.gz if compressed)
3. **`ground_truth_genomes.tsv`** - Genome-level ground truth metadata

Ground truth TSV contains:
- `genome_id` - Unique genome identifier
- `genome_type` - viral, host_dna, rrna, reagent_bacteria, phix
- `taxonomy` - Taxonomic lineage or organism name
- `length` - Genome length (bp)
- `gc_content` - GC content (%)
- `abundance` - Relative abundance (0-1)
- `source` - viral_community or contamination

### Validation

The module automatically validates FASTQ files to prevent common errors:

- ✅ Sequence length == quality length (prevents mismatches)
- ✅ Valid DNA characters (ATCGN only)
- ✅ Valid quality scores (Phred33: ASCII 33-126)
- ✅ File completeness (no truncation)

Validation can be disabled with `validate_output=False` (not recommended).

### Ground Truth Tracking

All generated reads can be traced back to their source genome using InSilicoSeq's read naming convention:

```
@genome_id_position_strand_readnumber
```

Example:
```
@genome123_1547_F_1
```

This allows you to:
- Validate read classification accuracy
- Calculate true positive/false positive rates
- Evaluate genome assembly completeness
- Benchmark analysis pipelines

### Examples

See `examples/generate_reads_example.py` for complete examples including:
- Quick generation
- Custom compositions
- Contamination level comparison
- Platform comparison (NovaSeq vs MiSeq vs HiSeq)
- File size estimation

### Performance

Approximate generation times (single CPU):

| Reads | NovaSeq | MiSeq | Notes |
|-------|---------|-------|-------|
| 100K | ~30 sec | ~45 sec | Quick testing |
| 1M | ~5 min | ~7 min | Standard dataset |
| 10M | ~50 min | ~70 min | Large dataset |
| 50M | ~4 hours | ~6 hours | Production scale |

Use `cpus` parameter to parallelize and speed up generation.

### Troubleshooting

**InSilicoSeq not found:**
```python
from viroforge.simulators import check_insilicoseq_installed

if not check_insilicoseq_installed():
    print("Install InSilicoSeq: conda install -c bioconda insilicoseq")
```

**Validation fails (compressed files):**

Validation doesn't support gzipped files yet. Either:
- Use `compress=False` during generation, or
- Decompress files before validation:
  ```bash
  gunzip *.fastq.gz
  ```

**Out of memory:**

For very large datasets (>50M reads), InSilicoSeq may use substantial memory. Consider:
- Reducing `n_reads`
- Generating multiple smaller batches
- Using a machine with more RAM

### Future Enhancements

Planned for future versions:
- Read-level ground truth export (exact read-to-genome mapping)
- Support for custom error models
- Integration with other simulators (ART, NEAT)
- Long-read support (PacBio, Oxford Nanopore)

---

For more information, see:
- [InSilicoSeq documentation](https://insilicoseq.readthedocs.io/)
- [ViroForge examples](../../examples/README.md)
- [ViroForge validation framework](../../docs/VALIDATION.md)
