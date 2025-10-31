# ViroForge Examples

This directory contains example scripts demonstrating how to use ViroForge modules.

## Available Examples

### create_community_example.py

Demonstrates the `core.community` module functionality:

- **Example 1**: Creating a custom viral community from scratch
- **Example 2**: Creating body-site specific communities (gut, oral, skin)
- **Example 3**: Comparing different abundance distributions (log-normal, power-law, even)
- **Example 4**: Exporting community data to files (TSV and FASTA)

**Run the example:**

```bash
python examples/create_community_example.py
```

**Output:**

The script will create an `example_output/` directory containing:
- `gut_virome_abundance.tsv` - Abundance table with taxonomy and metadata
- `gut_virome_genomes.fasta` - Genome sequences in FASTA format

### create_contamination_example.py

Demonstrates the `core.contamination` module functionality:

- **Example 1**: Pre-defined contamination profiles (clean, realistic, heavy, failed)
- **Example 2**: Building custom contamination profiles
- **Example 3**: Host DNA from different organisms (human, mouse, rat)
- **Example 4**: Creating complete mock virome compositions (viral + contamination)
- **Example 5**: Comparing VLP enrichment success vs failure
- **Example 6**: Exporting contamination data to files

**Run the example:**

```bash
python examples/create_contamination_example.py
```

**Output:**

The script will create an `example_output/` directory containing:
- `contamination_table.tsv` - Contamination profile with all contaminants
- `contaminants.fasta` - Contaminant sequences in FASTA format
- `complete_composition.tsv` - Complete virome composition (viral + contamination)

## Requirements

Make sure ViroForge is installed with its dependencies:

```bash
pip install -e .
```

Or with a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e .
```

### generate_reads_example.py

Demonstrates the `simulators.illumina` module for generating realistic Illumina reads:

- **Example 1**: Quick generation - one-liner for complete datasets
- **Example 2**: Custom composition with full control over parameters
- **Example 3**: Generating datasets with different contamination levels
- **Example 4**: Platform comparison (NovaSeq, MiSeq, HiSeq)
- **Example 5**: File size estimation before generation

**Requirements:**

This example requires InSilicoSeq to be installed:

```bash
conda install -c bioconda insilicoseq
# OR
pip install InSilicoSeq
```

**Run the example:**

```bash
python examples/generate_reads_example.py
```

**Output:**

The script will create FASTQ files and ground truth metadata:
- `*_R1.fastq` - Forward reads
- `*_R2.fastq` - Reverse reads
- `ground_truth_genomes.tsv` - Ground truth genome metadata
- `input_genomes.fasta` - Input genomes (if keep_temp_files=True)
- `input_abundances.txt` - Input abundances (if keep_temp_files=True)

**Note:** This example generates small datasets (50,000-100,000 reads) for demonstration. For production use, increase `n_reads` to 1M-10M+.

## VLP Enrichment Examples

### vlp_enrichment_basic.py

Demonstrates basic VLP (Virus-Like Particle) enrichment workflow:

- Creating a viral community and contamination profile
- Combining into a bulk metagenome composition (50% viral)
- Applying standard VLP enrichment protocol
- Calculating enrichment metrics and sequencing impact

**Run the example:**

```bash
python examples/vlp_enrichment_basic.py
```

**Key concepts:**
- VLP enrichment is the defining feature of viromics
- Increases viral fraction from ~50% to >95%
- Removes host and bacterial DNA contamination
- Dramatically improves viral genome coverage

### vlp_vs_bulk_comparison.py

Side-by-side comparison of VLP-enriched vs bulk metagenome sequencing:

- Creates paired samples from identical starting material
- Sample A: VLP enrichment applied
- Sample B: No enrichment (bulk metagenome)
- Compares viral fraction, sequencing yield, and cost efficiency
- Discusses research applications for each approach

**Run the example:**

```bash
python examples/vlp_vs_bulk_comparison.py
```

**Key insights:**
- VLP enrichment: ~2x increase in viral fraction
- Bulk metagenome: preserves virus-host ratios
- Cost efficiency: VLP saves ~50% sequencing costs for same viral reads
- Different use cases for virome vs metagenome studies

### vlp_protocol_comparison.py

Compares different VLP enrichment protocols used in the literature:

- **Standard VLP**: 0.2 μm TFF, 95% nuclease (most common)
- **Iron Chloride VLP**: FeCl3 precipitation, 98% nuclease (Conceição-Neto et al.)
- **Ultracentrifugation VLP**: Density gradient, 90% nuclease
- **Syringe Filter VLP**: Sharp cutoff, field-friendly
- **Bulk Metagenome**: No enrichment for comparison

**Run the example:**

```bash
python examples/vlp_protocol_comparison.py
```

**Guidance:**
- Protocol selection guide based on research goals
- Maximum purity vs standardization vs sample constraints
- Impact of methodology on final composition

### complete_vlp_workflow.py

Complete end-to-end workflow for creating VLP-enriched virome datasets:

- **Step 1**: Create viral community
- **Step 2**: Add realistic contamination
- **Step 3**: Create mock composition
- **Step 4**: Apply VLP enrichment
- **Step 5**: Generate sequencing reads
- **Step 6**: Export ground truth metadata

**Run the example:**

```bash
python examples/complete_vlp_workflow.py
```

**Output:**

The script creates an `output/vlp_virome_example/` directory containing:
- `ground_truth_composition.tsv` - Complete composition with abundances
- `summary_stats.txt` - Dataset summary statistics

**Use case:**
- Recommended workflow for benchmarking virome analysis pipelines
- Ideal for testing taxonomy assignment, assembly, and abundance estimation
- Includes ground truth for validation

### amplification_comparison.py

Demonstrates the `amplification` module for modeling library preparation biases:

- **No Amplification**: Control for high-biomass samples (no bias)
- **RdAB Amplification**: Random RT + dsDNA + PCR (most common method)
  - Length bias: exponential advantage for short genomes
  - GC bias: quadratic penalty at extreme GC values
  - Typical: 40 cycles, moderate bias
- **Linker Amplification**: Adapter ligation + PCR (modern protocols)
  - No length bias (all fragments have adapters)
  - Weak GC bias
  - Typical: 20 cycles, minimal bias
- **MDA Amplification**: Multiple Displacement Amplification (low biomass)
  - No length bias
  - Extreme GC bias (10-1000x range)
  - High stochasticity (random variations)
  - Chimera formation

**Run the example:**

```bash
python examples/amplification_comparison.py
```

**Key insights:**
- RdAB shows strongest bias (3x coefficient of variation)
- Linker has minimal bias (similar to no amplification)
- MDA has extreme GC bias but less predictable due to stochasticity
- Method choice depends on input DNA quantity and bias tolerance

### platform_comparison.py

Demonstrates the `artifacts` module for modeling platform-specific sequencing artifacts:

- **PolyG Tails**: Patterned flow cell artifact (NovaSeq, NextSeq)
  - Incomplete fluorophore quenching
  - R2 more affected than R1
  - Typical: 1-3% of reads, 15-45bp length
- **Optical Duplicates**: All platforms (rate varies)
  - Adjacent cluster signal bleed
  - NovaSeq: ~9%, MiSeq: ~2.5%
  - Can inflate coverage estimates
- **Index Hopping**: Barcode misassignment
  - Higher in patterned flow cells (1-2%)
  - Lower in cluster flow cells (~0.1%)
  - Critical for multiplexed studies

**Platforms compared:**
1. NovaSeq 6000 (patterned, high throughput)
2. NextSeq 2000 (patterned, mid throughput)
3. MiSeq (cluster, long reads, lowest artifacts)
4. HiSeq 2500 (cluster, legacy platform)
5. Ideal (no artifacts - control)

**Run the example:**

```bash
python examples/platform_comparison.py
```

**Key insights:**
- Patterned flow cells have polyG tails, cluster flow cells do not
- NovaSeq has highest optical duplicate rate (~9%)
- MiSeq has lowest artifacts, best for high-accuracy studies
- Index hopping higher in NovaSeq, requires mitigation strategies

**Use cases:**
- Testing cross-platform reproducibility
- Validating artifact removal pipelines
- Optimizing analysis for specific platforms
- Benchmarking with realistic platform biases

## Next Steps

As more ViroForge modules are implemented, additional examples will be added here for:

- Long-read sequencing (PacBio, Oxford Nanopore)
- Alternative platforms (MGI/DNBSEQ, Element Biosciences)
- Complete end-to-end workflows combining all modules

## Questions or Issues?

Please open an issue on the GitHub repository: https://github.com/shandley/viroforge/issues
