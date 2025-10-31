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

## Next Steps

As more ViroForge modules are implemented, additional examples will be added here for:

- VLP enrichment modeling
- Sequencing artifact simulation (polyG tails, optical duplicates)
- Long-read sequencing (PacBio, Oxford Nanopore)
- Complete benchmarking workflows

## Questions or Issues?

Please open an issue on the GitHub repository: https://github.com/shandley/viroforge/issues
