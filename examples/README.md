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

## Next Steps

As more ViroForge modules are implemented, additional examples will be added here for:

- Contamination profile generation
- VLP enrichment modeling
- Sequencing artifact simulation
- Complete end-to-end virome generation

## Questions or Issues?

Please open an issue on the GitHub repository: https://github.com/shandley/viroforge/issues
