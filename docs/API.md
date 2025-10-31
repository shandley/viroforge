# ViroForge API Reference

**Python API Documentation**

---

## Status

🚧 **Under Development**

The API reference is being prepared and will be auto-generated using Sphinx.

For now, please refer to:
- **[User Guide](USER_GUIDE.md)** - Complete parameter documentation
- **[Tutorial](TUTORIAL.md)** - Working code examples
- **[Examples](../examples/README.md)** - Real-world usage patterns

---

## Quick API Overview

### Core Modules

#### Community Module

```python
from viroforge.core.community import create_body_site_profile, ViralCommunity

# Create pre-built community
community = create_body_site_profile(
    body_site='gut',      # 'gut', 'oral', 'skin', 'respiratory', 'diverse'
    n_genomes=50,         # Number of viral genomes
    random_seed=42        # For reproducibility
)

# Or create custom community
community = ViralCommunity(name='my_virome')
```

#### Contamination Module

```python
from viroforge.core.contamination import create_contamination_profile

# Create pre-built contamination
contamination = create_contamination_profile(
    contamination_level='realistic',  # 'clean', 'realistic', 'heavy', 'failed'
    random_seed=42
)
```

#### Composition Module

```python
from viroforge.utils.composition import MockViromeComposition

# Combine viral + contamination
composition = MockViromeComposition(
    name='my_mock_virome',
    viral_community=community,
    contamination_profile=contamination,
    viral_fraction=0.5  # 50% viral, 50% contamination
)
```

#### Enrichment Module

```python
from viroforge.enrichment import standard_vlp, VLPEnrichment

# Use pre-built protocol
vlp = standard_vlp()
vlp.apply(composition)

# Or create custom protocol
vlp = VLPEnrichment(
    filtration_cutoff_um=0.2,
    nuclease_efficiency=0.95,
    random_seed=42
)
vlp.apply(composition)
```

#### Amplification Module

```python
from viroforge.amplification import rdab_40_cycles, RdABAmplification

# Use pre-built protocol
amp = rdab_40_cycles()
amp.apply(composition)

# Or create custom protocol
amp = RdABAmplification(
    cycles=40,
    efficiency=0.90,
    random_seed=42
)
amp.apply(composition)
```

#### Artifacts Module

```python
from viroforge.artifacts import novaseq_6000, ReadPair

# Use pre-built platform
platform = novaseq_6000()
reads_with_artifacts = platform.apply(reads, random_seed=42)

# Create ReadPair
read = ReadPair(
    read_id="read_001",
    forward_seq="ATCG...",
    reverse_seq="GCTA...",
    forward_qual="IIII...",
    reverse_qual="IIII...",
    genome_id="NC_001416"
)
```

---

## Complete Workflow Example

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000

# 1. Create viral community
community = create_body_site_profile('gut', n_genomes=50, random_seed=42)

# 2. Add contamination
contamination = create_contamination_profile('realistic', random_seed=42)

# 3. Create composition
composition = MockViromeComposition(
    name='complete_workflow',
    viral_community=community,
    contamination_profile=contamination,
    viral_fraction=0.5
)

# 4. Apply VLP enrichment
vlp = standard_vlp()
vlp.apply(composition)

# 5. Apply amplification
amplification = rdab_40_cycles()
amplification.apply(composition)

# 6. Apply platform artifacts
platform = novaseq_6000()
# reads_with_artifacts = platform.apply(reads, random_seed=42)

print(f"Final viral fraction: {composition.viral_fraction*100:.1f}%")
```

---

## Module Reference

### viroforge.core.community

**Functions:**
- `create_body_site_profile(body_site, n_genomes, random_seed)` - Create pre-built community
- `ViralCommunity(name)` - Custom community container
- `ViralGenome(...)` - Individual viral genome

**Body Sites:**
- `'gut'` - Gut-associated viruses
- `'oral'` - Oral cavity viruses
- `'skin'` - Skin-associated viruses
- `'respiratory'` - Respiratory tract viruses
- `'diverse'` - Maximum diversity

### viroforge.core.contamination

**Functions:**
- `create_contamination_profile(level, random_seed)` - Create pre-built contamination
- `ContaminationProfile(name)` - Custom contamination container

**Contamination Levels:**
- `'clean'` - ~2% total contamination
- `'realistic'` - ~5% total contamination
- `'heavy'` - ~10% total contamination
- `'failed'` - ~35% total contamination

### viroforge.enrichment

**Pre-built Protocols:**
- `standard_vlp()` - 0.2 μm TFF, 95% nuclease
- `iron_chloride_vlp()` - FeCl3 precipitation, 98% nuclease
- `ultracentrifugation_vlp()` - Density gradient, 90% nuclease
- `syringe_filter_vlp()` - 0.45 μm syringe, 80% nuclease

**Custom Protocol:**
- `VLPEnrichment(...)` - Full parameter control

### viroforge.amplification

**Pre-built Methods:**
- `no_amplification()` - Control (no bias)
- `rdab_30_cycles()` - RdAB, 30 cycles
- `rdab_40_cycles()` - RdAB, 40 cycles (standard)
- `mda_standard()` - MDA, 8 hours
- `mda_overnight()` - MDA, 16 hours
- `linker_standard()` - Linker, 20 cycles

**Custom Methods:**
- `RdABAmplification(...)` - Custom RdAB protocol
- `MDAAmplification(...)` - Custom MDA protocol
- `LinkerAmplification(...)` - Custom linker protocol

### viroforge.artifacts

**Pre-built Platforms:**
- `novaseq_6000()` - NovaSeq 6000 (patterned flow cell)
- `nextseq_2000()` - NextSeq 2000 (patterned flow cell)
- `miseq()` - MiSeq (cluster flow cell)
- `hiseq_2500()` - HiSeq 2500 (cluster flow cell)
- `no_artifacts()` - Ideal platform (control)

**Classes:**
- `ReadPair` - Paired-end read with metadata
- `PlatformProfile` - Custom platform configuration
- `PolyGTailArtifact` - PolyG tail artifact
- `OpticalDuplicateArtifact` - Optical duplicate artifact
- `IndexHoppingArtifact` - Index hopping artifact

---

## Return Types

### ViralCommunity

**Attributes:**
- `name` - Community name
- `genomes` - List of ViralGenome objects
- `metadata` - Additional metadata

**Methods:**
- `add_genome(genome)` - Add a viral genome
- `get_summary_stats()` - Get summary statistics
- `export_fasta(filename)` - Export to FASTA
- `export_abundance_table(filename)` - Export abundance table

### MockViromeComposition

**Attributes:**
- `name` - Composition name
- `viral_community` - ViralCommunity object
- `contamination_profile` - ContaminationProfile object
- `viral_fraction` - Viral fraction (0-1)

**Methods:**
- `get_total_abundance()` - Get total abundance (should be 1.0)
- `renormalize_abundances()` - Renormalize to 1.0
- `export_composition(filename)` - Export composition table

### ReadPair

**Attributes:**
- `read_id` - Read identifier
- `forward_seq` - Forward read sequence
- `reverse_seq` - Reverse read sequence
- `forward_qual` - Forward quality scores
- `reverse_qual` - Reverse quality scores
- `genome_id` - Source genome (ground truth)
- `tile_x`, `tile_y` - Flow cell coordinates
- `sample_index` - Sample barcode

---

## For More Information

- **[User Guide](USER_GUIDE.md)** - Complete parameter documentation
- **[Tutorial](TUTORIAL.md)** - Step-by-step learning
- **[Examples](../examples/)** - Working examples
- **Source Code** - See module docstrings

---

## Coming Soon

- Sphinx-generated API documentation
- Parameter type annotations
- Return value specifications
- Exception documentation
- Usage examples for each method

---

**Last Updated**: October 31, 2025
