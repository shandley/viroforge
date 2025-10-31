# ViroForge Tutorial

**A Step-by-Step Guide to Generating Synthetic Virome Data**

This tutorial will walk you through creating your first synthetic virome dataset with ViroForge, from basic viral communities to complete end-to-end pipelines.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation](#installation)
3. [Tutorial 1: Creating a Viral Community](#tutorial-1-creating-a-viral-community)
4. [Tutorial 2: Adding Contamination](#tutorial-2-adding-contamination)
5. [Tutorial 3: Applying VLP Enrichment](#tutorial-3-applying-vlp-enrichment)
6. [Tutorial 4: Modeling Amplification Bias](#tutorial-4-modeling-amplification-bias)
7. [Tutorial 5: Simulating Platform Artifacts](#tutorial-5-simulating-platform-artifacts)
8. [Tutorial 6: Complete End-to-End Workflow](#tutorial-6-complete-end-to-end-workflow)
9. [Tutorial 7: Cross-Platform Comparison](#tutorial-7-cross-platform-comparison)
10. [Next Steps](#next-steps)

---

## Prerequisites

- Python 3.8 or higher
- Basic familiarity with Python
- Understanding of virome sequencing concepts (helpful but not required)

---

## Installation

```bash
# Clone the repository
git clone https://github.com/shandley/viroforge.git
cd viroforge

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install ViroForge
pip install -e .
```

Verify installation:
```python
import viroforge
print(viroforge.__version__)
```

---

## Tutorial 1: Creating a Viral Community

Let's start by creating a simple viral community representing a gut virome.

### Step 1.1: Create a Basic Viral Community

```python
from viroforge.core.community import ViralCommunity, create_body_site_profile

# Create a gut-specific viral community with 20 genomes
community = create_body_site_profile(
    body_site='gut',
    n_genomes=20,
    random_seed=42  # For reproducibility
)

# View summary statistics
stats = community.get_summary_stats()
print(f"Total genomes: {stats['total_genomes']}")
print(f"Total families: {stats['families']}")
print(f"Mean genome length: {stats['mean_length']:,.0f} bp")
print(f"Mean GC content: {stats['mean_gc']*100:.1f}%")
```

**Expected Output:**
```
Total genomes: 20
Total families: 15
Mean genome length: 45,234 bp
Mean GC content: 48.3%
```

### Step 1.2: Explore the Community

```python
# View top 5 most abundant viruses
print("\nTop 5 Most Abundant Viruses:")
sorted_genomes = sorted(community.genomes, key=lambda g: g.abundance, reverse=True)
for i, genome in enumerate(sorted_genomes[:5], 1):
    print(f"{i}. {genome.organism} ({genome.family})")
    print(f"   Abundance: {genome.abundance:.4f}")
    print(f"   Length: {len(genome.sequence):,} bp")
    print()
```

### Step 1.3: Try Different Body Sites

```python
# Compare different body sites
body_sites = ['gut', 'oral', 'skin', 'respiratory']

for site in body_sites:
    community = create_body_site_profile(site, n_genomes=10, random_seed=42)
    stats = community.get_summary_stats()
    print(f"\n{site.upper()} Virome:")
    print(f"  Families: {stats['families']}")
    print(f"  Mean length: {stats['mean_length']:,.0f} bp")
```

### Key Concepts

- **Body Site Profiles**: Pre-configured viral compositions reflecting real body site diversity
- **Random Seed**: Ensures reproducibility - same seed produces same community
- **Abundance**: Relative abundance values sum to 1.0 across all genomes

---

## Tutorial 2: Adding Contamination

Real virome datasets always contain contamination. Let's model realistic contamination sources.

### Step 2.1: Create Contamination Profile

```python
from viroforge.core.contamination import create_contamination_profile

# Create realistic contamination profile
contamination = create_contamination_profile(
    contamination_level='realistic',
    random_seed=42
)

# View contamination sources
print("Contamination Sources:")
for source, abundance in contamination.get_source_summary().items():
    print(f"  {source}: {abundance*100:.2f}%")
```

**Expected Output:**
```
Contamination Sources:
  Host DNA (Human): 60.00%
  Bacterial DNA: 30.00%
  Fungal DNA: 8.00%
  PhiX Control: 2.00%
```

### Step 2.2: Compare Contamination Levels

```python
# Compare different contamination levels
levels = ['clean', 'realistic', 'heavy', 'failed']

for level in levels:
    contam = create_contamination_profile(level, random_seed=42)
    total_contam = sum(contam.get_source_summary().values())

    print(f"\n{level.upper()} Contamination:")
    print(f"  Total: {total_contam*100:.1f}%")
    for source, abund in contam.get_source_summary().items():
        print(f"  - {source}: {abund*100:.1f}%")
```

### Step 2.3: Create a Mock Virome Composition

```python
from viroforge.utils.composition import MockViromeComposition

# Create viral community
viral_community = create_body_site_profile('gut', n_genomes=20, random_seed=42)

# Create contamination
contamination = create_contamination_profile('realistic', random_seed=42)

# Combine into mock virome (50% viral, 50% contamination)
composition = MockViromeComposition(
    name='tutorial_virome',
    viral_community=viral_community,
    contamination_profile=contamination,
    viral_fraction=0.5  # 50% viral, 50% contamination
)

# Check composition
print(f"\nViral Fraction: {composition.viral_fraction*100:.1f}%")
print(f"Contamination Fraction: {(1-composition.viral_fraction)*100:.1f}%")
print(f"Total Abundance: {composition.get_total_abundance():.4f}")
```

### Key Concepts

- **Contamination Levels**: Pre-defined profiles (clean, realistic, heavy, failed)
- **Viral Fraction**: Proportion of reads from true viral genomes (0-1)
- **Mock Composition**: Combines viral community + contamination

---

## Tutorial 3: Applying VLP Enrichment

VLP (Virus-Like Particle) enrichment is the defining feature of viromics. Let's model it!

### Step 3.1: Apply Standard VLP Enrichment

```python
from viroforge.enrichment import standard_vlp

# Create composition (from Tutorial 2)
composition = MockViromeComposition(
    name='vlp_tutorial',
    viral_community=create_body_site_profile('gut', n_genomes=20, random_seed=42),
    contamination_profile=create_contamination_profile('realistic', random_seed=42),
    viral_fraction=0.5
)

print(f"Before VLP Enrichment:")
print(f"  Viral Fraction: {composition.viral_fraction*100:.1f}%")

# Apply standard VLP enrichment
vlp = standard_vlp()  # 0.2 μm TFF, 95% nuclease
vlp.apply(composition)

print(f"\nAfter VLP Enrichment:")
print(f"  Viral Fraction: {composition.viral_fraction*100:.1f}%")
print(f"  Enrichment Factor: {composition.viral_fraction / 0.5:.2f}x")
```

**Expected Output:**
```
Before VLP Enrichment:
  Viral Fraction: 50.0%

After VLP Enrichment:
  Viral Fraction: 97.1%
  Enrichment Factor: 1.94x
```

### Step 3.2: Compare VLP Protocols

```python
from viroforge.enrichment import (
    standard_vlp,
    iron_chloride_vlp,
    ultracentrifugation_vlp,
    syringe_filter_vlp
)

protocols = [
    ('Standard VLP', standard_vlp()),
    ('Iron Chloride', iron_chloride_vlp()),
    ('Ultracentrifugation', ultracentrifugation_vlp()),
    ('Syringe Filter', syringe_filter_vlp())
]

print("VLP Protocol Comparison:")
print("-" * 50)

for name, protocol in protocols:
    # Create fresh composition
    comp = MockViromeComposition(
        name=f'protocol_{name}',
        viral_community=create_body_site_profile('gut', n_genomes=20, random_seed=42),
        contamination_profile=create_contamination_profile('realistic', random_seed=42),
        viral_fraction=0.5
    )

    # Apply protocol
    protocol.apply(comp)

    print(f"\n{name}:")
    print(f"  Final Viral Fraction: {comp.viral_fraction*100:.1f}%")
    print(f"  Enrichment: {comp.viral_fraction / 0.5:.2f}x")
```

### Step 3.3: Model Failed VLP Enrichment

```python
from viroforge.enrichment import VLPEnrichment

# Successful enrichment
success = MockViromeComposition(
    name='success',
    viral_community=create_body_site_profile('gut', n_genomes=20, random_seed=42),
    contamination_profile=create_contamination_profile('realistic', random_seed=42),
    viral_fraction=0.5
)

# Failed enrichment
failed = MockViromeComposition(
    name='failed',
    viral_community=create_body_site_profile('gut', n_genomes=20, random_seed=42),
    contamination_profile=create_contamination_profile('realistic', random_seed=42),
    viral_fraction=0.5
)

# Apply protocols
standard_vlp().apply(success)

failed_protocol = VLPEnrichment(
    nuclease_efficiency=0.3,    # Poor nuclease treatment
    stochastic_variation=0.5    # High variability
)
failed_protocol.apply(failed)

print(f"Successful VLP: {success.viral_fraction*100:.1f}% viral")
print(f"Failed VLP: {failed.viral_fraction*100:.1f}% viral")
print(f"\nThis failed sample would be flagged by ViromeQC!")
```

### Key Concepts

- **VLP Enrichment**: Removes host/bacterial DNA, enriches viral particles
- **Filtration**: Size-based separation (0.2 μm typical for viruses)
- **Nuclease Treatment**: Enzymatic removal of free DNA/RNA
- **Protocol Variation**: Different methods give different purity levels

---

## Tutorial 4: Modeling Amplification Bias

Most virome protocols require PCR amplification, which introduces bias.

### Step 4.1: Apply RdAB Amplification

```python
from viroforge.amplification import rdab_40_cycles

# Create VLP-enriched composition
composition = MockViromeComposition(
    name='amplification_tutorial',
    viral_community=create_body_site_profile('gut', n_genomes=30, random_seed=42),
    contamination_profile=create_contamination_profile('clean', random_seed=42),
    viral_fraction=0.95  # Post-VLP enrichment level
)

# Get initial abundance distribution
initial_abundances = [g.abundance for g in composition.viral_community.genomes]
initial_cv = (
    sum((a - sum(initial_abundances)/len(initial_abundances))**2
        for a in initial_abundances) / len(initial_abundances)
) ** 0.5

print(f"Before Amplification:")
print(f"  Coefficient of Variation: {initial_cv:.3f}")

# Apply RdAB amplification (40 cycles)
amplification = rdab_40_cycles()
amplification.apply(composition)

# Get final abundance distribution
final_abundances = [g.abundance for g in composition.viral_community.genomes]
final_cv = (
    sum((a - sum(final_abundances)/len(final_abundances))**2
        for a in final_abundances) / len(final_abundances)
) ** 0.5

print(f"\nAfter Amplification:")
print(f"  Coefficient of Variation: {final_cv:.3f}")
print(f"  Bias Introduced: {(final_cv - initial_cv) / initial_cv * 100:.1f}% increase")
```

### Step 4.2: Compare Amplification Methods

```python
from viroforge.amplification import (
    no_amplification,
    rdab_30_cycles,
    rdab_40_cycles,
    mda_standard,
    linker_standard
)

methods = [
    ('No Amplification', no_amplification()),
    ('RdAB 30 cycles', rdab_30_cycles()),
    ('RdAB 40 cycles', rdab_40_cycles()),
    ('MDA', mda_standard()),
    ('Linker Amplification', linker_standard())
]

print("Amplification Method Comparison:")
print("-" * 60)

for name, method in methods:
    # Create fresh composition
    comp = MockViromeComposition(
        name=f'method_{name}',
        viral_community=create_body_site_profile('gut', n_genomes=30, random_seed=42),
        contamination_profile=create_contamination_profile('clean', random_seed=42),
        viral_fraction=0.95
    )

    # Apply method
    method.apply(comp)

    # Calculate coefficient of variation
    abundances = [g.abundance for g in comp.viral_community.genomes]
    mean_abund = sum(abundances) / len(abundances)
    variance = sum((a - mean_abund)**2 for a in abundances) / len(abundances)
    cv = (variance ** 0.5) / mean_abund if mean_abund > 0 else 0

    print(f"\n{name}:")
    print(f"  Coefficient of Variation: {cv:.3f}")
    print(f"  Bias Level: {'High' if cv > 2.0 else 'Moderate' if cv > 1.0 else 'Low'}")
```

### Key Concepts

- **Amplification Bias**: PCR preferentially amplifies certain genomes
- **RdAB**: Random RT + dsDNA + PCR (most common, moderate bias)
- **MDA**: Multiple Displacement Amplification (extreme GC bias, high stochasticity)
- **Linker**: Adapter ligation + PCR (minimal bias, modern method)
- **Coefficient of Variation**: Measures spread in abundance distribution

---

## Tutorial 5: Simulating Platform Artifacts

Different sequencing platforms introduce different artifacts.

### Step 5.1: Apply NovaSeq Artifacts

```python
from viroforge.artifacts import novaseq_6000, ReadPair

# Create mock reads (normally from FASTQ generation)
reads = []
for i in range(100):
    read = ReadPair(
        read_id=f"read_{i}",
        forward_seq="A" * 80,
        reverse_seq="T" * 80,
        forward_qual="I" * 80,
        reverse_qual="I" * 80,
        genome_id=f"genome_{i % 10}",
        tile_x=i * 100,
        tile_y=i * 100
    )
    reads.append(read)

print(f"Initial reads: {len(reads)}")

# Apply NovaSeq 6000 artifacts
platform = novaseq_6000()
reads_with_artifacts = platform.apply(reads, random_seed=42)

print(f"After artifacts: {len(reads_with_artifacts)}")
print(f"Optical duplicates added: {len(reads_with_artifacts) - len(reads)}")

# Count polyG tails
polyg_count = sum(1 for r in reads_with_artifacts
                  if len(r.forward_seq) > 80 or (r.reverse_seq and len(r.reverse_seq) > 80))
print(f"Reads with polyG tails: {polyg_count} ({polyg_count/len(reads)*100:.1f}%)")
```

### Step 5.2: Compare Sequencing Platforms

```python
from viroforge.artifacts import novaseq_6000, miseq, nextseq_2000, hiseq_2500

# Create test reads
reads = []
for i in range(1000):
    read = ReadPair(
        read_id=f"read_{i}",
        forward_seq="A" * 80,
        reverse_seq="T" * 80,
        forward_qual="I" * 80,
        reverse_qual="I" * 80,
        genome_id="genome_1",
        tile_x=i * 100,
        tile_y=i * 100
    )
    reads.append(read)

platforms = [
    ('NovaSeq 6000', novaseq_6000()),
    ('NextSeq 2000', nextseq_2000()),
    ('MiSeq', miseq()),
    ('HiSeq 2500', hiseq_2500())
]

print("Platform Artifact Comparison:")
print("-" * 70)
print(f"{'Platform':<15} {'Flow Cell':<12} {'Dup Rate':<10} {'PolyG':<10} {'Index Hop':<10}")
print("-" * 70)

for name, platform in platforms:
    # Apply artifacts
    artifact_reads = platform.apply(reads.copy(), random_seed=42)

    # Calculate metrics
    dup_rate = (len(artifact_reads) - 1000) / 1000 * 100
    polyg_count = sum(1 for r in artifact_reads
                      if len(r.forward_seq) > 80 or (r.reverse_seq and len(r.reverse_seq) > 80))
    polyg_pct = polyg_count / 1000 * 100 if len(artifact_reads) > 0 else 0

    flow_cell = platform.flow_cell_type if hasattr(platform, 'flow_cell_type') else 'N/A'

    print(f"{name:<15} {flow_cell:<12} {dup_rate:>6.1f}%    {polyg_pct:>5.1f}%     {'~1.5%' if 'NovaSeq' in name else '~0.1%'}")
```

### Key Concepts

- **Patterned Flow Cells**: NovaSeq, NextSeq (have polyG tails)
- **Cluster Flow Cells**: MiSeq, HiSeq (NO polyG tails)
- **Optical Duplicates**: Adjacent cluster signal bleed (varies by platform)
- **Index Hopping**: Barcode misassignment in multiplexed libraries

---

## Tutorial 6: Complete End-to-End Workflow

Now let's put it all together into a complete workflow!

### Step 6.1: Complete Pipeline

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000

# Step 1: Create viral community
print("Step 1: Creating viral community...")
viral_community = create_body_site_profile('gut', n_genomes=50, random_seed=42)
print(f"  ✓ Created community with {len(viral_community.genomes)} genomes")

# Step 2: Add contamination
print("\nStep 2: Adding contamination...")
contamination = create_contamination_profile('realistic', random_seed=42)
print(f"  ✓ Added realistic contamination sources")

# Step 3: Create initial composition
print("\nStep 3: Creating initial composition...")
composition = MockViromeComposition(
    name='tutorial_complete',
    viral_community=viral_community,
    contamination_profile=contamination,
    viral_fraction=0.5
)
print(f"  ✓ Initial viral fraction: {composition.viral_fraction*100:.1f}%")

# Step 4: Apply VLP enrichment
print("\nStep 4: Applying VLP enrichment...")
vlp = standard_vlp()
vlp.apply(composition)
print(f"  ✓ After VLP: {composition.viral_fraction*100:.1f}% viral (enrichment: {composition.viral_fraction/0.5:.2f}x)")

# Step 5: Apply amplification bias
print("\nStep 5: Applying amplification bias...")
amplification = rdab_40_cycles()
amplification.apply(composition)
print(f"  ✓ After amplification: {composition.viral_fraction*100:.1f}% viral")

# Step 6: Generate reads (simplified - see complete_workflow_integrated.py for full version)
print("\nStep 6: Generating reads...")
print(f"  ✓ Ready to generate FASTQ files")

# Step 7: Apply platform artifacts
print("\nStep 7: Applying platform artifacts...")
platform = novaseq_6000()
print(f"  ✓ NovaSeq 6000 artifacts configured")

print("\n✅ Complete workflow pipeline ready!")
print(f"   Final composition: {composition.viral_fraction*100:.1f}% viral")
```

### Step 6.2: Export Ground Truth

```python
import pandas as pd

# Export composition table
composition_data = []
for genome in composition.viral_community.genomes:
    composition_data.append({
        'genome_id': genome.genome_id,
        'organism': genome.organism,
        'family': genome.family,
        'abundance': genome.abundance,
        'length': len(genome.sequence),
        'gc_content': genome.gc
    })

df = pd.DataFrame(composition_data)
df.to_csv('ground_truth_composition.tsv', sep='\t', index=False)
print("\n✅ Ground truth exported to ground_truth_composition.tsv")
```

---

## Tutorial 7: Cross-Platform Comparison

Let's compare how the same viral community performs on different platforms.

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles
from viroforge.artifacts import novaseq_6000, miseq

# Create identical starting composition
def create_test_composition(name):
    return MockViromeComposition(
        name=name,
        viral_community=create_body_site_profile('gut', n_genomes=30, random_seed=42),
        contamination_profile=create_contamination_profile('realistic', random_seed=42),
        viral_fraction=0.5
    )

# Platform A: NovaSeq 6000
print("Platform A: NovaSeq 6000")
print("-" * 40)
comp_novaseq = create_test_composition('novaseq')
vlp = standard_vlp()
vlp.apply(comp_novaseq)
amp = rdab_40_cycles()
amp.apply(comp_novaseq)
print(f"  Viral fraction after processing: {comp_novaseq.viral_fraction*100:.1f}%")
print(f"  Expected artifacts:")
print(f"    - PolyG tails: ~2.5%")
print(f"    - Optical duplicates: ~9%")
print(f"    - Index hopping: ~1.5%")

# Platform B: MiSeq
print("\nPlatform B: MiSeq")
print("-" * 40)
comp_miseq = create_test_composition('miseq')
vlp = standard_vlp()
vlp.apply(comp_miseq)
amp = rdab_40_cycles()
amp.apply(comp_miseq)
print(f"  Viral fraction after processing: {comp_miseq.viral_fraction*100:.1f}%")
print(f"  Expected artifacts:")
print(f"    - PolyG tails: 0% (cluster flow cell)")
print(f"    - Optical duplicates: ~2.5%")
print(f"    - Index hopping: ~0.1%")

print("\nKey Differences:")
print("  • NovaSeq has polyG tails, MiSeq does not")
print("  • NovaSeq has higher optical duplicate rate")
print("  • NovaSeq has higher index hopping rate")
print("  • Both platforms give same biological composition")
```

---

## Next Steps

Congratulations! You've learned the basics of ViroForge. Here are some next steps:

### 1. Explore Complete Examples

Check out the `examples/` directory for full working examples:
```bash
python examples/complete_workflow_integrated.py
python examples/cross_platform_workflow.py
python examples/vlp_protocol_comparison.py
```

### 2. Generate Full FASTQ Files

See `examples/complete_workflow_integrated.py` for how to:
- Generate paired-end FASTQ files
- Export complete ground truth
- Create pipeline summary reports

### 3. Read the User Guide

See `docs/USER_GUIDE.md` for:
- Detailed parameter explanations
- Advanced configuration options
- Best practices
- Troubleshooting

### 4. Explore the API

See `docs/API.md` for complete API documentation of all modules.

### 5. Run the Test Suite

```bash
pytest tests/ -v
```

Explore the tests to see more usage examples!

### 6. Join the Community

- Open an issue on GitHub
- Contribute new features
- Share your datasets
- Provide feedback

---

## Common Questions

**Q: How do I generate actual FASTQ files?**
A: See `examples/complete_workflow_integrated.py` for the complete workflow including FASTQ generation.

**Q: Can I use my own viral genomes?**
A: Yes! Create a custom `ViralCommunity` and add your own `ViralGenome` objects.

**Q: How do I validate my analysis pipeline?**
A: Generate a dataset, run your pipeline, compare results to the ground truth files.

**Q: What's the difference between VLP and bulk metagenome?**
A: VLP enrichment removes host/bacterial DNA. See `examples/vlp_vs_bulk_comparison.py` for a direct comparison.

**Q: Which platform should I use?**
A: It depends on your study. See `examples/platform_comparison.py` for detailed comparisons.

---

## Additional Resources

- **Examples**: `examples/` directory
- **User Guide**: `docs/USER_GUIDE.md`
- **API Reference**: `docs/API.md`
- **Design Rationale**: `docs/DESIGN_RATIONALE.md`
- **VLP Biology**: `docs/VLP_ENRICHMENT_BIOLOGY.md`
- **GitHub Issues**: https://github.com/shandley/viroforge/issues

---

**Happy Forging! 🔨🦠**
