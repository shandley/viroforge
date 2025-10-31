# ViroForge Quickstart

**Get your first synthetic virome dataset in 5 minutes!**

---

## Install

```bash
git clone https://github.com/shandley/viroforge.git
cd viroforge
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -e .
```

---

## Generate Your First Dataset

Create a file `my_first_virome.py`:

```python
from viroforge.core.community import create_body_site_profile
from viroforge.core.contamination import create_contamination_profile
from viroforge.utils.composition import MockViromeComposition
from viroforge.enrichment import standard_vlp
from viroforge.amplification import rdab_40_cycles

# 1. Create a gut virome with 50 viral genomes
print("Creating viral community...")
community = create_body_site_profile('gut', n_genomes=50, random_seed=42)

# 2. Add realistic contamination (host DNA, bacteria)
print("Adding contamination...")
contamination = create_contamination_profile('realistic', random_seed=42)

# 3. Create initial composition (50% viral, 50% contamination)
print("Creating composition...")
composition = MockViromeComposition(
    name='my_first_virome',
    viral_community=community,
    contamination_profile=contamination,
    viral_fraction=0.5
)

print(f"  Initial viral fraction: {composition.viral_fraction*100:.1f}%")

# 4. Apply VLP enrichment (→ 97% viral)
print("Applying VLP enrichment...")
vlp = standard_vlp()
vlp.apply(composition)
print(f"  After VLP: {composition.viral_fraction*100:.1f}%")

# 5. Apply amplification bias (RdAB, 40 cycles)
print("Applying amplification...")
amplification = rdab_40_cycles()
amplification.apply(composition)
print(f"  After amplification: {composition.viral_fraction*100:.1f}%")

# 6. Export ground truth
import pandas as pd
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
print("\n✅ Success! Ground truth saved to: ground_truth_composition.tsv")
print(f"   Final viral fraction: {composition.viral_fraction*100:.1f}%")
print(f"   Total genomes: {len(composition_data)}")
```

Run it:
```bash
python my_first_virome.py
```

**Output:**
```
Creating viral community...
Adding contamination...
Creating composition...
  Initial viral fraction: 50.0%
Applying VLP enrichment...
  After VLP: 97.1%
Applying amplification...
  After amplification: 100.0%

✅ Success! Ground truth saved to: ground_truth_composition.tsv
   Final viral fraction: 100.0%
   Total genomes: 50
```

---

## What Just Happened?

1. **Viral Community** - Created 50 gut-specific viral genomes
2. **Contamination** - Added realistic host + bacterial DNA
3. **Initial Mix** - Started at 50% viral, 50% contamination
4. **VLP Enrichment** - Filtered + nuclease treatment → 97% viral
5. **Amplification** - RdAB PCR → 100% viral (removed remaining contam)
6. **Ground Truth** - Exported complete composition table

---

## Next Steps

### Run Complete Examples

```bash
# Full end-to-end workflow with FASTQ generation
python examples/complete_workflow_integrated.py

# Compare NovaSeq vs MiSeq platforms
python examples/cross_platform_workflow.py

# Compare VLP protocols
python examples/vlp_protocol_comparison.py
```

### Add Platform Artifacts

```python
from viroforge.artifacts import novaseq_6000, ReadPair

# Create mock reads
reads = []
for i in range(1000):
    read = ReadPair(
        read_id=f"read_{i}",
        forward_seq="A" * 80,
        reverse_seq="T" * 80,
        forward_qual="I" * 80,
        reverse_qual="I" * 80,
        genome_id=f"genome_{i % 50}",
        tile_x=i * 100,
        tile_y=i * 100
    )
    reads.append(read)

# Apply NovaSeq artifacts (polyG tails, optical dups, index hopping)
platform = novaseq_6000()
reads_with_artifacts = platform.apply(reads, random_seed=42)

print(f"Original reads: {len(reads)}")
print(f"With artifacts: {len(reads_with_artifacts)}")
print(f"Optical duplicates added: {len(reads_with_artifacts) - len(reads)}")
```

### Try Different Body Sites

```python
# Oral virome
oral = create_body_site_profile('oral', n_genomes=50, random_seed=42)

# Skin virome
skin = create_body_site_profile('skin', n_genomes=50, random_seed=42)

# Respiratory virome
respiratory = create_body_site_profile('respiratory', n_genomes=50, random_seed=42)
```

### Test Failed VLP Enrichment

```python
from viroforge.enrichment import VLPEnrichment

# Failed enrichment (should flag QC)
failed_vlp = VLPEnrichment(
    nuclease_efficiency=0.3,  # Poor nuclease
    stochastic_variation=0.5   # High variability
)
failed_vlp.apply(composition)
# Result: <80% viral (should be flagged by ViromeQC)
```

---

## Learn More

- **[Tutorial](TUTORIAL.md)** - Step-by-step learning (30 min)
- **[User Guide](USER_GUIDE.md)** - Complete documentation (2 hours)
- **[Examples](../examples/README.md)** - All example scripts explained
- **[API Reference](API.md)** - Full API documentation

---

## Common Use Cases

### 1. Benchmark Your Pipeline

```python
# Generate test dataset
community = create_body_site_profile('gut', n_genomes=100, random_seed=42)
# ... complete workflow ...
# Compare your pipeline results to ground_truth_composition.tsv
```

### 2. Test QC Tools

```python
# Successful VLP
success = create_contamination_profile('clean')  # 2% contam

# Failed VLP
failed = create_contamination_profile('failed')  # 35% contam

# Your QC tool should flag the failed sample!
```

### 3. Cross-Platform Testing

```python
from viroforge.artifacts import novaseq_6000, miseq

# Same community, different platforms
novaseq = novaseq_6000()  # Has polyG tails, high optical dups
miseq_platform = miseq()  # No polyG, low optical dups

# Test if your pipeline results are platform-independent
```

---

## Help & Support

- **Issues**: https://github.com/shandley/viroforge/issues
- **Discussions**: https://github.com/shandley/viroforge/discussions
- **Email**: scott.handley@wustl.edu

---

**Happy Forging! 🔨🦠**

*Generated datasets are ready for immediate use in benchmarking studies.*
