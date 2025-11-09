# VLP Protocol Comparison Tutorial

**ViroForge Phase 5 Feature**

This tutorial demonstrates how to compare different VLP enrichment protocols and understand their effects on viral recovery and contamination reduction.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Understanding VLP Protocols](#understanding-vlp-protocols)
4. [Step-by-Step Tutorial](#step-by-step-tutorial)
5. [Interpreting Results](#interpreting-results)
6. [Common Use Cases](#common-use-cases)
7. [Troubleshooting](#troubleshooting)

---

## Introduction

### Why Compare VLP Protocols?

Different VLP (Virus-Like Particle) enrichment protocols have distinct characteristics:

- **Viral recovery rates**: How much viral material is retained?
- **Contamination reduction**: How effectively are non-viral sequences removed?
- **Size bias**: Do protocols favor certain virus sizes?
- **Cost and complexity**: What are the practical trade-offs?

ViroForge enables systematic comparison with known ground truth, helping you:

1. **Optimize your experimental design**: Choose the best protocol for your study
2. **Validate QC pipelines**: Ensure your analysis detects protocol effects
3. **Understand limitations**: Know what biases to expect in real data

---

## Quick Start

### Generate All 5 VLP Protocols in One Command

```bash
# Generate datasets for all VLP protocols
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output data/vlp_comparison
```

This generates 5 datasets with identical starting material but different VLP protocols:
1. Tangential Flow Filtration (TFF)
2. Syringe Filter
3. Ultracentrifugation
4. Norgen Kit (column-based)
5. Bulk Metagenome (no VLP)

**Runtime**: ~20-30 minutes (depends on collection size)

**Output**:
```
data/vlp_comparison/
├── tangential_flow/
│   ├── fastq/
│   └── metadata/
├── syringe/
│   ├── fastq/
│   └── metadata/
├── ultracentrifugation/
│   ├── fastq/
│   └── metadata/
├── norgen/
│   ├── fastq/
│   └── metadata/
└── bulk/
    ├── fastq/
    └── metadata/
```

---

## Understanding VLP Protocols

### Protocol Characteristics

| Protocol | Method | Viral Recovery | Contamination Removal | Size Bias | Best For |
|----------|--------|---------------|----------------------|-----------|----------|
| **Tangential Flow** | 0.2 μm TFF | 85% | 91% | Moderate | Highest purity |
| **Syringe Filter** | 0.22 μm syringe | 60% | 86% | Strong | Field-friendly |
| **Ultracentrifugation** | Density gradient | 90% | 88% | Minimal | Highest recovery |
| **Norgen Kit** | Column-based | 70% | 87% | Moderate | Convenience |
| **Bulk** | None | 100% | 0% | None | Control |

### Key Differences

**Tangential Flow Filtration (TFF)**
- Most common in modern viromics
- Continuous flow allows processing large volumes
- Highest contamination reduction (91%)
- Moderate size bias (smaller viruses partially lost)

**Syringe Filter**
- Simple, field-deployable
- Lower recovery (60%) due to filter clogging
- Strong size bias (sharper retention curve)
- Cost-effective for small samples

**Ultracentrifugation**
- Traditional "gold standard" method
- Highest viral recovery (90%)
- Density-based separation (minimal size bias)
- Time-consuming, requires specialized equipment

**Norgen Kit**
- Commercial column-based extraction
- Moderate efficiency on all metrics
- Most convenient (no filtration equipment)
- Good for standardization across labs

**Bulk Metagenome**
- No VLP enrichment (control)
- All contamination retained
- Essential comparison for understanding VLP effects

---

## Step-by-Step Tutorial

### Step 1: Generate Individual Protocol Datasets

```bash
# Create output directory
mkdir -p tutorial_vlp_comparison

# Tangential flow (highest purity)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output tutorial_vlp_comparison/tangential_flow \
    --coverage 10 \
    --vlp-protocol tangential_flow \
    --contamination-level realistic \
    --seed 42

# Ultracentrifugation (highest recovery)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output tutorial_vlp_comparison/ultracentrifugation \
    --coverage 10 \
    --vlp-protocol ultracentrifugation \
    --contamination-level realistic \
    --seed 42

# Bulk metagenome (control)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output tutorial_vlp_comparison/bulk \
    --coverage 10 \
    --no-vlp \
    --contamination-level realistic \
    --seed 42
```

**Important**: Use the same `--seed` for reproducible comparisons!

### Step 2: Analyze Ground Truth Metadata

```bash
# View enrichment statistics for each protocol
cat tutorial_vlp_comparison/tangential_flow/metadata/*_metadata.json | jq '.enrichment_statistics'
cat tutorial_vlp_comparison/ultracentrifugation/metadata/*_metadata.json | jq '.enrichment_statistics'
cat tutorial_vlp_comparison/bulk/metadata/*_metadata.json | jq '.enrichment_statistics'
```

### Step 3: Extract Key Metrics

Create a Python script `compare_protocols.py`:

```python
#!/usr/bin/env python3
import json
import glob
from pathlib import Path
import pandas as pd

def extract_metrics(base_dir):
    """Extract comparison metrics from all protocol metadata files."""
    results = []

    protocols = ['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen', 'bulk']

    for protocol in protocols:
        metadata_file = glob.glob(f"{base_dir}/{protocol}/metadata/*_metadata.json")

        if not metadata_file:
            continue

        with open(metadata_file[0]) as f:
            data = json.load(f)

        # Extract metrics
        enrichment = data.get('enrichment_statistics', {})

        results.append({
            'Protocol': protocol.replace('_', ' ').title(),
            'Viral Fraction (%)': enrichment.get('viral_fraction_after', 0) * 100,
            'Contamination (%)': (1 - enrichment.get('viral_fraction_after', 1)) * 100,
            'Host DNA Removal (%)': enrichment.get('contamination_reduction', {}).get('host_dna_removal', 0) * 100,
            'rRNA Removal (%)': enrichment.get('contamination_reduction', {}).get('rrna_removal', 0) * 100,
            'Bacteria Removal (%)': enrichment.get('contamination_reduction', {}).get('bacteria_removal', 0) * 100,
            'Overall Reduction (%)': enrichment.get('contamination_reduction', {}).get('overall_reduction', 0) * 100
        })

    df = pd.DataFrame(results)
    return df

# Run comparison
df = extract_metrics('tutorial_vlp_comparison')
print(df.to_string(index=False))
df.to_csv('protocol_comparison.csv', index=False)
```

Run it:
```bash
python compare_protocols.py
```

**Expected Output**:
```
              Protocol  Viral Fraction (%)  Contamination (%)  Host DNA Removal (%)  rRNA Removal (%)  Bacteria Removal (%)  Overall Reduction (%)
  Tangential Flow                   99.35               0.65                 96.2             89.4                 98.5                   91.2
          Syringe                    98.95               1.05                 89.6             84.1                 94.5                   85.7
Ultracentrifugation                  99.15               0.85                 94.0             87.3                 93.6                   88.4
             Norgen                  99.06               0.94                 91.4             85.3                 95.5                   87.1
               Bulk                  92.60               7.40                  0.0              0.0                  0.0                    0.0
```

### Step 4: Visualize Protocol Comparison

Create `visualize_comparison.py`:

```python
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read comparison data
df = pd.read_csv('protocol_comparison.csv')

# Set style
sns.set_style("whitegrid")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Viral Fraction Comparison
ax1 = axes[0, 0]
df.plot(x='Protocol', y='Viral Fraction (%)', kind='bar', ax=ax1, legend=False, color='steelblue')
ax1.set_ylabel('Viral Fraction (%)')
ax1.set_title('Viral Fraction After VLP Enrichment')
ax1.axhline(y=95, color='red', linestyle='--', label='ViromeQC threshold (95%)')
ax1.legend()
ax1.set_ylim(90, 100)

# 2. Overall Contamination Reduction
ax2 = axes[0, 1]
df_vlp = df[df['Protocol'] != 'Bulk']
df_vlp.plot(x='Protocol', y='Overall Reduction (%)', kind='bar', ax=ax2, legend=False, color='forestgreen')
ax2.set_ylabel('Contamination Removed (%)')
ax2.set_title('Overall Contamination Reduction')
ax2.set_ylim(80, 95)

# 3. Type-Specific Contamination Removal
ax3 = axes[1, 0]
df_vlp[['Protocol', 'Host DNA Removal (%)', 'rRNA Removal (%)', 'Bacteria Removal (%)']].set_index('Protocol').plot(
    kind='bar', ax=ax3
)
ax3.set_ylabel('Removal (%)')
ax3.set_title('Type-Specific Contamination Removal')
ax3.legend(title='Contaminant Type', bbox_to_anchor=(1.05, 1), loc='upper left')
ax3.set_ylim(75, 100)

# 4. VLP vs Bulk Comparison
ax4 = axes[1, 1]
comparison_data = df[['Protocol', 'Contamination (%)']].sort_values('Contamination (%)')
comparison_data.plot(x='Protocol', y='Contamination (%)', kind='bar', ax=ax4, legend=False, color='coral')
ax4.set_ylabel('Contamination Remaining (%)')
ax4.set_title('Final Contamination Levels')
ax4.set_yscale('log')

plt.tight_layout()
plt.savefig('vlp_protocol_comparison.png', dpi=300, bbox_inches='tight')
print("Saved visualization to vlp_protocol_comparison.png")
```

Run it:
```bash
python visualize_comparison.py
```

---

## Interpreting Results

### Expected Patterns

**1. Viral Fraction**
- All VLP protocols should achieve >98% viral fraction
- Bulk metagenome typically 85-95% (depends on starting contamination)
- Fold difference: 5-20x between VLP and bulk

**2. Contamination Reduction**
- **Host DNA**: Highest removal (90-98%) due to DNase efficiency
- **rRNA**: Moderate removal (85-95%), slightly lower than DNA
- **Bacteria**: Variable (85-99%) depending on filtration method
- **PhiX**: Low removal (10-40%), treated as small virus

**3. Protocol Rankings**

**By Viral Purity** (highest to lowest):
1. Tangential Flow (99.35%)
2. Ultracentrifugation (99.15%)
3. Norgen (99.06%)
4. Syringe (98.95%)
5. Bulk (92.60%)

**By Viral Recovery** (if you had pre-VLP abundance data):
1. Ultracentrifugation (90%)
2. Tangential Flow (85%)
3. Norgen (70%)
4. Syringe (60%)

**By Contamination Reduction**:
1. Tangential Flow (91.2%)
2. Ultracentrifugation (88.4%)
3. Norgen (87.1%)
4. Syringe (85.7%)

### Red Flags

**Warning Signs**:
- Viral fraction <95% with VLP → Poor protocol execution
- Bulk metagenome >98% viral → Unrealistic starting material
- No difference between protocols → Check that VLP was actually applied

---

## Common Use Cases

### Use Case 1: Protocol Selection for New Study

**Goal**: Decide which VLP protocol to use for a gut virome study

**Approach**:
```bash
# Generate all protocols with gut collection
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output gut_protocol_test

# Run your pipeline on each dataset
for protocol in tangential_flow syringe ultracentrifugation norgen bulk; do
    your_pipeline gut_protocol_test/${protocol}/fastq/*_R{1,2}.fastq \
        --output results/${protocol}
done

# Compare results
python compare_protocols.py
```

**Decision Criteria**:
- **Maximum purity needed?** → Tangential Flow
- **Maximum recovery needed?** → Ultracentrifugation
- **Field work/portability?** → Syringe
- **Standardization across labs?** → Norgen
- **Budget constrained?** → Consider syringe or Norgen

### Use Case 2: Validate ViromeQC Implementation

**Goal**: Ensure your QC pipeline correctly flags poor VLP enrichment

**Approach**:
```bash
# Generate high-quality VLP (should PASS QC)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output qc_test/pass_vlp \
    --vlp-protocol tangential_flow \
    --contamination-level clean

# Generate poor VLP (should FAIL QC)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output qc_test/fail_vlp \
    --vlp-protocol tangential_flow \
    --contamination-level heavy

# Generate bulk (should FAIL QC)
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output qc_test/fail_bulk \
    --no-vlp \
    --contamination-level realistic

# Run ViromeQC
viromeqc --input qc_test/pass_vlp/fastq/*_R1.fastq --output qc_pass.txt
viromeqc --input qc_test/fail_vlp/fastq/*_R1.fastq --output qc_fail_vlp.txt
viromeqc --input qc_test/fail_bulk/fastq/*_R1.fastq --output qc_fail_bulk.txt
```

**Expected Results**:
- `qc_pass.txt`: Enrichment score >0.95 (PASS)
- `qc_fail_vlp.txt`: Enrichment score 0.80-0.90 (WARNING)
- `qc_fail_bulk.txt`: Enrichment score <0.80 (FAIL)

### Use Case 3: Size Bias Analysis

**Goal**: Understand how different protocols affect small vs large viruses

**Approach**:
```bash
# Generate with different protocols
python scripts/batch_generate_fastq.py \
    --preset vlp-protocol-comparison \
    --output size_bias_test

# Analyze size distribution in each dataset
python analyze_size_bias.py
```

Create `analyze_size_bias.py`:
```python
import json
import glob
import pandas as pd
import matplotlib.pyplot as plt

def analyze_size_bias(base_dir):
    protocols = ['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen', 'bulk']

    for protocol in protocols:
        metadata_file = glob.glob(f"{base_dir}/{protocol}/metadata/*_metadata.json")[0]

        with open(metadata_file) as f:
            data = json.load(f)

        # Extract viral genomes and their properties
        viral_seqs = [s for s in data['sequences'] if s['type'] == 'viral']

        df = pd.DataFrame(viral_seqs)
        df['size_category'] = pd.cut(df['length'],
                                     bins=[0, 10000, 50000, 100000, 1000000],
                                     labels=['Small (<10kb)', 'Medium (10-50kb)',
                                            'Large (50-100kb)', 'Jumbo (>100kb)'])

        # Plot size distribution
        size_dist = df.groupby('size_category')['relative_abundance'].sum()

        plt.figure(figsize=(8, 5))
        size_dist.plot(kind='bar')
        plt.title(f'Viral Size Distribution - {protocol.replace("_", " ").title()}')
        plt.ylabel('Cumulative Relative Abundance')
        plt.xlabel('Genome Size Category')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'size_dist_{protocol}.png', dpi=300)
        plt.close()

analyze_size_bias('size_bias_test')
```

**Expected Patterns**:
- **Tangential Flow**: Slight depletion of <5kb viruses
- **Syringe**: More pronounced depletion of <5kb viruses
- **Ultracentrifugation**: Minimal size bias
- **Norgen**: Moderate depletion of very small viruses
- **Bulk**: No size bias

---

## Troubleshooting

### Issue: All Protocols Show Same Results

**Possible Causes**:
1. VLP enrichment not actually applied
2. Starting contamination too low to see differences
3. Using same random seed without resetting

**Solutions**:
```bash
# Check that VLP was applied
cat output/metadata/*_metadata.json | jq '.configuration.vlp_protocol'

# Use higher contamination level
--contamination-level heavy

# Use different seeds for each run
--seed 42  # for protocol 1
--seed 43  # for protocol 2
```

### Issue: Viral Fraction Too Low (<95%)

**Possible Causes**:
1. Heavy contamination level with poor VLP protocol
2. Bug in contamination reduction (check metadata)

**Solutions**:
```bash
# Check enrichment statistics
cat output/metadata/*_metadata.json | jq '.enrichment_statistics'

# Try clean contamination level
--contamination-level clean

# Verify protocol was applied
--vlp-protocol tangential_flow  # explicit protocol
```

### Issue: No Difference Between VLP and Bulk

**Possible Causes**:
1. Forgot to use `--no-vlp` for bulk
2. Very low initial contamination

**Solutions**:
```bash
# Verify bulk mode
python scripts/generate_fastq_dataset.py \
    --no-vlp \
    --contamination-level realistic  # ensure some contamination

# Check metadata
cat output/metadata/*_metadata.json | jq '.configuration'
```

---

## Summary

This tutorial demonstrated:

✅ **Generating datasets** for all 5 VLP protocols
✅ **Extracting metrics** from ground truth metadata
✅ **Comparing protocols** across multiple dimensions
✅ **Visualizing results** for interpretation
✅ **Common use cases** for protocol comparison
✅ **Troubleshooting** common issues

### Key Takeaways

1. **Different protocols have distinct trade-offs**: Purity vs recovery vs convenience
2. **ViroForge provides complete ground truth**: Enables rigorous validation
3. **Use realistic contamination levels**: Too clean = no differences to observe
4. **Consistent random seeds**: Essential for reproducible comparisons
5. **Validate your QC pipeline**: Ensure it detects protocol effects correctly

### Next Steps

- Run your own pipeline on generated datasets
- Compare taxonomic classification across protocols
- Analyze assembly quality differences
- Validate detection limits for rare viruses

---

## References

**Literature on VLP Protocols**:
- Lim ES et al. (2020) mSystems - "Protocol comparison for viral metagenomics"
- Thurber RV et al. (2009) AME - "Laboratory procedures for viral metagenomes"
- Shkoporov AN et al. (2018) Nat Protoc - "Viral metagenomics protocol"
- Roux S et al. (2016) Microbiome - "ViromeQC assessment tool"

**ViroForge Documentation**:
- [VLP Integration Guide](VLP_CONTAMINATION_INTEGRATION.md)
- [FASTQ Generation Guide](PHASE4_FASTQ_GENERATION.md)
- [Phase 5 Validation Report](PHASE5_TASK3_VALIDATION_REPORT.md)

---

**Questions?** See [GitHub Issues](https://github.com/shandley/viroforge/issues) or email scott.handley@wustl.edu
