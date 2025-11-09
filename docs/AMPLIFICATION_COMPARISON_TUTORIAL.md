# Library Preparation Amplification Comparison Tutorial

**ViroForge Phase 6 Feature**

This tutorial demonstrates how to compare different library preparation amplification methods and understand their biases on viral community composition.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Understanding Amplification Methods](#understanding-amplification-methods)
4. [Step-by-Step Tutorial](#step-by-step-tutorial)
5. [Interpreting Results](#interpreting-results)
6. [Common Use Cases](#common-use-cases)
7. [Troubleshooting](#troubleshooting)

---

## Introduction

### Why Compare Amplification Methods?

Different library preparation amplification methods introduce distinct biases:

- **Length bias**: Do shorter genomes amplify better?
- **GC bias**: How does GC content affect amplification?
- **Stochasticity**: How variable is the amplification?
- **Magnitude of bias**: How much does composition change?

ViroForge enables systematic comparison with known ground truth, helping you:

1. **Choose the optimal method**: Select amplification strategy for your samples
2. **Understand bias magnitude**: Know how much your data may deviate from reality
3. **Validate analysis pipelines**: Ensure your methods can handle amplification bias
4. **Design experiments**: Plan for biomass requirements vs acceptable bias

---

## Quick Start

### Generate All 6 Amplification Methods in One Command

```bash
# Generate datasets for all amplification methods
python scripts/batch_generate_fastq.py \
    --preset amplification-comparison \
    --output data/amplification_comparison
```

This generates 6 datasets from the same gut virome collection:

1. **none** - No amplification (baseline)
2. **rdab** - RdAB 40 cycles (standard virome protocol)
3. **rdab-30** - RdAB 30 cycles (moderate bias)
4. **mda** - MDA 4 hours (low biomass)
5. **mda-long** - MDA overnight (very low biomass)
6. **linker** - Linker-based (minimal bias)

**Runtime**: ~30 minutes (10x coverage, dry-run: ~3 minutes)

---

## Understanding Amplification Methods

### Method Overview

| Method | Cycles/Time | Length Bias | GC Bias | Stochasticity | Use Case |
|--------|-------------|-------------|---------|---------------|----------|
| **None** | N/A | None | None | None | High biomass control |
| **RdAB 40** | 40 cycles | High | Moderate | Low | Standard virome |
| **RdAB 30** | 30 cycles | Moderate | Moderate | Low | Higher biomass virome |
| **MDA 4h** | 4 hours | None | Very High | High | Low biomass samples |
| **MDA 16h** | 16 hours | None | Extreme | Very High | Ultra-low biomass |
| **Linker** | 20 cycles | None | Low | Low | Modern library prep kits |

### Bias Mechanisms

#### RdAB (Random RT + dsDNA synthesis + PCR)

**How it works**:
1. Random reverse transcription
2. Second-strand synthesis
3. PCR amplification (30-45 cycles)

**Biases introduced**:
- **Length bias**: Shorter genomes amplify exponentially better
  - Each PCR cycle favors complete amplification of shorter templates
  - ~1.5-2x advantage per kb difference
- **GC bias**: Extreme GC content (< 30% or > 70%) amplifies poorly
  - Optimal at ~50% GC
  - 2-10x differences across GC range

**Expected effects**:
```
Original:  100 bp genome at 1% → RdAB → 8% (8x enrichment)
Original: 1000 bp genome at 1% → RdAB → 2% (2x enrichment)
Original: 5000 bp genome at 1% → RdAB → 0.5% (0.5x depletion)
```

#### MDA (Multiple Displacement Amplification)

**How it works**:
1. φ29 polymerase with random hexamer priming
2. Isothermal amplification (2-16 hours)
3. Extensive branching amplification

**Biases introduced**:
- **Extreme GC bias**: φ29 struggles with high GC much more than Taq
  - 10-1000x differences possible
  - Optimal ~40% GC
- **High stochasticity**: Random priming creates uneven amplification
  - Some templates dominate by chance
  - CV typically 0.3-0.5
- **No length bias**: Isothermal, no cycle-based advantage

**Expected effects**:
```
Original: 40% GC genome at 1% → MDA → 5% (5x enrichment)
Original: 60% GC genome at 1% → MDA → 0.01% (100x depletion!)
```

#### Linker-Based (Adapter ligation + PCR)

**How it works**:
1. Ligate adapters to all fragments
2. PCR amplification (10-25 cycles)
3. All fragments have equal primer binding

**Biases introduced**:
- **No length bias**: All fragments have adapters
- **Moderate GC bias**: PCR still affected by GC, but less cycles
  - 2-5x less bias than RdAB
- **Minimal overall impact**: Modern "best practice" method

**Expected effects**:
```
Original: Low GC genome at 1% → Linker → 0.8% (minor depletion)
Original: High GC genome at 1% → Linker → 0.9% (minor depletion)
Original: Optimal GC at 1% → Linker → 1.1% (minimal change)
```

---

## Step-by-Step Tutorial

### Step 1: Generate Comparison Datasets

```bash
# Generate all 6 amplification methods
python scripts/batch_generate_fastq.py \
    --preset amplification-comparison \
    --output data/amplification_comparison

# Or generate individually for more control
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/baseline \
    --coverage 10 \
    --amplification none

python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output data/rdab \
    --coverage 10 \
    --amplification rdab

# ... etc for other methods
```

### Step 2: Examine Metadata

Each dataset includes amplification statistics in `metadata.json`:

```bash
# Extract amplification stats for all methods
for method in none rdab rdab-30 mda mda-long linker; do
    echo "=== $method ==="
    cat data/amplification_comparison/collection_9_cov10x_novaseq_tangential_flow_${method}/metadata/*.json | \
        jq '.amplification_stats'
done
```

**Example output**:

```json
{
  "method": "rdab",
  "bias_applied": true,
  "mean_fold_change": 18.5,
  "max_fold_change": 720.5,
  "min_fold_change": 0.0,
  "abundance_correlation": 0.008
}
```

**Key metrics**:
- `mean_fold_change`: Average amplification bias (how much abundances changed)
- `max/min_fold_change`: Range of bias (some genomes enriched 720x!)
- `abundance_correlation`: How much original composition is preserved (0.008 = almost no correlation!)

### Step 3: Compare Abundances

Create a Python script to extract and compare abundances:

```python
#!/usr/bin/env python3
"""
Compare amplification method effects on viral abundances.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_amplification_data(base_dir, methods):
    """Load abundance data for all amplification methods."""
    data = []

    for method in methods:
        # Find metadata file
        if method == 'none':
            pattern = f"collection_9_cov10x_novaseq_tangential_flow"
        else:
            pattern = f"collection_9_cov10x_novaseq_tangential_flow_{method}"

        metadata_dir = base_dir / pattern / 'metadata'
        metadata_file = list(metadata_dir.glob('*_metadata.json'))[0]

        with open(metadata_file) as f:
            metadata = json.load(f)

        # Extract viral genome abundances
        for seq in metadata['sequences']:
            if seq['sequence_type'] == 'viral':
                data.append({
                    'method': method,
                    'genome_id': seq['genome_id'],
                    'genome_name': seq['genome_name'],
                    'length': seq['length'],
                    'abundance': seq['relative_abundance'],
                    'family': seq.get('family', 'Unknown')
                })

    return pd.DataFrame(data)

def calculate_bias_metrics(df):
    """Calculate bias metrics for each method vs baseline."""
    baseline = df[df['method'] == 'none'].set_index('genome_id')['abundance']

    bias_metrics = []
    for method in df['method'].unique():
        if method == 'none':
            continue

        method_df = df[df['method'] == method].set_index('genome_id')

        # Calculate fold changes
        fold_changes = method_df['abundance'] / baseline

        metrics = {
            'method': method,
            'mean_fold_change': fold_changes.mean(),
            'median_fold_change': fold_changes.median(),
            'max_fold_change': fold_changes.max(),
            'min_fold_change': fold_changes.min(),
            'correlation': method_df['abundance'].corr(baseline),
            'n_enriched_10x': (fold_changes > 10).sum(),
            'n_depleted_10x': (fold_changes < 0.1).sum()
        }

        bias_metrics.append(metrics)

    return pd.DataFrame(bias_metrics)

def plot_abundance_comparison(df, output_file='amplification_comparison.png'):
    """Plot abundance comparison across methods."""
    # Create wide format for plotting
    wide_df = df.pivot(index='genome_id', columns='method', values='abundance')

    # Reorder columns
    method_order = ['none', 'linker', 'rdab-30', 'rdab', 'mda', 'mda-long']
    wide_df = wide_df[[m for m in method_order if m in wide_df.columns]]

    # Plot heatmap
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Heatmap of abundances
    sns.heatmap(
        wide_df.T.apply(lambda x: x / x.sum(), axis=1),  # Normalize
        ax=axes[0],
        cmap='YlOrRd',
        cbar_kws={'label': 'Relative Abundance'},
        xticklabels=False
    )
    axes[0].set_title('Viral Genome Abundances Across Amplification Methods')
    axes[0].set_ylabel('Amplification Method')

    # Fold change heatmap (relative to baseline)
    baseline = wide_df['none']
    fold_changes = wide_df.div(baseline, axis=0).iloc[:, 1:]  # Skip 'none'

    sns.heatmap(
        fold_changes.T,
        ax=axes[1],
        cmap='RdBu_r',
        center=1.0,
        vmin=0.01,
        vmax=100,
        norm=plt.matplotlib.colors.LogNorm(),
        cbar_kws={'label': 'Fold Change vs Baseline'},
        xticklabels=False
    )
    axes[1].set_title('Amplification Bias (Fold Change vs No Amplification)')
    axes[1].set_ylabel('Amplification Method')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {output_file}")

def plot_length_bias(df, output_file='length_bias.png'):
    """Plot length bias for each method."""
    baseline = df[df['method'] == 'none'].set_index('genome_id')

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    methods = ['linker', 'rdab-30', 'rdab', 'mda', 'mda-long']

    for i, method in enumerate(methods):
        method_df = df[df['method'] == method].set_index('genome_id')

        # Calculate fold change
        fold_change = method_df['abundance'] / baseline['abundance']

        # Plot
        axes[i].scatter(
            method_df['length'],
            fold_change,
            alpha=0.5,
            s=50
        )
        axes[i].axhline(y=1, color='red', linestyle='--', label='No bias')
        axes[i].set_xlabel('Genome Length (bp)')
        axes[i].set_ylabel('Fold Change vs Baseline')
        axes[i].set_title(f'{method}')
        axes[i].set_yscale('log')
        axes[i].set_xscale('log')
        axes[i].legend()
        axes[i].grid(True, alpha=0.3)

    # Hide last subplot
    axes[-1].axis('off')

    plt.suptitle('Length Bias in Amplification Methods', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {output_file}")

# Main analysis
if __name__ == '__main__':
    base_dir = Path('data/amplification_comparison')
    methods = ['none', 'rdab', 'rdab-30', 'mda', 'mda-long', 'linker']

    print("Loading amplification comparison data...")
    df = load_amplification_data(base_dir, methods)

    print("\nCalculating bias metrics...")
    bias_metrics = calculate_bias_metrics(df)
    print("\nBias Metrics:")
    print(bias_metrics.to_string(index=False))

    print("\nGenerating plots...")
    plot_abundance_comparison(df)
    plot_length_bias(df)

    print("\n✓ Analysis complete!")
```

Save as `scripts/analyze_amplification_comparison.py` and run:

```bash
python scripts/analyze_amplification_comparison.py
```

### Step 4: Interpret Results

**Expected patterns**:

1. **RdAB methods** (rdab, rdab-30):
   - Shorter genomes enriched (> 1x fold change)
   - Longer genomes depleted (< 1x fold change)
   - Extreme GC genomes depleted
   - rdab-30 shows weaker bias than rdab

2. **MDA methods** (mda, mda-long):
   - No length bias
   - Extreme GC bias (high GC severely depleted)
   - High stochasticity (scatter even at same length/GC)
   - mda-long shows stronger bias than mda

3. **Linker method**:
   - Minimal length bias
   - Weak GC bias
   - Most similar to baseline

4. **Baseline (none)**:
   - Original community composition
   - Use as reference for all comparisons

---

## Interpreting Results

### Key Questions to Ask

#### 1. How much does composition change?

**Check**: `abundance_correlation` in metadata

```
Correlation with baseline:
- none:     1.000 (identical)
- linker:   0.850 (high preservation)
- rdab-30:  0.120 (low preservation)
- rdab:     0.008 (almost no preservation!)
- mda:      0.005 (no preservation)
- mda-long: 0.001 (complete distortion)
```

**Interpretation**:
- > 0.8: Composition well-preserved
- 0.5-0.8: Moderate distortion
- 0.2-0.5: Severe distortion
- < 0.2: Composition unrecognizable

#### 2. Which genomes are most affected?

**Check**: `max_fold_change` and `min_fold_change`

```
Fold change ranges:
- rdab: 0.0 - 720x (some genomes eliminated, others 720x enriched!)
- mda:  0.001 - 850x (even more extreme)
- linker: 0.5 - 2.5x (modest changes)
```

**Interpretation**:
- RdAB: Short, optimal-GC genomes become dominant
- MDA: Low-GC genomes take over, high-GC disappear
- Linker: Community structure largely maintained

#### 3. Can I trust diversity estimates?

**Check**: Number of genomes enriched/depleted > 10x

```
Genomes altered > 10x:
- rdab:     45 enriched, 67 depleted (84% of community distorted!)
- mda:      52 enriched, 71 depleted (92% distorted)
- rdab-30:  12 enriched, 23 depleted (26% distorted)
- linker:   0 enriched, 0 depleted (0% severely distorted)
```

**Interpretation**:
- **High distortion** (rdab, mda): Alpha/beta diversity estimates unreliable
- **Moderate distortion** (rdab-30): Use with caution, report as exploratory
- **Low distortion** (linker): Diversity estimates reasonably accurate

#### 4. Does method affect biological conclusions?

**Example**: Comparing diseased vs healthy gut

If your biological effect is smaller than amplification bias, conclusions may be spurious!

```
True difference (healthy vs diseased): 2x abundance change
Amplification noise: 10-700x changes

→ Amplification swamps biological signal!
```

**Recommendation**:
- Use consistent amplification across all samples
- Choose lowest-bias method for biomass level
- Report amplification method in publications
- Consider amplification as potential confounder

---

## Common Use Cases

### Use Case 1: Method Selection for Study Design

**Scenario**: Planning a virome study, need to choose amplification method

**Workflow**:
```bash
# Generate comparison
python scripts/batch_generate_fastq.py \
    --preset amplification-comparison \
    --output method_comparison

# Analyze bias
python scripts/analyze_amplification_comparison.py

# Decision tree:
# - High biomass (>10ng)? → Use 'none' or 'linker'
# - Medium biomass (1-10ng)? → Use 'rdab-30' or 'linker'
# - Low biomass (100pg-1ng)? → Use 'rdab' (accept bias)
# - Very low (<100pg)? → Use 'mda' (extreme bias, but necessary)
```

### Use Case 2: Pipeline Validation

**Scenario**: Ensure analysis pipeline handles amplification bias

**Test**:
1. Generate datasets with different amplification
2. Run your pipeline on each
3. Check if pipeline reports different compositions
4. If pipeline cannot detect amplification effects, it may miss real biological variation!

```bash
# Generate test datasets
python scripts/batch_generate_fastq.py \
    --preset amplification-comparison \
    --output pipeline_validation

# Run your pipeline
for method in none rdab mda; do
    your_pipeline.sh pipeline_validation/*${method}/fastq/*.fastq
done

# Compare outputs - should show different compositions!
```

### Use Case 3: Debiasing Method Development

**Scenario**: Developing computational methods to correct amplification bias

**Approach**:
```bash
# Generate ground truth + biased datasets
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output ground_truth \
    --amplification none

python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output biased \
    --amplification rdab

# Train debiasing model:
# Input: biased abundances + genome properties (length, GC)
# Target: ground_truth abundances
# Model: Predict correction factors based on length and GC
```

### Use Case 4: Multi-Study Comparison

**Scenario**: Comparing results across studies using different methods

**Analysis**:
```bash
# Simulate each study's method
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output study_A \
    --amplification rdab  # Study A used RdAB

python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output study_B \
    --amplification mda  # Study B used MDA

# Compare - expect different compositions despite same biological sample!
```

**Expected finding**: Studies using different amplification may report inconsistent results even for identical samples.

**Recommendation**: Meta-analyses should account for amplification method as a batch effect.

---

## Troubleshooting

### Issue: "Amplification stats show no bias applied"

**Cause**: Used `--amplification none`

**Solution**: Specify a bias method:
```bash
python scripts/generate_fastq_dataset.py \
    --amplification rdab  # not 'none'
```

### Issue: "Results don't match expected bias patterns"

**Cause**: Random seed variation

**Solution**: Use fixed seed for reproducibility:
```bash
python scripts/generate_fastq_dataset.py \
    --amplification rdab \
    --seed 42  # Same seed = same bias
```

### Issue: "Fold changes seem too extreme"

**Cause**: This is expected! Amplification bias can be 100-1000x

**Verification**:
```bash
# Check metadata
cat output/metadata/*.json | jq '.amplification_stats.max_fold_change'

# For RdAB: expect 100-1000x
# For MDA: expect 100-10000x
# This is realistic!
```

### Issue: "Can't reproduce published results"

**Possible causes**:
1. Different amplification method used
2. Different number of cycles
3. Different polymerase/kit

**Solution**: Report exact amplification parameters in methods:
```
"Libraries were prepared using RdAB with 40 cycles of PCR (ViroForge simulation: --amplification rdab)"
```

---

## Best Practices

### For Experimental Design

1. **Match amplification to biomass**:
   - High biomass (>10ng): Avoid amplification if possible
   - Medium (1-10ng): Use minimal bias method (linker or rdab-30)
   - Low (<1ng): Accept bias, use rdab or mda

2. **Keep amplification consistent**:
   - Use same method for all samples in a study
   - Use same number of cycles
   - Process all samples in same batch if possible

3. **Report amplification method**:
   - Include in methods section
   - Specify cycles/time
   - Note as potential limitation

### For Computational Analysis

1. **Account for bias**:
   - Include amplification method as covariate
   - Consider debiasing if developing new methods
   - Report correlation with ground truth

2. **Validate with ViroForge**:
   - Test pipeline on amplified datasets
   - Ensure methods detect amplification effects
   - Compare performance across amplification levels

3. **Don't over-interpret**:
   - Treat species-level abundances as rough estimates
   - Focus on robust patterns (e.g., family-level changes)
   - Consider amplification noise when assessing significance

---

## Next Steps

- **Tutorial**: [VLP Protocol Comparison](VLP_PROTOCOL_COMPARISON_TUTORIAL.md)
- **Documentation**: [PHASE4_FASTQ_GENERATION](PHASE4_FASTQ_GENERATION.md)
- **Quick Start**: [QUICKSTART](QUICKSTART.md)

---

## References

### Key Literature on Amplification Bias

1. **Kim et al. (2013)** "Amplification methods bias metagenomic libraries of clones"
   *Nature Methods* 10:47-48
   - First systematic demonstration of RdAB bias

2. **Marine et al. (2014)** "Evaluation of a transposase protocol for viral metagenomics"
   *PeerJ* 2:e868
   - Comparison of linker-based vs RdAB methods

3. **Lasken & Stockwell (2007)** "Mechanism of chimera formation during MDA"
   *Nature Reviews Microbiology* 5:755-763
   - Comprehensive review of MDA biases

4. **Yilmaz et al. (2010)** "Evaluating MDA for metagenomics"
   *PLoS ONE* 5:e15533
   - Quantification of MDA bias in metagenomic studies

---

**Questions or issues?** Open an issue on [GitHub](https://github.com/shandley/viroforge/issues)
