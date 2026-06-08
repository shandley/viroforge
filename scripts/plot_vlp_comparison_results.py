#!/usr/bin/env python3
"""
Generate comparison plots and tables for VLP 0.2μm vs 0.45μm datasets.

Outputs:
  1. Stacked bar: Pre-VLP composition (all 4 datasets)
  2. Stacked bar: Post-VLP composition (all 4 datasets)
  3. Grouped bar: virome-qc module removal rates
  4. Heatmap: VLP contamination reduction by type
  5. Summary tables as PNG
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import json
import os

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

OUTDIR = "data/analysis_plots"
os.makedirs(OUTDIR, exist_ok=True)

COLORS = {
    'viral': '#2ca02c',
    'dark_matter': '#98df8a',
    'bacterial_background': '#1f77b4',
    'host_dna': '#d62728',
    'rrna': '#ff7f0e',
    'fungal_background': '#9467bd',
    'phix': '#e377c2',
    'reagent_bacteria': '#8c564b',
    'artifact_low_complexity': '#7f7f7f',
}

LABELS = {
    'viral': 'Known Viral',
    'dark_matter': 'Dark Matter',
    'bacterial_background': 'Bacteria',
    'host_dna': 'Host DNA',
    'rrna': 'rRNA',
    'fungal_background': 'Fungi',
    'phix': 'PhiX',
    'reagent_bacteria': 'Reagent',
    'artifact_low_complexity': 'Low Complexity',
}

CAT_ORDER = ['viral', 'dark_matter', 'bacterial_background', 'host_dna',
             'rrna', 'fungal_background', 'phix', 'artifact_low_complexity', 'reagent_bacteria']

# Load metadata
datasets = {
    'IBD 0.2μm': 'data/ibd_vlp_novaseq/metadata/ibd_gut_virome_inflammatory_bowel_disease_metadata.json',
    'IBD 0.45μm': 'data/ibd_vlp045_novaseq/metadata/ibd_gut_virome_inflammatory_bowel_disease_metadata.json',
    'Vag 0.2μm': 'data/vaginal_vlp_novaseq/metadata/vaginal_virome_healthy_metadata.json',
    'Vag 0.45μm': 'data/vaginal_vlp045_novaseq/metadata/vaginal_virome_healthy_metadata.json',
}

# Load virome-qc passports
passports = {
    'IBD 0.2μm': 'data/ibd_vlp_novaseq/analysis/virome_qc/passport.json',
    'IBD 0.45μm': 'data/ibd_vlp045_novaseq/analysis/virome_qc/passport.json',
    'Vag 0.2μm': 'data/vaginal_vlp_novaseq/analysis/virome_qc/passport.json',
    'Vag 0.45μm': 'data/vaginal_vlp045_novaseq/analysis/virome_qc/passport.json',
}

# Load FASTQ source counts
fastq_sources = {}
for name, meta_path in datasets.items():
    ds_dir = os.path.dirname(os.path.dirname(meta_path))
    import glob
    r1_files = glob.glob(f"{ds_dir}/fastq/*_R1.fastq")
    if r1_files:
        counts = {}
        with open(r1_files[0]) as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    for part in line.strip().split():
                        if part.startswith('source='):
                            src = part.split('=')[1]
                            counts[src] = counts.get(src, 0) + 1
        fastq_sources[name] = counts

meta = {}
for name, path in datasets.items():
    with open(path) as f:
        meta[name] = json.load(f)

qc = {}
for name, path in passports.items():
    with open(path) as f:
        qc[name] = json.load(f)


# ============================================================
# PLOT 1: Pre-VLP vs Post-VLP Composition (stacked bars)
# ============================================================
def plot_composition():
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    names = list(datasets.keys())

    contam_cats = ['host_dna', 'rrna', 'bacterial_background', 'fungal_background', 'phix', 'reagent_bacteria']

    # Pre-VLP
    ax = axes[0]
    ax.set_title('Pre-VLP Composition', fontweight='bold')
    x = np.arange(len(names))
    bottom = np.zeros(len(names))

    # Viral fraction
    for name_idx, name in enumerate(names):
        cr = meta[name]['enrichment_stats']['contamination_reduction']
        total_contam = cr['original_total_contamination']
        viral_frac = 1.0 - total_contam

    viral_vals = []
    for name in names:
        cr = meta[name]['enrichment_stats']['contamination_reduction']
        viral_vals.append((1.0 - cr['original_total_contamination']) * 100)
    ax.bar(x, viral_vals, 0.6, label='Viral (known + dark matter)', color=COLORS['viral'], bottom=bottom)
    bottom = np.array(viral_vals)

    for cat in contam_cats:
        vals = []
        for name in names:
            cr = meta[name]['enrichment_stats']['contamination_reduction']
            rt = cr.get('reduction_by_type', {})
            vals.append(rt.get(cat, {}).get('original_abundance', 0) * 100)
        ax.bar(x, vals, 0.6, label=LABELS.get(cat, cat), color=COLORS.get(cat, '#999'), bottom=bottom)
        bottom += np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha='right')
    ax.set_ylabel('Fraction (%)')
    ax.set_ylim(0, 105)

    # Post-VLP
    ax = axes[1]
    ax.set_title('Post-VLP Composition (Final FASTQ)', fontweight='bold')
    bottom = np.zeros(len(names))

    for cat in ['viral', 'dark_matter'] + contam_cats + ['artifact_low_complexity']:
        vals = []
        for name in names:
            total = sum(fastq_sources[name].values())
            vals.append(fastq_sources[name].get(cat, 0) * 100 / total)
        ax.bar(x, vals, 0.6, label=LABELS.get(cat, cat), color=COLORS.get(cat, '#999'), bottom=bottom)
        bottom += np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha='right')
    ax.set_ylabel('Fraction (%)')
    ax.set_ylim(0, 105)

    # Single legend below both panels
    handles1, labels1 = axes[0].get_legend_handles_labels()
    handles2, labels2 = axes[1].get_legend_handles_labels()
    # Combine unique labels
    all_handles = handles1 + [h for h, l in zip(handles2, labels2) if l not in labels1]
    all_labels = labels1 + [l for l in labels2 if l not in labels1]
    fig.legend(all_handles, all_labels, loc='lower center', ncol=5, fontsize=9,
              bbox_to_anchor=(0.5, -0.02))

    plt.suptitle('ViroForge Synthetic Dataset Composition: 0.2μm vs 0.45μm VLP', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(f'{OUTDIR}/plot1_composition_comparison.png')
    print(f'Saved: {OUTDIR}/plot1_composition_comparison.png')
    plt.close()


# ============================================================
# PLOT 2: virome-qc Module Removal Rates
# ============================================================
def plot_qc_modules():
    fig, ax = plt.subplots(figsize=(12, 6))
    names = list(qc.keys())
    modules = ['adapter', 'quality', 'complexity', 'dedup', 'contaminant', 'rrna', 'host']

    x = np.arange(len(modules))
    width = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5]

    colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']

    for i, name in enumerate(names):
        vals = []
        for mod in modules:
            ri = qc[name]['reads_input']
            for m in qc[name]['modules']:
                if m['name'] == mod:
                    vals.append(m['reads_removed'] * 100 / ri)
                    break
            else:
                vals.append(0)
        ax.bar(x + offsets[i] * width, vals, width, label=name, color=colors[i])

    ax.set_xticks(x)
    ax.set_xticklabels([m.title() for m in modules])
    ax.set_ylabel('Reads Removed (%)')
    ax.set_title('virome-qc Module Removal Rates by Dataset', fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 30)

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/plot2_qc_module_comparison.png')
    print(f'Saved: {OUTDIR}/plot2_qc_module_comparison.png')
    plt.close()


# ============================================================
# PLOT 3: VLP Reduction Heatmap
# ============================================================
def plot_vlp_heatmap():
    fig, ax = plt.subplots(figsize=(10, 5))

    names = list(datasets.keys())
    contam_cats = ['host_dna', 'rrna', 'bacterial_background', 'fungal_background', 'phix', 'reagent_bacteria']

    data = []
    for name in names:
        cr = meta[name]['enrichment_stats']['contamination_reduction']
        row = []
        for cat in contam_cats:
            rt = cr.get('reduction_by_type', {})
            row.append(rt.get(cat, {}).get('reduction_pct', 0))
        data.append(row)

    data = np.array(data)
    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)

    ax.set_xticks(range(len(contam_cats)))
    ax.set_xticklabels([LABELS.get(c, c) for c in contam_cats], rotation=30, ha='right')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)

    # Add text
    for i in range(len(names)):
        for j in range(len(contam_cats)):
            ax.text(j, i, f'{data[i, j]:.0f}%', ha='center', va='center', fontsize=10,
                    color='white' if data[i, j] > 60 else 'black')

    plt.colorbar(im, label='Removal Rate (%)')
    ax.set_title('VLP Contamination Reduction by Type and Pore Size', fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/plot3_vlp_reduction_heatmap.png')
    print(f'Saved: {OUTDIR}/plot3_vlp_reduction_heatmap.png')
    plt.close()


# ============================================================
# PLOT 4: Survival Rate Comparison
# ============================================================
def plot_survival():
    fig, ax = plt.subplots(figsize=(8, 5))
    names = list(qc.keys())
    survival = [qc[n]['survival_rate'] * 100 for n in names]

    colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']
    bars = ax.bar(range(len(names)), survival, color=colors)

    for bar, val in zip(bars, survival):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{val:.1f}%', ha='center', fontsize=11, fontweight='bold')

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_ylabel('Survival Rate (%)')
    ax.set_title('virome-qc Survival Rate by Dataset', fontweight='bold')
    ax.set_ylim(0, 100)
    ax.axhline(y=80, color='gray', linestyle='--', alpha=0.5, label='80% threshold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/plot4_survival_comparison.png')
    print(f'Saved: {OUTDIR}/plot4_survival_comparison.png')
    plt.close()


# ============================================================
# TABLE: Summary as PNG
# ============================================================
def table_summary():
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')

    names = list(datasets.keys())
    headers = ['Metric'] + names

    rows = [
        ['Viral genomes'] + ['89' if 'IBD' in n else '26' for n in names],
        ['Dark matter genomes'] + ['44' if 'IBD' in n else '13' for n in names],
        ['VLP pore size'] + ['0.2 μm' if '0.2' in n else '0.45 μm' for n in names],
        ['Total reads (R1)'] + [f"{sum(fastq_sources[n].values()):,}" for n in names],
    ]

    # Source composition
    for cat in ['viral', 'dark_matter', 'bacterial_background', 'host_dna', 'rrna', 'fungal_background']:
        row = [LABELS.get(cat, cat)]
        for n in names:
            total = sum(fastq_sources[n].values())
            count = fastq_sources[n].get(cat, 0)
            row.append(f'{count:,} ({count*100/total:.1f}%)')
        rows.append(row)

    # VLP stats
    rows.append(['VLP reduction'] + [
        f"{meta[n]['enrichment_stats']['contamination_reduction']['overall_reduction_factor']*100:.1f}%"
        for n in names])

    # virome-qc stats
    rows.append(['virome-qc survival'] + [f"{qc[n]['survival_rate']*100:.1f}%" for n in names])

    for mod_name in ['dedup', 'host', 'rrna', 'contaminant']:
        row = [f'QC: {mod_name}']
        for n in names:
            for m in qc[n]['modules']:
                if m['name'] == mod_name:
                    ri = qc[n]['reads_input']
                    row.append(f"{m['reads_removed']:,} ({m['reads_removed']*100/ri:.1f}%)")
                    break
        rows.append(row)

    table = ax.table(cellText=rows, colLabels=headers, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)

    # Style header
    for j in range(len(headers)):
        table[(0, j)].set_facecolor('#4C72B0')
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    ax.set_title('ViroForge + virome-qc Results Summary', fontsize=14, fontweight='bold', pad=20)
    plt.savefig(f'{OUTDIR}/table_summary.png')
    print(f'Saved: {OUTDIR}/table_summary.png')
    plt.close()


if __name__ == '__main__':
    print("Generating plots and tables...")
    plot_composition()
    plot_qc_modules()
    plot_vlp_heatmap()
    plot_survival()
    table_summary()
    print(f"\nAll outputs saved to {OUTDIR}/")
