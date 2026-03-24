#!/usr/bin/env python3
"""
Create key figures for ViroForge grant proposal (HVP R03).

Generates 2 publication-quality figures:
1. ViroForge Integrated Benchmarking Ecosystem
2. Sparse Data Challenge & Validation Strategy
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("colorblind")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2


def create_figure1_ecosystem():
    """
    Figure 1: ViroForge Integrated Benchmarking Ecosystem

    Shows the complete workflow: synthetic data → tools → validation → HVP
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Color scheme
    viroforge_color = '#2E86AB'  # Blue
    viromedata_color = '#A23B72'  # Purple
    validation_color = '#F18F01'  # Orange
    hvp_color = '#C73E1D'  # Red

    # Title
    ax.text(5, 9.5, 'ViroForge Integrated Benchmarking Ecosystem',
            ha='center', va='top', fontsize=14, fontweight='bold')

    # ========== AIM 1: ViroForge (Top Left) ==========
    # ViroForge box
    vf_box = FancyBboxPatch((0.3, 6.5), 3.0, 2.3,
                            boxstyle="round,pad=0.1",
                            edgecolor=viroforge_color,
                            facecolor=viroforge_color, alpha=0.15, linewidth=2)
    ax.add_patch(vf_box)

    ax.text(1.8, 8.4, 'Aim 1: ViroForge v2.0', ha='center', fontsize=11,
            fontweight='bold', color=viroforge_color)
    ax.text(1.8, 8.0, 'Synthetic Ground-Truth Data', ha='center', fontsize=9)

    # ViroForge capabilities (bullets)
    capabilities = [
        '• 28 curated collections',
        '• 8 HVP body sites',
        '• 5 sequencing platforms',
        '• DNA + RNA workflows',
        '• Longitudinal & batch effects',
        '• Complete ground truth'
    ]
    y_start = 7.5
    for i, cap in enumerate(capabilities):
        ax.text(0.5, y_start - i*0.25, cap, fontsize=7, va='top')

    # ========== AIM 2: ViromeData (Top Right) ==========
    # ViromeData box
    vd_box = FancyBboxPatch((5.5, 6.5), 3.8, 2.3,
                            boxstyle="round,pad=0.1",
                            edgecolor=viromedata_color,
                            facecolor=viromedata_color, alpha=0.15, linewidth=2)
    ax.add_patch(vd_box)

    ax.text(7.4, 8.4, 'Aim 2: ViromeData Infrastructure', ha='center',
            fontsize=11, fontweight='bold', color=viromedata_color)
    ax.text(7.4, 8.0, 'Standardized Benchmarking Platform', ha='center', fontsize=9)

    # ViromeData features
    features = [
        '• Extends SummarizedExperiment',
        '• Baltimore classification',
        '• RdRP classifier integration',
        '• Virus-host associations',
        '• Multi-tool interoperability',
        '• R + Python versions'
    ]
    y_start = 7.5
    for i, feat in enumerate(features):
        ax.text(5.7, y_start - i*0.25, feat, fontsize=7, va='top')

    # ========== AIM 3: Statistical Validation (Middle) ==========
    # Validation box
    val_box = FancyBboxPatch((1.5, 3.5), 6.5, 2.3,
                             boxstyle="round,pad=0.1",
                             edgecolor=validation_color,
                             facecolor=validation_color, alpha=0.15, linewidth=2)
    ax.add_patch(val_box)

    ax.text(4.75, 5.4, 'Aim 3: Statistical Method Validation', ha='center',
            fontsize=11, fontweight='bold', color=validation_color)
    ax.text(4.75, 5.0, 'Systematic Testing on Sparse Virome Data', ha='center', fontsize=9)

    # Three validation components (side by side)
    components = [
        ('Differential\nAbundance', ['DESeq2', 'ANCOM-BC', 'ALDEx2']),
        ('Diversity\nMetrics', ['Shannon', 'Bray-Curtis', 'Aitchison']),
        ('Normalization\nMethods', ['TSS', 'CSS', 'GMPR'])
    ]

    x_positions = [2.5, 4.75, 7.0]
    for x, (title, methods) in zip(x_positions, components):
        ax.text(x, 4.6, title, ha='center', fontsize=8, fontweight='bold')
        for i, method in enumerate(methods):
            ax.text(x, 4.3 - i*0.2, method, ha='center', fontsize=7)

    # Sparsity gradient note
    ax.text(4.75, 3.7, 'Test across sparsity gradient: 50% → 95% zeros',
            ha='center', fontsize=7, style='italic',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # ========== HVP Integration (Bottom) ==========
    # HVP box
    hvp_box = FancyBboxPatch((1.5, 0.5), 6.5, 2.3,
                             boxstyle="round,pad=0.1",
                             edgecolor=hvp_color,
                             facecolor=hvp_color, alpha=0.15, linewidth=2)
    ax.add_patch(hvp_box)

    ax.text(4.75, 2.4, 'HVP Consortium Integration', ha='center',
            fontsize=11, fontweight='bold', color=hvp_color)

    # HVP deliverables
    deliverables = [
        '✓ Tool validation for consortium workflows',
        '✓ Quality control standards for cross-site comparisons',
        '✓ Statistical best practices for HVP analyses',
        '✓ Power analysis & sample size recommendations',
        '✓ Integrated with CODCC data portal'
    ]
    y_start = 2.0
    for i, deliv in enumerate(deliverables):
        ax.text(2.0, y_start - i*0.25, deliv, fontsize=7, va='top')

    # ========== Arrows showing workflow ==========
    # ViroForge → ViromeData
    arrow1 = FancyArrowPatch((3.3, 7.6), (5.5, 7.6),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5)
    ax.add_patch(arrow1)
    ax.text(4.4, 7.8, 'feeds', ha='center', fontsize=7, style='italic')

    # Both → Validation
    arrow2 = FancyArrowPatch((1.8, 6.5), (3.0, 5.8),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5)
    ax.add_patch(arrow2)

    arrow3 = FancyArrowPatch((7.4, 6.5), (6.5, 5.8),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5)
    ax.add_patch(arrow3)

    # Validation → HVP
    arrow4 = FancyArrowPatch((4.75, 3.5), (4.75, 2.8),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5)
    ax.add_patch(arrow4)
    ax.text(5.2, 3.15, 'enables', ha='left', fontsize=7, style='italic')

    plt.tight_layout()
    return fig


def create_figure2_sparsity():
    """
    Figure 2: Sparse Data Challenge & Validation Strategy

    Shows the zero-inflation problem and ViroForge's solution.
    """
    fig = plt.figure(figsize=(10, 6))

    # Create 2x2 grid
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], width_ratios=[1.2, 1],
                          hspace=0.35, wspace=0.3)

    # ========== Panel A: Zero-Inflation Comparison ==========
    ax1 = fig.add_subplot(gs[0, 0])

    # Data types and their sparsity ranges
    data_types = ['Bulk\nRNA-seq', 'Bacterial\n16S rRNA', 'Bacterial\nWGS', 'Virome\nData']
    sparsity_min = [10, 55, 35, 70]
    sparsity_max = [40, 83, 89, 95]
    colors = ['#4A90E2', '#7B68EE', '#9370DB', '#C73E1D']

    y_pos = np.arange(len(data_types))

    # Plot ranges as horizontal bars
    for i, (dtype, smin, smax, color) in enumerate(zip(data_types, sparsity_min, sparsity_max, colors)):
        ax1.barh(i, smax - smin, left=smin, height=0.6,
                color=color, alpha=0.7, edgecolor='black', linewidth=1)
        # Add range text
        ax1.text(smax + 2, i, f'{smin}-{smax}%', va='center', fontsize=8)

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(data_types, fontsize=9)
    ax1.set_xlabel('Percentage of Zeros (%)', fontsize=10, fontweight='bold')
    ax1.set_xlim(0, 105)
    ax1.set_title('A. The Sparse Data Challenge', fontsize=11, fontweight='bold', pad=10)
    ax1.axvspan(90, 95, alpha=0.2, color='red', zorder=0)
    ax1.text(92.5, 3.5, 'Extreme\nsparsity', ha='center', fontsize=7,
            style='italic', color='darkred')
    ax1.grid(axis='x', alpha=0.3, linestyle='--')

    # ========== Panel B: Method Performance Question ==========
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis('off')

    # Question box
    question_box = FancyBboxPatch((0.05, 0.3), 0.9, 0.55,
                                  boxstyle="round,pad=0.05",
                                  edgecolor='#F18F01', facecolor='#F18F01',
                                  alpha=0.1, linewidth=2)
    ax2.add_patch(question_box)

    ax2.text(0.5, 0.8, 'B. Critical Question', ha='center', fontsize=11,
            fontweight='bold', transform=ax2.transAxes)

    question_text = (
        "Do compositionally-aware methods\n"
        "(ANCOM-BC, ALDEx2, Aitchison)\n\n"
        "maintain superior performance\n"
        "vs. borrowed methods\n"
        "(DESeq2, Bray-Curtis)\n\n"
        "at extreme virome sparsity\n"
        "(90-95% zeros)?"
    )
    ax2.text(0.5, 0.5, question_text, ha='center', va='center', fontsize=9,
            transform=ax2.transAxes, linespacing=1.5)

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)

    # ========== Panel C: ViroForge Validation Strategy ==========
    ax3 = fig.add_subplot(gs[1, :])

    ax3.text(0.5, 0.95, 'C. ViroForge Systematic Validation Strategy',
            ha='center', fontsize=11, fontweight='bold', transform=ax3.transAxes)

    # Three-stage workflow
    stages = [
        {
            'title': '1. Generate Synthetic Data',
            'x': 0.15,
            'items': [
                'Known composition',
                'Vary sparsity:\n50%, 70%, 90%, 95%',
                'Vary effect size:\n1.5×, 2×, 5×, 10×',
                'Vary sample size:\nn=5 to n=100'
            ],
            'color': '#2E86AB'
        },
        {
            'title': '2. Test Methods',
            'x': 0.5,
            'items': [
                'Differential abundance:\nDESeq2 vs ANCOM-BC\nvs ALDEx2',
                'Diversity metrics:\nBray-Curtis vs Aitchison',
                'Normalization:\nTSS vs CSS vs GMPR'
            ],
            'color': '#F18F01'
        },
        {
            'title': '3. Deliver Guidance',
            'x': 0.85,
            'items': [
                'Performance curves\nacross sparsity',
                'Method-switching\nthresholds',
                'Power analysis &\nsample size tables',
                'HVP best practices'
            ],
            'color': '#C73E1D'
        }
    ]

    for stage in stages:
        # Box for each stage
        box = FancyBboxPatch((stage['x']-0.12, 0.15), 0.24, 0.65,
                            boxstyle="round,pad=0.02",
                            edgecolor=stage['color'], facecolor=stage['color'],
                            alpha=0.15, linewidth=2, transform=ax3.transAxes)
        ax3.add_patch(box)

        # Title
        ax3.text(stage['x'], 0.75, stage['title'], ha='center', fontsize=9,
                fontweight='bold', transform=ax3.transAxes, color=stage['color'])

        # Items
        y_start = 0.65
        for i, item in enumerate(stage['items']):
            ax3.text(stage['x'], y_start - i*0.12, item, ha='center', fontsize=7,
                    transform=ax3.transAxes, linespacing=1.3)

    # Arrows between stages
    arrow_y = 0.47
    arrow1 = FancyArrowPatch((0.28, arrow_y), (0.37, arrow_y),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5,
                            transform=ax3.transAxes)
    ax3.add_patch(arrow1)

    arrow2 = FancyArrowPatch((0.63, arrow_y), (0.72, arrow_y),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='black', alpha=0.5,
                            transform=ax3.transAxes)
    ax3.add_patch(arrow2)

    # Bottom note
    ax3.text(0.5, 0.05, 'Ground truth enables definitive performance evaluation impossible with real data',
            ha='center', fontsize=8, style='italic', transform=ax3.transAxes,
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.axis('off')

    plt.tight_layout()
    return fig


def main():
    """Generate both figures."""
    print("Creating ViroForge grant proposal figures...")

    # Figure 1: Ecosystem
    print("\n[1/2] Creating Figure 1: ViroForge Integrated Benchmarking Ecosystem...")
    fig1 = create_figure1_ecosystem()
    fig1.savefig('Figure1_ViroForge_Ecosystem.png', dpi=300, bbox_inches='tight')
    fig1.savefig('Figure1_ViroForge_Ecosystem.pdf', bbox_inches='tight')
    print("  ✓ Saved: Figure1_ViroForge_Ecosystem.png (300 DPI)")
    print("  ✓ Saved: Figure1_ViroForge_Ecosystem.pdf")

    # Figure 2: Sparsity Challenge
    print("\n[2/2] Creating Figure 2: Sparse Data Challenge & Validation Strategy...")
    fig2 = create_figure2_sparsity()
    fig2.savefig('Figure2_Sparsity_Validation.png', dpi=300, bbox_inches='tight')
    fig2.savefig('Figure2_Sparsity_Validation.pdf', bbox_inches='tight')
    print("  ✓ Saved: Figure2_Sparsity_Validation.png (300 DPI)")
    print("  ✓ Saved: Figure2_Sparsity_Validation.pdf")

    print("\n" + "="*60)
    print("SUCCESS: Both figures generated!")
    print("="*60)
    print("\nFigure 1: Shows complete 3-aim integration workflow")
    print("  - Aim 1: ViroForge synthetic data generation")
    print("  - Aim 2: ViromeData infrastructure")
    print("  - Aim 3: Statistical validation")
    print("  - HVP integration at bottom")
    print("\nFigure 2: Highlights the sparse data challenge")
    print("  - Panel A: Virome sparsity vs other data types")
    print("  - Panel B: Central hypothesis question")
    print("  - Panel C: Three-stage validation strategy")
    print("\nBoth figures ready for grant proposal!")
    print("="*60)

    plt.show()


if __name__ == '__main__':
    main()
