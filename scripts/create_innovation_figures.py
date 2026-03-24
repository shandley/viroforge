#!/usr/bin/env python3
"""
Create ViroForge innovation figures for grant proposal.

Generates 2 figures:
1. ViroForge Workflow & Innovative Capabilities Overview
2. Contamination Modeling Details
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import numpy as np
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("colorblind")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2


def create_workflow_figure():
    """
    Figure 1: ViroForge Complete Workflow & Innovative Capabilities

    Shows the end-to-end virome data generation workflow with
    ViroForge's key innovations highlighted.
    """
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Title
    ax.text(7, 7.5, 'ViroForge: Complete Virome Workflow Modeling',
            ha='center', va='top', fontsize=16, fontweight='bold')

    # ========== WORKFLOW STAGES (Top) ==========
    stages = [
        {
            'x': 1.5,
            'title': '1. VLP\nEnrichment',
            'details': [
                '5 Protocols:',
                '• Tangential flow',
                '• Ultracentrifuge',
                '• Density gradient',
                '• Size exclusion',
                '• PEG precipitation'
            ],
            'color': '#2E86AB',
            'innovation': 'Size biases\nRecovery rates'
        },
        {
            'x': 4.5,
            'title': '2. Contamination',
            'details': [
                'Realistic sources:',
                '• Host DNA (5-40%)',
                '• rRNA (80-95%)',
                '• Bacterial reagents',
                '• PhiX control'
            ],
            'color': '#C73E1D',
            'innovation': 'Complete\nmanifest'
        },
        {
            'x': 7.5,
            'title': '3. Library Prep',
            'details': [
                'Amplification:',
                '• MDA',
                '• SISPA',
                '• TruSeq',
                'GC bias',
                'Chimera formation'
            ],
            'color': '#F18F01',
            'innovation': 'Sequence-\ndependent bias'
        },
        {
            'x': 10.5,
            'title': '4. Sequencing',
            'details': [
                '5 Platforms:',
                '• NovaSeq',
                '• MiSeq/HiSeq',
                '• PacBio HiFi',
                '• ONT',
                'Error models'
            ],
            'color': '#A23B72',
            'innovation': 'Long-read\nsupport'
        },
        {
            'x': 13,
            'title': '5. Ground\nTruth',
            'details': [
                'Complete metadata:',
                '• Read provenance',
                '• Contamination',
                '• Expected coverage',
                '• Benchmarking'
            ],
            'color': '#06A77D',
            'innovation': 'Full\ntraceability'
        }
    ]

    for i, stage in enumerate(stages):
        # Stage box
        box_width = 2.0 if i < 4 else 1.5
        box = FancyBboxPatch((stage['x'] - 0.75, 4.5), box_width, 2.3,
                            boxstyle="round,pad=0.1",
                            edgecolor=stage['color'],
                            facecolor=stage['color'], alpha=0.15, linewidth=2)
        ax.add_patch(box)

        # Title
        ax.text(stage['x'], 6.6, stage['title'], ha='center', fontsize=10,
                fontweight='bold', color=stage['color'])

        # Details
        y_start = 6.3
        for j, detail in enumerate(stage['details']):
            ax.text(stage['x'], y_start - j*0.25, detail, ha='center', fontsize=7)

        # Innovation badge
        innovation_box = FancyBboxPatch((stage['x'] - 0.5, 4.6), 1.0, 0.5,
                                       boxstyle="round,pad=0.05",
                                       edgecolor=stage['color'],
                                       facecolor='yellow', alpha=0.3, linewidth=1.5)
        ax.add_patch(innovation_box)
        ax.text(stage['x'], 4.85, stage['innovation'], ha='center', fontsize=6,
                fontweight='bold', color=stage['color'], linespacing=1.2)

        # Arrows between stages
        if i < len(stages) - 1:
            arrow = FancyArrowPatch((stage['x'] + 0.85, 5.6),
                                   (stages[i+1]['x'] - 0.85, 5.6),
                                   arrowstyle='->', mutation_scale=20,
                                   linewidth=2.5, color='black', alpha=0.6)
            ax.add_patch(arrow)

    # ========== RNA VIROME BRANCH (Middle) ==========
    rna_y = 3.5
    rna_box = FancyBboxPatch((0.5, rna_y - 0.8), 5.5, 1.3,
                            boxstyle="round,pad=0.1",
                            edgecolor='#8B4513', facecolor='#8B4513',
                            alpha=0.1, linewidth=2, linestyle='--')
    ax.add_patch(rna_box)

    ax.text(3.25, rna_y + 0.3, 'RNA Virome Workflow (Optional Branch)',
            ha='center', fontsize=11, fontweight='bold', color='#8B4513')

    rna_steps = [
        ('Reverse\nTranscription', 'Virus-type\nspecific\nefficiency'),
        ('rRNA\nDepletion', 'Ribo-Zero\nRiboMinus\n90-95% removal'),
        ('RNA\nDegradation', 'Fragment\npatterns\n5\'/3\' bias')
    ]

    for i, (step, detail) in enumerate(rna_steps):
        x_pos = 1.5 + i * 1.6
        # Step box
        step_box = FancyBboxPatch((x_pos - 0.4, rna_y - 0.6), 0.8, 0.5,
                                 boxstyle="round,pad=0.05",
                                 edgecolor='#8B4513', facecolor='white',
                                 linewidth=1.5)
        ax.add_patch(step_box)
        ax.text(x_pos, rna_y - 0.35, step, ha='center', fontsize=7,
                fontweight='bold', color='#8B4513', linespacing=1.2)

        # Detail text
        ax.text(x_pos, rna_y - 0.8, detail, ha='center', fontsize=6,
                color='#8B4513', linespacing=1.1)

        # Arrow
        if i < len(rna_steps) - 1:
            arrow = FancyArrowPatch((x_pos + 0.45, rna_y - 0.35),
                                   (x_pos + 1.15, rna_y - 0.35),
                                   arrowstyle='->', mutation_scale=15,
                                   linewidth=1.5, color='#8B4513', alpha=0.6)
            ax.add_patch(arrow)

    # ========== KEY INNOVATIONS (Bottom) ==========
    ax.text(7, 2.2, 'ViroForge Innovations vs. Other Simulators',
            ha='center', fontsize=12, fontweight='bold')

    innovations = [
        {
            'feature': 'VLP Enrichment',
            'viroforge': 'Yes (5 protocols)',
            'others': 'No',
            'impact': 'Protocol batch effects'
        },
        {
            'feature': 'Contamination',
            'viroforge': 'Yes (complete manifest)',
            'others': 'No',
            'impact': 'QC validation'
        },
        {
            'feature': 'RNA Workflow',
            'viroforge': 'Yes (RT, rRNA depletion)',
            'others': 'No',
            'impact': 'RNA virome support'
        },
        {
            'feature': 'Long-read',
            'viroforge': 'Yes (PacBio, ONT)',
            'others': 'Limited',
            'impact': 'Hybrid assembly'
        },
        {
            'feature': 'Ground Truth',
            'viroforge': 'Complete (every stage)',
            'others': 'Partial',
            'impact': 'Full traceability'
        }
    ]

    # Table header
    headers = ['Feature', 'ViroForge', 'InSilicoSeq/CAMISIM', 'HVP Impact']
    x_positions = [1.5, 5.0, 8.5, 11.5]
    for x, header in zip(x_positions, headers):
        ax.text(x, 1.7, header, ha='center', fontsize=9, fontweight='bold')

    # Table rows
    for i, innov in enumerate(innovations):
        y_pos = 1.4 - i * 0.25

        # Feature
        ax.text(x_positions[0], y_pos, innov['feature'], ha='left', fontsize=8)

        # ViroForge (green checkmark style)
        ax.text(x_positions[1], y_pos, innov['viroforge'], ha='center', fontsize=8,
                color='darkgreen', fontweight='bold')

        # Others (red/orange)
        color = 'darkred' if innov['others'] == 'No' else 'orange'
        ax.text(x_positions[2], y_pos, innov['others'], ha='center', fontsize=8,
                color=color, fontweight='bold')

        # Impact
        ax.text(x_positions[3], y_pos, innov['impact'], ha='center', fontsize=8,
                style='italic')

    # Bottom banner
    fig.text(0.5, 0.02,
            'ViroForge v0.11.0: Only simulator modeling complete virome experimental workflow | '
            'Production-ready | 20 curated collections | 14,423 genomes',
            ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    return fig


def create_contamination_figure():
    """
    Figure 2: ViroForge Contamination Modeling Details

    Shows contamination sources, abundances, and properties for
    DNA vs RNA virome workflows.
    """
    fig = plt.figure(figsize=(14, 8))

    # Create 2x2 grid
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1, 1],
                          hspace=0.35, wspace=0.3)

    # ========== Title ==========
    fig.suptitle('ViroForge Contamination Modeling: Realistic Non-Viral Background',
                fontsize=14, fontweight='bold', y=0.98)

    # ========== Panel A: DNA Virome Composition ==========
    ax1 = fig.add_subplot(gs[0, 0])

    dna_labels = ['Viral\nReads', 'Host\nDNA', 'Bacterial\nReagents', 'PhiX\nControl', 'Other']
    dna_sizes = [30, 40, 15, 10, 5]
    dna_colors = ['#2E86AB', '#C73E1D', '#F18F01', '#A23B72', '#999999']

    wedges, texts, autotexts = ax1.pie(dna_sizes, labels=dna_labels,
                                        autopct='%d%%', colors=dna_colors,
                                        startangle=90, textprops={'fontsize': 9})

    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)

    ax1.set_title('A. DNA Virome Composition\n(After VLP enrichment)',
                 fontsize=11, fontweight='bold', pad=10)

    # ========== Panel B: RNA Virome Before/After Depletion ==========
    ax2 = fig.add_subplot(gs[0, 1])

    categories = ['Before\nrRNA Depletion', 'After\nrRNA Depletion']
    viral = [1, 20]
    rrna = [90, 10]
    host = [5, 40]
    bacterial = [3, 25]
    other = [1, 5]

    x = np.arange(len(categories))
    width = 0.5

    p1 = ax2.bar(x, viral, width, label='Viral Reads', color='#2E86AB')
    p2 = ax2.bar(x, rrna, width, bottom=viral, label='rRNA', color='#8B4513')
    p3 = ax2.bar(x, host, width, bottom=np.array(viral)+np.array(rrna),
                label='Host DNA', color='#C73E1D')
    p4 = ax2.bar(x, bacterial, width,
                bottom=np.array(viral)+np.array(rrna)+np.array(host),
                label='Bacterial', color='#F18F01')
    p5 = ax2.bar(x, other, width,
                bottom=np.array(viral)+np.array(rrna)+np.array(host)+np.array(bacterial),
                label='Other', color='#999999')

    ax2.set_ylabel('Percentage (%)', fontsize=10, fontweight='bold')
    ax2.set_title('B. RNA Virome: Impact of rRNA Depletion',
                 fontsize=11, fontweight='bold', pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories, fontsize=9)
    ax2.legend(loc='upper left', fontsize=8)
    ax2.set_ylim(0, 105)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    # Add annotations
    ax2.text(0, 0.5, '1%', ha='center', va='center', fontsize=9,
            fontweight='bold', color='white')
    ax2.text(1, 10, '20%', ha='center', va='center', fontsize=9,
            fontweight='bold', color='white')

    # ========== Panel C: Contamination Sources & Properties ==========
    ax3 = fig.add_subplot(gs[1, :])
    ax3.axis('off')

    ax3.text(0.5, 0.95, 'C. Contamination Sources & Properties',
            ha='center', fontsize=11, fontweight='bold', transform=ax3.transAxes)

    # Table of contamination sources
    contam_data = [
        {
            'source': 'Host DNA',
            'abundance': '5-40%',
            'origin': 'Tissue/cell lysis',
            'properties': 'High molecular weight, GC bias',
            'hvp_impact': 'Body site specific (gut: high, plasma: low)'
        },
        {
            'source': 'Ribosomal RNA',
            'abundance': '80-95% (RNA)',
            'origin': 'Abundant cellular RNA',
            'properties': 'Dominant in RNA preps, depleted 90-95%',
            'hvp_impact': 'Critical for respiratory/enteric RNA viromes'
        },
        {
            'source': 'Bacterial Reagents',
            'abundance': '0.1-5%',
            'origin': 'Kit contamination',
            'properties': 'Low biomass samples most affected',
            'hvp_impact': 'Batch effects across sites/protocols'
        },
        {
            'source': 'PhiX Control',
            'abundance': '5-10%',
            'origin': 'Sequencing spike-in',
            'properties': 'Known sequence, easy to remove',
            'hvp_impact': 'QC metric, cluster density control'
        }
    ]

    # Table headers
    headers = ['Source', 'Abundance', 'Origin', 'Properties', 'HVP Impact']
    x_positions = [0.05, 0.2, 0.35, 0.5, 0.7]

    for x, header in zip(x_positions, headers):
        ax3.text(x, 0.8, header, ha='left', fontsize=9, fontweight='bold',
                transform=ax3.transAxes)

    # Horizontal line under headers
    ax3.plot([0.02, 0.98], [0.77, 0.77], 'k-', linewidth=1,
            transform=ax3.transAxes)

    # Table rows
    colors_row = ['#C73E1D', '#8B4513', '#F18F01', '#A23B72']
    for i, (data, color) in enumerate(zip(contam_data, colors_row)):
        y_pos = 0.7 - i * 0.15

        # Colored box for source
        source_box = Rectangle((0.03, y_pos - 0.03), 0.12, 0.08,
                               transform=ax3.transAxes,
                               facecolor=color, alpha=0.2, edgecolor=color,
                               linewidth=1.5)
        ax3.add_patch(source_box)

        ax3.text(x_positions[0], y_pos, data['source'], ha='left', fontsize=8,
                fontweight='bold', color=color, transform=ax3.transAxes)
        ax3.text(x_positions[1], y_pos, data['abundance'], ha='left', fontsize=7,
                transform=ax3.transAxes)
        ax3.text(x_positions[2], y_pos, data['origin'], ha='left', fontsize=7,
                transform=ax3.transAxes)
        ax3.text(x_positions[3], y_pos, data['properties'], ha='left', fontsize=7,
                transform=ax3.transAxes, wrap=True)
        ax3.text(x_positions[4], y_pos, data['hvp_impact'], ha='left', fontsize=7,
                style='italic', transform=ax3.transAxes)

    # Bottom note
    ax3.text(0.5, 0.02,
            'ViroForge exports complete contamination manifest: every read tracked to source | '
            'Enables systematic QC validation',
            ha='center', fontsize=9, style='italic', transform=ax3.transAxes,
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def main():
    """Generate both innovation figures."""
    print("Creating ViroForge Innovation Figures...")

    # Figure 1: Workflow & Capabilities
    print("\n[1/2] Creating Figure: ViroForge Workflow & Capabilities...")
    fig1 = create_workflow_figure()
    fig1.savefig('Figure_ViroForge_Workflow.png', dpi=300, bbox_inches='tight')
    fig1.savefig('Figure_ViroForge_Workflow.pdf', bbox_inches='tight')
    print("  ✓ Saved: Figure_ViroForge_Workflow.png (300 DPI)")
    print("  ✓ Saved: Figure_ViroForge_Workflow.pdf (vector)")

    # Figure 2: Contamination Details
    print("\n[2/2] Creating Figure: Contamination Modeling Details...")
    fig2 = create_contamination_figure()
    fig2.savefig('Figure_ViroForge_Contamination.png', dpi=300, bbox_inches='tight')
    fig2.savefig('Figure_ViroForge_Contamination.pdf', bbox_inches='tight')
    print("  ✓ Saved: Figure_ViroForge_Contamination.png (300 DPI)")
    print("  ✓ Saved: Figure_ViroForge_Contamination.pdf (vector)")

    print("\n" + "="*70)
    print("SUCCESS: ViroForge Innovation Figures Generated!")
    print("="*70)
    print("\nFigure 1: ViroForge Workflow & Capabilities")
    print("  - Complete workflow: VLP → Contamination → Library → Sequencing")
    print("  - RNA virome branch highlighted")
    print("  - Innovation comparison table")
    print("\nFigure 2: Contamination Modeling Details")
    print("  - DNA virome composition pie chart")
    print("  - RNA virome before/after rRNA depletion")
    print("  - Contamination sources table with HVP impact")
    print("\nReady for grant proposal Innovation section!")
    print("="*70)

    plt.show()


if __name__ == '__main__':
    main()
