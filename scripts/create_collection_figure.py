#!/usr/bin/env python3
"""
Create ViroForge collection overview figure for grant proposal.

Shows the 28 curated collections with viral family composition,
diversity metrics, and coverage across body sites.
"""

import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("colorblind")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 1.2


def load_database():
    """Load ViroForge database."""
    db_path = 'viroforge/data/viral_genomes.db'
    conn = sqlite3.connect(db_path)
    return conn


def get_collection_data(conn):
    """Get all collections with genome counts."""
    query = """
    SELECT
        c.collection_id,
        c.collection_name,
        c.n_genomes,
        c.description
    FROM body_site_collections c
    ORDER BY c.collection_id
    """
    df = pd.read_sql_query(query, conn)
    return df


def get_family_composition(conn):
    """Get viral family composition for each collection."""
    query = """
    SELECT
        cg.collection_id,
        t.family,
        COUNT(*) as n_genomes,
        SUM(cg.relative_abundance) as total_abundance
    FROM collection_genomes cg
    JOIN taxonomy t ON cg.genome_id = t.genome_id
    WHERE t.family != 'Unknown'
    GROUP BY cg.collection_id, t.family
    ORDER BY cg.collection_id, total_abundance DESC
    """
    df = pd.read_sql_query(query, conn)
    return df


def get_top_families(family_df, n=15):
    """Get top N most abundant viral families across all collections."""
    family_totals = family_df.groupby('family')['total_abundance'].sum()
    top_families = family_totals.nlargest(n).index.tolist()
    return top_families


def create_collection_figure():
    """
    Create comprehensive ViroForge collection figure.

    3-panel layout:
    A. Collection overview (bar chart by category)
    B. Viral family heatmap across collections
    C. Key statistics (diversity, coverage, etc.)
    """
    conn = load_database()

    # Load data
    collections = get_collection_data(conn)
    family_comp = get_family_composition(conn)

    # Get top families
    top_families = get_top_families(family_comp, n=15)

    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 2, 1], width_ratios=[2, 1],
                          hspace=0.4, wspace=0.3)

    # ========== Panel A: Collection Categories ==========
    ax1 = fig.add_subplot(gs[0, :])

    # Categorize collections
    categories = {
        'Host-Associated\n(Original)': list(range(1, 9)),
        'VLP Comparison': list(range(9, 16)),
        'Amplification': [16],
        'Disease/Environment': list(range(17, 21)),
        'RNA Viromes': list(range(21, 24)),
        'Additional Host Sites': list(range(24, 29))
    }

    colors_cat = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#06A77D', '#8B4513']

    x_pos = 0
    x_ticks = []
    x_labels = []

    for (cat_name, ids), color in zip(categories.items(), colors_cat):
        cat_collections = collections[collections['collection_id'].isin(ids)]
        widths = cat_collections['n_genomes'].values

        for i, (idx, row) in enumerate(cat_collections.iterrows()):
            ax1.bar(x_pos, row['n_genomes'], color=color, alpha=0.8,
                   edgecolor='black', linewidth=0.5)
            x_ticks.append(x_pos)
            # Shortened labels
            label = row['collection_name'].replace('Human ', '').replace(' Virome', '')
            if len(label) > 15:
                label = label[:12] + '...'
            x_labels.append(f"{row['collection_id']}")
            x_pos += 1

        x_pos += 0.5  # Gap between categories

    ax1.set_ylabel('Number of Genomes', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Collection ID', fontsize=11, fontweight='bold')
    ax1.set_title('A. ViroForge Collection Landscape: 28 Curated Viral Communities',
                 fontsize=12, fontweight='bold', pad=10)
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')

    # Legend for categories
    legend_elements = [mpatches.Patch(facecolor=color, edgecolor='black',
                                     label=cat.replace('\n', ' '))
                      for cat, color in zip(categories.keys(), colors_cat)]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=3)

    # Add summary stats
    total_genomes = collections['n_genomes'].sum()
    ax1.text(0.02, 0.95, f'28 Collections | {len(collections)} unique viral communities',
            transform=ax1.transAxes, fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # ========== Panel B: Viral Family Heatmap ==========
    ax2 = fig.add_subplot(gs[1, :])

    # Create pivot table for heatmap
    heatmap_data = family_comp[family_comp['family'].isin(top_families)].pivot_table(
        index='family', columns='collection_id', values='total_abundance', fill_value=0
    )

    # Reorder by total abundance
    heatmap_data = heatmap_data.loc[top_families]

    # Log transform for better visualization
    heatmap_log = np.log10(heatmap_data + 0.001)

    # Plot heatmap
    sns.heatmap(heatmap_log, ax=ax2, cmap='YlOrRd',
               cbar_kws={'label': 'log₁₀(Relative Abundance)'},
               linewidths=0.5, linecolor='lightgray')

    # Create collection ID to name mapping
    id_to_name = dict(zip(collections['collection_id'], collections['collection_name']))

    # Replace collection IDs with names on x-axis
    collection_names = []
    for cid in heatmap_data.columns:
        if cid in id_to_name:
            name = id_to_name[cid]
            # Shorten names for readability
            name = name.replace('Human ', '').replace(' Virome', '')
            name = name.replace('Healthy ', '').replace('Collection', 'Coll.')
            if len(name) > 20:
                name = name[:17] + '...'
            collection_names.append(name)
        else:
            collection_names.append(f'Coll. {cid}')

    ax2.set_xticklabels(collection_names, rotation=90, ha='center', fontsize=7)
    ax2.set_xlabel('Collection', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Viral Family', fontsize=11, fontweight='bold')
    ax2.set_title('B. Viral Family Composition Across Collections (Top 15 Families)',
                 fontsize=12, fontweight='bold', pad=10)
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0, fontsize=9)

    # ========== Panel C: Key Statistics ==========
    ax3 = fig.add_subplot(gs[2, 0])

    # Calculate diversity metrics
    stats_data = []
    for cid in collections['collection_id']:
        coll_families = family_comp[family_comp['collection_id'] == cid]
        n_families = len(coll_families)
        n_genomes = collections[collections['collection_id'] == cid]['n_genomes'].values[0]

        # Shannon diversity
        if len(coll_families) > 0:
            abundances = coll_families['total_abundance'].values
            abundances = abundances / abundances.sum()
            shannon = -np.sum(abundances * np.log(abundances + 1e-10))
        else:
            shannon = 0

        stats_data.append({
            'collection_id': cid,
            'n_families': n_families,
            'n_genomes': n_genomes,
            'shannon': shannon
        })

    stats_df = pd.DataFrame(stats_data)

    # Scatter plot: genomes vs families
    scatter = ax3.scatter(stats_df['n_genomes'], stats_df['n_families'],
                         c=stats_df['shannon'], s=100, cmap='viridis',
                         alpha=0.7, edgecolors='black', linewidth=1)

    ax3.set_xlabel('Number of Genomes', fontsize=10, fontweight='bold')
    ax3.set_ylabel('Number of Viral Families', fontsize=10, fontweight='bold')
    ax3.set_title('C. Collection Complexity', fontsize=11, fontweight='bold')
    ax3.grid(alpha=0.3, linestyle='--')

    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Shannon Diversity', fontsize=9)

    # ========== Panel D: Body Site Coverage ==========
    ax4 = fig.add_subplot(gs[2, 1])

    # Define body sites and their collection counts
    body_sites = {
        'Gut': [1, 9, 10, 11, 12, 13, 14, 15, 18, 19, 23],
        'Respiratory': [5, 20, 21],
        'Skin': [2, 25],
        'Oral': [3, 26],
        'Urogenital': [4, 24],
        'Blood': [27],
        'Ocular': [28],
        'Environmental': [6, 7, 8, 17, 22]
    }

    site_counts = {site: len(ids) for site, ids in body_sites.items()}

    # Pie chart
    colors_pie = sns.color_palette("Set2", len(site_counts))
    wedges, texts, autotexts = ax4.pie(site_counts.values(), labels=site_counts.keys(),
                                       autopct='%d', colors=colors_pie,
                                       startangle=90, textprops={'fontsize': 8})

    ax4.set_title('D. Body Site Coverage', fontsize=11, fontweight='bold')

    # Add summary text at bottom
    fig.text(0.5, 0.02,
            'ViroForge v0.11.0: Production-ready synthetic virome generator | '
            '14,423 RefSeq genomes | 5 sequencing platforms | '
            'DNA + RNA workflows | Complete ground truth',
            ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout(rect=[0, 0.03, 1, 1])

    conn.close()
    return fig


def main():
    """Generate ViroForge collection figure."""
    print("Creating ViroForge Collection Overview Figure...")
    print("\nQuerying database for collection data...")

    fig = create_collection_figure()

    print("\nSaving figure...")
    fig.savefig('Figure_ViroForge_Collections.png', dpi=300, bbox_inches='tight')
    fig.savefig('Figure_ViroForge_Collections.pdf', bbox_inches='tight')

    print("\n" + "="*60)
    print("SUCCESS: ViroForge Collection Figure Generated!")
    print("="*60)
    print("\n✓ Saved: Figure_ViroForge_Collections.png (300 DPI)")
    print("✓ Saved: Figure_ViroForge_Collections.pdf (vector)")
    print("\nFigure panels:")
    print("  A. Collection landscape (28 curated communities)")
    print("  B. Viral family heatmap (top 15 families)")
    print("  C. Collection complexity (genomes vs families)")
    print("  D. Body site coverage (pie chart)")
    print("\nReady for grant proposal!")
    print("="*60)

    plt.show()


if __name__ == '__main__':
    main()
