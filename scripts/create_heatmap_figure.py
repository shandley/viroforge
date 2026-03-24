#!/usr/bin/env python3
"""
Create standalone viral family composition heatmap for grant proposal.

Shows viral family composition across all 28 ViroForge collections.
"""

import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("colorblind")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
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
        AND cg.collection_id >= 9
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


def shorten_collection_name(name):
    """
    Create clear, specific collection names for heatmap.

    Strategy:
    - Keep body site/disease state clear
    - Add context (RNA, disease, etc.)
    - Keep concise and readable
    """
    # Accurate mapping for collections 9-28 (20 total)
    name_map = {
        'Gut Virome - Adult Healthy (Western Diet)': 'Healthy Gut',
        'Oral Virome - Saliva (Healthy)': 'Healthy Oral',
        'Skin Virome - Sebaceous Sites (Healthy)': 'Healthy Skin',
        'Respiratory Virome - Nasopharynx (Healthy)': 'Healthy Respiratory',
        'Marine Virome - Coastal Surface Water': 'Marine',
        'Soil Virome - Agricultural': 'Soil',
        'Freshwater Virome - Lake Surface Water': 'Freshwater',
        'Mouse Gut Virome - Laboratory (C57BL/6)': 'Mouse Gut',
        'Wastewater Virome - Urban Treatment Plant': 'Wastewater',
        'IBD Gut Virome (Inflammatory Bowel Disease)': 'IBD Gut',
        'HIV+ Gut Virome': 'HIV+ Gut',
        'Cystic Fibrosis (CF) Respiratory Virome': 'CF Respiratory',
        'Human Respiratory RNA Virome': 'Respiratory RNA',
        'Arbovirus Environmental (Mosquito Virome)': 'Arbovirus (Env.)',
        'Fecal RNA Virome': 'Fecal RNA',
        'Vaginal Virome (Healthy)': 'Vaginal',
        'Blood/Plasma Virome (Healthy)': 'Blood/Plasma',
        'Ocular Surface Virome (Healthy)': 'Ocular Surface',
        'Lower Respiratory (Lung) Virome (Healthy)': 'Lung',
        'Urinary Virome (Healthy)': 'Urinary'
    }

    if name in name_map:
        return name_map[name]
    else:
        # Fallback: remove common words
        short = name.replace('Human ', '').replace(' Virome', '')
        short = short.replace('Healthy ', '').replace('Collection', '')
        return short[:25]


def create_heatmap_figure():
    """Create standalone viral family composition heatmap."""
    conn = load_database()

    # Load data
    collections = get_collection_data(conn)
    family_comp = get_family_composition(conn)

    # Get top families
    top_families = get_top_families(family_comp, n=15)

    # Create pivot table for heatmap
    heatmap_data = family_comp[family_comp['family'].isin(top_families)].pivot_table(
        index='family', columns='collection_id', values='total_abundance', fill_value=0
    )

    # Reorder by total abundance
    heatmap_data = heatmap_data.loc[top_families]

    # Log transform for better visualization
    heatmap_log = np.log10(heatmap_data + 0.001)

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Plot heatmap
    sns.heatmap(heatmap_log, ax=ax, cmap='YlOrRd',
               cbar_kws={'label': 'log₁₀(Relative Abundance)', 'shrink': 0.8},
               linewidths=0.5, linecolor='lightgray',
               xticklabels=True, yticklabels=True)

    # Create collection ID to name mapping
    id_to_name = dict(zip(collections['collection_id'], collections['collection_name']))

    # Replace collection IDs with clear names on x-axis
    collection_names = []
    for cid in heatmap_data.columns:
        if cid in id_to_name:
            name = shorten_collection_name(id_to_name[cid])
            collection_names.append(name)
        else:
            collection_names.append(f'Collection {cid}')

    ax.set_xticklabels(collection_names, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    ax.set_xlabel('ViroForge Collection', fontsize=12, fontweight='bold')
    ax.set_ylabel('Viral Family', fontsize=12, fontweight='bold')
    ax.set_title('ViroForge Viral Family Composition Across 20 Curated Collections',
                 fontsize=14, fontweight='bold', pad=15)

    # Add summary stats at bottom
    fig.text(0.5, 0.02,
            f'Top 15 viral families | 20 curated collections | '
            f'14,423 RefSeq genomes | ViroForge v0.11.0',
            ha='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout(rect=[0, 0.04, 1, 1])

    conn.close()
    return fig


def main():
    """Generate standalone heatmap figure."""
    print("Creating ViroForge Viral Family Heatmap...")
    print("\nQuerying database for family composition data...")

    fig = create_heatmap_figure()

    print("\nSaving figure...")
    fig.savefig('Figure_ViroForge_Heatmap.png', dpi=300, bbox_inches='tight')
    fig.savefig('Figure_ViroForge_Heatmap.pdf', bbox_inches='tight')

    print("\n" + "="*60)
    print("SUCCESS: ViroForge Heatmap Figure Generated!")
    print("="*60)
    print("\n✓ Saved: Figure_ViroForge_Heatmap.png (300 DPI)")
    print("✓ Saved: Figure_ViroForge_Heatmap.pdf (vector)")
    print("\nStandalone heatmap showing:")
    print("  - Top 15 viral families (y-axis)")
    print("  - 20 curated collections (x-axis, IDs 9-28)")
    print("  - Clear, specific collection names")
    print("  - Log-transformed relative abundances")
    print("  - Collections 1-8 excluded (deprecated data)")
    print("\nReady for grant proposal!")
    print("="*60)

    plt.show()


if __name__ == '__main__':
    main()
