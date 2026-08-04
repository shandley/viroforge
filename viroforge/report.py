"""
ViroForge Generation Report

Produces an HTML report with diagnostic plots after dataset generation.
The report helps users evaluate whether their synthetic dataset looks
realistic compared to real virome data.

Usage (called automatically by run_generation):
    from viroforge.report import generate_report
    generate_report(metadata_path, output_dir)
"""

import base64
import json
import logging
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _fig_to_base64(fig) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    return encoded


def _composition_pie(metadata: Dict) -> str:
    """Create a pie chart of dataset composition (viral vs contamination types).

    Uses a legend instead of inline labels to avoid overlap.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    enrichment = metadata.get("enrichment_stats", {})
    reduction = enrichment.get("contamination_reduction", {})
    reduction_by_type = reduction.get("reduction_by_type", {})

    # Color scheme: Okabe-Ito palette (colorblind-safe)
    color_map = {
        "Viral": "#0072B2",
        "Host DNA": "#D55E00",
        "rRNA": "#CC79A7",
        "Bacterial": "#009E73",
        "Fungal": "#F0E442",
        "Archaeal": "#56B4E9",
        "Reagent": "#E69F00",
        "PhiX": "#999999",
    }

    type_name_map = {
        "host_dna": "Host DNA",
        "rrna": "rRNA",
        "bacterial_background": "Bacterial",
        "fungal_background": "Fungal",
        "archaeal_background": "Archaeal",
        "reagent_bacteria": "Reagent",
        "phix": "PhiX",
    }

    labels = []
    sizes = []
    colors = []

    viral_frac = enrichment.get("viral_fraction", 0)
    if viral_frac > 0:
        labels.append("Viral")
        sizes.append(viral_frac)
        colors.append(color_map["Viral"])

    for key, display_name in type_name_map.items():
        if key in reduction_by_type:
            val = reduction_by_type[key].get("reduced_abundance", 0)
            if val > 0.0001:
                labels.append(display_name)
                sizes.append(val)
                colors.append(color_map.get(display_name, "#AAAAAA"))

    if not sizes:
        return ""

    fig, ax = plt.subplots(figsize=(7, 5))
    wedges, _ = ax.pie(
        sizes,
        colors=colors,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1},
    )
    # Use a legend instead of inline labels to avoid overlap
    legend_labels = [f"{name} ({val*100:.1f}%)" for name, val in zip(labels, sizes)]
    ax.legend(wedges, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=10, frameon=False)
    ax.set_title("Dataset Composition (Post-VLP Enrichment)", fontsize=13)

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _before_after_vlp_bar(metadata: Dict) -> str:
    """Create side-by-side bar chart showing composition before and after VLP."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    enrichment = metadata.get("enrichment_stats", {})
    reduction = enrichment.get("contamination_reduction", {})
    reduction_by_type = reduction.get("reduction_by_type", {})

    if not reduction_by_type:
        return ""

    color_map = {
        "Viral": "#0072B2",
        "Host DNA": "#D55E00",
        "rRNA": "#CC79A7",
        "Bacterial": "#009E73",
        "Fungal": "#F0E442",
        "Archaeal": "#56B4E9",
        "Reagent": "#E69F00",
        "PhiX": "#999999",
    }

    type_name_map = {
        "host_dna": "Host DNA",
        "rrna": "rRNA",
        "bacterial_background": "Bacterial",
        "fungal_background": "Fungal",
        "archaeal_background": "Archaeal",
        "reagent_bacteria": "Reagent",
        "phix": "PhiX",
    }

    # Compute before-VLP viral fraction
    # Before VLP: viral fraction = 1 - total contamination
    total_contam_before = reduction.get("original_total_contamination", 0)
    viral_before = 1.0 - total_contam_before
    viral_after = enrichment.get("viral_fraction", 0)

    categories = ["Viral"]
    before_vals = [viral_before * 100]
    after_vals = [viral_after * 100]
    bar_colors = [color_map["Viral"]]

    for key, display_name in type_name_map.items():
        if key in reduction_by_type:
            data = reduction_by_type[key]
            orig = data.get("original_abundance", 0)
            reduced = data.get("reduced_abundance", 0)
            if orig > 0.0001 or reduced > 0.0001:
                categories.append(display_name)
                before_vals.append(orig * 100)
                after_vals.append(reduced * 100)
                bar_colors.append(color_map.get(display_name, "#AAAAAA"))

    if len(categories) < 2:
        return ""

    x = np.arange(len(categories))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))
    bars1 = ax.bar(x - width / 2, before_vals, width, label="Before VLP",
                   color=bar_colors, alpha=0.5, edgecolor="gray", linewidth=0.5)
    bars2 = ax.bar(x + width / 2, after_vals, width, label="After VLP",
                   color=bar_colors, alpha=1.0, edgecolor="gray", linewidth=0.5)

    ax.set_ylabel("Percentage (%)")
    ax.set_title("Composition Before vs After VLP Enrichment", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=30, ha="right", fontsize=9)
    ax.legend()
    ax.set_ylim(bottom=0)

    # Add value labels on bars
    for bar in bars1:
        h = bar.get_height()
        if h > 1:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                    f"{h:.1f}%", ha="center", va="bottom", fontsize=7)
    for bar in bars2:
        h = bar.get_height()
        if h > 1:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                    f"{h:.1f}%", ha="center", va="bottom", fontsize=7)

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _viral_families_bar(metadata: Dict) -> str:
    """Create a horizontal bar chart of top viral families by abundance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sequences = metadata.get("sequences", [])
    if not sequences:
        return ""

    # Aggregate abundance by family
    family_abundance = {}
    for seq in sequences:
        if seq.get("sequence_type") == "viral":
            family = seq.get("family", "Unknown")
            family_abundance[family] = family_abundance.get(family, 0) + seq.get(
                "relative_abundance", 0
            )

    if not family_abundance:
        return ""

    # Sort and take top 15
    sorted_families = sorted(family_abundance.items(), key=lambda x: x[1], reverse=True)
    top = sorted_families[:15]
    if len(sorted_families) > 15:
        other = sum(v for _, v in sorted_families[15:])
        top.append(("Other", other))

    names = [f for f, _ in reversed(top)]
    values = [v for _, v in reversed(top)]

    fig, ax = plt.subplots(figsize=(8, max(4, len(names) * 0.35)))
    bars = ax.barh(names, values, color="#0072B2")
    ax.set_xlabel("Relative Abundance")
    ax.set_title("Top Viral Families", fontsize=13)
    ax.tick_params(axis="y", labelsize=9)

    # Add value labels
    for bar, val in zip(bars, values):
        if val > 0.001:
            ax.text(
                bar.get_width() + max(values) * 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}",
                va="center",
                fontsize=8,
            )

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _rank_abundance_curve(metadata: Dict) -> str:
    """Create a rank-abundance curve (log scale) for all genomes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    sequences = metadata.get("sequences", [])
    viral_seqs = [s for s in sequences if s.get("sequence_type") == "viral"]
    if not viral_seqs:
        return ""

    abundances = sorted(
        [s.get("relative_abundance", 0) for s in viral_seqs], reverse=True
    )
    abundances = [a for a in abundances if a > 0]

    if not abundances:
        return ""

    ranks = list(range(1, len(abundances) + 1))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(ranks, abundances, "o-", markersize=3, color="#0072B2", linewidth=1)
    ax.set_xlabel("Genome Rank")
    ax.set_ylabel("Relative Abundance (log scale)")
    ax.set_title("Rank-Abundance Curve (Viral Genomes)", fontsize=13)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _gc_content_histogram(metadata: Dict) -> str:
    """Create a GC content distribution histogram."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sequences = metadata.get("sequences", [])
    if not sequences:
        return ""

    # Try to get GC from benchmarking data or compute from sequence info
    gc_values = []
    for seq in sequences:
        gc = seq.get("gc_content")
        if gc is not None:
            gc_values.append(gc)

    if not gc_values:
        return ""

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(gc_values, bins=30, color="#0072B2", edgecolor="white", alpha=0.8)
    ax.set_xlabel("GC Content")
    ax.set_ylabel("Number of Genomes")
    ax.set_title("GC Content Distribution", fontsize=13)
    ax.axvline(x=0.5, color="#D55E00", linestyle="--", alpha=0.5, label="50% GC")
    ax.legend()

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _genome_length_vs_abundance(metadata: Dict) -> str:
    """Create a scatter plot of genome length vs abundance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    sequences = metadata.get("sequences", [])
    viral_seqs = [s for s in sequences if s.get("sequence_type") == "viral"]
    if not viral_seqs:
        return ""

    lengths = [s.get("length", 0) for s in viral_seqs]
    abundances = [s.get("relative_abundance", 0) for s in viral_seqs]

    if not lengths or not abundances:
        return ""

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(
        [l / 1000 for l in lengths],
        abundances,
        alpha=0.6,
        s=20,
        color="#0072B2",
        edgecolors="white",
        linewidth=0.3,
    )
    ax.set_xlabel("Genome Length (kb)")
    ax.set_ylabel("Relative Abundance")
    ax.set_title("Genome Length vs Abundance", fontsize=13)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    encoded = _fig_to_base64(fig)
    plt.close(fig)
    return encoded


def _vlp_reduction_table(metadata: Dict) -> str:
    """Generate HTML table for VLP contamination reduction."""
    enrichment = metadata.get("enrichment_stats", {})
    reduction = enrichment.get("contamination_reduction", {})
    reduction_by_type = reduction.get("reduction_by_type", {})

    if not reduction_by_type:
        return ""

    display_names = {
        "host_dna": "Host DNA",
        "rrna": "rRNA",
        "bacterial_background": "Bacterial",
        "fungal_background": "Fungal",
        "archaeal_background": "Archaeal",
        "reagent_bacteria": "Reagent Bacteria",
        "phix": "PhiX",
    }

    rows = ""
    for key, data in reduction_by_type.items():
        name = display_names.get(key, key)
        orig = data.get("original_abundance", 0) * 100
        reduced = data.get("reduced_abundance", 0) * 100
        removal = data.get("reduction_pct", 0)
        rows += f"""
        <tr>
            <td>{name}</td>
            <td>{orig:.2f}%</td>
            <td>{reduced:.2f}%</td>
            <td>{removal:.1f}%</td>
        </tr>"""

    return f"""
    <table>
        <tr>
            <th>Contaminant</th>
            <th>Before VLP</th>
            <th>After VLP</th>
            <th>Removal</th>
        </tr>
        {rows}
    </table>"""


def _artifact_summary_table(metadata: Dict) -> str:
    """Generate HTML table for artifact injection summary."""
    rows = ""

    # Adapter stats
    adapter = metadata.get("adapter_stats", {})
    if adapter:
        mode = adapter.get("mode", "N/A")
        if mode == "fixed_rate":
            rate = adapter.get("adapter_rate", 0)
            modified = adapter.get("reads_modified", 0)
            total = adapter.get("reads_total", 0)
            rows += f"""
            <tr>
                <td>Adapter Read-through</td>
                <td>{rate:.1%}</td>
                <td>{modified:,} / {total:,} reads</td>
            </tr>"""
        elif mode == "insert_size_driven":
            emergent = adapter.get("emergent_adapter_rate", 0)
            readthrough = adapter.get("reads_readthrough", 0)
            total = adapter.get("reads_total", 0)
            rows += f"""
            <tr>
                <td>Adapter Read-through (insert-size)</td>
                <td>{emergent:.1%}</td>
                <td>{readthrough:,} / {total:,} reads</td>
            </tr>"""

    # Duplicate stats
    dup = metadata.get("duplicate_stats", {})
    if dup:
        rate = dup.get("duplicate_rate", 0)
        copies = dup.get("copies_generated", 0)
        templates = dup.get("templates", 0)
        method = dup.get("amplification_method", "pcr")
        rows += f"""
        <tr>
            <td>PCR Duplicates ({method.upper()})</td>
            <td>{rate:.1%}</td>
            <td>{copies:,} copies from {templates:,} templates</td>
        </tr>"""

    # Low complexity stats
    lc = metadata.get("low_complexity_stats", {})
    if lc:
        rate = lc.get("rate", 0)
        modified = lc.get("reads_modified", 0)
        rows += f"""
        <tr>
            <td>Low-Complexity Artifacts</td>
            <td>{rate:.1%}</td>
            <td>{modified:,} reads modified</td>
        </tr>"""

    # Chimera stats
    chimera = metadata.get("chimera_stats", {})
    if chimera:
        rate = chimera.get("chimera_rate", 0)
        created = chimera.get("chimeras_created", 0)
        rows += f"""
        <tr>
            <td>MDA Chimeras</td>
            <td>{rate:.1%}</td>
            <td>{created:,} chimeric reads</td>
        </tr>"""

    if not rows:
        return ""

    return f"""
    <table>
        <tr>
            <th>Artifact</th>
            <th>Rate</th>
            <th>Details</th>
        </tr>
        {rows}
    </table>"""


def generate_report(metadata_path: Path, output_dir: Path) -> Path:
    """Generate an HTML report with diagnostic plots for a generated dataset.

    Args:
        metadata_path: Path to the metadata JSON file.
        output_dir: Directory to write the report HTML file.

    Returns:
        Path to the generated HTML report.
    """
    with open(metadata_path) as f:
        metadata = json.load(f)

    collection = metadata.get("collection", {})
    config = metadata.get("configuration", {})
    enrichment = metadata.get("enrichment_stats", {})
    amplification = metadata.get("amplification_stats", {})

    # Generate plots
    logger.info("Generating report plots...")
    pie_b64 = _composition_pie(metadata)
    before_after_b64 = _before_after_vlp_bar(metadata)
    families_b64 = _viral_families_bar(metadata)
    rank_b64 = _rank_abundance_curve(metadata)
    gc_b64 = _gc_content_histogram(metadata)
    scatter_b64 = _genome_length_vs_abundance(metadata)

    # Generate tables
    vlp_table = _vlp_reduction_table(metadata)
    artifact_table = _artifact_summary_table(metadata)

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ViroForge Generation Report — {collection.get('name', 'Dataset')}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
            background: #fafafa;
            color: #333;
        }}
        h1 {{
            color: #0072B2;
            border-bottom: 2px solid #0072B2;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #444;
            margin-top: 30px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 5px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 15px;
            margin: 15px 0;
        }}
        .summary-card {{
            background: white;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 15px;
        }}
        .summary-card h3 {{
            margin: 0 0 10px 0;
            color: #555;
            font-size: 14px;
        }}
        .summary-card .value {{
            font-size: 22px;
            font-weight: bold;
            color: #0072B2;
        }}
        .summary-card .detail {{
            font-size: 12px;
            color: #888;
            margin-top: 5px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 10px 0;
            background: white;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px 12px;
            text-align: left;
        }}
        th {{
            background: #f0f0f0;
            font-weight: 600;
        }}
        tr:nth-child(even) {{
            background: #f9f9f9;
        }}
        .plot {{
            text-align: center;
            margin: 15px 0;
        }}
        .plot img {{
            max-width: 100%;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 15px;
            border-top: 1px solid #ddd;
            color: #888;
            font-size: 12px;
        }}
    </style>
</head>
<body>

<h1>ViroForge Generation Report</h1>

<div class="summary-grid">
    <div class="summary-card">
        <h3>Collection</h3>
        <div class="value">{collection.get('name', 'N/A')}</div>
        <div class="detail">ID: {collection.get('id', 'N/A')}</div>
    </div>
    <div class="summary-card">
        <h3>Platform</h3>
        <div class="value">{config.get('platform', 'N/A').upper()}</div>
        <div class="detail">Coverage: {config.get('coverage', 'N/A')}x</div>
    </div>
    <div class="summary-card">
        <h3>Viral Genomes</h3>
        <div class="value">{collection.get('n_viral_genomes', 0)}</div>
        <div class="detail">+ {collection.get('n_contaminants', 0)} contaminants</div>
    </div>
    <div class="summary-card">
        <h3>Viral Fraction</h3>
        <div class="value">{enrichment.get('viral_fraction', 0) * 100:.1f}%</div>
        <div class="detail">VLP: {config.get('vlp_protocol', 'none')} | {config.get('contamination_level', 'N/A')}</div>
    </div>
</div>

<h2>1. Dataset Composition</h2>
<p>Breakdown of read types in the generated dataset after VLP enrichment and contamination injection.</p>
{"<div class='plot'><img src='data:image/png;base64," + pie_b64 + "' alt='Composition Pie Chart'></div>" if pie_b64 else "<p><em>No composition data available.</em></p>"}

<h2>2. Before vs After VLP Enrichment</h2>
<p>Comparison of component fractions before and after VLP enrichment. VLP filtration removes most bacterial and host DNA contamination while retaining viral particles.</p>
{"<div class='plot'><img src='data:image/png;base64," + before_after_b64 + "' alt='Before vs After VLP'></div>" if before_after_b64 else "<p><em>No VLP enrichment data available.</em></p>"}

<h2>3. Viral Community Structure</h2>
<p>Most abundant viral families in the dataset. The community composition reflects the curated collection for this body site.</p>
{"<div class='plot'><img src='data:image/png;base64," + families_b64 + "' alt='Viral Families Bar Chart'></div>" if families_b64 else "<p><em>No viral family data available.</em></p>"}

<h2>4. Rank-Abundance Curve</h2>
<p>Log-scale rank-abundance plot showing community evenness. A steep curve indicates a few dominant genomes; a flat curve indicates even abundance.</p>
{"<div class='plot'><img src='data:image/png;base64," + rank_b64 + "' alt='Rank-Abundance Curve'></div>" if rank_b64 else "<p><em>No abundance data available.</em></p>"}

<h2>5. VLP Contamination Reduction</h2>
<p>How effectively VLP enrichment removed each type of contamination.</p>
{vlp_table if vlp_table else "<p><em>No VLP enrichment applied.</em></p>"}

<h2>6. GC Content Distribution</h2>
<p>Distribution of GC content across all genomes in the dataset. Amplification bias can skew this distribution.</p>
{"<div class='plot'><img src='data:image/png;base64," + gc_b64 + "' alt='GC Content Histogram'></div>" if gc_b64 else "<p><em>No GC content data available.</em></p>"}

<h2>7. Genome Length vs Abundance</h2>
<p>Scatter plot showing whether genome size correlates with abundance. VLP size-based filtration can introduce length-dependent bias.</p>
{"<div class='plot'><img src='data:image/png;base64," + scatter_b64 + "' alt='Genome Length vs Abundance'></div>" if scatter_b64 else "<p><em>No genome length data available.</em></p>"}

<h2>8. Sequencing Artifacts</h2>
<p>Summary of injected sequencing artifacts for QC tool validation.</p>
{artifact_table if artifact_table else "<p><em>No artifacts injected.</em></p>"}

<div class="footer">
    Generated by ViroForge v{metadata.get('generation_info', {}).get('viroforge_version', 'unknown')}
    | Seed: {metadata.get('generation_info', {}).get('random_seed', 'N/A')}
    | {metadata.get('generation_info', {}).get('timestamp', '')}
</div>

</body>
</html>"""

    # Write report
    report_path = output_dir / "generation_report.html"
    with open(report_path, "w") as f:
        f.write(html)

    logger.info(f"Report saved to {report_path}")
    return report_path
