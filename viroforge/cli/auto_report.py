"""
Auto-generated HTML diagnostic report for ViroForge datasets.

Produces a self-contained HTML report immediately after dataset generation,
summarizing the generation parameters, community composition, enrichment,
contamination, artifact injection, and output files.
"""

import json
import html as html_lib
from pathlib import Path
from typing import Dict, List, Optional


def generate_auto_report(dataset_path: Path) -> Optional[Path]:
    """Generate an HTML diagnostic report for a dataset directory.

    Looks for the metadata JSON, composition TSV, and output files under
    *dataset_path*, builds a self-contained HTML page, and writes it to
    ``<dataset_path>/<name>_report.html``.

    Returns the Path to the report file, or None on failure.
    """
    dataset_path = Path(dataset_path)
    metadata = _load_metadata(dataset_path)
    if metadata is None:
        return None

    composition = _load_composition(dataset_path)
    file_inventory = _inventory_files(dataset_path)
    html = _render_html(metadata, composition, file_inventory, dataset_path)

    report_path = dataset_path / f"{dataset_path.name}_report.html"
    report_path.write_text(html, encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _load_metadata(dataset_path: Path) -> Optional[Dict]:
    metadata_dir = dataset_path / "metadata"
    candidates = []
    if metadata_dir.is_dir():
        candidates.extend(sorted(metadata_dir.glob("*_metadata.json")))
    candidates.append(dataset_path / "metadata.json")
    for p in candidates:
        if p.is_file():
            try:
                return json.loads(p.read_text(encoding="utf-8"))
            except Exception:
                continue
    return None


def _load_composition(dataset_path: Path) -> List[Dict]:
    metadata_dir = dataset_path / "metadata"
    if not metadata_dir.is_dir():
        return []
    for tsv in sorted(metadata_dir.glob("*_composition.tsv")):
        try:
            rows = []
            lines = tsv.read_text(encoding="utf-8").strip().splitlines()
            if not lines:
                continue
            headers = lines[0].split("\t")
            for line in lines[1:]:
                vals = line.split("\t")
                row = dict(zip(headers, vals))
                for num_key in ("relative_abundance", "gc_content", "length"):
                    if num_key in row:
                        try:
                            row[num_key] = float(row[num_key])
                        except (ValueError, TypeError):
                            pass
                rows.append(row)
            return rows
        except Exception:
            continue
    return []


def _inventory_files(dataset_path: Path) -> List[Dict]:
    inventory = []
    for subdir in ("fastq", "fasta", "metadata"):
        d = dataset_path / subdir
        if not d.is_dir():
            continue
        for f in sorted(d.iterdir()):
            if not f.is_file():
                continue
            size = f.stat().st_size
            inventory.append({
                "name": f.name,
                "category": subdir.upper(),
                "size": size,
                "size_str": _fmt_size(size),
            })
    return inventory


def _fmt_size(nbytes: int) -> str:
    if nbytes >= 1_073_741_824:
        return f"{nbytes / 1_073_741_824:.1f} GB"
    if nbytes >= 1_048_576:
        return f"{nbytes / 1_048_576:.1f} MB"
    if nbytes >= 1024:
        return f"{nbytes / 1024:.1f} KB"
    return f"{nbytes} B"


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

_ESC = html_lib.escape


def _render_html(
    metadata: Dict,
    composition: List[Dict],
    file_inventory: List[Dict],
    dataset_path: Path,
) -> str:
    collection = metadata.get("collection", {})
    gen_info = metadata.get("generation_info", {})
    config = metadata.get("configuration", {})
    enrichment = metadata.get("enrichment_stats", {})
    amplification = metadata.get("amplification_stats", {})

    dataset_name = _ESC(dataset_path.name)
    coll_name_raw = str(collection.get("name", "Unknown"))
    coll_name = _ESC(coll_name_raw)
    coll_id = collection.get("id", "N/A")
    timestamp_raw = str(gen_info.get("timestamp", "N/A"))
    timestamp = _ESC(timestamp_raw)
    seed = gen_info.get("random_seed", "N/A")
    version_raw = str(gen_info.get("viroforge_version", "N/A"))
    version = _ESC(version_raw)

    platform_raw = str(config.get("platform", "unknown")).upper()
    platform_name = _ESC(platform_raw)
    coverage = config.get("coverage") or config.get("target_coverage")
    depth = config.get("depth") or config.get("target_depth")
    cov_str = f"{coverage}x" if coverage else (f"{depth}x" if depth else "N/A")
    read_length = config.get("read_length", "N/A")
    insert_size = config.get("insert_size", "N/A")
    molecule_raw = str(config.get("molecule_type", "dna")).upper()
    amplification_method = str(config.get("amplification", "N/A"))
    contamination_level = str(config.get("contamination_level", "N/A"))
    vlp_protocol = config.get("vlp_protocol", "none")

    n_viral = collection.get("n_viral_genomes", 0)
    n_contam = collection.get("n_contaminants", 0)
    n_total = collection.get("total_sequences", n_viral + n_contam)
    viral_frac = enrichment.get("viral_fraction")
    contam_frac = enrichment.get("contamination_fraction")

    # _kv_rows escapes its values, so pass raw strings here
    summary_rows = _kv_rows([
        ("Collection", f"{coll_name_raw} (ID: {coll_id})"),
        ("Platform", platform_raw),
        ("Molecule Type", molecule_raw),
        ("Coverage / Depth", cov_str),
        ("Read Length", f"{read_length} bp" if read_length != "N/A" else "N/A"),
        ("Insert Size", f"{insert_size} bp" if insert_size != "N/A" else "N/A"),
        ("Amplification", amplification_method),
        ("Contamination Level", contamination_level),
        ("VLP Protocol", vlp_protocol),
        ("Random Seed", str(seed)),
        ("ViroForge Version", version_raw),
        ("Generated", timestamp_raw),
    ])

    composition_section = _render_composition(
        composition, n_viral, n_contam, n_total, viral_frac, contam_frac
    )

    rank_abundance_chart = _render_rank_abundance(composition)
    size_abundance_chart = _render_size_vs_abundance(composition)
    coverage_chart = _render_coverage_distribution(metadata)
    vlp_section = _render_vlp(enrichment, vlp_protocol)
    amp_section = _render_amplification(amplification)
    artifact_section = _render_artifacts(metadata)
    files_section = _render_files(file_inventory)
    dark_matter_section = _render_dark_matter(config)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ViroForge Report - {dataset_name}</title>
<style>
:root {{
  --bg: #f5f6fa; --card: #fff; --text: #2d3436; --muted: #636e72;
  --accent: #6c5ce7; --accent2: #00b894; --border: #dfe6e9;
  --green: #00b894; --blue: #0984e3; --orange: #e17055; --red: #d63031;
}}
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
  Oxygen, Ubuntu, sans-serif; background: var(--bg); color: var(--text);
  line-height: 1.6; }}
.header {{ background: linear-gradient(135deg, #6c5ce7 0%, #a29bfe 100%);
  color: #fff; padding: 2rem 1.5rem; }}
.header h1 {{ font-size: 1.6rem; font-weight: 700; }}
.header p {{ opacity: .85; margin-top: .25rem; }}
.container {{ max-width: 960px; margin: 0 auto; padding: 1.5rem; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px,1fr));
  gap: 1rem; margin-bottom: 1.5rem; }}
.card-stat {{ background: var(--card); border-radius: 8px; padding: 1rem 1.25rem;
  text-align: center; border: 1px solid var(--border); }}
.card-stat .num {{ font-size: 1.8rem; font-weight: 700; }}
.card-stat .label {{ font-size: .85rem; color: var(--muted); }}
.section {{ background: var(--card); border-radius: 8px; border: 1px solid var(--border);
  margin-bottom: 1.5rem; overflow: hidden; }}
.section-header {{ padding: .75rem 1.25rem; font-weight: 600; font-size: 1rem;
  border-bottom: 1px solid var(--border); background: #fafafa; }}
.section-body {{ padding: 1rem 1.25rem; }}
table {{ width: 100%; border-collapse: collapse; }}
th, td {{ text-align: left; padding: .4rem .6rem; font-size: .9rem; }}
th {{ color: var(--muted); font-weight: 600; border-bottom: 2px solid var(--border); }}
td {{ border-bottom: 1px solid #f0f0f0; }}
.kv th {{ width: 40%; font-weight: 500; color: var(--text); }}
.bar {{ display: inline-block; height: 14px; border-radius: 3px; }}
code {{ background: #f0f0f0; padding: .15rem .35rem; border-radius: 3px;
  font-size: .85rem; }}
.footer {{ text-align: center; color: var(--muted); font-size: .8rem;
  padding: 1.5rem 0; }}
.badge {{ display: inline-block; padding: .15rem .5rem; border-radius: 4px;
  font-size: .8rem; font-weight: 600; }}
.badge-green {{ background: #dfe6dc; color: #27ae60; }}
.badge-blue {{ background: #d6eaf8; color: #2980b9; }}
.badge-orange {{ background: #fdebd0; color: #e67e22; }}
@media (prefers-color-scheme: dark) {{
  :root {{ --bg: #1a1a2e; --card: #16213e; --text: #e0e0e0; --muted: #a0a0a0;
    --border: #2a2a4a; }}
  .section-header {{ background: #1e2d4a; }}
  th {{ border-bottom-color: var(--border); }}
  td {{ border-bottom-color: #1e1e3a; }}
  code {{ background: #2a2a4a; }}
  .badge-green {{ background: #1e3a2a; }} .badge-blue {{ background: #1e2a3a; }}
  .badge-orange {{ background: #3a2a1e; }}
}}
</style>
</head>
<body>
<div class="header">
  <h1>ViroForge Dataset Report</h1>
  <p>{dataset_name} &mdash; generated {timestamp}</p>
</div>
<div class="container">

<!-- summary cards -->
<div class="cards">
  <div class="card-stat"><div class="num">{n_total}</div><div class="label">Total Genomes</div></div>
  <div class="card-stat"><div class="num">{n_viral}</div><div class="label">Viral Genomes</div></div>
  <div class="card-stat"><div class="num">{platform_name}</div><div class="label">Platform</div></div>
  <div class="card-stat"><div class="num">{cov_str}</div><div class="label">Coverage</div></div>
</div>

<!-- generation summary -->
<div class="section">
  <div class="section-header">Generation Summary</div>
  <div class="section-body">
    <table class="kv">{summary_rows}</table>
  </div>
</div>

{composition_section}
{rank_abundance_chart}
{size_abundance_chart}
{coverage_chart}
{vlp_section}
{amp_section}
{dark_matter_section}
{artifact_section}
{files_section}

<div class="footer">
  Generated by ViroForge {version}
</div>
</div>
</body>
</html>"""


def _kv_rows(pairs):
    return "\n".join(
        f"<tr><th>{_ESC(str(k))}</th><td>{_ESC(str(v))}</td></tr>"
        for k, v in pairs if v and str(v) != "N/A" and str(v) != "None"
    )


def _render_composition(
    composition, n_viral, n_contam, n_total, viral_frac, contam_frac
):
    frac_html = ""
    if viral_frac is not None:
        vp = viral_frac * 100
        cp = (contam_frac or 0) * 100
        frac_html = f"""
        <p style="margin-top:.75rem">
          <span class="badge badge-green">Viral: {vp:.1f}%</span>
          <span class="badge badge-orange">Contaminant: {cp:.1f}%</span>
        </p>"""

    top_html = ""
    if composition:
        sorted_comp = sorted(
            composition, key=lambda g: g.get("relative_abundance", 0), reverse=True
        )[:10]
        rows = ""
        for i, g in enumerate(sorted_comp, 1):
            name = _ESC(str(g.get("genome_name", g.get("genome_id", "?")))[:60])
            gid = _ESC(str(g.get("genome_id", "")))
            abund = g.get("relative_abundance", 0)
            seq_type = g.get("sequence_type", g.get("genome_type", ""))
            badge_cls = "badge-green" if seq_type == "viral" else "badge-orange"
            bar_w = max(1, int(abund * 300))
            bar_color = "#00b894" if seq_type == "viral" else "#e17055"
            rows += (
                f"<tr><td>{i}</td><td title=\"{gid}\">{name}</td>"
                f"<td><span class=\"badge {badge_cls}\">{seq_type}</span></td>"
                f"<td>{abund:.4%}</td>"
                f"<td><span class=\"bar\" style=\"width:{bar_w}px;"
                f"background:{bar_color}\"></span></td></tr>"
            )
        top_html = f"""
<div class="section">
  <div class="section-header">Top 10 Most Abundant Genomes</div>
  <div class="section-body" style="overflow-x:auto">
    <table>
      <thead><tr><th>#</th><th>Genome</th><th>Type</th>
        <th>Abundance</th><th></th></tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
</div>"""

    return f"""
<div class="section">
  <div class="section-header">Genome Composition</div>
  <div class="section-body">
    <table class="kv">
      <tr><th>Viral genomes</th><td>{n_viral}</td></tr>
      <tr><th>Contaminant sequences</th><td>{n_contam}</td></tr>
      <tr><th>Total sequences</th><td>{n_total}</td></tr>
    </table>
    {frac_html}
  </div>
</div>
{top_html}"""


def _render_rank_abundance(composition):
    """Log-scale rank-abundance curve."""
    if not composition or len(composition) < 2:
        return ""
    import math

    sorted_comp = sorted(
        composition, key=lambda g: g.get("relative_abundance", 0), reverse=True
    )
    abundances = [g.get("relative_abundance", 0) for g in sorted_comp if g.get("relative_abundance", 0) > 0]
    if len(abundances) < 2:
        return ""

    chart_w, chart_h = 560, 220
    pad_l, pad_r, pad_t, pad_b = 70, 20, 20, 40
    plot_w = chart_w - pad_l - pad_r
    plot_h = chart_h - pad_t - pad_b

    log_min = math.floor(math.log10(min(abundances)))
    log_max = math.ceil(math.log10(max(abundances)))
    if log_max == log_min:
        log_max = log_min + 1

    def y_pos(val):
        if val <= 0:
            return pad_t + plot_h
        lv = math.log10(val)
        frac = (lv - log_min) / (log_max - log_min)
        return pad_t + plot_h * (1 - frac)

    def x_pos(rank):
        return pad_l + (rank / max(len(abundances) - 1, 1)) * plot_w

    # gridlines
    grid = []
    for exp in range(log_min, log_max + 1):
        y = y_pos(10 ** exp)
        grid.append(
            f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l + plot_w}" y2="{y:.1f}" '
            f'stroke="var(--border)" stroke-width="0.5"/>'
        )
        grid.append(
            f'<text x="{pad_l - 8}" y="{y + 4:.1f}" text-anchor="end" '
            f'fill="var(--muted)" font-size="10">1e{exp}</text>'
        )

    # rank axis ticks
    n = len(abundances)
    tick_step = max(1, n // 5)
    for i in range(0, n, tick_step):
        x = x_pos(i)
        grid.append(
            f'<text x="{x:.1f}" y="{pad_t + plot_h + 16}" text-anchor="middle" '
            f'fill="var(--muted)" font-size="10">{i + 1}</text>'
        )

    # polyline for the curve
    points = " ".join(f"{x_pos(i):.1f},{y_pos(a):.1f}" for i, a in enumerate(abundances))

    # dots colored by type
    dots = []
    for i, g in enumerate(sorted_comp):
        a = g.get("relative_abundance", 0)
        if a <= 0:
            continue
        stype = g.get("sequence_type", g.get("genome_type", "viral"))
        color = "#00b894" if stype == "viral" else "#e17055"
        cx, cy = x_pos(i), y_pos(a)
        name = _ESC(str(g.get("genome_name", g.get("genome_id", "")))[:50])
        dots.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3" fill="{color}" '
            f'opacity="0.8"><title>#{i+1} {name} ({a:.4%})</title></circle>'
        )

    return f"""
<div class="section">
  <div class="section-header">Rank-Abundance Curve</div>
  <div class="section-body" style="overflow-x:auto">
    <svg viewBox="0 0 {chart_w} {chart_h}"
         style="width:100%;max-width:{chart_w}px;height:auto"
         xmlns="http://www.w3.org/2000/svg">
      {chr(10).join(grid)}
      <polyline points="{points}" fill="none" stroke="#6c5ce7"
                stroke-width="2" stroke-linejoin="round"/>
      {chr(10).join(dots)}
      <text x="{pad_l + plot_w // 2}" y="{chart_h - 4}" text-anchor="middle"
            fill="var(--muted)" font-size="11">Rank</text>
      <text x="12" y="{pad_t + plot_h // 2}" text-anchor="middle"
            fill="var(--muted)" font-size="11"
            transform="rotate(-90, 12, {pad_t + plot_h // 2})">Relative Abundance</text>
    </svg>
    <p style="margin-top:.5rem;font-size:.8rem;color:var(--muted)">
      {len(abundances)} genomes ranked by abundance (log scale).
      <span class="badge badge-green">Viral</span>
      <span class="badge badge-orange">Contaminant</span>
    </p>
  </div>
</div>"""


def _render_size_vs_abundance(composition):
    """Scatter plot of genome length vs relative abundance."""
    if not composition or len(composition) < 2:
        return ""
    import math

    valid = [
        g for g in composition
        if g.get("relative_abundance", 0) > 0 and g.get("length", 0) > 0
    ]
    if len(valid) < 2:
        return ""

    chart_w, chart_h = 560, 260
    pad_l, pad_r, pad_t, pad_b = 70, 20, 20, 50
    plot_w = chart_w - pad_l - pad_r
    plot_h = chart_h - pad_t - pad_b

    lengths = [g["length"] for g in valid]
    abundances = [g["relative_abundance"] for g in valid]

    log_len_min = math.floor(math.log10(min(lengths)))
    log_len_max = math.ceil(math.log10(max(lengths)))
    log_ab_min = math.floor(math.log10(min(abundances)))
    log_ab_max = math.ceil(math.log10(max(abundances)))
    if log_len_max == log_len_min:
        log_len_max += 1
    if log_ab_max == log_ab_min:
        log_ab_max += 1

    def x_pos(length):
        lv = math.log10(length)
        return pad_l + (lv - log_len_min) / (log_len_max - log_len_min) * plot_w

    def y_pos(abund):
        lv = math.log10(abund)
        return pad_t + plot_h * (1 - (lv - log_ab_min) / (log_ab_max - log_ab_min))

    # gridlines
    grid = []
    for exp in range(log_ab_min, log_ab_max + 1):
        y = y_pos(10 ** exp)
        grid.append(
            f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l + plot_w}" y2="{y:.1f}" '
            f'stroke="var(--border)" stroke-width="0.5"/>'
        )
        grid.append(
            f'<text x="{pad_l - 8}" y="{y + 4:.1f}" text-anchor="end" '
            f'fill="var(--muted)" font-size="10">1e{exp}</text>'
        )
    for exp in range(log_len_min, log_len_max + 1):
        x = x_pos(10 ** exp)
        grid.append(
            f'<line x1="{x:.1f}" y1="{pad_t}" x2="{x:.1f}" y2="{pad_t + plot_h}" '
            f'stroke="var(--border)" stroke-width="0.5"/>'
        )
        label = f"{10**exp:,}" if exp < 6 else f"{10**exp/1e6:.0f}M"
        grid.append(
            f'<text x="{x:.1f}" y="{pad_t + plot_h + 16}" text-anchor="middle" '
            f'fill="var(--muted)" font-size="10">{label} bp</text>'
        )

    # dots
    dots = []
    for g in valid:
        stype = g.get("sequence_type", g.get("genome_type", "viral"))
        color = "#00b894" if stype == "viral" else "#e17055"
        cx, cy = x_pos(g["length"]), y_pos(g["relative_abundance"])
        name = _ESC(str(g.get("genome_name", g.get("genome_id", "")))[:50])
        dots.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" fill="{color}" '
            f'opacity="0.6"><title>{name} ({g["length"]:,} bp, '
            f'{g["relative_abundance"]:.4%})</title></circle>'
        )

    return f"""
<div class="section">
  <div class="section-header">Genome Size vs Abundance</div>
  <div class="section-body" style="overflow-x:auto">
    <svg viewBox="0 0 {chart_w} {chart_h}"
         style="width:100%;max-width:{chart_w}px;height:auto"
         xmlns="http://www.w3.org/2000/svg">
      {chr(10).join(grid)}
      {chr(10).join(dots)}
      <text x="{pad_l + plot_w // 2}" y="{chart_h - 4}" text-anchor="middle"
            fill="var(--muted)" font-size="11">Genome Length</text>
      <text x="12" y="{pad_t + plot_h // 2}" text-anchor="middle"
            fill="var(--muted)" font-size="11"
            transform="rotate(-90, 12, {pad_t + plot_h // 2})">Relative Abundance</text>
    </svg>
    <p style="margin-top:.5rem;font-size:.8rem;color:var(--muted)">
      Log-log scale. Each dot is one genome.
      <span class="badge badge-green">Viral</span>
      <span class="badge badge-orange">Contaminant</span>
    </p>
  </div>
</div>"""


def _render_coverage_distribution(metadata):
    """Bar chart of expected coverage categories."""
    ec = metadata.get("benchmarking", {}).get("expected_coverage", {})
    per_genome = ec.get("per_genome", []) if isinstance(ec, dict) else []
    if not per_genome:
        return ""

    from collections import Counter
    cats_order = ["complete", "high", "partial", "fragmented", "missing"]
    cat_colors = {
        "complete": "#00b894", "high": "#0984e3", "partial": "#fdcb6e",
        "fragmented": "#e17055", "missing": "#d63031",
    }
    counts = Counter(g.get("coverage_category", "unknown") for g in per_genome)
    total = sum(counts.values())
    if total == 0:
        return ""

    chart_w, chart_h = 480, 160
    pad_l, pad_r, pad_t, pad_b = 10, 10, 30, 40
    plot_w = chart_w - pad_l - pad_r
    plot_h = chart_h - pad_t - pad_b

    bar_items = [(cat, counts.get(cat, 0)) for cat in cats_order if counts.get(cat, 0) > 0]
    if not bar_items:
        return ""
    max_count = max(c for _, c in bar_items)
    n_bars = len(bar_items)
    bar_gap = 12
    bar_w = min(60, (plot_w - bar_gap * (n_bars - 1)) / n_bars)
    total_bar_w = n_bars * bar_w + (n_bars - 1) * bar_gap
    start_x = pad_l + (plot_w - total_bar_w) / 2

    bars = []
    for i, (cat, cnt) in enumerate(bar_items):
        x = start_x + i * (bar_w + bar_gap)
        h = max(2, (cnt / max_count) * plot_h)
        y = pad_t + plot_h - h
        color = cat_colors.get(cat, "#b2bec3")
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" '
            f'fill="{color}" rx="3"/>'
        )
        bars.append(
            f'<text x="{x + bar_w/2:.1f}" y="{y - 6:.1f}" text-anchor="middle" '
            f'fill="var(--text)" font-size="12" font-weight="600">{cnt}</text>'
        )
        bars.append(
            f'<text x="{x + bar_w/2:.1f}" y="{pad_t + plot_h + 18}" '
            f'text-anchor="middle" fill="var(--muted)" font-size="10">'
            f'{cat.title()}</text>'
        )

    pct_complete = counts.get("complete", 0) / total * 100

    return f"""
<div class="section">
  <div class="section-header">Expected Genome Recovery</div>
  <div class="section-body" style="overflow-x:auto">
    <svg viewBox="0 0 {chart_w} {chart_h}"
         style="width:100%;max-width:{chart_w}px;height:auto"
         xmlns="http://www.w3.org/2000/svg">
      <line x1="{pad_l}" y1="{pad_t + plot_h}" x2="{pad_l + plot_w}"
            y2="{pad_t + plot_h}" stroke="var(--border)" stroke-width="1"/>
      {chr(10).join(bars)}
    </svg>
    <p style="margin-top:.5rem;font-size:.8rem;color:var(--muted)">
      Coverage categories from Lander-Waterman expected coverage.
      {total} genomes total, {pct_complete:.0f}% at complete coverage.
    </p>
  </div>
</div>"""


def _render_vlp(enrichment, vlp_protocol):
    if not enrichment or vlp_protocol in (None, "none", ""):
        return ""
    viral_enrich = enrichment.get("viral_enrichment", {})
    contam_reduc = enrichment.get("contamination_reduction", {})
    mean_enrich = viral_enrich.get("mean_enrichment_factor")
    reduction_by_type = contam_reduc.get("reduction_by_type", {})

    rows = [("Protocol", vlp_protocol)]
    if mean_enrich is not None:
        rows.append(("Mean Viral Enrichment", f"{mean_enrich:.2f}x"))
    for rtype, rdata in sorted(reduction_by_type.items()):
        rf = rdata.get("reduction_factor")
        if rf is not None:
            rows.append((f"{_label(rtype)} Reduction", f"{rf*100:.1f}%"))

    chart_html = _render_vlp_composition_chart(enrichment, reduction_by_type)

    return f"""
<div class="section">
  <div class="section-header">VLP Enrichment</div>
  <div class="section-body">
    {chart_html}
    <table class="kv">{_kv_rows(rows)}</table>
  </div>
</div>"""


_PARTICLE_COLORS = {
    "viral": "#00b894",
    "host_dna": "#e17055",
    "rrna": "#fdcb6e",
    "bacterial_background": "#6c5ce7",
    "fungal_background": "#00cec9",
    "archaeal_background": "#a29bfe",
    "reagent_bacteria": "#fab1a0",
    "phix": "#636e72",
    "dietary_background": "#e84393",
}


def _render_vlp_composition_chart(enrichment, reduction_by_type):
    if not reduction_by_type:
        return ""

    viral_frac = enrichment.get("viral_fraction")
    contam_frac = enrichment.get("contamination_fraction")
    orig_total_contam = enrichment.get("contamination_reduction", {}).get(
        "original_total_contamination"
    )
    if viral_frac is None or orig_total_contam is None:
        return ""

    pre_viral = 1.0 - orig_total_contam
    post_viral = viral_frac

    pre_segments = [("viral", pre_viral)]
    post_segments = [("viral", post_viral)]
    for rtype, rdata in sorted(reduction_by_type.items()):
        orig = rdata.get("original_abundance", 0)
        reduced = rdata.get("reduced_abundance", 0)
        if orig > 0:
            pre_segments.append((rtype, orig))
        if reduced > 0:
            post_segments.append((rtype, reduced))

    bar_w = 520
    bar_h = 32
    svg_h = 160
    label_x = 130
    bar_x = label_x + 10

    def make_bar(segments, y):
        rects = []
        x = 0
        for stype, frac in segments:
            w = max(0, frac * bar_w)
            if w < 0.5:
                continue
            color = _PARTICLE_COLORS.get(stype, "#b2bec3")
            rects.append(
                f'<rect x="{bar_x + x:.1f}" y="{y}" width="{w:.1f}" '
                f'height="{bar_h}" fill="{color}" rx="2"/>'
            )
            if w > 28:
                tx = bar_x + x + w / 2
                rects.append(
                    f'<text x="{tx:.1f}" y="{y + bar_h/2 + 4:.0f}" '
                    f'text-anchor="middle" fill="#fff" font-size="11" '
                    f'font-weight="600">{frac*100:.1f}%</text>'
                )
            x += w
        return "\n    ".join(rects)

    pre_rects = make_bar(pre_segments, 10)
    post_rects = make_bar(post_segments, 56)

    legend_items = []
    seen = set()
    for segments in (pre_segments, post_segments):
        for stype, frac in segments:
            if stype not in seen and frac > 0.001:
                seen.add(stype)
                color = _PARTICLE_COLORS.get(stype, "#b2bec3")
                legend_items.append((stype, color))

    max_legend_w = bar_x + bar_w
    legend_parts = []
    lx = bar_x
    ly = 0
    legend_row_h = 18
    legend_rows = 1
    for stype, color in legend_items:
        label = _label(stype)
        item_w = len(label) * 7 + 24
        if lx + item_w > max_legend_w and lx > bar_x:
            legend_rows += 1
            lx = bar_x
            ly += legend_row_h
        legend_parts.append(
            f'<rect x="{lx}" y="{ly}" width="10" height="10" '
            f'fill="{color}" rx="2"/>'
            f'<text x="{lx + 14}" y="{ly + 9}" fill="var(--text)" '
            f'font-size="11">{_ESC(label)}</text>'
        )
        lx += item_w

    legend_h = legend_rows * legend_row_h + 8
    svg_h = 100 + legend_h
    legend_y_offset = 100

    return f"""
    <div style="margin-bottom:1rem;overflow-x:auto">
      <svg viewBox="0 0 {bar_x + bar_w + 20} {svg_h}"
           style="width:100%;max-width:{bar_x + bar_w + 20}px;height:auto"
           xmlns="http://www.w3.org/2000/svg">
        <text x="0" y="32" fill="var(--muted)" font-size="13"
              font-weight="600">Pre-VLP</text>
        {pre_rects}
        <text x="0" y="78" fill="var(--muted)" font-size="13"
              font-weight="600">Post-VLP</text>
        {post_rects}
        <line x1="{bar_x}" y1="48" x2="{bar_x + bar_w}" y2="48"
              stroke="var(--border)" stroke-width="1" stroke-dasharray="4,3"/>
        <g transform="translate(0,{legend_y_offset})">
          {chr(10).join(legend_parts)}
        </g>
      </svg>
    </div>"""


def _render_amplification(amplification):
    if not amplification:
        return ""
    rows = []
    for k, v in amplification.items():
        if isinstance(v, dict):
            continue
        rows.append((_label(k), str(v)))
    if not rows:
        return ""
    return f"""
<div class="section">
  <div class="section-header">Amplification</div>
  <div class="section-body">
    <table class="kv">{_kv_rows(rows)}</table>
  </div>
</div>"""


def _render_dark_matter(config):
    dm = config.get("dark_matter")
    if not dm:
        return ""
    rows = [
        ("Dark Matter Fraction", f"{dm.get('dark_matter_fraction', 0)*100:.0f}%"),
        ("Dark Matter Genomes", str(dm.get("dark_matter_genomes", 0))),
        ("Known Viral Genomes", str(dm.get("known_viral_genomes", 0))),
    ]
    return f"""
<div class="section">
  <div class="section-header">Dark Matter</div>
  <div class="section-body">
    <table class="kv">{_kv_rows(rows)}</table>
  </div>
</div>"""


def _render_artifacts(metadata):
    sections = []

    # adapter stats
    adapter = metadata.get("adapter_stats") or metadata.get("insert_size_adapter_stats")
    if adapter:
        rows = []
        for k in ("emergent_adapter_rate", "reads_readthrough", "reads_chimera",
                   "adapter_type", "mean_insert_size"):
            v = adapter.get(k)
            if v is not None:
                if isinstance(v, float) and v < 1:
                    rows.append((_label(k), f"{v:.2%}"))
                else:
                    rows.append((_label(k), str(v)))
        if rows:
            sections.append(("Adapter Read-Through", _kv_rows(rows)))

    # low-complexity
    lc = metadata.get("low_complexity_stats")
    if lc:
        rows = [
            ("Rate", f"{lc.get('rate', 0):.2%}"),
            ("Reads Modified", str(lc.get("reads_modified", 0))),
        ]
        counts = lc.get("artifact_counts", {})
        for atype, cnt in sorted(counts.items()):
            rows.append((_label(atype), str(cnt)))
        sections.append(("Low-Complexity Artifacts", _kv_rows(rows)))

    # duplicates
    dup = metadata.get("duplicate_stats")
    if dup:
        rows = [
            ("Rate", f"{dup.get('rate', 0):.2%}"),
            ("Duplicates Created", str(dup.get("duplicates_created", 0))),
        ]
        sections.append(("PCR Duplicates", _kv_rows(rows)))

    # chimeras
    chim = metadata.get("chimera_stats")
    if chim:
        rows = [
            ("Rate", f"{chim.get('chimera_rate', 0):.2%}"),
            ("Chimeras Created", str(chim.get("chimeras_created", 0))),
            ("Inter-Genomic", str(chim.get("inter_genomic", 0))),
            ("Intra-Genomic", str(chim.get("intra_genomic", 0))),
        ]
        sections.append(("Chimeric Reads", _kv_rows(rows)))

    # ERV
    for erv_key, erv_label in [
        ("erv_endogenous_stats", "Endogenous Retroviral (HERV)"),
        ("erv_exogenous_stats", "Exogenous Retroviral"),
    ]:
        erv = metadata.get(erv_key)
        if erv:
            rows = [(k.replace("_", " ").title(), str(v)) for k, v in erv.items()
                    if not isinstance(v, (dict, list))]
            if rows:
                sections.append((erv_label, _kv_rows(rows)))

    if not sections:
        return ""

    parts = []
    for title, content in sections:
        parts.append(
            f'<h4 style="margin-top:1rem;font-size:.95rem">{_ESC(title)}</h4>'
            f'<table class="kv">{content}</table>'
        )

    return f"""
<div class="section">
  <div class="section-header">Artifact Injection</div>
  <div class="section-body">{"".join(parts)}</div>
</div>"""


def _render_files(file_inventory):
    if not file_inventory:
        return ""
    rows = ""
    for f in file_inventory:
        rows += (
            f'<tr><td><span class="badge badge-green">OK</span></td>'
            f'<td><code>{_ESC(f["name"])}</code></td>'
            f'<td>{_ESC(f["category"])}</td>'
            f'<td>{_ESC(f["size_str"])}</td></tr>'
        )
    return f"""
<div class="section">
  <div class="section-header">Output Files</div>
  <div class="section-body" style="overflow-x:auto">
    <table>
      <thead><tr><th></th><th>File</th><th>Type</th><th>Size</th></tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
</div>"""


_CASE_FIXES = {"Gc": "GC", "Dna": "DNA", "Rna": "RNA", "Rrna": "rRNA", "Phix": "PhiX",
               "Vlp": "VLP", "Pcr": "PCR", "Erv": "ERV", "Herv": "HERV", "Mda": "MDA"}


def _label(key: str) -> str:
    label = key.replace("_", " ").replace("-", " ").title()
    for wrong, right in _CASE_FIXES.items():
        label = label.replace(wrong, right)
    return label
