"""HTML report generator for ViroForge benchmark results.

Produces a self-contained HTML file with inline SVG charts for each
benchmark module (QC, Assembly, Taxonomy).
"""

from __future__ import annotations

import html as _html
import json
from pathlib import Path
from typing import Any


def generate_benchmark_report(
    metrics: dict[str, Any],
    module: str,
    output_path: Path | str | None = None,
) -> str:
    """Generate a self-contained HTML benchmark report.

    Parameters
    ----------
    metrics : dict
        The metrics dict from a benchmark run (QC, assembly, or taxonomy).
    module : str
        One of 'qc', 'assembly', 'taxonomy'.
    output_path : Path or str, optional
        If given, write the HTML to this file and return the path.
        Otherwise return the HTML string.

    Returns
    -------
    str
        The HTML string, or the output path as a string if output_path was given.
    """
    if module == 'qc':
        body = _render_qc(metrics)
        title = 'QC Benchmark Report'
    elif module == 'assembly':
        body = _render_assembly(metrics)
        title = 'Assembly Benchmark Report'
    elif module == 'taxonomy':
        body = _render_taxonomy(metrics)
        title = 'Taxonomy Benchmark Report'
    else:
        raise ValueError(f'Unknown module: {module}')

    page = _wrap_page(title, body)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.write_text(page, encoding='utf-8')
        return str(output_path)
    return page


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _esc(text: Any) -> str:
    return _html.escape(str(text))


def _pct(v: float | None, decimals: int = 1) -> str:
    if v is None:
        return 'N/A'
    return f'{v * 100:.{decimals}f}%'


def _fmt(v: Any) -> str:
    if v is None:
        return 'N/A'
    if isinstance(v, float):
        if 0 <= v <= 1 and v not in (0.0, 1.0):
            return f'{v * 100:.1f}%'
        return f'{v:.4f}'
    if isinstance(v, int):
        return f'{v:,}'
    return str(v)


_COLORS = [
    '#00b894', '#0984e3', '#e17055', '#fdcb6e', '#6c5ce7',
    '#00cec9', '#d63031', '#e84393', '#2d3436', '#636e72',
]


def _bar_color(i: int) -> str:
    return _COLORS[i % len(_COLORS)]


def _svg_bar_chart(
    labels: list[str],
    values: list[float],
    title: str = '',
    value_fmt: str = 'pct',
    max_val: float | None = None,
    colors: list[str] | None = None,
    width: int = 600,
    bar_h: int = 28,
    gap: int = 6,
) -> str:
    if not labels:
        return ''
    n = len(labels)
    label_w = 160
    chart_w = width - label_w - 80
    top_pad = 40 if title else 10
    svg_h = top_pad + n * (bar_h + gap) + 10

    if max_val is None:
        max_val = max(values) if values else 1
    if max_val == 0:
        max_val = 1

    lines = [f'<svg viewBox="0 0 {width} {svg_h}" xmlns="http://www.w3.org/2000/svg"'
             f' style="width:100%;max-width:{width}px;font-family:system-ui,sans-serif;font-size:12px">']
    if title:
        lines.append(f'<text x="{width // 2}" y="22" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="14" font-weight="600">{_esc(title)}</text>')

    for i, (label, val) in enumerate(zip(labels, values)):
        y = top_pad + i * (bar_h + gap)
        bar_w = (val / max_val) * chart_w if max_val else 0
        c = (colors[i] if colors else _bar_color(i))
        lines.append(f'<text x="{label_w - 8}" y="{y + bar_h * 0.65}" text-anchor="end" '
                      f'fill="var(--fg)" font-size="11">{_esc(label)}</text>')
        lines.append(f'<rect x="{label_w}" y="{y}" width="{max(bar_w, 1)}" height="{bar_h}" '
                      f'rx="3" fill="{c}" opacity="0.85"/>')
        if value_fmt == 'pct':
            vtext = _pct(val)
        elif value_fmt == 'int':
            vtext = f'{int(val):,}'
        else:
            vtext = f'{val:.2f}'
        lines.append(f'<text x="{label_w + bar_w + 6}" y="{y + bar_h * 0.65}" '
                      f'fill="var(--fg)" font-size="11">{vtext}</text>')
    lines.append('</svg>')
    return '\n'.join(lines)


def _svg_grouped_bar_chart(
    groups: list[str],
    series: dict[str, list[float]],
    title: str = '',
    width: int = 600,
    bar_h: int = 22,
) -> str:
    if not groups:
        return ''
    series_names = list(series.keys())
    n_series = len(series_names)
    n_groups = len(groups)
    label_w = 100
    chart_w = width - label_w - 100
    group_h = n_series * (bar_h + 2) + 10
    top_pad = 40 if title else 10
    svg_h = top_pad + n_groups * group_h + 40

    all_vals = [v for vals in series.values() for v in vals if v is not None]
    max_val = max(all_vals) if all_vals else 1
    if max_val == 0:
        max_val = 1

    series_colors = {name: _COLORS[i % len(_COLORS)] for i, name in enumerate(series_names)}

    lines = [f'<svg viewBox="0 0 {width} {svg_h}" xmlns="http://www.w3.org/2000/svg"'
             f' style="width:100%;max-width:{width}px;font-family:system-ui,sans-serif;font-size:12px">']
    if title:
        lines.append(f'<text x="{width // 2}" y="22" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="14" font-weight="600">{_esc(title)}</text>')

    for gi, group in enumerate(groups):
        gy = top_pad + gi * group_h
        lines.append(f'<text x="{label_w - 8}" y="{gy + group_h * 0.4}" text-anchor="end" '
                      f'fill="var(--fg)" font-size="12" font-weight="500">{_esc(group)}</text>')
        for si, sname in enumerate(series_names):
            val = series[sname][gi]
            by = gy + si * (bar_h + 2)
            if val is None:
                continue
            bar_w = (val / max_val) * chart_w
            c = series_colors[sname]
            lines.append(f'<rect x="{label_w}" y="{by}" width="{max(bar_w, 1)}" height="{bar_h}" '
                          f'rx="3" fill="{c}" opacity="0.85"/>')
            lines.append(f'<text x="{label_w + bar_w + 6}" y="{by + bar_h * 0.7}" '
                          f'fill="var(--fg)" font-size="10">{_pct(val)}</text>')

    # Legend
    ly = svg_h - 20
    lx = label_w
    for sname in series_names:
        c = series_colors[sname]
        lines.append(f'<rect x="{lx}" y="{ly}" width="12" height="12" rx="2" fill="{c}"/>')
        lines.append(f'<text x="{lx + 16}" y="{ly + 10}" fill="var(--fg)" font-size="11">{_esc(sname)}</text>')
        lx += len(sname) * 7 + 30

    lines.append('</svg>')
    return '\n'.join(lines)


def _svg_stacked_bar(
    labels: list[str],
    segments: dict[str, list[float]],
    title: str = '',
    width: int = 600,
    bar_h: int = 32,
    gap: int = 8,
) -> str:
    if not labels:
        return ''
    seg_names = list(segments.keys())
    n = len(labels)
    label_w = 160
    chart_w = width - label_w - 20
    top_pad = 40 if title else 10
    legend_h = 30
    svg_h = top_pad + n * (bar_h + gap) + legend_h + 10

    seg_colors = {name: _COLORS[i % len(_COLORS)] for i, name in enumerate(seg_names)}

    lines = [f'<svg viewBox="0 0 {width} {svg_h}" xmlns="http://www.w3.org/2000/svg"'
             f' style="width:100%;max-width:{width}px;font-family:system-ui,sans-serif;font-size:12px">']
    if title:
        lines.append(f'<text x="{width // 2}" y="22" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="14" font-weight="600">{_esc(title)}</text>')

    for i, label in enumerate(labels):
        y = top_pad + i * (bar_h + gap)
        total = sum(segments[s][i] for s in seg_names)
        if total == 0:
            total = 1
        lines.append(f'<text x="{label_w - 8}" y="{y + bar_h * 0.65}" text-anchor="end" '
                      f'fill="var(--fg)" font-size="11">{_esc(label)}</text>')
        x = label_w
        for sname in seg_names:
            val = segments[sname][i]
            w = (val / total) * chart_w
            if w < 0.5:
                x += w
                continue
            lines.append(f'<rect x="{x}" y="{y}" width="{w}" height="{bar_h}" '
                          f'fill="{seg_colors[sname]}"/>')
            if w > 30:
                lines.append(f'<text x="{x + w / 2}" y="{y + bar_h * 0.65}" text-anchor="middle" '
                              f'fill="white" font-size="10" font-weight="500">{int(val)}</text>')
            x += w

    ly = svg_h - legend_h + 5
    lx = label_w
    for sname in seg_names:
        c = seg_colors[sname]
        lines.append(f'<rect x="{lx}" y="{ly}" width="12" height="12" rx="2" fill="{c}"/>')
        lines.append(f'<text x="{lx + 16}" y="{ly + 10}" fill="var(--fg)" font-size="11">{_esc(sname)}</text>')
        lx += len(sname) * 7 + 30

    lines.append('</svg>')
    return '\n'.join(lines)


def _svg_scatter(
    x_vals: list[float],
    y_vals: list[float],
    title: str = '',
    x_label: str = '',
    y_label: str = '',
    width: int = 500,
    height: int = 400,
) -> str:
    if not x_vals or not y_vals:
        return ''
    pad_l, pad_r, pad_t, pad_b = 60, 20, 40 if title else 15, 40
    cw = width - pad_l - pad_r
    ch = height - pad_t - pad_b

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    if x_max == x_min:
        x_max = x_min + 1
    if y_max == y_min:
        y_max = y_min + 1

    def tx(v):
        return pad_l + (v - x_min) / (x_max - x_min) * cw

    def ty(v):
        return pad_t + ch - (v - y_min) / (y_max - y_min) * ch

    lines = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg"'
             f' style="width:100%;max-width:{width}px;font-family:system-ui,sans-serif;font-size:11px">']
    if title:
        lines.append(f'<text x="{width // 2}" y="20" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="14" font-weight="600">{_esc(title)}</text>')

    # Diagonal reference line (y=x)
    diag_min = max(x_min, y_min)
    diag_max = min(x_max, y_max)
    if diag_max > diag_min:
        lines.append(f'<line x1="{tx(diag_min)}" y1="{ty(diag_min)}" '
                      f'x2="{tx(diag_max)}" y2="{ty(diag_max)}" '
                      f'stroke="var(--fg)" stroke-opacity="0.2" stroke-dasharray="4"/>')

    # Axes
    lines.append(f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t + ch}" '
                  f'stroke="var(--fg)" stroke-opacity="0.3"/>')
    lines.append(f'<line x1="{pad_l}" y1="{pad_t + ch}" x2="{pad_l + cw}" y2="{pad_t + ch}" '
                  f'stroke="var(--fg)" stroke-opacity="0.3"/>')

    # Axis labels
    if x_label:
        lines.append(f'<text x="{pad_l + cw / 2}" y="{height - 5}" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="11">{_esc(x_label)}</text>')
    if y_label:
        lines.append(f'<text x="14" y="{pad_t + ch / 2}" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="11" '
                      f'transform="rotate(-90,14,{pad_t + ch / 2})">{_esc(y_label)}</text>')

    # Tick labels
    for frac in (0, 0.5, 1.0):
        xv = x_min + frac * (x_max - x_min)
        yv = y_min + frac * (y_max - y_min)
        lines.append(f'<text x="{tx(xv)}" y="{pad_t + ch + 15}" text-anchor="middle" '
                      f'fill="var(--fg)" font-size="10">{xv:.2f}</text>')
        lines.append(f'<text x="{pad_l - 5}" y="{ty(yv) + 4}" text-anchor="end" '
                      f'fill="var(--fg)" font-size="10">{yv:.2f}</text>')

    # Points
    for xv, yv in zip(x_vals, y_vals):
        lines.append(f'<circle cx="{tx(xv)}" cy="{ty(yv)}" r="4" '
                      f'fill="#0984e3" opacity="0.7"/>')

    lines.append('</svg>')
    return '\n'.join(lines)


def _metric_card(label: str, value: str, color: str = '#0984e3') -> str:
    return (f'<div style="display:inline-block;min-width:140px;margin:6px;padding:12px 18px;'
            f'border-radius:8px;border:1px solid var(--border);text-align:center">'
            f'<div style="font-size:24px;font-weight:700;color:{color}">{_esc(value)}</div>'
            f'<div style="font-size:12px;color:var(--fg-muted)">{_esc(label)}</div>'
            f'</div>')


def _section(title: str, content: str) -> str:
    return (f'<div style="margin-bottom:28px">'
            f'<h3 style="border-bottom:2px solid var(--border);padding-bottom:6px;'
            f'margin-bottom:12px;color:var(--fg)">{_esc(title)}</h3>'
            f'{content}</div>')


def _table(headers: list[str], rows: list[list[str]]) -> str:
    parts = ['<div style="overflow-x:auto"><table style="width:100%;border-collapse:collapse;font-size:13px">']
    parts.append('<thead><tr>')
    for h in headers:
        parts.append(f'<th style="text-align:left;padding:6px 10px;border-bottom:2px solid var(--border);'
                     f'color:var(--fg)">{_esc(h)}</th>')
    parts.append('</tr></thead><tbody>')
    for row in rows:
        parts.append('<tr>')
        for cell in row:
            parts.append(f'<td style="padding:5px 10px;border-bottom:1px solid var(--border);'
                         f'color:var(--fg)">{cell}</td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    return '\n'.join(parts)


# ---------------------------------------------------------------------------
# QC report
# ---------------------------------------------------------------------------

def _render_qc(m: dict) -> str:
    parts = []

    # Headline metrics
    cards = ''
    cards += _metric_card('Viral Retention', _pct(m.get('viral_retention')), '#00b894')
    cards += _metric_card('Over-filtering', _pct(m.get('over_filtering_rate')), '#e17055')
    cont = m.get('contamination', {})
    cards += _metric_card('Precision', _pct(cont.get('precision')), '#0984e3')
    cards += _metric_card('Recall', _pct(cont.get('recall')), '#0984e3')
    cards += _metric_card('F1 Score', _pct(cont.get('f1')), '#6c5ce7')
    dedup = m.get('dedup', {})
    if dedup.get('n_duplicates', 0) > 0:
        cards += _metric_card('Dedup Rate', _pct(dedup.get('dedup_rate')), '#fdcb6e')
    parts.append(_section('Summary', f'<div style="display:flex;flex-wrap:wrap">{cards}</div>'))

    # Read counts
    read_info = (f'<p style="color:var(--fg)">Raw reads: <strong>{m.get("n_raw_reads", 0):,}</strong> | '
                 f'Cleaned reads: <strong>{m.get("n_cleaned_reads", 0):,}</strong> | '
                 f'Match rate: <strong>{_pct(m.get("match_rate"))}</strong></p>')
    parts.append(read_info)

    # Removal rate by source type
    rbt = m.get('removal_rate_by_type', {})
    if rbt:
        labels = list(rbt.keys())
        values = [rbt[k].get('removal_rate', 0) or 0 for k in labels]
        display_labels = [k.replace('_', ' ').title() for k in labels]
        chart = _svg_bar_chart(display_labels, values,
                               title='Contamination Removal Rate by Source Type',
                               max_val=1.0,
                               colors=['#d63031'] * len(labels))
        counts_rows = [[k.replace('_', ' ').title(), f'{rbt[k].get("n", 0):,}', _pct(rbt[k].get("removal_rate"))]
                       for k in labels]
        tbl = _table(['Source Type', 'Read Count', 'Removal Rate'], counts_rows)
        parts.append(_section('Contamination Removal by Source Type', chart + tbl))

    # Confusion matrix
    if cont:
        cm_html = ('<div style="display:inline-grid;grid-template-columns:auto auto auto;gap:2px;'
                   'text-align:center;font-size:13px;color:var(--fg)">')
        cm_html += '<div></div><div style="font-weight:600;padding:6px">Predicted Remove</div><div style="font-weight:600;padding:6px">Predicted Keep</div>'
        cm_html += (f'<div style="font-weight:600;padding:6px">Should Remove</div>'
                    f'<div style="background:#00b894;color:white;padding:10px;border-radius:4px">TP={cont.get("tp", 0):,}</div>'
                    f'<div style="background:#e17055;color:white;padding:10px;border-radius:4px">FN={cont.get("fn", 0):,}</div>')
        cm_html += (f'<div style="font-weight:600;padding:6px">Should Keep</div>'
                    f'<div style="background:#fdcb6e;color:#333;padding:10px;border-radius:4px">FP={cont.get("fp", 0):,}</div>'
                    f'<div style="background:#0984e3;color:white;padding:10px;border-radius:4px">TN={cont.get("tn", 0):,}</div>')
        cm_html += '</div>'
        parts.append(_section('Confusion Matrix', cm_html))

    return '\n'.join(parts)


# ---------------------------------------------------------------------------
# Assembly report
# ---------------------------------------------------------------------------

def _render_assembly(m: dict) -> str:
    parts = []

    # Assembly stats
    astats = m.get('assembly_stats', {})
    cards = ''
    cards += _metric_card('Contigs', f'{astats.get("n_contigs", 0):,}', '#0984e3')
    cards += _metric_card('Total bp', f'{astats.get("total_bp", 0):,}', '#0984e3')
    cards += _metric_card('N50', f'{astats.get("n50", 0):,}', '#00b894')
    cards += _metric_card('L50', f'{astats.get("l50", 0):,}', '#00b894')
    cards += _metric_card('Longest', f'{astats.get("longest", 0):,}', '#6c5ce7')
    parts.append(_section('Assembly Statistics', f'<div style="display:flex;flex-wrap:wrap">{cards}</div>'))

    # Genome recovery categories bar chart
    gr = m.get('genome_recovery', {})
    cats = gr.get('categories', {})
    if cats:
        cat_labels = ['Complete (≥95%)', 'High Quality (≥75%)', 'Partial (≥50%)', 'Fragmented (>0%)', 'Missing (0%)']
        cat_values = [cats.get('complete', 0), cats.get('high_quality', 0),
                      cats.get('partial', 0), cats.get('fragmented', 0), cats.get('missing', 0)]
        cat_colors = ['#00b894', '#0984e3', '#fdcb6e', '#e17055', '#d63031']
        chart = _svg_bar_chart(cat_labels, cat_values,
                               title='Genome Recovery',
                               value_fmt='int', colors=cat_colors)

        recovery_cards = ''
        recovery_cards += _metric_card('Genomes', f'{gr.get("n_viral_genomes", 0):,}', '#0984e3')
        recovery_cards += _metric_card('Recovered', f'{gr.get("recovered", 0):,}', '#00b894')
        recovery_cards += _metric_card('Recovery Rate', _pct(gr.get('recovery_rate')), '#00b894')
        recovery_cards += _metric_card('Mean Completeness', _pct(gr.get('mean_completeness')), '#6c5ce7')
        parts.append(_section('Genome Recovery',
                              f'<div style="display:flex;flex-wrap:wrap;margin-bottom:12px">{recovery_cards}</div>{chart}'))

    # Chimeras
    chimeras = m.get('chimeras', {})
    if chimeras.get('n', 0) > 0:
        chi_cards = ''
        chi_cards += _metric_card('Chimeric Contigs', f'{chimeras.get("n", 0):,}', '#e17055')
        chi_cards += _metric_card('Chimera Rate', _pct(chimeras.get('rate')), '#e17055')
        examples = chimeras.get('examples', [])
        tbl = ''
        if examples:
            rows = [[_esc(e.get('contig', '')), f'{e.get("length", 0):,}', str(len(e.get('segments', [])))]
                    for e in examples[:5]]
            tbl = _table(['Contig', 'Length', 'Segments'], rows)
        parts.append(_section('Chimeras', f'<div style="display:flex;flex-wrap:wrap">{chi_cards}</div>{tbl}'))

    # Observed vs expected completeness scatter
    per_genome = m.get('per_genome', {})
    if per_genome:
        x_vals, y_vals = [], []
        for gid, gdata in per_genome.items():
            exp = gdata.get('expected_completeness')
            obs = gdata.get('completeness')
            if exp is not None and obs is not None:
                x_vals.append(exp)
                y_vals.append(obs)
        if len(x_vals) >= 3:
            scatter = _svg_scatter(x_vals, y_vals,
                                   title='Observed vs Expected Completeness',
                                   x_label='Expected Completeness',
                                   y_label='Observed Completeness')
            parts.append(_section('Completeness: Observed vs Expected', scatter))

    # Abundance accuracy
    abund = m.get('abundance_accuracy', {})
    if abund.get('available'):
        ab_cards = ''
        ab_cards += _metric_card('Spearman ρ', f'{abund.get("spearman", 0):.3f}', '#0984e3')
        ab_cards += _metric_card('Mean Abs Error', f'{abund.get("mean_abs_error", 0):.4f}', '#e17055')
        parts.append(_section('Abundance Accuracy', f'<div style="display:flex;flex-wrap:wrap">{ab_cards}</div>'))

    return '\n'.join(parts)


# ---------------------------------------------------------------------------
# Taxonomy report
# ---------------------------------------------------------------------------

def _render_taxonomy(m: dict) -> str:
    parts = []

    # Headline
    cards = ''
    n_assign = m.get('n_assignments', m.get('n_contigs', 0))
    cards += _metric_card('Assignments', f'{n_assign:,}', '#0984e3')

    kv = m.get('known_viruses', {})
    cards += _metric_card('Known Sensitivity', _pct(kv.get('sensitivity')), '#00b894')
    cards += _metric_card('Known Precision', _pct(kv.get('precision')), '#00b894')

    dm = m.get('dark_matter', {})
    cards += _metric_card('Dark Matter Unclassified', _pct(dm.get('unclassified_rate')), '#6c5ce7')
    parts.append(_section('Summary', f'<div style="display:flex;flex-wrap:wrap">{cards}</div>'))

    # Known viruses vs dark matter stacked bars
    strata = {}
    for label, data in [('Known Viruses', kv), ('Dark Matter', dm)]:
        if data and data.get('n', 0) > 0:
            strata[label] = data
    if strata:
        bar_labels = list(strata.keys())
        segments = {
            'Correct': [strata[k].get('correct', 0) for k in bar_labels],
            'Unclassified': [strata[k].get('unclassified', 0) for k in bar_labels],
            'Misclassified': [strata[k].get('misclassified', 0) for k in bar_labels],
        }
        chart = _svg_stacked_bar(bar_labels, segments,
                                  title='Classification Breakdown')
        rows = []
        for label in bar_labels:
            d = strata[label]
            rows.append([label, f'{d.get("n", 0):,}', f'{d.get("correct", 0):,}',
                        f'{d.get("unclassified", 0):,}', f'{d.get("misclassified", 0):,}',
                        _pct(d.get('sensitivity')), _pct(d.get('precision'))])
        tbl = _table(['Stratum', 'Total', 'Correct', 'Unclassified', 'Misclassified',
                       'Sensitivity', 'Precision'], rows)
        parts.append(_section('Classification Breakdown', chart + tbl))

    # Per-rank accuracy
    per_rank = m.get('per_rank')
    if per_rank:
        ranks = ['species', 'genus', 'family']
        available_ranks = [r for r in ranks if r in per_rank]
        if available_ranks:
            series = {
                'Precision': [per_rank[r].get('precision') for r in available_ranks],
                'Recall': [per_rank[r].get('recall') for r in available_ranks],
                'F1': [per_rank[r].get('f1') for r in available_ranks],
            }
            display_ranks = [r.title() for r in available_ranks]
            chart = _svg_grouped_bar_chart(display_ranks, series,
                                           title='Classification Accuracy by Taxonomic Rank')
            rows = []
            for r in available_ranks:
                rd = per_rank[r]
                rows.append([r.title(), f'{rd.get("n", 0):,}', f'{rd.get("correct", 0):,}',
                            _pct(rd.get('precision')), _pct(rd.get('recall')), _pct(rd.get('f1'))])
            tbl = _table(['Rank', 'Total', 'Correct', 'Precision', 'Recall', 'F1'], rows)
            parts.append(_section('Accuracy by Taxonomic Rank', chart + tbl))

    # Abundance profile
    abund = m.get('abundance_profile', {})
    if abund.get('n_taxa', 0) > 0:
        ab_cards = ''
        ab_cards += _metric_card('Taxa', f'{abund.get("n_taxa", 0):,}', '#0984e3')
        if abund.get('bray_curtis') is not None:
            ab_cards += _metric_card('Bray-Curtis', f'{abund["bray_curtis"]:.3f}', '#e17055')
        if abund.get('pearson') is not None:
            ab_cards += _metric_card('Pearson r', f'{abund["pearson"]:.3f}', '#00b894')
        if abund.get('spearman') is not None:
            ab_cards += _metric_card('Spearman ρ', f'{abund["spearman"]:.3f}', '#00b894')
        if abund.get('mean_abs_error') is not None:
            ab_cards += _metric_card('Mean Abs Error', f'{abund["mean_abs_error"]:.4f}', '#e17055')
        parts.append(_section('Abundance Profile', f'<div style="display:flex;flex-wrap:wrap">{ab_cards}</div>'))

    return '\n'.join(parts)


# ---------------------------------------------------------------------------
# Page wrapper
# ---------------------------------------------------------------------------

def _wrap_page(title: str, body: str) -> str:
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_esc(title)} | ViroForge</title>
<style>
:root {{
  --bg: #ffffff;
  --fg: #1a1a2e;
  --fg-muted: #6c757d;
  --border: #dee2e6;
  --card-bg: #f8f9fa;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #1a1a2e;
    --fg: #e8e8e8;
    --fg-muted: #a0a0a0;
    --border: #3a3a5c;
    --card-bg: #24243e;
  }}
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
  font-family: system-ui, -apple-system, sans-serif;
  background: var(--bg);
  color: var(--fg);
  padding: 24px;
  max-width: 900px;
  margin: 0 auto;
  line-height: 1.5;
}}
h1 {{
  font-size: 22px;
  margin-bottom: 4px;
}}
h3 {{
  font-size: 16px;
}}
p {{
  margin-bottom: 8px;
  font-size: 13px;
}}
a {{ color: #0984e3; }}
.header {{
  margin-bottom: 24px;
  padding-bottom: 12px;
  border-bottom: 2px solid var(--border);
}}
.header small {{
  color: var(--fg-muted);
  font-size: 13px;
}}
</style>
</head>
<body>
<div class="header">
  <h1>{_esc(title)}</h1>
  <small>Generated by ViroForge Benchmarking Framework</small>
</div>
{body}
</body>
</html>'''
