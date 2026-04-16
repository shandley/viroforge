#!/usr/bin/env python3
"""
Compare Command

Compare multiple datasets side-by-side.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import json
from pathlib import Path
from typing import List, Dict, Optional

try:
    from rich.console import Console
    from rich.table import Table
    from rich import box
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

from .report import load_dataset_metadata, load_composition_file

console = Console()


def run_compare(args):
    """Compare multiple datasets."""
    datasets = [Path(d) for d in args.datasets]

    # Validate datasets exist
    missing = [d for d in datasets if not d.exists()]
    if missing:
        console.print("[red]Error: Datasets not found:[/red]")
        for d in missing:
            console.print(f"  • {d}")
        return 1

    # Load metadata for all datasets
    metadata_list = []
    for dataset in datasets:
        metadata = load_dataset_metadata(dataset)
        if not metadata:
            console.print(f"[yellow]Warning: Could not load metadata for {dataset.name}[/yellow]")
            metadata = {'_path': dataset, '_name': dataset.name}
        else:
            metadata['_path'] = dataset
            metadata['_name'] = dataset.name
        metadata_list.append(metadata)

    # Display comparison
    if args.format == 'terminal':
        show_terminal_comparison(metadata_list)
    elif args.format == 'json':
        show_json_comparison(metadata_list, args.export)
    elif args.format == 'html':
        output_path = show_html_comparison(metadata_list, args.export)
        if output_path:
            console.print(f"[green]✓ HTML comparison generated: {output_path}[/green]")
            # Try to open in browser
            try:
                import webbrowser
                webbrowser.open(f'file://{output_path.absolute()}')
                console.print("[dim]Opening in browser...[/dim]")
            except:
                console.print(f"[dim]Open in browser: file://{output_path.absolute()}[/dim]")
        else:
            console.print("[red]Failed to generate HTML comparison[/red]")
            return 1

    return 0


def show_terminal_comparison(metadata_list: List[Dict]):
    """Display dataset comparison in terminal."""
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(" Dataset Comparison", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Basic comparison table
    show_basic_comparison(metadata_list)
    console.print()

    # Composition consistency check
    show_composition_consistency(metadata_list)
    console.print()

    # Platform comparison
    show_platform_comparison(metadata_list)
    console.print()

    # Recommendations
    show_recommendations(metadata_list)
    console.print()


def show_basic_comparison(metadata_list: List[Dict]):
    """Show basic comparison table."""
    console.print("[bold]Dataset Summary:[/bold]")
    console.print()

    table = Table(box=box.ROUNDED)
    table.add_column("Dataset", style="cyan")
    table.add_column("Collection")
    table.add_column("Platform")
    table.add_column("Coverage/Depth")
    table.add_column("Genomes")
    table.add_column("Viral%")

    for metadata in metadata_list:
        name = metadata.get('_name', 'Unknown')

        # Collection
        collection = metadata.get('collection', {})
        coll_name = collection.get('name', 'N/A')
        coll_id = collection.get('id', '')
        coll_str = f"{coll_name} ({coll_id})" if coll_id else coll_name

        # Platform (v1.1: configuration.platform, v1.0: platform.name)
        config = metadata.get('configuration', {})
        plat_name = config.get('platform', metadata.get('platform', {}).get('name', 'N/A')).upper()

        # Coverage/Depth
        coverage = config.get('coverage') or config.get('target_coverage')
        depth = config.get('depth') or config.get('target_depth')
        cov_str = f"{coverage}x" if coverage else f"{depth}x" if depth else 'N/A'

        # Genomes (v1.1: enrichment_stats, v1.0: vlp_enrichment)
        enrichment = metadata.get('enrichment_stats', metadata.get('vlp_enrichment', {}))
        n_genomes = enrichment.get('n_viral_genomes', 'N/A')

        # Viral fraction
        viral_frac = enrichment.get('viral_fraction')
        viral_str = f"{viral_frac*100:.1f}%" if viral_frac is not None else 'N/A'

        table.add_row(name, coll_str, plat_name, cov_str, str(n_genomes), viral_str)

    console.print(table)


def show_composition_consistency(metadata_list: List[Dict]):
    """Check if datasets have consistent compositions."""
    console.print("[bold]Composition Consistency:[/bold]")

    # Load compositions
    compositions = []
    for metadata in metadata_list:
        dataset_path = metadata.get('_path')
        comp = load_composition_file(dataset_path) if dataset_path else None
        compositions.append(comp)

    # Check if all have compositions
    if all(c is None for c in compositions):
        console.print("  [dim]No composition files found[/dim]")
        return

    # Check collection IDs
    collection_ids = set()
    for metadata in metadata_list:
        coll_id = metadata.get('collection', {}).get('id')
        if coll_id:
            collection_ids.add(coll_id)

    if len(collection_ids) == 1:
        console.print("  [green]✓[/green] All datasets from same collection")
    elif len(collection_ids) > 1:
        console.print(f"  [yellow]⚠[/yellow] Datasets from {len(collection_ids)} different collections")

    # Check random seeds
    seeds = set()
    for metadata in metadata_list:
        seed = metadata.get('generation_info', {}).get('random_seed')
        if seed is not None:
            seeds.add(seed)

    if len(seeds) == 1:
        console.print(f"  [green]✓[/green] Same random seed ({list(seeds)[0]})")
    elif len(seeds) > 1:
        console.print(f"  [yellow]⚠[/yellow] Different random seeds: {sorted(seeds)}")

    # Check genome counts
    genome_counts = set()
    for comp in compositions:
        if comp:
            viral_count = len([g for g in comp if g.get('sequence_type', g.get('genome_type')) == 'viral'])
            genome_counts.add(viral_count)

    if len(genome_counts) == 1:
        console.print(f"  [green]✓[/green] Same number of viral genomes ({list(genome_counts)[0]})")
    elif len(genome_counts) > 1:
        console.print(f"  [yellow]⚠[/yellow] Different viral genome counts: {sorted(genome_counts)}")


def show_platform_comparison(metadata_list: List[Dict]):
    """Show platform-specific comparison."""
    console.print("[bold]Platform Comparison:[/bold]")

    platforms = {}
    for metadata in metadata_list:
        config = metadata.get('configuration', {})
        plat_name = config.get('platform', metadata.get('platform', {}).get('name', 'unknown'))

        if plat_name not in platforms:
            platforms[plat_name] = []

        platforms[plat_name].append(metadata)

    for plat_name, datasets in platforms.items():
        count = len(datasets)
        console.print(f"  • {plat_name.upper()}: {count} dataset(s)")


def _get_platform_name(metadata: Dict) -> Optional[str]:
    """Extract platform name from metadata (supports v1.0 and v1.1 schema)."""
    config = metadata.get('configuration', {})
    return config.get('platform', metadata.get('platform', {}).get('name'))


def show_recommendations(metadata_list: List[Dict]):
    """Show recommendations based on comparison."""
    console.print("[bold]Recommendations:[/bold]")

    # Check if suitable for technology comparison
    collection_ids = set(m.get('collection', {}).get('id') for m in metadata_list)
    seeds = set(m.get('generation_info', {}).get('random_seed') for m in metadata_list if m.get('generation_info'))
    platforms = set(_get_platform_name(m) for m in metadata_list if _get_platform_name(m))

    if len(collection_ids) == 1 and len(seeds) == 1 and len(platforms) > 1:
        console.print("  [green]✓[/green] Suitable for technology/platform comparison")
        console.print("    • Same collection and seed ensures identical genome composition")
        console.print("    • Multiple platforms enable direct performance comparison")

    elif len(collection_ids) > 1:
        console.print("  [cyan]ℹ[/cyan] Multi-collection comparison")
        console.print("    • Compare virome characteristics across different environments")
        console.print("    • Note: Different genomes, so assembly metrics not directly comparable")

    elif len(platforms) == 1:
        console.print("  [cyan]ℹ[/cyan] Same-platform comparison")
        console.print("    • Useful for testing different parameters (coverage, VLP protocol, etc.)")

    # Check if suitable for hybrid assembly
    has_short = any(_get_platform_name(m) in ['novaseq', 'miseq', 'hiseq'] for m in metadata_list)
    has_long = any(_get_platform_name(m) in ['pacbio-hifi', 'nanopore'] for m in metadata_list)

    if has_short and has_long and len(collection_ids) == 1 and len(seeds) == 1:
        console.print()
        console.print("  [green]✓[/green] Suitable for hybrid assembly!")
        console.print("    • Short + long reads with matched compositions")
        console.print("    • Try: Unicycler, SPAdes hybrid mode, or MaSuRCA")


def show_json_comparison(metadata_list: List[Dict], export_path: Optional[str]):
    """Show comparison in JSON format."""
    comparison = {
        'datasets': [
            {
                'name': m.get('_name'),
                'collection_id': m.get('collection', {}).get('id'),
                'platform': m.get('platform', {}).get('name'),
                'random_seed': m.get('generation_info', {}).get('random_seed'),
            }
            for m in metadata_list
        ]
    }

    if export_path:
        with open(export_path, 'w') as f:
            json.dump(comparison, f, indent=2)
        console.print(f"[green]Comparison exported to {export_path}[/green]")
    else:
        console.print(json.dumps(comparison, indent=2))


def show_html_comparison(metadata_list: List[Dict], export_path: Optional[str]) -> Optional[Path]:
    """Generate HTML comparison report."""
    # Determine output path
    if export_path:
        output_path = Path(export_path)
    else:
        output_path = Path("viroforge_comparison.html")

    # Generate HTML
    html_content = generate_html_comparison_content(metadata_list)

    try:
        with open(output_path, 'w') as f:
            f.write(html_content)
        return output_path
    except Exception as e:
        console.print(f"[red]Error writing HTML comparison: {e}[/red]")
        return None


def generate_html_comparison_content(metadata_list: List[Dict]) -> str:
    """Generate HTML comparison content."""
    # Build comparison table
    table_rows = ""
    for metadata in metadata_list:
        name = metadata.get('_name', 'Unknown')
        collection = metadata.get('collection', {})
        coll_name = collection.get('name', 'N/A')
        coll_id = collection.get('id', '')

        # Get platform from configuration (v1.1 structure)
        config = metadata.get('configuration', {})
        plat_name = config.get('platform', 'N/A').upper()

        coverage = config.get('coverage') or config.get('target_coverage')
        depth = config.get('depth') or config.get('target_depth')
        cov_str = f"{coverage}x" if coverage else f"{depth}x" if depth else 'N/A'

        # Get stats from enrichment_stats (v1.1 structure)
        enrichment = metadata.get('enrichment_stats', {})
        n_genomes = enrichment.get('n_viral_genomes', collection.get('n_viral_genomes', 'N/A'))

        viral_frac = enrichment.get('viral_fraction')
        viral_str = f"{viral_frac*100:.1f}%" if viral_frac is not None else 'N/A'

        seed = metadata.get('generation_info', {}).get('random_seed', 'N/A')

        table_rows += f"""
        <tr>
            <td><strong>{name}</strong></td>
            <td>{coll_name}<br><small class="text-muted">ID: {coll_id}</small></td>
            <td><span class="badge bg-info">{plat_name}</span></td>
            <td>{cov_str}</td>
            <td>{n_genomes}</td>
            <td>{viral_str}</td>
            <td><code>{seed}</code></td>
        </tr>
        """

    # Check consistency
    collection_ids = set(m.get('collection', {}).get('id') for m in metadata_list if m.get('collection', {}).get('id'))
    seeds = set(m.get('generation_info', {}).get('random_seed') for m in metadata_list if m.get('generation_info'))
    # Get platforms from configuration (v1.1 structure)
    platforms = set(m.get('configuration', {}).get('platform') for m in metadata_list if m.get('configuration', {}).get('platform'))

    # Remove None values
    collection_ids.discard(None)
    seeds.discard(None)
    platforms.discard(None)

    # Consistency checks
    consistency_html = ""

    if len(collection_ids) == 1:
        consistency_html += '<div class="alert alert-success"><strong>✓</strong> All datasets from same collection</div>'
    elif len(collection_ids) > 1:
        consistency_html += f'<div class="alert alert-warning"><strong>⚠</strong> Datasets from {len(collection_ids)} different collections</div>'

    if len(seeds) == 1:
        consistency_html += f'<div class="alert alert-success"><strong>✓</strong> Same random seed ({list(seeds)[0]})</div>'
    elif len(seeds) > 1:
        consistency_html += f'<div class="alert alert-warning"><strong>⚠</strong> Different random seeds: {sorted(seeds)}</div>'

    # Recommendations
    recommendations_html = ""

    # Technology comparison
    if len(collection_ids) == 1 and len(seeds) == 1 and len(platforms) > 1:
        recommendations_html += """
        <div class="alert alert-success">
            <h5>✓ Suitable for Technology/Platform Comparison</h5>
            <ul>
                <li>Same collection and seed ensures identical genome composition</li>
                <li>Multiple platforms enable direct performance comparison</li>
                <li>Ideal for benchmarking sequencing technologies</li>
            </ul>
        </div>
        """

    # Multi-collection comparison
    elif len(collection_ids) > 1:
        recommendations_html += """
        <div class="alert alert-info">
            <h5>ℹ Multi-Collection Comparison</h5>
            <ul>
                <li>Compare virome characteristics across different environments</li>
                <li>Note: Different genomes, so assembly metrics not directly comparable</li>
                <li>Useful for ecosystem or disease state comparisons</li>
            </ul>
        </div>
        """

    # Same platform comparison
    elif len(platforms) == 1:
        recommendations_html += """
        <div class="alert alert-info">
            <h5>ℹ Same-Platform Comparison</h5>
            <ul>
                <li>Useful for testing different parameters</li>
                <li>Compare VLP protocols, coverage levels, or amplification methods</li>
            </ul>
        </div>
        """

    # Hybrid assembly check
    has_short = any(m.get('platform', {}).get('name', m.get('configuration', {}).get('platform')) in ['novaseq', 'miseq', 'hiseq'] for m in metadata_list)
    has_long = any(m.get('platform', {}).get('name', m.get('configuration', {}).get('platform')) in ['pacbio-hifi', 'nanopore'] for m in metadata_list)

    if has_short and has_long and len(collection_ids) == 1 and len(seeds) == 1:
        recommendations_html += """
        <div class="alert alert-success">
            <h5>✓ Suitable for Hybrid Assembly!</h5>
            <ul>
                <li>Short + long reads with matched compositions</li>
                <li>Recommended assemblers: Unicycler, SPAdes hybrid mode, MaSuRCA</li>
                <li>Enables validation of hybrid assembly strategies</li>
            </ul>
        </div>
        """

    # Platform comparison chart
    platform_chart_html = ""
    if platforms:
        platform_counts = {}
        for metadata in metadata_list:
            plat = metadata.get('platform', {}).get('name', metadata.get('configuration', {}).get('platform', 'unknown'))
            platform_counts[plat] = platform_counts.get(plat, 0) + 1

        platform_rows = ""
        for plat, count in sorted(platform_counts.items()):
            bar_width = (count / len(metadata_list)) * 100
            platform_rows += f"""
            <tr>
                <td><span class="badge bg-info">{plat.upper()}</span></td>
                <td>{count} dataset(s)</td>
                <td>
                    <div class="progress" style="height: 20px;">
                        <div class="progress-bar bg-info" role="progressbar"
                             style="width: {bar_width}%">{count}</div>
                    </div>
                </td>
            </tr>
            """

        platform_chart_html = f"""
        <div class="card mb-4">
            <div class="card-header bg-secondary text-white">
                <h5 class="mb-0">Platform Distribution</h5>
            </div>
            <div class="card-body">
                <table class="table">
                    <thead>
                        <tr>
                            <th style="width: 30%;">Platform</th>
                            <th style="width: 20%;">Count</th>
                            <th style="width: 50%;">Visual</th>
                        </tr>
                    </thead>
                    <tbody>
                        {platform_rows}
                    </tbody>
                </table>
            </div>
        </div>
        """

    # Build HTML
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ViroForge Dataset Comparison</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: #f8f9fa; }}
        .comparison-header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem;
            margin-bottom: 2rem;
        }}
        .card {{ margin-bottom: 1.5rem; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .card-header {{ font-weight: bold; }}
        .progress {{ background-color: #e9ecef; }}
    </style>
</head>
<body>
    <div class="comparison-header">
        <div class="container">
            <h1><i class="bi bi-file-earmark-diff"></i> ViroForge Dataset Comparison</h1>
            <h3>Comparing {len(metadata_list)} Datasets</h3>
        </div>
    </div>

    <div class="container">
        <!-- Comparison Table -->
        <div class="card mb-4">
            <div class="card-header bg-primary text-white">
                <h5 class="mb-0">Dataset Summary</h5>
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-hover">
                        <thead>
                            <tr>
                                <th>Dataset</th>
                                <th>Collection</th>
                                <th>Platform</th>
                                <th>Coverage/Depth</th>
                                <th>Genomes</th>
                                <th>Viral%</th>
                                <th>Seed</th>
                            </tr>
                        </thead>
                        <tbody>
                            {table_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Consistency Checks -->
        <div class="card mb-4">
            <div class="card-header bg-success text-white">
                <h5 class="mb-0">Consistency Checks</h5>
            </div>
            <div class="card-body">
                {consistency_html if consistency_html else '<p class="text-muted">No consistency issues detected</p>'}
            </div>
        </div>

        <!-- Platform Chart -->
        {platform_chart_html}

        <!-- Recommendations -->
        {f'''
        <div class="card mb-4">
            <div class="card-header bg-warning text-dark">
                <h5 class="mb-0">Recommendations</h5>
            </div>
            <div class="card-body">
                {recommendations_html}
            </div>
        </div>
        ''' if recommendations_html else ''}

        <!-- Footer -->
        <div class="text-center text-muted mt-4 mb-4">
            <p>Generated by ViroForge 0.11.0</p>
            <p><small>Synthetic Virome Data Generator | <a href="https://github.com/hecatomb/viroforge">github.com/hecatomb/viroforge</a></small></p>
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
    """

    return html
