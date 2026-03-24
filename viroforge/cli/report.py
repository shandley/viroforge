#!/usr/bin/env python3
"""
Report Command

Generate quality reports for datasets.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import json
from pathlib import Path
from typing import Dict, Optional
import subprocess

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich import box
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()


def run_report(args):
    """Generate dataset report."""
    dataset_path = Path(args.dataset)

    if not dataset_path.exists():
        console.print(f"[red]Error: Dataset not found: {dataset_path}[/red]")
        return 1

    # Load metadata
    metadata = load_dataset_metadata(dataset_path)

    if not metadata:
        console.print(f"[red]Error: Could not load dataset metadata from {dataset_path}[/red]")
        console.print("[dim]Make sure the dataset was generated with ViroForge[/dim]")
        return 1

    # Display report
    if args.format == 'terminal':
        show_terminal_report(metadata, dataset_path)
    elif args.format == 'json':
        show_json_report(metadata, args.export)
    elif args.format == 'html':
        output_path = show_html_report(metadata, dataset_path, args.export)
        if output_path:
            console.print(f"[green]✓ HTML report generated: {output_path}[/green]")
            # Try to open in browser
            try:
                import webbrowser
                webbrowser.open(f'file://{output_path.absolute()}')
                console.print("[dim]Opening in browser...[/dim]")
            except:
                console.print(f"[dim]Open in browser: file://{output_path.absolute()}[/dim]")
        else:
            console.print("[red]Failed to generate HTML report[/red]")
            return 1

    return 0


def load_dataset_metadata(dataset_path: Path) -> Optional[Dict]:
    """Load dataset metadata from standard ViroForge locations."""
    # Try multiple possible metadata locations
    possible_locations = [
        dataset_path / "metadata" / f"{dataset_path.name}_metadata.json",
        dataset_path / f"{dataset_path.name}_metadata.json",
        dataset_path / "metadata.json",
    ]

    # Also try glob pattern
    metadata_dir = dataset_path / "metadata"
    if metadata_dir.exists():
        metadata_files = list(metadata_dir.glob("*_metadata.json"))
        if metadata_files:
            possible_locations.insert(0, metadata_files[0])

    for metadata_path in possible_locations:
        if metadata_path.exists():
            try:
                with open(metadata_path) as f:
                    return json.load(f)
            except Exception as e:
                console.print(f"[dim]Warning: Could not load {metadata_path}: {e}[/dim]")

    return None


def show_terminal_report(metadata: Dict, dataset_path: Path):
    """Display report in terminal."""
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(f" Dataset Report: {dataset_path.name}", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Generation Summary
    show_generation_summary(metadata)
    console.print()

    # Platform Information
    show_platform_info(metadata)
    console.print()

    # Composition Summary
    show_composition_summary(metadata, dataset_path)
    console.print()

    # VLP Enrichment
    if metadata.get('vlp_enrichment'):
        show_vlp_summary(metadata)
        console.print()

    # File Summary
    show_file_summary(dataset_path)
    console.print()


def show_generation_summary(metadata: Dict):
    """Show generation summary."""
    console.print("[bold]Generation Summary:[/bold]")

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Property", style="cyan")
    table.add_column("Value")

    # Collection info
    if 'collection' in metadata:
        collection = metadata['collection']
        table.add_row("Collection", f"{collection.get('name', 'Unknown')} (ID: {collection.get('id', 'N/A')})")

    # Generation info
    if 'generation_info' in metadata:
        gen_info = metadata['generation_info']
        if 'timestamp' in gen_info:
            table.add_row("Generated", gen_info['timestamp'])
        if 'random_seed' in gen_info:
            table.add_row("Random Seed", str(gen_info['random_seed']))
        if 'viroforge_version' in gen_info:
            table.add_row("ViroForge Version", gen_info['viroforge_version'])

    console.print(table)


def show_platform_info(metadata: Dict):
    """Show platform and sequencing information."""
    console.print("[bold]Platform Information:[/bold]")

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Property", style="cyan")
    table.add_column("Value")

    if 'platform' in metadata:
        platform = metadata['platform']
        table.add_row("Platform", platform.get('name', 'Unknown').upper())

        if 'read_type' in platform:
            table.add_row("Read Type", platform['read_type'])

    # Configuration
    if 'configuration' in metadata:
        config = metadata['configuration']

        if 'target_coverage' in config:
            table.add_row("Target Coverage", f"{config['target_coverage']}x")
        elif 'target_depth' in config:
            table.add_row("Target Depth", f"{config['target_depth']}x")

        if 'read_length' in config:
            table.add_row("Read Length", f"{config['read_length']} bp")

        if 'insert_size' in config:
            table.add_row("Insert Size", f"{config['insert_size']} bp")

    console.print(table)


def show_composition_summary(metadata: Dict, dataset_path: Path):
    """Show genome composition summary."""
    console.print("[bold]Genome Composition:[/bold]")

    # Try to load composition file
    composition = load_composition_file(dataset_path)

    if composition is not None:
        viral_count = len([g for g in composition if g.get('genome_type') == 'viral'])
        total_count = len(composition)

        console.print(f"  Total genomes: [green]{total_count}[/green]")
        console.print(f"  Viral genomes: [green]{viral_count}[/green]")

        if 'vlp_enrichment' in metadata:
            vlp = metadata['vlp_enrichment']
            if 'viral_fraction' in vlp:
                console.print(f"  Viral fraction: [green]{vlp['viral_fraction']*100:.1f}%[/green]")
            if 'contamination_fraction' in vlp:
                console.print(f"  Contamination: [yellow]{vlp['contamination_fraction']*100:.1f}%[/yellow]")

        # Show top genomes
        if composition and len(composition) > 0:
            console.print()
            console.print("[bold]Top 5 Most Abundant:[/bold]")
            sorted_comp = sorted(composition, key=lambda x: x.get('relative_abundance', 0), reverse=True)
            for i, genome in enumerate(sorted_comp[:5], 1):
                name = genome.get('genome_name', genome.get('genome_id', 'Unknown'))[:50]
                abundance = genome.get('relative_abundance', 0)
                bar = "█" * int(abundance * 40)
                console.print(f"  {i}. {name:50s} {abundance:6.2%} {bar}")
    else:
        console.print("[dim]Composition file not found[/dim]")


def show_vlp_summary(metadata: Dict):
    """Show VLP enrichment summary."""
    console.print("[bold]VLP Enrichment:[/bold]")

    vlp = metadata['vlp_enrichment']

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Property", style="cyan")
    table.add_column("Value")

    if 'protocol' in vlp:
        table.add_row("Protocol", vlp['protocol'])

    if 'viral_enrichment' in vlp:
        enrich = vlp['viral_enrichment']
        if 'mean_enrichment' in enrich:
            table.add_row("Mean Viral Enrichment", f"{enrich['mean_enrichment']:.2f}x")

    if 'contamination_reduction' in vlp:
        contam = vlp['contamination_reduction']
        if 'bacterial_reduction' in contam:
            table.add_row("Bacterial Reduction", f"{contam['bacterial_reduction']*100:.1f}%")
        if 'host_reduction' in contam:
            table.add_row("Host Reduction", f"{contam['host_reduction']*100:.1f}%")

    console.print(table)


def show_file_summary(dataset_path: Path):
    """Show file summary."""
    console.print("[bold]Output Files:[/bold]")

    # Look for FASTQ files
    fastq_dir = dataset_path / "fastq"
    if fastq_dir.exists():
        fastq_files = list(fastq_dir.glob("*.fastq*"))
        if fastq_files:
            for fq_file in sorted(fastq_files):
                size_mb = fq_file.stat().st_size / (1024 * 1024)
                console.print(f"  [green]✓[/green] {fq_file.name:40s} {size_mb:8.1f} MB")

    # Look for FASTA files
    fasta_dir = dataset_path / "fasta"
    if fasta_dir.exists():
        fasta_files = list(fasta_dir.glob("*.fasta"))
        if fasta_files:
            for fa_file in sorted(fasta_files):
                size_mb = fa_file.stat().st_size / (1024 * 1024)
                console.print(f"  [green]✓[/green] {fa_file.name:40s} {size_mb:8.1f} MB")

    # Look for metadata files
    metadata_dir = dataset_path / "metadata"
    if metadata_dir.exists():
        metadata_files = list(metadata_dir.glob("*"))
        if metadata_files:
            for md_file in sorted(metadata_files):
                size_kb = md_file.stat().st_size / 1024
                console.print(f"  [green]✓[/green] {md_file.name:40s} {size_kb:8.1f} KB")


def load_composition_file(dataset_path: Path) -> Optional[list]:
    """Load composition TSV file."""
    # Try multiple locations
    possible_locations = [
        dataset_path / "metadata" / f"{dataset_path.name}_composition.tsv",
    ]

    # Also try glob
    metadata_dir = dataset_path / "metadata"
    if metadata_dir.exists():
        comp_files = list(metadata_dir.glob("*_composition.tsv"))
        if comp_files:
            possible_locations.insert(0, comp_files[0])

    for comp_path in possible_locations:
        if comp_path.exists():
            try:
                import pandas as pd
                df = pd.read_csv(comp_path, sep='\t')
                return df.to_dict('records')
            except Exception:
                # If pandas not available, try basic parsing
                try:
                    composition = []
                    with open(comp_path) as f:
                        headers = f.readline().strip().split('\t')
                        for line in f:
                            values = line.strip().split('\t')
                            entry = dict(zip(headers, values))
                            # Convert abundance to float
                            if 'relative_abundance' in entry:
                                try:
                                    entry['relative_abundance'] = float(entry['relative_abundance'])
                                except:
                                    pass
                            composition.append(entry)
                    return composition
                except Exception:
                    pass

    return None


def show_json_report(metadata: Dict, export_path: Optional[str]):
    """Show report in JSON format."""
    if export_path:
        with open(export_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        console.print(f"[green]Report exported to {export_path}[/green]")
    else:
        console.print(json.dumps(metadata, indent=2))


def show_html_report(metadata: Dict, dataset_path: Path, export_path: Optional[str]) -> Optional[Path]:
    """Generate HTML report."""
    # Determine output path
    if export_path:
        output_path = Path(export_path)
    else:
        output_path = dataset_path / f"{dataset_path.name}_report.html"

    # Load composition
    composition = load_composition_file(dataset_path)

    # Generate HTML
    html_content = generate_html_report_content(metadata, dataset_path, composition)

    try:
        with open(output_path, 'w') as f:
            f.write(html_content)
        return output_path
    except Exception as e:
        console.print(f"[red]Error writing HTML report: {e}[/red]")
        return None


def generate_html_report_content(metadata: Dict, dataset_path: Path, composition: Optional[list]) -> str:
    """Generate HTML report content."""
    # Extract metadata fields (v1.1 structure)
    collection = metadata.get('collection', {})
    gen_info = metadata.get('generation_info', {})
    config = metadata.get('configuration', {})
    enrichment = metadata.get('enrichment_stats', {})

    # Collection info
    coll_name = collection.get('name', 'Unknown')
    coll_id = collection.get('id', 'N/A')

    # Generation info
    timestamp = gen_info.get('timestamp', 'N/A')
    seed = gen_info.get('random_seed', 'N/A')
    version = gen_info.get('viroforge_version', 'N/A')

    # Platform info (from configuration)
    platform_name = config.get('platform', 'Unknown').upper()
    # Determine read type from platform
    read_type = 'long' if platform_name.lower() in ['pacbio-hifi', 'nanopore'] else 'short'
    coverage = config.get('coverage') or config.get('target_coverage')
    depth = config.get('depth') or config.get('target_depth')
    cov_str = f"{coverage}x" if coverage else f"{depth}x" if depth else 'N/A'
    read_length = config.get('read_length', 'N/A')
    insert_size = config.get('insert_size', 'N/A')

    # Composition info
    viral_count = 0
    total_count = 0
    viral_fraction = None
    contam_fraction = None
    top_genomes_html = ""

    if composition:
        viral_count = len([g for g in composition if g.get('genome_type') == 'viral'])
        total_count = len(composition)

        # Get viral fraction from enrichment_stats
        if enrichment:
            viral_fraction = enrichment.get('viral_fraction')
            contam_fraction = enrichment.get('contamination_fraction')

        # Top 5 genomes
        sorted_comp = sorted(composition, key=lambda x: x.get('relative_abundance', 0), reverse=True)
        for i, genome in enumerate(sorted_comp[:5], 1):
            name = genome.get('genome_name', genome.get('genome_id', 'Unknown'))
            abundance = genome.get('relative_abundance', 0)
            bar_width = int(abundance * 100)
            top_genomes_html += f"""
            <tr>
                <td>{i}</td>
                <td class="text-truncate" style="max-width: 300px;" title="{name}">{name}</td>
                <td>{abundance:.4%}</td>
                <td>
                    <div class="progress" style="height: 20px;">
                        <div class="progress-bar bg-success" role="progressbar"
                             style="width: {bar_width}%" aria-valuenow="{bar_width}"
                             aria-valuemin="0" aria-valuemax="100">{abundance:.2%}</div>
                    </div>
                </td>
            </tr>
            """

    # VLP info (from enrichment_stats)
    vlp_html = ""
    vlp_protocol = config.get('vlp_protocol')
    if enrichment and vlp_protocol and vlp_protocol != 'none':
        viral_enrich = enrichment.get('viral_enrichment', {})
        contam_reduc = enrichment.get('contamination_reduction', {})

        # Get mean enrichment factor
        mean_enrichment = viral_enrich.get('mean_enrichment_factor', 'N/A')

        # Get reduction by type for host and bacterial
        reduction_by_type = contam_reduc.get('reduction_by_type', {})
        host_reduction = reduction_by_type.get('host_dna', {}).get('reduction_factor')
        bacterial_reduction = reduction_by_type.get('reagent_bacteria', {}).get('reduction_factor')

        vlp_html = f"""
        <div class="col-md-6">
            <div class="card h-100">
                <div class="card-header bg-info text-white">
                    <h5 class="mb-0">VLP Enrichment</h5>
                </div>
                <div class="card-body">
                    <table class="table table-sm">
                        <tr><th>Protocol:</th><td>{vlp_protocol}</td></tr>
                        <tr><th>Mean Enrichment:</th><td>{mean_enrichment if mean_enrichment != 'N/A' else mean_enrichment}{"x" if mean_enrichment != 'N/A' else ""}</td></tr>
                        {f'<tr><th>Bacterial Reduction:</th><td>{bacterial_reduction*100:.1f}%</td></tr>' if bacterial_reduction else ''}
                        {f'<tr><th>Host Reduction:</th><td>{host_reduction*100:.1f}%</td></tr>' if host_reduction else ''}
                    </table>
                </div>
            </div>
        </div>
        """

    # File summary
    file_list_html = ""
    for subdir in ['fastq', 'fasta', 'metadata']:
        dir_path = dataset_path / subdir
        if dir_path.exists():
            files = sorted(dir_path.glob("*"))
            for file in files:
                if file.is_file():
                    size_mb = file.stat().st_size / (1024 * 1024)
                    size_str = f"{size_mb:.1f} MB" if size_mb >= 1 else f"{file.stat().st_size / 1024:.1f} KB"
                    file_list_html += f"""
                    <tr>
                        <td><span class="badge bg-success">✓</span></td>
                        <td><code>{file.name}</code></td>
                        <td>{subdir.upper()}</td>
                        <td>{size_str}</td>
                    </tr>
                    """

    # Build HTML
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ViroForge Dataset Report - {dataset_path.name}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: #f8f9fa; }}
        .report-header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem;
            margin-bottom: 2rem;
        }}
        .card {{ margin-bottom: 1.5rem; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .card-header {{ font-weight: bold; }}
        th {{ width: 40%; }}
        .progress {{ background-color: #e9ecef; }}
    </style>
</head>
<body>
    <div class="report-header">
        <div class="container">
            <h1><i class="bi bi-file-earmark-bar-graph"></i> ViroForge Dataset Report</h1>
            <h3>{dataset_path.name}</h3>
            <p class="mb-0">Generated: {timestamp}</p>
        </div>
    </div>

    <div class="container">
        <!-- Summary Cards -->
        <div class="row mb-4">
            <div class="col-md-3">
                <div class="card text-center bg-primary text-white">
                    <div class="card-body">
                        <h2 class="mb-0">{total_count}</h2>
                        <p class="mb-0">Total Genomes</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-center bg-success text-white">
                    <div class="card-body">
                        <h2 class="mb-0">{viral_count}</h2>
                        <p class="mb-0">Viral Genomes</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-center bg-info text-white">
                    <div class="card-body">
                        <h2 class="mb-0">{platform_name}</h2>
                        <p class="mb-0">Platform</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card text-center bg-warning text-dark">
                    <div class="card-body">
                        <h2 class="mb-0">{cov_str}</h2>
                        <p class="mb-0">Coverage/Depth</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- Main Content -->
        <div class="row">
            <!-- Generation Summary -->
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h5 class="mb-0">Generation Summary</h5>
                    </div>
                    <div class="card-body">
                        <table class="table table-sm">
                            <tr><th>Collection:</th><td>{coll_name} (ID: {coll_id})</td></tr>
                            <tr><th>Generated:</th><td>{timestamp}</td></tr>
                            <tr><th>Random Seed:</th><td>{seed}</td></tr>
                            <tr><th>ViroForge Version:</th><td>{version}</td></tr>
                        </table>
                    </div>
                </div>
            </div>

            <!-- Platform Information -->
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header bg-secondary text-white">
                        <h5 class="mb-0">Platform Information</h5>
                    </div>
                    <div class="card-body">
                        <table class="table table-sm">
                            <tr><th>Platform:</th><td>{platform_name}</td></tr>
                            <tr><th>Read Type:</th><td>{read_type}</td></tr>
                            <tr><th>Coverage/Depth:</th><td>{cov_str}</td></tr>
                            <tr><th>Read Length:</th><td>{read_length} bp</td></tr>
                            {f'<tr><th>Insert Size:</th><td>{insert_size} bp</td></tr>' if insert_size != 'N/A' else ''}
                        </table>
                    </div>
                </div>
            </div>
        </div>

        <!-- Composition & VLP -->
        <div class="row">
            <!-- Genome Composition -->
            <div class="col-md-{6 if vlp_html else 12}">
                <div class="card">
                    <div class="card-header bg-success text-white">
                        <h5 class="mb-0">Genome Composition</h5>
                    </div>
                    <div class="card-body">
                        <table class="table table-sm">
                            <tr><th>Total Genomes:</th><td>{total_count}</td></tr>
                            <tr><th>Viral Genomes:</th><td>{viral_count}</td></tr>
                            {f'<tr><th>Viral Fraction:</th><td>{viral_fraction*100:.1f}%</td></tr>' if viral_fraction else ''}
                            {f'<tr><th>Contamination:</th><td>{contam_fraction*100:.1f}%</td></tr>' if contam_fraction else ''}
                        </table>
                    </div>
                </div>
            </div>

            {vlp_html}
        </div>

        <!-- Top Genomes -->
        {f'''
        <div class="row">
            <div class="col-12">
                <div class="card">
                    <div class="card-header bg-dark text-white">
                        <h5 class="mb-0">Top 5 Most Abundant Genomes</h5>
                    </div>
                    <div class="card-body">
                        <table class="table table-hover">
                            <thead>
                                <tr>
                                    <th style="width: 5%;">#</th>
                                    <th style="width: 40%;">Genome</th>
                                    <th style="width: 15%;">Abundance</th>
                                    <th style="width: 40%;">Visual</th>
                                </tr>
                            </thead>
                            <tbody>
                                {top_genomes_html}
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>
        ''' if top_genomes_html else ''}

        <!-- Output Files -->
        <div class="row">
            <div class="col-12">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h5 class="mb-0">Output Files</h5>
                    </div>
                    <div class="card-body">
                        <table class="table table-sm table-hover">
                            <thead>
                                <tr>
                                    <th style="width: 5%;">Status</th>
                                    <th style="width: 50%;">Filename</th>
                                    <th style="width: 20%;">Type</th>
                                    <th style="width: 25%;">Size</th>
                                </tr>
                            </thead>
                            <tbody>
                                {file_list_html if file_list_html else '<tr><td colspan="4" class="text-center text-muted">No files found</td></tr>'}
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>

        <!-- Footer -->
        <div class="text-center text-muted mt-4 mb-4">
            <p>Generated by ViroForge {version}</p>
            <p><small>Synthetic Virome Data Generator | <a href="https://github.com/hecatomb/viroforge">github.com/hecatomb/viroforge</a></small></p>
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
    """

    return html
