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
        console.print("[yellow]HTML format not yet implemented[/yellow]")
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
