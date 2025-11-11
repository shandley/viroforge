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
        console.print("[yellow]HTML format not yet implemented[/yellow]")
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

        # Platform
        platform = metadata.get('platform', {})
        plat_name = platform.get('name', 'N/A').upper()

        # Coverage/Depth
        config = metadata.get('configuration', {})
        coverage = config.get('target_coverage')
        depth = config.get('target_depth')
        cov_str = f"{coverage}x" if coverage else f"{depth}x" if depth else 'N/A'

        # Genomes
        vlp = metadata.get('vlp_enrichment', {})
        n_genomes = vlp.get('n_viral_genomes', 'N/A')

        # Viral fraction
        viral_frac = vlp.get('viral_fraction')
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
            viral_count = len([g for g in comp if g.get('genome_type') == 'viral'])
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
        platform = metadata.get('platform', {})
        plat_name = platform.get('name', 'unknown')

        if plat_name not in platforms:
            platforms[plat_name] = []

        platforms[plat_name].append(metadata)

    for plat_name, datasets in platforms.items():
        count = len(datasets)
        console.print(f"  • {plat_name.upper()}: {count} dataset(s)")


def show_recommendations(metadata_list: List[Dict]):
    """Show recommendations based on comparison."""
    console.print("[bold]Recommendations:[/bold]")

    # Check if suitable for technology comparison
    collection_ids = set(m.get('collection', {}).get('id') for m in metadata_list)
    seeds = set(m.get('generation_info', {}).get('random_seed') for m in metadata_list if m.get('generation_info'))
    platforms = set(m.get('platform', {}).get('name') for m in metadata_list if m.get('platform'))

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
    has_short = any(m.get('platform', {}).get('name') in ['novaseq', 'miseq', 'hiseq'] for m in metadata_list)
    has_long = any(m.get('platform', {}).get('name') in ['pacbio-hifi', 'nanopore'] for m in metadata_list)

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
