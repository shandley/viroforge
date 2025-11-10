#!/usr/bin/env python3
"""
Presets Command

Manage configuration presets.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import json
from pathlib import Path

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.columns import Columns
    from rich.markdown import Markdown
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

from .preset_loader import (
    list_all_presets,
    load_preset,
    get_preset_names,
    create_preset_from_dict,
    validate_preset
)

console = Console()


# Category display names
CATEGORY_NAMES = {
    'virome_type': 'Virome Type Presets',
    'assembly_benchmarking': 'Assembly Benchmarking',
    'testing': 'Quick Testing',
    'technology_comparison': 'Technology Comparison',
    'method_comparison': 'Method Comparison',
    'other': 'Other'
}


def run_presets(args):
    """Manage presets."""
    if not args.presets_command or args.presets_command == 'list':
        list_presets()
    elif args.presets_command == 'show':
        show_preset(args.name)
    elif args.presets_command == 'create':
        create_preset(args.name, args.from_dataset)
    else:
        console.print(f"[red]Unknown presets command: {args.presets_command}[/red]")
        return 1

    return 0


def list_presets():
    """List available presets."""
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(" " * 25 + "ViroForge Presets", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    presets_by_category = list_all_presets()

    if not presets_by_category:
        console.print("[yellow]No presets found.[/yellow]")
        console.print()
        return

    for category_key, presets in sorted(presets_by_category.items()):
        category_name = CATEGORY_NAMES.get(category_key, category_key.replace('_', ' ').title())

        console.print(f"[bold blue]{category_name}:[/bold blue]")

        for preset in presets:
            name = preset['preset_name']
            desc = preset.get('description', 'No description')
            source_icon = "📦" if preset['source'] == 'built-in' else "👤"

            console.print(f"  {source_icon} [cyan]{name:30s}[/cyan] {desc}")

        console.print()

    console.print("[dim]Use 'viroforge presets show <name>' for details[/dim]")
    console.print()


def show_preset(name: str):
    """Show preset details."""
    preset = load_preset(name)

    if not preset:
        console.print(f"[red]Preset '{name}' not found.[/red]")
        console.print()
        console.print("Available presets:")
        for preset_name in get_preset_names():
            console.print(f"  • {preset_name}")
        console.print()
        return

    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(f" Preset: {preset['name']}", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Description
    if preset.get('description'):
        console.print("[bold]Description:[/bold]")
        console.print(f"  {preset['description']}")
        console.print()

    # Parameters
    console.print("[bold]Parameters:[/bold]")
    params = preset.get('parameters', {})

    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Parameter", style="cyan")
    table.add_column("Value")

    for key, value in params.items():
        if isinstance(value, bool):
            value_str = "✓" if value else "✗"
        else:
            value_str = str(value)
        table.add_row(key, value_str)

    console.print(table)
    console.print()

    # Metadata
    if preset.get('metadata'):
        metadata = preset['metadata']

        if metadata.get('recommended_for'):
            console.print("[bold]Recommended For:[/bold]")
            for use_case in metadata['recommended_for']:
                console.print(f"  • {use_case}")
            console.print()

        if metadata.get('estimated_time'):
            console.print(f"[bold]Estimated Time:[/bold] {metadata['estimated_time']}")

        if metadata.get('estimated_size'):
            console.print(f"[bold]Estimated Size:[/bold] {metadata['estimated_size']}")

        if metadata.get('tags'):
            tags = ' '.join([f"[dim]#{tag}[/dim]" for tag in metadata['tags']])
            console.print(f"[bold]Tags:[/bold] {tags}")

        console.print()

    # Usage
    console.print("[bold cyan]Usage:[/bold cyan]")
    console.print(f"  viroforge generate --preset {name}")
    console.print()
    console.print("[dim]Override parameters:[/dim]")
    console.print(f"  viroforge generate --preset {name} --seed 123 --output my_data")
    console.print()


def create_preset(name: str, from_dataset: str = None):
    """Create a new preset."""
    console.print()
    console.print(f"[bold cyan]Creating preset: {name}[/bold cyan]")
    console.print()

    if from_dataset:
        # Load metadata from existing dataset
        dataset_path = Path(from_dataset)
        metadata_file = dataset_path / "metadata" / f"{dataset_path.name}_metadata.json"

        if not metadata_file.exists():
            console.print(f"[red]Metadata file not found: {metadata_file}[/red]")
            return

        try:
            with open(metadata_file) as f:
                metadata = json.load(f)

            # Extract parameters from metadata
            config = {
                'name': name,
                'description': f"Preset created from {from_dataset}",
                'category': 'custom',
                'parameters': {
                    'collection_id': metadata.get('collection', {}).get('id'),
                    'platform': metadata.get('platform', {}).get('name'),
                    'coverage': metadata.get('configuration', {}).get('target_coverage'),
                }
            }

            # Add VLP protocol if present
            if metadata.get('configuration', {}).get('vlp_protocol'):
                config['parameters']['vlp_protocol'] = metadata['configuration']['vlp_protocol']

            # Validate and save
            is_valid, errors = validate_preset(config)
            if not is_valid:
                console.print("[red]Preset validation failed:[/red]")
                for error in errors:
                    console.print(f"  • {error}")
                return

            preset_file = create_preset_from_dict(name, config, user_preset=True)
            console.print(f"[green]Preset created: {preset_file}[/green]")

        except Exception as e:
            console.print(f"[red]Error creating preset: {e}[/red]")

    else:
        console.print("[yellow]Interactive preset creation not yet implemented.[/yellow]")
        console.print()
        console.print("Use --from-dataset to create from existing dataset:")
        console.print(f"  viroforge presets create {name} --from-dataset data/my_dataset")
        console.print()

    console.print()
