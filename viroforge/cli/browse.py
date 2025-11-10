#!/usr/bin/env python3
"""
Interactive Collection Browser

Terminal-based interactive browser for exploring ViroForge collections.
Uses the rich library for beautiful terminal output.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
from typing import Optional, List, Dict

try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.layout import Layout
    from rich.text import Text
    from rich.prompt import Prompt, Confirm
    from rich import box
    from rich.columns import Columns
    from rich.markdown import Markdown
except ImportError:
    print("Error: 'rich' library is required for the interactive browser.", file=sys.stderr)
    print("Install it with: pip install rich", file=sys.stderr)
    sys.exit(1)

from .db_utils import (
    load_all_collections,
    load_collection_details,
    search_collections,
    get_collection_categories,
    get_database_stats
)


console = Console()


def format_collection_row(collection: Dict, icons: bool = True) -> str:
    """Format a collection as a table row."""
    icon = "" if icons else ""
    popular_tag = "" if collection['n_genomes'] > 50 else ""
    name = collection['name']
    n_genomes = collection['n_genomes']

    return f"{icon} [{collection['id']:2d}] {name:40s} {n_genomes:3d} genomes {popular_tag}"


def show_collection_list(categories: Dict[str, List[Dict]], search_query: str = "", icons: bool = True):
    """Display list of collections organized by category."""
    console.clear()

    # Header
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(" " * 20 + "ViroForge Collection Browser", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Search box
    if search_query:
        console.print(f"[yellow]Search:[/yellow] [bold]{search_query}[/bold]", style="dim")
    else:
        console.print("[dim]Type a collection ID to view details, or 'q' to quit[/dim]")
    console.print()

    # Database stats
    try:
        stats = get_database_stats()
        stats_text = f"[dim]Database: {stats['total_genomes']} genomes, {stats['total_collections']} collections, {stats['total_families']} families[/dim]"
        console.print(stats_text)
        console.print()
    except:
        pass

    # Collections by category
    for category_name, collections in categories.items():
        if not collections:
            continue

        # Category header
        console.print(f"─ {category_name} ({len(collections)} collections) ".ljust(80, "─"), style="bold blue")
        console.print()

        # Collection rows
        for collection in collections[:10]:  # Limit to 10 per category
            console.print(f"  {format_collection_row(collection, icons)}")

        if len(collections) > 10:
            console.print(f"  [dim]... and {len(collections) - 10} more[/dim]")

        console.print()


def show_collection_details(collection_id: int, icons: bool = True):
    """Show detailed view of a specific collection."""
    try:
        collection = load_collection_details(collection_id)
    except Exception as e:
        console.print(f"[red]Error loading collection {collection_id}: {e}[/red]")
        console.input("\n[dim]Press Enter to continue...[/dim]")
        return

    console.clear()

    # Header
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(f" Collection {collection['collection_id']}: {collection['collection_name']}", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Overview
    console.print("[bold]Overview:[/bold]")
    console.print(f"  Genomes: [green]{collection['actual_genome_count']}[/green] viral genomes")
    if collection.get('curation_date'):
        console.print(f"  Curated: {collection['curation_date']}")
    console.print()

    # Description
    if collection.get('description'):
        console.print("[bold]Description:[/bold]")
        # Wrap text at 76 chars
        desc = collection['description']
        words = desc.split()
        line = "  "
        for word in words:
            if len(line) + len(word) + 1 > 78:
                console.print(line)
                line = "  " + word
            else:
                line += " " + word if len(line) > 2 else word
        if line.strip():
            console.print(line)
        console.print()

    # Genome types
    if collection.get('genome_types'):
        console.print("[bold]Genome Types:[/bold]")
        for genome_type, count in collection['genome_types'].items():
            percentage = (count / collection['actual_genome_count']) * 100
            bar = "█" * int(percentage / 5)  # Scale to ~20 chars
            console.print(f"  {genome_type:20s} {count:3d} ({percentage:5.1f}%) {bar}")
        console.print()

    # Top families
    if collection.get('top_families'):
        console.print("[bold]Top Viral Families:[/bold]")
        for family in collection['top_families'][:5]:
            console.print(f"  • {family['family']:30s} {family['count']:3d} genomes")
        console.print()

    # Top genomes
    if collection.get('genomes'):
        console.print("[bold]Most Abundant Genomes:[/bold]")
        for genome in collection['genomes'][:5]:
            name = genome['genome_name'][:50]
            abundance = genome.get('relative_abundance', 0) or 0
            if abundance > 0:
                bar = "█" * int(abundance * 40)  # Scale for display
                console.print(f"  {name:50s} {abundance:6.2%} {bar}")
            else:
                console.print(f"  {name:50s}")

        if len(collection['genomes']) > 5:
            console.print(f"  [dim]... and {len(collection['genomes']) - 5} more genomes[/dim]")
        console.print()

    # Selection criteria
    if collection.get('selection_criteria'):
        console.print("[bold]Selection Criteria:[/bold]")
        console.print(f"  {collection['selection_criteria']}")
        console.print()

    # Literature
    if collection.get('literature_references'):
        console.print("[bold]Literature:[/bold]")
        refs = collection['literature_references']
        for ref in refs.split(';')[:3]:
            console.print(f"  • {ref.strip()}")
        console.print()

    # Quick actions
    console.print("─" * 80, style="dim")
    console.print()
    console.print("[bold cyan]Quick Actions:[/bold cyan]")
    console.print("  [G] Generate dataset with this collection")
    console.print("  [E] Export genome list")
    console.print("  [B] Back to browser")
    console.print()


def run_browser(args):
    """Main browser loop."""
    icons = not args.no_icons
    search_query = ""

    while True:
        try:
            # Get collections (filtered if searching)
            if search_query:
                collections = search_collections(search_query)
                categories = {'Search Results': collections}
            else:
                categories = get_collection_categories()

            # Show collection list
            show_collection_list(categories, search_query, icons)

            # Get user input
            console.print()
            choice = Prompt.ask(
                "[bold]Enter collection ID, 's' to search, or 'q' to quit[/bold]",
                default="q"
            )

            if choice.lower() == 'q':
                console.print("\n[green]Goodbye![/green]")
                return 0
            elif choice.lower() == 's':
                search_query = Prompt.ask("[bold]Search query[/bold]")
                continue
            elif choice.lower() == 'c':
                search_query = ""
                continue
            else:
                try:
                    collection_id = int(choice)
                    show_collection_details(collection_id, icons)

                    # Action menu
                    action = Prompt.ask(
                        "[bold]Action[/bold]",
                        choices=['g', 'e', 'b'],
                        default='b'
                    )

                    if action == 'g':
                        # Generate dataset
                        console.print("\n[yellow]Launching generator...[/yellow]")
                        if generate_from_collection(collection_id):
                            console.print("[green]Generation complete![/green]")
                        console.input("\n[dim]Press Enter to continue...[/dim]")
                    elif action == 'e':
                        # Export genome list
                        export_genome_list(collection_id)
                        console.input("\n[dim]Press Enter to continue...[/dim]")
                    # 'b' goes back to list

                except ValueError:
                    console.print(f"[red]Invalid collection ID: {choice}[/red]")
                    console.input("\n[dim]Press Enter to continue...[/dim]")

        except KeyboardInterrupt:
            console.print("\n\n[green]Goodbye![/green]")
            return 0
        except Exception as e:
            console.print(f"\n[red]Error: {e}[/red]")
            console.input("\n[dim]Press Enter to continue...[/dim]")

    return 0


def generate_from_collection(collection_id: int) -> bool:
    """
    Launch dataset generation for a collection.

    Parameters
    ----------
    collection_id : int
        Collection ID to generate

    Returns
    -------
    bool
        True if generation was launched, False otherwise
    """
    import subprocess
    from pathlib import Path

    console.print()
    console.print("─" * 80, style="bold cyan")
    console.print(" Generate Dataset", style="bold cyan")
    console.print("─" * 80, style="bold cyan")
    console.print()

    # Get parameters from user
    output = Prompt.ask(
        "[bold]Output directory[/bold]",
        default=f"data/collection_{collection_id}"
    )

    platform = Prompt.ask(
        "[bold]Platform[/bold]",
        choices=['novaseq', 'miseq', 'hiseq', 'pacbio-hifi', 'nanopore'],
        default='novaseq'
    )

    if platform in ['novaseq', 'miseq', 'hiseq']:
        coverage = Prompt.ask(
            "[bold]Coverage[/bold]",
            default="30"
        )
        param_flag = '--coverage'
        param_value = coverage
    else:
        depth = Prompt.ask(
            "[bold]Depth[/bold]",
            default="15"
        )
        param_flag = '--depth'
        param_value = depth

    seed = Prompt.ask(
        "[bold]Random seed[/bold]",
        default="42"
    )

    # Confirm
    console.print()
    console.print("[yellow]Generation parameters:[/yellow]")
    console.print(f"  Collection ID: {collection_id}")
    console.print(f"  Output: {output}")
    console.print(f"  Platform: {platform}")
    console.print(f"  {param_flag.replace('--', '').title()}: {param_value}")
    console.print(f"  Seed: {seed}")
    console.print()

    if not Confirm.ask("[bold]Proceed with generation?[/bold]", default=True):
        console.print("[yellow]Cancelled[/yellow]")
        return False

    # Find generate_fastq_dataset.py script
    script_dir = Path(__file__).parent.parent.parent / "scripts"
    script_path = script_dir / "generate_fastq_dataset.py"

    if not script_path.exists():
        console.print(f"[red]Error: Could not find {script_path}[/red]")
        return False

    # Build command
    cmd = [
        'python3',
        str(script_path),
        '--collection-id', str(collection_id),
        '--output', output,
        '--platform', platform,
        param_flag, param_value,
        '--seed', seed
    ]

    # Execute
    console.print()
    console.print("[cyan]Running generation...[/cyan]")
    console.print(f"[dim]Command: {' '.join(cmd)}[/dim]")
    console.print()

    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Generation failed: {e}[/red]")
        return False
    except KeyboardInterrupt:
        console.print("\n[yellow]Generation cancelled[/yellow]")
        return False


def export_genome_list(collection_id: int):
    """Export genome list for a collection."""
    try:
        collection = load_collection_details(collection_id)
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        return

    output_file = Prompt.ask(
        "[bold]Output file[/bold]",
        default=f"collection_{collection_id}_genomes.tsv"
    )

    try:
        with open(output_file, 'w') as f:
            # Header
            f.write("genome_id\tgenome_name\tfamily\tgenus\tspecies\trelative_abundance\n")

            # Genomes
            for genome in collection['genomes']:
                f.write(f"{genome['genome_id']}\t")
                f.write(f"{genome['genome_name']}\t")
                f.write(f"{genome.get('family', 'N/A')}\t")
                f.write(f"{genome.get('genus', 'N/A')}\t")
                f.write(f"{genome.get('species', 'N/A')}\t")
                f.write(f"{genome.get('relative_abundance', 0.0)}\n")

        console.print(f"[green]Exported {len(collection['genomes'])} genomes to {output_file}[/green]")

    except Exception as e:
        console.print(f"[red]Export failed: {e}[/red]")
