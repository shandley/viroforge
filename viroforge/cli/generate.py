#!/usr/bin/env python3
"""
Generate Command

Generate synthetic virome datasets using presets or custom parameters.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

try:
    from rich.console import Console
    from rich.progress import (
        Progress,
        SpinnerColumn,
        TextColumn,
        BarColumn,
        TaskProgressColumn,
        TimeRemainingColumn,
        TimeElapsedColumn
    )
    from rich.panel import Panel
    from rich.table import Table
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

from .preset_loader import load_preset, get_preset_names

console = Console()


def run_generate(args):
    """Generate datasets using presets or custom parameters."""

    # If using preset
    if args.preset:
        return generate_with_preset(args)
    elif args.collection_id and args.output:
        return generate_with_params(args)
    else:
        console.print("[red]Error: Either --preset or (--collection-id and --output) required[/red]")
        console.print()
        console.print("Examples:")
        console.print("  viroforge generate --preset gut-standard")
        console.print("  viroforge generate --collection-id 9 --output data/gut --platform novaseq")
        console.print()
        console.print("Available presets:")
        for name in get_preset_names():
            console.print(f"  • {name}")
        console.print()
        console.print("Or use the interactive browser:")
        console.print("  viroforge browse")
        console.print()
        return 1


def generate_with_preset(args):
    """Generate dataset using a preset."""
    preset = load_preset(args.preset)

    if not preset:
        console.print(f"[red]Error: Preset '{args.preset}' not found[/red]")
        console.print()
        console.print("Available presets:")
        for name in get_preset_names():
            console.print(f"  • {name}")
        console.print()
        return 1

    # Show preset info
    console.print()
    console.print(f"[bold cyan]Using preset: {preset['name']}[/bold cyan]")
    if preset.get('description'):
        console.print(f"[dim]{preset['description']}[/dim]")
    console.print()

    # Build parameters from preset and overrides
    params = preset['parameters'].copy()

    # Apply command-line overrides
    if args.output:
        params['output'] = args.output
    elif 'output' not in params:
        # Default output based on preset name
        params['output'] = f"data/{args.preset}"

    if args.seed is not None:
        params['seed'] = args.seed

    if args.coverage is not None:
        params['coverage'] = args.coverage

    if args.depth is not None:
        params['depth'] = args.depth

    # Show final parameters
    show_parameters(params)

    # Execute generation
    return execute_generation(params, args.verbose)


def generate_with_params(args):
    """Generate dataset using direct parameters."""
    params = {
        'collection_id': args.collection_id,
        'output': args.output,
        'platform': args.platform or 'novaseq',
    }

    # Add optional parameters
    if args.coverage:
        params['coverage'] = args.coverage
    if args.depth:
        params['depth'] = args.depth
    if args.seed is not None:
        params['seed'] = args.seed

    console.print()
    console.print("[bold cyan]Generating dataset with custom parameters[/bold cyan]")
    console.print()

    show_parameters(params)

    return execute_generation(params, args.verbose)


def show_parameters(params: Dict):
    """Display generation parameters in a table."""
    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="bold")

    for key, value in sorted(params.items()):
        table.add_row(key.replace('_', ' ').title(), str(value))

    console.print(table)
    console.print()


def execute_generation(params: Dict, verbose: bool = False):
    """
    Execute dataset generation with progress reporting.

    Parameters
    ----------
    params : Dict
        Generation parameters
    verbose : bool
        Show detailed progress

    Returns
    -------
    int
        Exit code (0 = success, 1 = failure)
    """
    # Find generate_fastq_dataset.py script
    script_dir = Path(__file__).parent.parent.parent / "scripts"
    script_path = script_dir / "generate_fastq_dataset.py"

    if not script_path.exists():
        console.print(f"[red]Error: Could not find {script_path}[/red]")
        return 1

    # Build command
    cmd = build_command(script_path, params)

    if verbose:
        console.print(f"[dim]Command: {' '.join(str(c) for c in cmd)}[/dim]")
        console.print()

    # Execute with progress monitoring
    console.print("[bold green]Starting generation...[/bold green]")
    console.print()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=False
    ) as progress:

        # Add main task
        task = progress.add_task("[cyan]Generating dataset...", total=None)

        try:
            # Run the command
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )

            # Monitor output
            for line in process.stdout:
                line = line.rstrip()

                if verbose and line:
                    console.print(f"[dim]{line}[/dim]")

                # Update progress based on log messages
                if "Loading collection" in line:
                    progress.update(task, description="[cyan]Loading collection...")
                elif "VLP enrichment" in line:
                    progress.update(task, description="[cyan]Applying VLP enrichment...")
                elif "contamination" in line:
                    progress.update(task, description="[cyan]Adding contamination...")
                elif "Writing FASTA" in line:
                    progress.update(task, description="[cyan]Writing FASTA files...")
                elif "Generating" in line or "simulating" in line.lower():
                    progress.update(task, description="[cyan]Generating FASTQ reads...")
                elif "Writing metadata" in line:
                    progress.update(task, description="[cyan]Writing metadata...")
                elif "Complete" in line or "✓" in line:
                    progress.update(task, description="[green]✓ Generation complete")

            # Wait for completion
            return_code = process.wait()

            if return_code == 0:
                progress.update(task, description="[green]✓ Generation complete", completed=True)
                console.print()
                console.print("[bold green]✓ Dataset generated successfully![/bold green]")
                console.print()
                console.print(f"Output: {params.get('output', 'data/')}")
                console.print()
                return 0
            else:
                progress.update(task, description="[red]✗ Generation failed", completed=True)
                console.print()
                console.print("[bold red]✗ Generation failed[/bold red]")
                return 1

        except KeyboardInterrupt:
            console.print()
            console.print("[yellow]Generation cancelled by user[/yellow]")
            return 130
        except Exception as e:
            console.print()
            console.print(f"[red]Error: {e}[/red]")
            return 1


def build_command(script_path: Path, params: Dict) -> List[str]:
    """
    Build command-line arguments from parameters.

    Parameters
    ----------
    script_path : Path
        Path to generate_fastq_dataset.py
    params : Dict
        Generation parameters

    Returns
    -------
    List[str]
        Command and arguments
    """
    cmd = ['python3', str(script_path)]

    # Required parameters
    if 'collection_id' in params:
        cmd.extend(['--collection-id', str(params['collection_id'])])

    if 'output' in params:
        cmd.extend(['--output', str(params['output'])])

    # Platform
    if 'platform' in params:
        cmd.extend(['--platform', params['platform']])

    # Coverage/depth
    platform = params.get('platform', 'novaseq')
    if platform in ['novaseq', 'miseq', 'hiseq']:
        if 'coverage' in params:
            cmd.extend(['--coverage', str(params['coverage'])])
    else:  # Long-read platforms
        if 'depth' in params:
            cmd.extend(['--depth', str(params['depth'])])

    # Read parameters
    if 'read_length' in params:
        cmd.extend(['--read-length', str(params['read_length'])])

    if 'insert_size' in params:
        cmd.extend(['--insert-size', str(params['insert_size'])])

    # VLP and contamination
    if params.get('no_vlp'):
        cmd.append('--no-vlp')
    elif 'vlp_protocol' in params:
        cmd.extend(['--vlp-protocol', params['vlp_protocol']])

    if 'contamination_level' in params:
        cmd.extend(['--contamination-level', params['contamination_level']])

    # Amplification
    if 'amplification' in params:
        cmd.extend(['--amplification', params['amplification']])

    # Molecule type
    if 'molecule_type' in params:
        cmd.extend(['--molecule-type', params['molecule_type']])

    # RNA-specific
    if 'rna_depletion' in params:
        cmd.extend(['--rna-depletion', params['rna_depletion']])

    # Long-read specific
    if 'pacbio_passes' in params:
        cmd.extend(['--pacbio-passes', str(params['pacbio_passes'])])

    if 'pacbio_read_length' in params:
        cmd.extend(['--pacbio-read-length', str(params['pacbio_read_length'])])

    if 'ont_chemistry' in params:
        cmd.extend(['--ont-chemistry', params['ont_chemistry']])

    if 'ont_read_length' in params:
        cmd.extend(['--ont-read-length', str(params['ont_read_length'])])

    # Random seed
    if 'seed' in params:
        cmd.extend(['--seed', str(params['seed'])])

    return cmd
