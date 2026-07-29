#!/usr/bin/env python3
"""
Generate Command

Generate synthetic virome datasets using presets or custom parameters.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
from typing import Dict, List, Optional

try:
    from rich.console import Console
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

    # Contamination and VLP
    if getattr(args, 'contamination_level', None):
        params['contamination_level'] = args.contamination_level
    if getattr(args, 'vlp_protocol', None):
        params['vlp_protocol'] = args.vlp_protocol
    if getattr(args, 'pore_size', None):
        params['pore_size'] = args.pore_size
    if getattr(args, 'no_vlp', False):
        params['no_vlp'] = True

    # Molecule type and RNA options
    if getattr(args, 'molecule_type', None):
        params['molecule_type'] = args.molecule_type
    if getattr(args, 'rna_depletion', None):
        params['rna_depletion'] = args.rna_depletion

    # Artifact injection. Use "is not None" so an explicit 0 disables the
    # artifact (the underlying script now defaults these to non-zero rates);
    # a truthy check would swallow 0 and leave the default in place.
    if getattr(args, 'adapter_rate', None) is not None:
        params['adapter_rate'] = args.adapter_rate
    if getattr(args, 'low_complexity_rate', None) is not None:
        params['low_complexity_rate'] = args.low_complexity_rate
    if getattr(args, 'duplicate_rate', None) is not None:
        params['duplicate_rate'] = args.duplicate_rate

    # ERV injection
    if getattr(args, 'erv_endogenous_rate', None):
        params['erv_endogenous_rate'] = args.erv_endogenous_rate
    if getattr(args, 'erv_exogenous_rate', None):
        params['erv_exogenous_rate'] = args.erv_exogenous_rate

    # Dark matter
    dm_frac = getattr(args, 'dark_matter_fraction', None)
    if dm_frac is not None:
        params['dark_matter_fraction'] = dm_frac

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
    Execute dataset generation by calling run_generation() directly.

    Converts the params dict into an argparse-style namespace and calls the
    generation engine in-process. Generation runs without a progress bar: the
    engine logs its own progress, and the previous bar was driven by parsing
    the output of a subprocess that no longer exists.

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
    import argparse
    import logging

    # Build an argparse-compatible namespace from the params dict
    args = _params_to_namespace(params)

    if verbose:
        logging.basicConfig(level=logging.DEBUG)
        console.print(f"[dim]Parameters: {params}[/dim]")
        console.print()

    console.print("[bold green]Starting generation...[/bold green]")
    console.print()

    try:
        # Import and call run_generation directly (no subprocess)
        from viroforge.generator import run_generation

        result = run_generation(args)

        if result is None or result == 0:
            console.print()
            console.print("[bold green]✓ Dataset generated successfully![/bold green]")
            console.print()
            console.print(f"Output: {params.get('output', 'data/')}")
            console.print()
            return 0
        else:
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
        if verbose:
            import traceback
            traceback.print_exc()
        return 1


def _params_to_namespace(params: Dict):
    """Convert a params dict to the Namespace run_generation() expects.

    Defaults come from the real parser via parse_args([]), not a second copy
    of them. An earlier version restated all 41 defaults here; it had drifted
    to 6 missing attributes and 11 wrong values, so `viroforge generate`
    raised AttributeError before reaching the generator.
    """
    from viroforge.generator import build_parser

    args = build_parser().parse_args([])

    unknown = sorted(k for k in params if not hasattr(args, k))
    if unknown:
        raise ValueError(
            f"unknown generation parameter(s): {', '.join(unknown)}. "
            "Add the argument to build_parser() in viroforge/generator.py."
        )

    for key, value in params.items():
        setattr(args, key, value)

    if params.get("no_vlp"):
        args.vlp_protocol = "none"

    return args


