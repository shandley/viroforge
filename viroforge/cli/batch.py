#!/usr/bin/env python3
"""
Batch Command

Generate multiple datasets from YAML configuration.

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import yaml
from pathlib import Path
from typing import Dict, List
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import json

try:
    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
    from rich.table import Table
    from rich.panel import Panel
    from rich.prompt import Confirm
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

from .generate import execute_generation

console = Console()


def run_batch(args):
    """Run batch generation from YAML configuration."""
    config_path = Path(args.config)

    if not config_path.exists():
        console.print(f"[red]Error: Config file not found: {config_path}[/red]")
        return 1

    # Load batch configuration
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        console.print(f"[red]Error loading config: {e}[/red]")
        return 1

    # Validate configuration
    if not validate_batch_config(config):
        return 1

    # Show batch summary
    show_batch_summary(config, args)

    # Confirm unless dry-run
    if args.dry_run:
        console.print("[yellow]DRY RUN - No datasets will be generated[/yellow]")
        return 0

    if not Confirm.ask("[bold]Proceed with batch generation?[/bold]", default=True):
        console.print("[yellow]Cancelled[/yellow]")
        return 0

    # Execute batch
    return execute_batch(config, args)


def validate_batch_config(config: Dict) -> bool:
    """Validate batch configuration."""
    errors = []

    if 'batch_name' not in config:
        errors.append("Missing 'batch_name' field")

    if 'datasets' not in config and 'parameter_sweep' not in config:
        errors.append("Must have either 'datasets' or 'parameter_sweep' field")

    if errors:
        console.print("[red]Configuration errors:[/red]")
        for error in errors:
            console.print(f"  • {error}")
        console.print()
        return False

    return True


def show_batch_summary(config: Dict, args):
    """Display batch configuration summary."""
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(f" Batch: {config.get('batch_name', 'Unnamed')}", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Count datasets
    datasets = config.get('datasets', [])
    sweep_datasets = []

    if 'parameter_sweep' in config:
        sweep_datasets = expand_parameter_sweep(config['parameter_sweep'])

    total = len(datasets) + len(sweep_datasets)

    console.print(f"[bold]Total datasets:[/bold] {total}")
    console.print(f"  • Individual datasets: {len(datasets)}")
    console.print(f"  • Parameter sweep: {len(sweep_datasets)}")
    console.print()

    if args.parallel > 1:
        console.print(f"[bold]Parallel execution:[/bold] {args.parallel} concurrent")
    else:
        console.print(f"[bold]Execution:[/bold] Sequential")

    console.print()

    # Estimate
    avg_time = 5  # minutes per dataset (rough estimate)
    if args.parallel > 1:
        total_time = (total / args.parallel) * avg_time
    else:
        total_time = total * avg_time

    console.print(f"[dim]Estimated time: ~{int(total_time)} minutes[/dim]")
    console.print()


def expand_parameter_sweep(sweep_config: Dict) -> List[Dict]:
    """Expand parameter sweep into individual dataset configs."""
    base_config = sweep_config.get('base_config', {})
    sweep_params = sweep_config.get('sweep_parameters', {})
    name_template = sweep_config.get('name_template', 'sweep_{index}')

    datasets = []

    # Get all parameter values to sweep
    param_names = list(sweep_params.keys())
    param_values = list(sweep_params.values())

    # Generate all combinations
    import itertools
    for combo in itertools.product(*param_values):
        dataset = base_config.copy()

        # Apply sweep values
        for param_name, value in zip(param_names, combo):
            dataset[param_name] = value

        # Generate name
        name_parts = {param: value for param, value in zip(param_names, combo)}
        dataset['name'] = name_template.format(**name_parts, index=len(datasets))

        datasets.append(dataset)

    return datasets


def execute_batch(config: Dict, args) -> int:
    """Execute batch generation."""
    # Collect all datasets
    datasets = config.get('datasets', []).copy()

    if 'parameter_sweep' in config:
        sweep_datasets = expand_parameter_sweep(config['parameter_sweep'])
        datasets.extend(sweep_datasets)

    # Setup output directory
    output_base = Path(config.get('output_base', 'data/batch'))
    output_base.mkdir(parents=True, exist_ok=True)

    # Track results
    results = {
        'batch_name': config.get('batch_name'),
        'start_time': datetime.now().isoformat(),
        'datasets': [],
        'successful': 0,
        'failed': 0
    }

    console.print("[bold green]Starting batch generation...[/bold green]")
    console.print()

    # Execute based on parallel setting
    if args.parallel > 1:
        exit_code = execute_parallel(datasets, output_base, args.parallel, results)
    else:
        exit_code = execute_sequential(datasets, output_base, results)

    # Save batch report
    results['end_time'] = datetime.now().isoformat()
    report_path = output_base / 'batch_report.json'
    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)

    # Show summary
    show_batch_results(results, output_base)

    return exit_code


def execute_sequential(datasets: List[Dict], output_base: Path, results: Dict) -> int:
    """Execute datasets sequentially with progress tracking."""
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:

        batch_task = progress.add_task(
            "[cyan]Batch progress...",
            total=len(datasets)
        )

        for i, dataset in enumerate(datasets, 1):
            name = dataset.get('name', f'dataset_{i}')

            progress.update(
                batch_task,
                description=f"[cyan]Generating {name} ({i}/{len(datasets)})..."
            )

            # Prepare parameters
            params = prepare_dataset_params(dataset, output_base)

            # Execute
            try:
                return_code = execute_generation(params, verbose=False)

                if return_code == 0:
                    results['successful'] += 1
                    status = 'success'
                else:
                    results['failed'] += 1
                    status = 'failed'

            except Exception as e:
                console.print(f"[red]Error generating {name}: {e}[/red]")
                results['failed'] += 1
                status = 'error'

            results['datasets'].append({
                'name': name,
                'status': status,
                'output': str(params.get('output'))
            })

            progress.update(batch_task, advance=1)

    return 0 if results['failed'] == 0 else 1


def execute_parallel(datasets: List[Dict], output_base: Path, max_workers: int, results: Dict) -> int:
    """Execute datasets in parallel."""
    console.print(f"[cyan]Running {max_workers} datasets in parallel...[/cyan]")
    console.print()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:

        batch_task = progress.add_task(
            "[cyan]Batch progress...",
            total=len(datasets)
        )

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_dataset = {}
            for i, dataset in enumerate(datasets, 1):
                params = prepare_dataset_params(dataset, output_base)
                future = executor.submit(execute_generation, params, False)
                future_to_dataset[future] = (dataset.get('name', f'dataset_{i}'), params)

            # Process completions
            for future in as_completed(future_to_dataset):
                name, params = future_to_dataset[future]

                try:
                    return_code = future.result()
                    if return_code == 0:
                        results['successful'] += 1
                        status = 'success'
                        console.print(f"[green]✓ {name} complete[/green]")
                    else:
                        results['failed'] += 1
                        status = 'failed'
                        console.print(f"[red]✗ {name} failed[/red]")

                except Exception as e:
                    results['failed'] += 1
                    status = 'error'
                    console.print(f"[red]✗ {name} error: {e}[/red]")

                results['datasets'].append({
                    'name': name,
                    'status': status,
                    'output': str(params.get('output'))
                })

                progress.update(batch_task, advance=1)

    return 0 if results['failed'] == 0 else 1


def prepare_dataset_params(dataset: Dict, output_base: Path) -> Dict:
    """Prepare parameters for dataset generation."""
    params = dataset.copy()

    # Set output path
    name = params.pop('name', 'unnamed')
    if 'output' not in params:
        params['output'] = str(output_base / name)

    return params


def show_batch_results(results: Dict, output_base: Path):
    """Display batch generation results."""
    console.print()
    console.print("═" * 80, style="bold cyan")
    console.print(" Batch Generation Complete", style="bold cyan")
    console.print("═" * 80, style="bold cyan")
    console.print()

    # Summary table
    table = Table(show_header=False, box=None, padding=(0, 2))
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="bold")

    table.add_row("Batch name", results['batch_name'])
    table.add_row("Total datasets", str(len(results['datasets'])))
    table.add_row("Successful", f"[green]{results['successful']}[/green]")
    table.add_row("Failed", f"[red]{results['failed']}[/red]")

    console.print(table)
    console.print()

    # Dataset list
    if results['datasets']:
        console.print("[bold]Datasets:[/bold]")
        for dataset in results['datasets']:
            status_icon = "✓" if dataset['status'] == 'success' else "✗"
            status_color = "green" if dataset['status'] == 'success' else "red"
            console.print(f"  [{status_color}]{status_icon}[/{status_color}] {dataset['name']}")
        console.print()

    console.print(f"Output: {output_base}")
    console.print(f"Report: {output_base / 'batch_report.json'}")
    console.print()
