"""
Summary Command

Compute expected QC metric ranges from a batch of generated datasets.
Reads metadata JSON files and outputs YAML ranges for virome-qc profiles.

Usage:
    viroforge summary data/virome_qc_references/stool_vlp_tagmentation/
    viroforge summary data/virome_qc_references/stool_vlp_tagmentation/ --format yaml
    viroforge summary data/virome_qc_references/stool_vlp_tagmentation/ --output ranges.yaml
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    from rich.console import Console
    from rich.table import Table
except ImportError:
    print("Error: 'rich' library required. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()


def _find_metadata_files(dataset_dirs: List[Path]) -> List[Path]:
    """Find all metadata JSON files in dataset directories."""
    metadata_files = []
    for d in dataset_dirs:
        if d.is_file() and d.suffix == ".json":
            metadata_files.append(d)
        elif d.is_dir():
            # Search for metadata JSON in standard locations
            for pattern in ["metadata/*_metadata.json", "*_metadata.json", "metadata.json"]:
                metadata_files.extend(d.glob(pattern))
            # Also check subdirectories (batch output has named subdirs)
            for subdir in d.iterdir():
                if subdir.is_dir():
                    for pattern in ["metadata/*_metadata.json", "*_metadata.json"]:
                        metadata_files.extend(subdir.glob(pattern))
    return sorted(set(metadata_files))


def _extract_metrics(metadata: dict) -> dict:
    """Extract QC-relevant metrics from a ViroForge metadata JSON."""
    metrics = {}

    # Enrichment stats
    enrichment = metadata.get("enrichment_stats", {})
    if enrichment:
        metrics["viral_fraction"] = enrichment.get("viral_fraction")
        metrics["contamination_fraction"] = enrichment.get("contamination_fraction")
        n_viral = enrichment.get("n_viral_genomes", 0)
        n_contam = enrichment.get("n_contaminants", 0)
        total = n_viral + n_contam
        if total > 0:
            metrics["host_fraction"] = None  # computed below
            metrics["rrna_fraction"] = None

    # Configuration
    config = metadata.get("configuration", {})
    metrics["platform"] = config.get("platform")
    metrics["vlp_protocol"] = config.get("vlp_protocol")
    metrics["contamination_level"] = config.get("contamination_level")

    # Contamination manifest for per-type fractions
    benchmarking = metadata.get("benchmarking", {})
    manifest = benchmarking.get("contamination_manifest", {})
    contaminants = manifest.get("contaminants", [])

    host_total = 0.0
    rrna_total = 0.0
    phix_total = 0.0
    reagent_total = 0.0
    for c in contaminants:
        ctype = c.get("contaminant_type", "")
        abundance = c.get("abundance", 0.0)
        if ctype == "host_dna":
            host_total += abundance
        elif ctype == "rrna":
            rrna_total += abundance
        elif ctype == "phix":
            phix_total += abundance
        elif ctype == "reagent_bacteria":
            reagent_total += abundance

    if contaminants:
        metrics["host_fraction"] = host_total
        metrics["rrna_fraction"] = rrna_total
        metrics["phix_fraction"] = phix_total
        metrics["reagent_fraction"] = reagent_total

    # Survival (viral fraction = 1 - total contamination)
    vf = metrics.get("viral_fraction")
    if vf is not None:
        metrics["survival"] = vf

    # Adapter stats
    adapter = metadata.get("adapter_stats", {})
    if adapter:
        metrics["adapter_rate"] = adapter.get(
            "emergent_adapter_rate", adapter.get("adapter_rate", 0.0)
        )
        metrics["mean_adapter_length"] = adapter.get("mean_adapter_length")

    # Duplicate stats
    dup = metadata.get("duplicate_stats", {})
    if dup:
        metrics["duplication_rate"] = dup.get("duplicate_fraction")

    # Low-complexity stats
    lc = metadata.get("low_complexity_stats", {})
    if lc:
        reads_mod = lc.get("reads_modified", 0)
        # Approximate rate from reads_modified / total
        # (total is from adapter_stats or other)
        total_reads = adapter.get("reads_total", dup.get("reads_original", 0))
        if total_reads > 0:
            metrics["low_complexity_rate"] = reads_mod / total_reads

    # Filter out None values
    return {k: v for k, v in metrics.items() if v is not None}


def _compute_ranges(
    all_metrics: List[dict], percentiles: tuple = (5, 25, 50, 75, 95)
) -> dict:
    """Compute per-metric ranges across multiple datasets."""
    import statistics

    # Collect numeric values per metric
    metric_values: Dict[str, List[float]] = {}
    for metrics in all_metrics:
        for key, value in metrics.items():
            if isinstance(value, (int, float)):
                metric_values.setdefault(key, []).append(float(value))

    ranges = {}
    for metric, values in sorted(metric_values.items()):
        if len(values) < 2:
            ranges[metric] = {
                "min": values[0],
                "max": values[0],
                "n": len(values),
            }
            continue

        values_sorted = sorted(values)
        n = len(values_sorted)
        ranges[metric] = {
            "min": values_sorted[0],
            "max": values_sorted[-1],
            "mean": statistics.mean(values),
            "median": statistics.median(values),
            "n": n,
        }
        # Percentiles
        for p in percentiles:
            idx = int(n * p / 100)
            idx = max(0, min(idx, n - 1))
            ranges[metric][f"p{p}"] = values_sorted[idx]

    return ranges


def _format_yaml_ranges(ranges: dict, profile_name: str) -> str:
    """Format ranges as YAML for virome-qc profile expected_ranges."""
    lines = [f"# Expected ranges for {profile_name}"]
    lines.append(f"# Derived from {ranges.get('survival', {}).get('n', '?')} ViroForge reference datasets")
    lines.append("expected_ranges:")

    # Key metrics for virome-qc
    key_metrics = [
        "survival", "host_fraction", "rrna_fraction", "phix_fraction",
        "adapter_rate", "duplication_rate", "low_complexity_rate",
    ]

    for metric in key_metrics:
        if metric in ranges:
            r = ranges[metric]
            lo = r.get("p5", r["min"])
            hi = r.get("p95", r["max"])
            lines.append(f"  {metric}: [{lo:.4f}, {hi:.4f}]")

    return "\n".join(lines)


def run_summary(args):
    """Compute expected QC metric ranges from generated datasets."""
    dataset_paths = [Path(p) for p in args.datasets]

    # Find metadata files
    metadata_files = _find_metadata_files(dataset_paths)

    if not metadata_files:
        console.print("[red]No metadata JSON files found in the specified paths.[/red]")
        console.print("Expected: metadata/*_metadata.json in each dataset directory.")
        return 1

    console.print(f"Found [bold]{len(metadata_files)}[/bold] metadata files")
    console.print()

    # Extract metrics from each dataset
    all_metrics = []
    for mf in metadata_files:
        try:
            with open(mf) as f:
                metadata = json.load(f)
            metrics = _extract_metrics(metadata)
            metrics["_source"] = str(mf)
            all_metrics.append(metrics)
        except (json.JSONDecodeError, KeyError) as e:
            console.print(f"[yellow]Warning: Could not parse {mf}: {e}[/yellow]")

    if not all_metrics:
        console.print("[red]No valid metadata found.[/red]")
        return 1

    # Compute ranges
    ranges = _compute_ranges(all_metrics)

    # Display table
    table = Table(title="QC Metric Ranges", show_header=True)
    table.add_column("Metric", style="cyan")
    table.add_column("Min", justify="right")
    table.add_column("p5", justify="right")
    table.add_column("Median", justify="right")
    table.add_column("p95", justify="right")
    table.add_column("Max", justify="right")
    table.add_column("N", justify="right")

    for metric in sorted(ranges.keys()):
        r = ranges[metric]
        if metric.startswith("_"):
            continue
        table.add_row(
            metric,
            f"{r['min']:.4f}",
            f"{r.get('p5', r['min']):.4f}",
            f"{r.get('median', r['min']):.4f}",
            f"{r.get('p95', r['max']):.4f}",
            f"{r['max']:.4f}",
            str(r["n"]),
        )

    console.print(table)
    console.print()

    # Output YAML
    profile_name = args.profile or "unknown"
    yaml_output = _format_yaml_ranges(ranges, profile_name)

    if args.format == "yaml" or args.output:
        if args.output:
            Path(args.output).write_text(yaml_output + "\n")
            console.print(f"Wrote ranges to [bold]{args.output}[/bold]")
        else:
            console.print(yaml_output)
    else:
        console.print("[bold]YAML output (use --format yaml or --output FILE):[/bold]")
        console.print(yaml_output)

    return 0
