#!/usr/bin/env python3
"""
ViroForge Batch FASTQ Generator

Generate FASTQ datasets for multiple collections with different configurations.
Useful for creating comprehensive benchmark datasets.

Usage:
    # Generate standard benchmark suite
    python batch_generate_fastq.py \\
        --config configs/benchmark_suite.yaml \\
        --output data/benchmark_datasets

    # Generate quick test datasets
    python batch_generate_fastq.py \\
        --preset quick-test \\
        --output data/test_datasets

Author: ViroForge Development Team
Date: 2025-11-01
"""

import argparse
import sys
import yaml
import logging
import subprocess
from pathlib import Path
from typing import Dict, List
from datetime import datetime
import json

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Preset configurations
PRESETS = {
    'quick-test': {
        'description': 'Quick test with small collections and low coverage',
        'collections': [16, 11, 10],  # Mouse gut, Skin, Oral (smallest collections)
        'coverage': 1,
        'platform': 'miseq',
        'vlp_protocol': 'tangential_flow',
        'contamination_level': 'clean'
    },
    'benchmark-standard': {
        'description': 'Standard benchmark suite - all collections, 10x coverage',
        'collections': [9, 10, 11, 12, 13, 14, 15, 16],  # All collections
        'coverage': 10,
        'platform': 'novaseq',
        'vlp_protocol': 'tangential_flow',
        'contamination_level': 'realistic'
    },
    'vlp-comparison': {
        'description': 'VLP vs bulk comparison for selected collections',
        'collections': [9, 13],  # Gut and Marine
        'coverage': 10,
        'platforms': ['novaseq'],
        'vlp_protocol': 'tangential_flow',
        'contamination_level': 'realistic',
        'generate_both': True  # Both VLP and non-VLP
    },
    'vlp-protocol-comparison': {
        'description': 'Compare different VLP protocols on gut virome',
        'collections': [9],  # Gut
        'coverage': 10,
        'platform': 'novaseq',
        'contamination_level': 'realistic',
        'vlp_protocols': ['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen', None]
    },
    'platform-comparison': {
        'description': 'Cross-platform comparison',
        'collections': [9],  # Gut
        'coverage': 10,
        'platforms': ['novaseq', 'miseq', 'hiseq'],
        'vlp_protocol': 'tangential_flow',
        'contamination_level': 'realistic'
    },
    'coverage-series': {
        'description': 'Coverage depth series for gut virome',
        'collections': [9],  # Gut
        'coverages': [1, 5, 10, 20, 50],
        'platform': 'novaseq',
        'vlp_protocol': 'tangential_flow',
        'contamination_level': 'realistic'
    }
}


class BatchGenerator:
    """Batch FASTQ generation manager."""

    def __init__(self, output_dir: Path, dry_run: bool = False):
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run
        self.generation_script = Path(__file__).parent / 'generate_fastq_dataset.py'

        if not self.generation_script.exists():
            raise FileNotFoundError(f"Generation script not found: {self.generation_script}")

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Track results
        self.results = []
        self.start_time = datetime.now()

    def generate_dataset(
        self,
        collection_id: int,
        name: str,
        coverage: float,
        platform: str = 'novaseq',
        vlp_protocol: str = None,
        contamination_level: str = 'realistic',
        seed: int = 42
    ) -> Dict:
        """Generate single dataset."""
        # Create output path
        protocol_suffix = vlp_protocol if vlp_protocol else 'bulk'
        dataset_name = f"{name}_cov{coverage}x_{platform}_{protocol_suffix}"
        output_path = self.output_dir / dataset_name

        # Build command
        cmd = [
            'python', str(self.generation_script),
            '--collection-id', str(collection_id),
            '--output', str(output_path),
            '--coverage', str(coverage),
            '--platform', platform,
            '--contamination-level', contamination_level,
            '--seed', str(seed)
        ]

        if vlp_protocol is None:
            cmd.append('--no-vlp')
        else:
            cmd.extend(['--vlp-protocol', vlp_protocol])

        if self.dry_run:
            cmd.append('--dry-run')

        logger.info(f"Generating: {dataset_name}")
        logger.info(f"  Collection: {collection_id}")
        logger.info(f"  Coverage: {coverage}x")
        logger.info(f"  Platform: {platform}")
        logger.info(f"  VLP protocol: {vlp_protocol if vlp_protocol else 'none (bulk)'}")
        logger.info(f"  Contamination: {contamination_level}")

        result = {
            'dataset_name': dataset_name,
            'collection_id': collection_id,
            'coverage': coverage,
            'platform': platform,
            'vlp_protocol': vlp_protocol,
            'contamination_level': contamination_level,
            'output_path': str(output_path),
            'command': ' '.join(cmd),
            'start_time': datetime.now().isoformat()
        }

        try:
            proc_result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            result['status'] = 'success'
            result['end_time'] = datetime.now().isoformat()
            logger.info(f"✓ Completed: {dataset_name}")

        except subprocess.CalledProcessError as e:
            result['status'] = 'failed'
            result['error'] = e.stderr
            result['end_time'] = datetime.now().isoformat()
            logger.error(f"✗ Failed: {dataset_name}")
            logger.error(f"  Error: {e.stderr}")

        self.results.append(result)
        return result

    def run_preset(self, preset_name: str):
        """Run a preset configuration."""
        if preset_name not in PRESETS:
            raise ValueError(f"Unknown preset: {preset_name}. Available: {list(PRESETS.keys())}")

        preset = PRESETS[preset_name]
        logger.info(f"Running preset: {preset_name}")
        logger.info(f"  Description: {preset['description']}")

        # Handle different preset types
        if 'generate_both' in preset and preset['generate_both']:
            # VLP comparison - generate both VLP and non-VLP
            for collection_id in preset['collections']:
                for vlp_protocol in [preset.get('vlp_protocol', 'tangential_flow'), None]:
                    self.generate_dataset(
                        collection_id=collection_id,
                        name=f"collection_{collection_id}",
                        coverage=preset['coverage'],
                        platform=preset.get('platform', 'novaseq'),
                        vlp_protocol=vlp_protocol,
                        contamination_level=preset.get('contamination_level', 'realistic')
                    )

        elif 'vlp_protocols' in preset:
            # VLP protocol comparison
            for collection_id in preset['collections']:
                for vlp_protocol in preset['vlp_protocols']:
                    self.generate_dataset(
                        collection_id=collection_id,
                        name=f"collection_{collection_id}",
                        coverage=preset['coverage'],
                        platform=preset.get('platform', 'novaseq'),
                        vlp_protocol=vlp_protocol,
                        contamination_level=preset.get('contamination_level', 'realistic')
                    )

        elif 'platforms' in preset:
            # Platform comparison
            for collection_id in preset['collections']:
                for platform in preset['platforms']:
                    self.generate_dataset(
                        collection_id=collection_id,
                        name=f"collection_{collection_id}",
                        coverage=preset['coverage'],
                        platform=platform,
                        vlp_protocol=preset.get('vlp_protocol', 'tangential_flow'),
                        contamination_level=preset.get('contamination_level', 'realistic')
                    )

        elif 'coverages' in preset:
            # Coverage series
            for collection_id in preset['collections']:
                for coverage in preset['coverages']:
                    self.generate_dataset(
                        collection_id=collection_id,
                        name=f"collection_{collection_id}",
                        coverage=coverage,
                        platform=preset['platform'],
                        vlp_protocol=preset.get('vlp_protocol', 'tangential_flow'),
                        contamination_level=preset.get('contamination_level', 'realistic')
                    )

        else:
            # Standard batch generation
            for collection_id in preset['collections']:
                self.generate_dataset(
                    collection_id=collection_id,
                    name=f"collection_{collection_id}",
                    coverage=preset['coverage'],
                    platform=preset.get('platform', 'novaseq'),
                    vlp_protocol=preset.get('vlp_protocol', 'tangential_flow'),
                    contamination_level=preset.get('contamination_level', 'realistic')
                )

    def run_config(self, config_path: Path):
        """Run from configuration file."""
        with open(config_path) as f:
            config = yaml.safe_load(f)

        logger.info(f"Running configuration: {config_path}")

        for dataset_config in config.get('datasets', []):
            self.generate_dataset(**dataset_config)

    def export_summary(self):
        """Export generation summary."""
        summary = {
            'batch_info': {
                'start_time': self.start_time.isoformat(),
                'end_time': datetime.now().isoformat(),
                'output_directory': str(self.output_dir),
                'dry_run': self.dry_run
            },
            'results': self.results,
            'summary_stats': {
                'total_datasets': len(self.results),
                'successful': sum(1 for r in self.results if r['status'] == 'success'),
                'failed': sum(1 for r in self.results if r['status'] == 'failed')
            }
        }

        summary_path = self.output_dir / 'batch_generation_summary.json'
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"\nBatch Summary:")
        logger.info(f"  Total datasets: {summary['summary_stats']['total_datasets']}")
        logger.info(f"  Successful: {summary['summary_stats']['successful']}")
        logger.info(f"  Failed: {summary['summary_stats']['failed']}")
        logger.info(f"  Summary saved to: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge Batch FASTQ Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Presets:
{chr(10).join(f"  {name}: {p['description']}" for name, p in PRESETS.items())}

Examples:
  # Quick test with small datasets
  python batch_generate_fastq.py --preset quick-test --output data/test

  # Full benchmark suite
  python batch_generate_fastq.py --preset benchmark-standard --output data/benchmark

  # From config file
  python batch_generate_fastq.py --config my_config.yaml --output data/custom
        """
    )

    parser.add_argument(
        '--preset',
        choices=list(PRESETS.keys()),
        help='Use a preset configuration'
    )

    parser.add_argument(
        '--config',
        type=Path,
        help='Path to YAML configuration file'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output directory for all datasets'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Dry run - show what would be generated'
    )

    args = parser.parse_args()

    if not args.preset and not args.config:
        parser.error("Either --preset or --config required")

    # Create batch generator
    generator = BatchGenerator(args.output, dry_run=args.dry_run)

    try:
        if args.preset:
            generator.run_preset(args.preset)
        elif args.config:
            generator.run_config(args.config)

        # Export summary
        generator.export_summary()

        logger.info("\n✓ Batch generation complete!")

    except Exception as e:
        logger.error(f"Batch generation encountered error: {e}")
        # Still export partial summary
        logger.info("Exporting partial results summary...")
        try:
            generator.export_summary()
        except Exception as summary_error:
            logger.error(f"Failed to export summary: {summary_error}")
        sys.exit(1)


if __name__ == '__main__':
    main()
