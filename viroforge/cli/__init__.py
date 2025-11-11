#!/usr/bin/env python3
"""
ViroForge CLI

Command-line interface for ViroForge synthetic virome data generator.
Provides interactive tools for browsing collections, generating datasets,
and analyzing results.

Usage:
    viroforge browse              # Interactive collection browser
    viroforge generate            # Generate datasets
    viroforge report <dataset>    # View dataset reports
    viroforge compare <datasets>  # Compare multiple datasets
    viroforge batch <config>      # Batch generation
    viroforge presets             # Manage presets
    viroforge web                 # Launch web interface

Author: ViroForge Development Team
Date: 2025-11-10
"""

import sys
import argparse
from pathlib import Path


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog='viroforge',
        description='ViroForge: Synthetic virome data generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  viroforge browse                    Browse available collections
  viroforge generate --preset gut     Generate using preset
  viroforge report data/gut_virome    View dataset report
  viroforge compare data/gut_*        Compare multiple datasets
  viroforge batch config.yaml         Generate from batch config
  viroforge web                       Launch web interface

For more information: https://github.com/hecatomb/viroforge
        """
    )

    parser.add_argument(
        '--version',
        action='version',
        version='ViroForge 0.10.0 (Phase 12: CLI Enhancements)'
    )

    # Create subparsers for subcommands
    subparsers = parser.add_subparsers(
        title='commands',
        description='Available ViroForge commands',
        dest='command',
        help='Command to run'
    )

    # ========================================================================
    # BROWSE command
    # ========================================================================
    browse_parser = subparsers.add_parser(
        'browse',
        help='Browse available collections interactively',
        description='Interactive terminal browser for exploring ViroForge collections'
    )
    browse_parser.add_argument(
        '--no-icons',
        action='store_true',
        help='Disable Unicode icons (for compatibility)'
    )

    # ========================================================================
    # GENERATE command
    # ========================================================================
    generate_parser = subparsers.add_parser(
        'generate',
        help='Generate synthetic virome datasets',
        description='Generate datasets using presets or custom parameters'
    )
    generate_parser.add_argument(
        '--preset',
        help='Use a named preset (e.g., gut-standard, marine-standard)'
    )
    generate_parser.add_argument(
        '--collection-id',
        type=int,
        help='Collection ID to generate from'
    )
    generate_parser.add_argument(
        '--output',
        help='Output directory'
    )
    generate_parser.add_argument(
        '--platform',
        choices=['novaseq', 'miseq', 'hiseq', 'pacbio-hifi', 'nanopore'],
        help='Sequencing platform'
    )
    generate_parser.add_argument(
        '--coverage',
        type=float,
        help='Short-read coverage (e.g., 30)'
    )
    generate_parser.add_argument(
        '--depth',
        type=float,
        help='Long-read depth (e.g., 15)'
    )
    generate_parser.add_argument(
        '--seed',
        type=int,
        help='Random seed for reproducibility'
    )
    generate_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Show detailed progress'
    )

    # ========================================================================
    # REPORT command
    # ========================================================================
    report_parser = subparsers.add_parser(
        'report',
        help='Generate quality report for a dataset',
        description='View detailed statistics and quality metrics for generated datasets'
    )
    report_parser.add_argument(
        'dataset',
        help='Path to dataset directory'
    )
    report_parser.add_argument(
        '--format',
        choices=['terminal', 'json', 'html'],
        default='terminal',
        help='Output format (default: terminal)'
    )
    report_parser.add_argument(
        '--export',
        help='Export report to file'
    )

    # ========================================================================
    # COMPARE command
    # ========================================================================
    compare_parser = subparsers.add_parser(
        'compare',
        help='Compare multiple datasets',
        description='Side-by-side comparison of multiple datasets'
    )
    compare_parser.add_argument(
        'datasets',
        nargs='+',
        help='Paths to dataset directories'
    )
    compare_parser.add_argument(
        '--format',
        choices=['terminal', 'json', 'html'],
        default='terminal',
        help='Output format (default: terminal)'
    )
    compare_parser.add_argument(
        '--export',
        help='Export comparison to file'
    )

    # ========================================================================
    # BATCH command
    # ========================================================================
    batch_parser = subparsers.add_parser(
        'batch',
        help='Generate multiple datasets from YAML config',
        description='Batch generation from YAML configuration file'
    )
    batch_parser.add_argument(
        'config',
        help='Path to batch configuration YAML file'
    )
    batch_parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be generated without executing'
    )
    batch_parser.add_argument(
        '--parallel',
        type=int,
        default=1,
        help='Number of datasets to generate in parallel (default: 1)'
    )

    # ========================================================================
    # PRESETS command
    # ========================================================================
    presets_parser = subparsers.add_parser(
        'presets',
        help='Manage configuration presets',
        description='List, create, and manage configuration presets'
    )
    presets_subparsers = presets_parser.add_subparsers(
        dest='presets_command',
        help='Preset command'
    )

    # List presets
    presets_subparsers.add_parser(
        'list',
        help='List available presets'
    )

    # Show preset details
    show_parser = presets_subparsers.add_parser(
        'show',
        help='Show preset details'
    )
    show_parser.add_argument(
        'name',
        help='Preset name'
    )

    # Create preset
    create_parser = presets_subparsers.add_parser(
        'create',
        help='Create a new preset'
    )
    create_parser.add_argument(
        'name',
        help='Preset name'
    )
    create_parser.add_argument(
        '--from-dataset',
        help='Create preset from existing dataset metadata'
    )

    # ========================================================================
    # WEB command
    # ========================================================================
    web_parser = subparsers.add_parser(
        'web',
        help='Launch web interface',
        description='Start ViroForge web interface for browser-based interactions'
    )
    web_parser.add_argument(
        '--host',
        default='127.0.0.1',
        help='Host to bind to (default: 127.0.0.1)'
    )
    web_parser.add_argument(
        '--port',
        type=int,
        default=5000,
        help='Port to bind to (default: 5000)'
    )
    web_parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug mode'
    )
    web_parser.add_argument(
        '--no-browser',
        action='store_true',
        help='Don\'t open browser automatically'
    )

    # Parse arguments
    args = parser.parse_args()

    # If no command specified, show help
    if not args.command:
        parser.print_help()
        return 0

    # Route to appropriate subcommand
    try:
        if args.command == 'browse':
            from .browse import run_browser
            return run_browser(args)
        elif args.command == 'generate':
            from .generate import run_generate
            return run_generate(args)
        elif args.command == 'report':
            from .report import run_report
            return run_report(args)
        elif args.command == 'compare':
            from .compare import run_compare
            return run_compare(args)
        elif args.command == 'batch':
            from .batch import run_batch
            return run_batch(args)
        elif args.command == 'presets':
            from .presets import run_presets
            return run_presets(args)
        elif args.command == 'web':
            from .web import run_web
            return run_web(args)
        else:
            parser.print_help()
            return 1

    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user.")
        return 130
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        if '--verbose' in sys.argv or '-v' in sys.argv:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
