#!/usr/bin/env python3
"""
ViroForge Hybrid Dataset Generator

Generate matched short-read and long-read datasets from the same collection
for hybrid assembly benchmarking. This script ensures identical genome
compositions by using the same random seed and collection for both platforms.

Hybrid assembly combines short-read accuracy with long-read length to produce
superior viral genome assemblies. Supported assemblers include:
- Unicycler (SPAdes + miniasm hybrid)
- SPAdes (--nanopore mode)
- MaSuRCA (hybrid OLC assembler)
- HybridSPAdes

Usage:
    # Generate NovaSeq + Nanopore hybrid dataset
    python generate_hybrid_dataset.py \\
        --collection-id 9 \\
        --output data/gut_hybrid \\
        --short-platform novaseq \\
        --coverage 30 \\
        --depth 15 \\
        --seed 42

    # Generate MiSeq + Nanopore hybrid dataset
    python generate_hybrid_dataset.py \\
        --collection-id 1 \\
        --output data/soil_hybrid \\
        --short-platform miseq \\
        --coverage 50 \\
        --depth 20 \\
        --seed 123

Author: ViroForge Development Team
Date: 2025-11-10
"""

import argparse
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n{'='*80}")
    print(f"{description}")
    print(f"{'='*80}")
    print(f"Command: {' '.join(cmd)}\n")

    try:
        result = subprocess.run(cmd, check=True, capture_output=False, text=True)
        print(f"\n✓ {description} complete")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ {description} failed: {e}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print(f"\n✗ Command not found: {cmd[0]}", file=sys.stderr)
        print(f"Make sure ViroForge is properly installed", file=sys.stderr)
        return False


def create_hybrid_metadata(args, short_output, long_output, output_dir):
    """Create hybrid metadata JSON linking the two datasets."""
    metadata = {
        "hybrid_dataset": True,
        "generation_timestamp": datetime.now().isoformat(),
        "viroforge_version": "0.13.0",
        "random_seed": args.seed,
        "collection": {
            "id": args.collection_id
        },
        "vlp_protocol": args.vlp_protocol if not args.no_vlp else "none",
        "contamination_level": args.contamination_level,
        "short_reads": {
            "platform": args.short_platform,
            "coverage": args.coverage,
            "read_length": args.read_length,
            "insert_size": args.insert_size,
            "output_dir": str(short_output.relative_to(output_dir)),
            "r1": str((short_output / "fastq" / f"*_R1.fastq").relative_to(output_dir)),
            "r2": str((short_output / "fastq" / f"*_R2.fastq").relative_to(output_dir))
        },
        "long_reads": {
            "platform": "nanopore",
            "depth": args.depth,
            "output_dir": str(long_output.relative_to(output_dir)),
            "chemistry": args.ont_chemistry,
            "read_length_mean": args.ont_read_length
        },
        "composition_consistency": {
            "same_seed": True,
            "same_collection": True,
            "same_vlp_protocol": True,
            "note": "Both datasets generated from identical genome composition"
        },
        "usage": {
            "unicycler_example": f"unicycler -1 {short_output}/fastq/*_R1.fastq -2 {short_output}/fastq/*_R2.fastq -l {long_output}/fastq/*.fastq* -o results/unicycler",
            "spades_example": f"spades.py --meta -1 {short_output}/fastq/*_R1.fastq -2 {short_output}/fastq/*_R2.fastq --nanopore {long_output}/fastq/*.fastq* -o results/spades"
        }
    }

    # Write metadata
    metadata_path = output_dir / "hybrid_metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"\n✓ Created hybrid metadata: {metadata_path}")
    return metadata_path


def main():
    parser = argparse.ArgumentParser(
        description='Generate matched short-read and long-read datasets for hybrid assembly',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # NovaSeq + Nanopore (recommended)
  python generate_hybrid_dataset.py \\
      --collection-id 9 \\
      --output data/gut_hybrid \\
      --short-platform novaseq \\
      --coverage 30 --depth 15

  # MiSeq + Nanopore (cost-effective)
  python generate_hybrid_dataset.py \\
      --collection-id 1 \\
      --output data/soil_hybrid \\
      --short-platform miseq \\
      --coverage 50 --depth 20

Supported Hybrid Assemblers:
  - Unicycler: Combines SPAdes (short) + miniasm (long)
  - SPAdes: --nanopore mode
  - MaSuRCA: Overlap-layout-consensus hybrid
  - HybridSPAdes: Specialized for metagenomes
        """
    )

    # Required arguments
    parser.add_argument(
        '--collection-id',
        type=int,
        required=True,
        help='Collection ID to generate from'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output directory (will create short_reads/ and long_reads/ subdirs)'
    )

    # Platform selection
    parser.add_argument(
        '--short-platform',
        choices=['novaseq', 'miseq', 'hiseq'],
        default='novaseq',
        help='Short-read platform (default: novaseq)'
    )

    # Short-read parameters
    parser.add_argument(
        '--coverage',
        type=float,
        default=30.0,
        help='Short-read coverage (default: 30x)'
    )

    parser.add_argument(
        '--read-length',
        type=int,
        default=150,
        help='Short-read length in bp (default: 150)'
    )

    parser.add_argument(
        '--insert-size',
        type=int,
        default=350,
        help='Insert size for paired-end (default: 350)'
    )

    # Long-read parameters (Nanopore only)
    parser.add_argument(
        '--depth',
        type=float,
        default=15.0,
        help='Nanopore long-read depth (default: 15x)'
    )

    parser.add_argument(
        '--ont-chemistry',
        choices=['R9.4', 'R10.4'],
        default='R10.4',
        help='Nanopore chemistry (default: R10.4)'
    )

    parser.add_argument(
        '--ont-read-length',
        type=int,
        default=20000,
        help='Nanopore mean read length (default: 20000)'
    )

    # VLP and contamination
    parser.add_argument(
        '--no-vlp',
        action='store_true',
        help='Skip VLP enrichment (bulk metagenome)'
    )

    parser.add_argument(
        '--vlp-protocol',
        choices=['tangential_flow', 'syringe', 'ultracentrifugation', 'norgen'],
        default='tangential_flow',
        help='VLP enrichment protocol (default: tangential_flow)'
    )

    parser.add_argument(
        '--contamination-level',
        choices=['clean', 'realistic', 'heavy'],
        default='realistic',
        help='Contamination level (default: realistic)'
    )

    # Other parameters
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show commands without executing'
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    short_output = output_dir / "short_reads"
    long_output = output_dir / "long_reads"

    print(f"""
{'='*80}
ViroForge Hybrid Dataset Generator
{'='*80}

Configuration:
  Collection ID:      {args.collection_id}
  Output directory:   {output_dir}
  Random seed:        {args.seed}

Short-read platform:  {args.short_platform}
  Coverage:           {args.coverage}x
  Read length:        {args.read_length} bp
  Output:             {short_output}

Long-read platform:   nanopore
  Depth:              {args.depth}x
  Chemistry:          {args.ont_chemistry}
  Read length:        {args.ont_read_length} bp
  Output:             {long_output}

VLP protocol:         {args.vlp_protocol if not args.no_vlp else 'none (bulk)'}
Contamination:        {args.contamination_level}
{'='*80}
""")

    if args.dry_run:
        print("DRY RUN MODE - Commands will be shown but not executed\n")

    # Path to generate_fastq_dataset.py
    generator_script = Path(__file__).parent / "generate_fastq_dataset.py"

    if not generator_script.exists():
        print(f"Error: Could not find {generator_script}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Step 1: Generate SHORT-READ dataset
    # ========================================================================
    short_cmd = [
        sys.executable, str(generator_script),
        '--collection-id', str(args.collection_id),
        '--output', str(short_output),
        '--platform', args.short_platform,
        '--coverage', str(args.coverage),
        '--read-length', str(args.read_length),
        '--insert-size', str(args.insert_size),
        '--contamination-level', args.contamination_level,
        '--seed', str(args.seed)
    ]

    if args.no_vlp:
        short_cmd.append('--no-vlp')
    else:
        short_cmd.extend(['--vlp-protocol', args.vlp_protocol])

    if args.dry_run:
        print(f"SHORT-READ COMMAND:\n{' '.join(short_cmd)}\n")
    else:
        success = run_command(short_cmd, f"Step 1/2: Generating short-read dataset ({args.short_platform})")
        if not success:
            print("\n✗ Failed to generate short-read dataset", file=sys.stderr)
            sys.exit(1)

    # ========================================================================
    # Step 2: Generate LONG-READ dataset (Nanopore)
    # ========================================================================
    long_cmd = [
        sys.executable, str(generator_script),
        '--collection-id', str(args.collection_id),
        '--output', str(long_output),
        '--platform', 'nanopore',
        '--depth', str(args.depth),
        '--ont-chemistry', args.ont_chemistry,
        '--ont-read-length', str(args.ont_read_length),
        '--contamination-level', args.contamination_level,
        '--seed', str(args.seed)
    ]

    if args.no_vlp:
        long_cmd.append('--no-vlp')
    else:
        long_cmd.extend(['--vlp-protocol', args.vlp_protocol])

    if args.dry_run:
        print(f"LONG-READ COMMAND:\n{' '.join(long_cmd)}\n")
    else:
        success = run_command(long_cmd, "Step 2/2: Generating long-read dataset (nanopore)")
        if not success:
            print("\n✗ Failed to generate long-read dataset", file=sys.stderr)
            sys.exit(1)

    # ========================================================================
    # Step 3: Create hybrid metadata
    # ========================================================================
    if not args.dry_run:
        metadata_path = create_hybrid_metadata(args, short_output, long_output, output_dir)

    # ========================================================================
    # Success message
    # ========================================================================
    print(f"""
{'='*80}
✓ HYBRID DATASET GENERATION COMPLETE
{'='*80}

Output structure:
  {output_dir}/
  ├── short_reads/        # {args.short_platform} reads
  │   ├── fastq/
  │   │   ├── *_R1.fastq
  │   │   └── *_R2.fastq
  │   └── metadata/
  ├── long_reads/         # nanopore reads
  │   ├── fastq/
  │   │   └── *.fastq*
  │   └── metadata/
  └── hybrid_metadata.json

Next steps - Hybrid Assembly:

1. Unicycler (recommended):
   unicycler \\
       -1 {short_output}/fastq/*_R1.fastq \\
       -2 {short_output}/fastq/*_R2.fastq \\
       -l {long_output}/fastq/*.fastq* \\
       -o results/unicycler_hybrid

2. SPAdes hybrid mode:
   spades.py --meta \\
       -1 {short_output}/fastq/*_R1.fastq \\
       -2 {short_output}/fastq/*_R2.fastq \\
       --nanopore {long_output}/fastq/*.fastq* \\
       -o results/spades_hybrid

3. Evaluate against ground truth:
   # Both datasets have identical genome compositions
   # Compare assembly to {short_output}/metadata/*_composition.tsv

See docs/HYBRID_ASSEMBLY_TUTORIAL.md for more examples.
{'='*80}
""")


if __name__ == '__main__':
    main()
