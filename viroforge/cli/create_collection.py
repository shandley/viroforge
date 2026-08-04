"""
Create Collection Command

Register a custom set of viral genomes as a new collection in the
ViroForge database. The collection can then be used with any
viroforge generate command.

Usage:
    viroforge create-collection \
        --genomes my_phages.fasta \
        --name "Penguin Gut Virome" \
        --host penguin \
        --body-site gut

    # With custom host genome for accurate contamination
    viroforge create-collection \
        --genomes my_phages.fasta \
        --name "Penguin Gut Virome" \
        --host penguin \
        --body-site gut \
        --host-genome penguin_genome.fasta
"""

import shutil
import sys
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def run_create_collection(args) -> int:
    """Create a custom collection from a FASTA file and register it in the database."""
    from viroforge.core.collection import CustomCollectionBuilder

    genomes_path = Path(args.genomes)
    if not genomes_path.exists():
        print(f"Error: FASTA file not found: {genomes_path}", file=sys.stderr)
        return 1

    db_path = Path(args.database)
    if not db_path.exists():
        print(f"Error: Database not found: {db_path}", file=sys.stderr)
        print("Run 'viroforge setup-db' first to create the database.", file=sys.stderr)
        return 1

    # Validate host genome if provided
    host_genome_path = None
    if args.host_genome:
        host_genome_path = Path(args.host_genome)
        if not host_genome_path.exists():
            print(f"Error: Host genome not found: {host_genome_path}", file=sys.stderr)
            return 1

    # Build the collection from FASTA
    builder = CustomCollectionBuilder(
        name=args.name,
        description=args.description,
        host=args.host,
        body_site=args.body_site,
    )

    abundances_path = Path(args.abundances) if args.abundances else None
    collection_meta, genomes = builder.from_fasta(
        genomes_path,
        abundances_tsv=abundances_path,
        random_seed=args.seed,
    )

    # Register in the database
    collection_id = builder.register_in_database(
        db_path=db_path,
        collection_meta=collection_meta,
        genomes=genomes,
    )

    # Install custom host genome if provided
    if host_genome_path:
        custom_refs_dir = Path("viroforge/data/references/custom")
        custom_refs_dir.mkdir(parents=True, exist_ok=True)
        dest = custom_refs_dir / f"host_{collection_id}_{args.host}.fasta"
        shutil.copy2(host_genome_path, dest)
        logger.info(f"Installed custom host genome: {dest}")

        # Update collection metadata in database to point to custom host
        import sqlite3
        conn = sqlite3.connect(db_path)
        conn.execute(
            "UPDATE body_site_collections SET host_organism = ? WHERE collection_id = ?",
            (f"custom:{dest}", collection_id),
        )
        conn.commit()
        conn.close()

    # Print summary
    print(f"\nCollection created successfully!")
    print(f"  Name: {args.name}")
    print(f"  Collection ID: {collection_id}")
    print(f"  Genomes: {len(genomes)}")
    print(f"  Host: {args.host}")
    print(f"  Body site: {args.body_site}")
    if host_genome_path:
        print(f"  Host genome: {host_genome_path} (installed)")
    else:
        print(f"  Host genome: default (human)")
        print(f"    Note: host DNA contamination will use human reference.")
        print(f"    For accurate {args.host} host contamination, re-run with:")
        print(f"    --host-genome <path_to_{args.host}_genome.fasta>")
    print(f"\nTo generate a dataset from this collection:")
    print(f"  viroforge generate --collection-id {collection_id} --output data/my_dataset")
    print()

    return 0
