#!/usr/bin/env python3
"""
ViroForge Database Inspector

Interactive CLI tool for exploring the ViroForge viral genome database.

Usage:
    viroforge-db stats                          # Database statistics
    viroforge-db search --family Siphoviridae   # Search by taxonomy
    viroforge-db search --host Escherichia      # Search by host
    viroforge-db collection gut_virome          # Collection details
    viroforge-db taxonomy --rank family         # Browse taxonomy
    viroforge-db export --family Crassvirales   # Export genomes

Author: ViroForge Development Team
Date: 2025-11-01
"""

import sqlite3
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
import json


class DatabaseInspector:
    """Inspector for ViroForge database."""

    def __init__(self, db_path: str):
        self.db_path = Path(db_path)

        if not self.db_path.exists():
            print(f"❌ Error: Database not found at {db_path}", file=sys.stderr)
            print(f"   Run the database creation pipeline first.", file=sys.stderr)
            sys.exit(1)

    def get_connection(self):
        """Get database connection with row factory."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def stats(self, verbose: bool = False) -> Dict:
        """Get comprehensive database statistics."""
        conn = self.get_connection()

        stats = {
            'genomes': self._genome_stats(conn),
            'taxonomy': self._taxonomy_stats(conn),
            'collections': self._collection_stats(conn),
            'hosts': self._host_stats(conn)
        }

        conn.close()
        return stats

    def _genome_stats(self, conn: sqlite3.Connection) -> Dict:
        """Get genome statistics."""
        # Total genomes
        cursor = conn.execute("SELECT COUNT(*) FROM genomes")
        total = cursor.fetchone()[0]

        # Genome types
        cursor = conn.execute("""
            SELECT genome_type, COUNT(*) as count
            FROM genomes
            GROUP BY genome_type
            ORDER BY count DESC
        """)
        types = dict(cursor.fetchall())

        # Length statistics
        cursor = conn.execute("""
            SELECT
                AVG(length) as mean,
                MIN(length) as min,
                MAX(length) as max
            FROM genomes
        """)
        row = cursor.fetchone()
        length_stats = {
            'mean': int(row['mean']) if row['mean'] else 0,
            'min': row['min'] or 0,
            'max': row['max'] or 0
        }

        # GC content statistics
        cursor = conn.execute("""
            SELECT
                AVG(gc_content) as mean,
                MIN(gc_content) as min,
                MAX(gc_content) as max
            FROM genomes
            WHERE gc_content IS NOT NULL
        """)
        row = cursor.fetchone()
        gc_stats = {
            'mean': float(row['mean']) if row['mean'] else 0,
            'min': float(row['min']) if row['min'] else 0,
            'max': float(row['max']) if row['max'] else 0
        }

        return {
            'total': total,
            'types': types,
            'length': length_stats,
            'gc_content': gc_stats
        }

    def _taxonomy_stats(self, conn: sqlite3.Connection) -> Dict:
        """Get taxonomy statistics."""
        # Count with ICTV taxonomy
        cursor = conn.execute("""
            SELECT COUNT(*) FROM taxonomy WHERE realm IS NOT NULL
        """)
        with_realm = cursor.fetchone()[0]

        cursor = conn.execute("""
            SELECT COUNT(*) FROM taxonomy WHERE species IS NOT NULL
        """)
        with_species = cursor.fetchone()[0]

        # Top realms
        cursor = conn.execute("""
            SELECT realm, COUNT(*) as count
            FROM taxonomy
            WHERE realm IS NOT NULL
            GROUP BY realm
            ORDER BY count DESC
            LIMIT 5
        """)
        top_realms = dict(cursor.fetchall())

        # Top families
        cursor = conn.execute("""
            SELECT family, COUNT(*) as count
            FROM taxonomy
            WHERE family IS NOT NULL
            GROUP BY family
            ORDER BY count DESC
            LIMIT 10
        """)
        top_families = dict(cursor.fetchall())

        return {
            'with_realm': with_realm,
            'with_species': with_species,
            'top_realms': top_realms,
            'top_families': top_families
        }

    def _collection_stats(self, conn: sqlite3.Connection) -> Dict:
        """Get collection statistics."""
        # Check if table exists
        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name='body_site_collections'
        """)
        if not cursor.fetchone():
            return {
                'total': 0,
                'collections': []
            }

        # Check schema version (old vs new)
        cursor = conn.execute("PRAGMA table_info(body_site_collections)")
        columns = [row[1] for row in cursor.fetchall()]

        cursor = conn.execute("SELECT COUNT(*) FROM body_site_collections")
        total = cursor.fetchone()[0]

        # Handle both old and new schema
        if 'collection_name' in columns:
            # Old schema
            cursor = conn.execute("""
                SELECT
                    collection_id,
                    collection_name as name,
                    '' as body_site,
                    n_genomes as genome_count
                FROM body_site_collections
                ORDER BY n_genomes DESC
            """)
        else:
            # New schema
            cursor = conn.execute("""
                SELECT
                    collection_id,
                    name,
                    body_site,
                    genome_count
                FROM body_site_collections
                ORDER BY genome_count DESC
            """)

        collections = [dict(row) for row in cursor.fetchall()]

        return {
            'total': total,
            'collections': collections
        }

    def _host_stats(self, conn: sqlite3.Connection) -> Dict:
        """Get host statistics."""
        # Check if table exists
        cursor = conn.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name='host_associations'
        """)
        if not cursor.fetchone():
            return {
                'unique_hosts': 0,
                'top_hosts': {}
            }

        cursor = conn.execute("""
            SELECT COUNT(DISTINCT host_name) FROM host_associations
        """)
        unique_hosts = cursor.fetchone()[0]

        cursor = conn.execute("""
            SELECT host_name, COUNT(*) as count
            FROM host_associations
            GROUP BY host_name
            ORDER BY count DESC
            LIMIT 10
        """)
        top_hosts = dict(cursor.fetchall())

        return {
            'unique_hosts': unique_hosts,
            'top_hosts': top_hosts
        }

    def print_stats(self, verbose: bool = False):
        """Print formatted statistics."""
        stats = self.stats(verbose)

        print("=" * 70)
        print("ViroForge Database Statistics")
        print("=" * 70)
        print()

        # Genome statistics
        print("📊 Genome Statistics:")
        print(f"  Total genomes:      {stats['genomes']['total']:,}")
        print(f"  Mean length:        {stats['genomes']['length']['mean']:,} bp")
        print(f"  Length range:       {stats['genomes']['length']['min']:,} - {stats['genomes']['length']['max']:,} bp")
        print(f"  Mean GC content:    {stats['genomes']['gc_content']['mean']:.1%}")
        print()

        print("  Genome types:")
        for gtype, count in stats['genomes']['types'].items():
            pct = count / stats['genomes']['total']
            print(f"    {gtype:15s} {count:>6,} ({pct:>5.1%})")
        print()

        # Taxonomy statistics
        print("🧬 Taxonomy Coverage:")
        total = stats['genomes']['total']
        with_realm = stats['taxonomy']['with_realm']
        with_species = stats['taxonomy']['with_species']

        print(f"  With ICTV realm:    {with_realm:,} ({with_realm/total:.1%})")
        print(f"  With species:       {with_species:,} ({with_species/total:.1%})")
        print()

        if stats['taxonomy']['top_realms']:
            print("  Top 5 Realms:")
            for realm, count in stats['taxonomy']['top_realms'].items():
                pct = count / total
                print(f"    {realm:30s} {count:>6,} ({pct:>5.1%})")
            print()

        if stats['taxonomy']['top_families']:
            print("  Top 10 Families:")
            for family, count in list(stats['taxonomy']['top_families'].items())[:10]:
                pct = count / total
                print(f"    {family:30s} {count:>6,} ({pct:>5.1%})")
            print()

        # Collection statistics
        if stats['collections']['total'] > 0:
            print("📚 Collections:")
            print(f"  Total collections:  {stats['collections']['total']}")
            print()
            print("  Available collections:")
            for coll in stats['collections']['collections']:
                print(f"    {coll['name'][:50]:50s} {coll['genome_count']:>4} genomes")
            print()

        # Host statistics
        if stats['hosts']['unique_hosts'] > 0:
            print("🦠 Host Associations:")
            print(f"  Unique hosts:       {stats['hosts']['unique_hosts']:,}")
            print()
            if stats['hosts']['top_hosts']:
                print("  Top 10 Hosts:")
                for host, count in list(stats['hosts']['top_hosts'].items())[:10]:
                    print(f"    {host[:40]:40s} {count:>6,} phages")
                print()

    def search(
        self,
        realm: Optional[str] = None,
        kingdom: Optional[str] = None,
        phylum: Optional[str] = None,
        tax_class: Optional[str] = None,
        order: Optional[str] = None,
        family: Optional[str] = None,
        genus: Optional[str] = None,
        species: Optional[str] = None,
        host: Optional[str] = None,
        genome_type: Optional[str] = None,
        length_min: Optional[int] = None,
        length_max: Optional[int] = None,
        gc_min: Optional[float] = None,
        gc_max: Optional[float] = None,
        limit: int = 100
    ) -> List[Dict]:
        """
        Search genomes with multiple filters.

        Returns list of matching genomes with metadata.
        """
        conn = self.get_connection()

        query = """
            SELECT
                g.genome_id,
                g.genome_name as species_name,
                g.genome_type,
                g.length,
                g.gc_content,
                t.realm,
                t.family,
                t.genus,
                t.species as tax_species
            FROM genomes g
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
        """

        where_clauses = []
        params = []

        # Taxonomy filters
        if realm:
            where_clauses.append("t.realm LIKE ?")
            params.append(f"%{realm}%")
        if kingdom:
            where_clauses.append("t.kingdom LIKE ?")
            params.append(f"%{kingdom}%")
        if phylum:
            where_clauses.append("t.phylum LIKE ?")
            params.append(f"%{phylum}%")
        if tax_class:
            where_clauses.append("t.class LIKE ?")
            params.append(f"%{tax_class}%")
        if order:
            where_clauses.append("t.tax_order LIKE ?")
            params.append(f"%{order}%")
        if family:
            where_clauses.append("t.family LIKE ?")
            params.append(f"%{family}%")
        if genus:
            where_clauses.append("t.genus LIKE ?")
            params.append(f"%{genus}%")
        if species:
            where_clauses.append("(t.species LIKE ? OR g.genome_name LIKE ?)")
            params.extend([f"%{species}%", f"%{species}%"])

        # Host filter
        if host:
            where_clauses.append("""
                EXISTS (
                    SELECT 1 FROM host_associations h
                    WHERE h.genome_id = g.genome_id
                    AND h.host_name LIKE ?
                )
            """)
            params.append(f"%{host}%")

        # Genome characteristics
        if genome_type:
            where_clauses.append("g.genome_type = ?")
            params.append(genome_type)
        if length_min:
            where_clauses.append("g.length >= ?")
            params.append(length_min)
        if length_max:
            where_clauses.append("g.length <= ?")
            params.append(length_max)
        if gc_min:
            where_clauses.append("g.gc_content >= ?")
            params.append(gc_min)
        if gc_max:
            where_clauses.append("g.gc_content <= ?")
            params.append(gc_max)

        if where_clauses:
            query += " WHERE " + " AND ".join(where_clauses)

        query += f" LIMIT {limit}"

        cursor = conn.execute(query, params)
        results = [dict(row) for row in cursor.fetchall()]

        conn.close()
        return results

    def collection_info(self, collection_id: str) -> Optional[Dict]:
        """Get detailed information about a collection."""
        conn = self.get_connection()

        # Get collection metadata
        cursor = conn.execute("""
            SELECT * FROM body_site_collections
            WHERE collection_id LIKE ?
        """, (f"%{collection_id}%",))

        meta = cursor.fetchone()
        if not meta:
            conn.close()
            return None

        meta = dict(meta)

        # Get genomes in collection
        cursor = conn.execute("""
            SELECT
                cg.genome_id,
                cg.abundance,
                g.genome_name as species_name,
                g.length,
                t.family
            FROM collection_genomes cg
            JOIN genomes g ON cg.genome_id = g.genome_id
            LEFT JOIN taxonomy t ON g.genome_id = t.genome_id
            WHERE cg.collection_id = ?
            ORDER BY cg.abundance DESC
        """, (meta['collection_id'],))

        genomes = [dict(row) for row in cursor.fetchall()]

        # Get taxonomic composition
        cursor = conn.execute("""
            SELECT t.family, COUNT(*) as count
            FROM collection_genomes cg
            JOIN taxonomy t ON cg.genome_id = t.genome_id
            WHERE cg.collection_id = ?
            AND t.family IS NOT NULL
            GROUP BY t.family
            ORDER BY count DESC
            LIMIT 10
        """, (meta['collection_id'],))

        family_composition = dict(cursor.fetchall())

        conn.close()

        return {
            'metadata': meta,
            'genomes': genomes,
            'family_composition': family_composition
        }

    def taxonomy_browse(self, rank: str = 'family', limit: int = 50) -> List[Tuple[str, int]]:
        """Browse taxonomy at specified rank."""
        conn = self.get_connection()

        rank_map = {
            'realm': 'realm',
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'class': 'class',
            'order': 'tax_order',
            'family': 'family',
            'genus': 'genus',
            'species': 'species'
        }

        if rank not in rank_map:
            print(f"❌ Error: Invalid rank '{rank}'", file=sys.stderr)
            print(f"   Valid ranks: {', '.join(rank_map.keys())}", file=sys.stderr)
            return []

        db_column = rank_map[rank]

        cursor = conn.execute(f"""
            SELECT {db_column}, COUNT(*) as count
            FROM taxonomy
            WHERE {db_column} IS NOT NULL
            GROUP BY {db_column}
            ORDER BY count DESC
            LIMIT ?
        """, (limit,))

        results = cursor.fetchall()
        conn.close()

        return results

    def export_genomes(
        self,
        output_file: Path,
        format: str = 'fasta',
        **search_kwargs
    ) -> int:
        """Export genomes matching search criteria."""
        # Search for genomes
        genomes = self.search(**search_kwargs)

        if not genomes:
            print(f"⚠️  No genomes match the search criteria.", file=sys.stderr)
            return 0

        conn = self.get_connection()

        if format == 'fasta':
            with open(output_file, 'w') as f:
                for genome in genomes:
                    cursor = conn.execute(
                        "SELECT sequence FROM genomes WHERE genome_id = ?",
                        (genome['genome_id'],)
                    )
                    row = cursor.fetchone()
                    if row and row['sequence']:
                        f.write(f">{genome['genome_id']} {genome['species_name']}\n")
                        seq = row['sequence']
                        # Write in 80-character lines
                        for i in range(0, len(seq), 80):
                            f.write(seq[i:i+80] + '\n')

        elif format == 'tsv':
            with open(output_file, 'w') as f:
                # Header
                f.write("genome_id\tspecies_name\tgenome_type\tlength\tgc_content\trealm\tfamily\tgenus\n")
                for genome in genomes:
                    f.write(f"{genome['genome_id']}\t")
                    f.write(f"{genome['species_name']}\t")
                    f.write(f"{genome['genome_type']}\t")
                    f.write(f"{genome['length']}\t")
                    f.write(f"{genome.get('gc_content', '')}\t")
                    f.write(f"{genome.get('realm', '')}\t")
                    f.write(f"{genome.get('family', '')}\t")
                    f.write(f"{genome.get('genus', '')}\n")

        elif format == 'json':
            with open(output_file, 'w') as f:
                json.dump(genomes, f, indent=2)

        conn.close()
        return len(genomes)


def main():
    parser = argparse.ArgumentParser(
        description='ViroForge Database Inspector',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Database statistics
  viroforge-db stats

  # Search by family
  viroforge-db search --family Siphoviridae --limit 20

  # Search by host
  viroforge-db search --host "Escherichia coli"

  # Complex search
  viroforge-db search --family Crassvirales --length-min 100000 --gc-min 0.4

  # Collection details
  viroforge-db collection gut_virome

  # Browse taxonomy
  viroforge-db taxonomy --rank family

  # Export genomes
  viroforge-db export --family Crassvirales --output crassvirales.fasta
        """
    )

    parser.add_argument(
        '--database',
        default='viroforge/data/viral_genomes.db',
        help='Path to ViroForge database (default: viroforge/data/viral_genomes.db)'
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Database statistics')
    stats_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    # Search command
    search_parser = subparsers.add_parser('search', help='Search genomes')
    search_parser.add_argument('--realm', help='Filter by realm')
    search_parser.add_argument('--kingdom', help='Filter by kingdom')
    search_parser.add_argument('--phylum', help='Filter by phylum')
    search_parser.add_argument('--class', dest='tax_class', help='Filter by class')
    search_parser.add_argument('--order', help='Filter by order')
    search_parser.add_argument('--family', help='Filter by family')
    search_parser.add_argument('--genus', help='Filter by genus')
    search_parser.add_argument('--species', help='Filter by species')
    search_parser.add_argument('--host', help='Filter by host')
    search_parser.add_argument('--genome-type', help='Filter by genome type (dsDNA, ssRNA, etc)')
    search_parser.add_argument('--length-min', type=int, help='Minimum genome length (bp)')
    search_parser.add_argument('--length-max', type=int, help='Maximum genome length (bp)')
    search_parser.add_argument('--gc-min', type=float, help='Minimum GC content (0-1)')
    search_parser.add_argument('--gc-max', type=float, help='Maximum GC content (0-1)')
    search_parser.add_argument('--limit', type=int, default=100, help='Maximum results (default: 100)')
    search_parser.add_argument('--format', choices=['table', 'tsv', 'json'], default='table', help='Output format')

    # Collection command
    collection_parser = subparsers.add_parser('collection', help='Collection details')
    collection_parser.add_argument('collection_id', help='Collection ID (partial match supported)')
    collection_parser.add_argument('--top', type=int, default=20, help='Show top N genomes (default: 20)')

    # Taxonomy command
    taxonomy_parser = subparsers.add_parser('taxonomy', help='Browse taxonomy')
    taxonomy_parser.add_argument(
        '--rank',
        choices=['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
        default='family',
        help='Taxonomic rank to browse (default: family)'
    )
    taxonomy_parser.add_argument('--limit', type=int, default=50, help='Maximum results (default: 50)')

    # Export command
    export_parser = subparsers.add_parser('export', help='Export genomes')
    export_parser.add_argument('--realm', help='Filter by realm')
    export_parser.add_argument('--family', help='Filter by family')
    export_parser.add_argument('--genus', help='Filter by genus')
    export_parser.add_argument('--species', help='Filter by species')
    export_parser.add_argument('--host', help='Filter by host')
    export_parser.add_argument('--genome-type', help='Filter by genome type')
    export_parser.add_argument('--length-min', type=int, help='Minimum genome length (bp)')
    export_parser.add_argument('--length-max', type=int, help='Maximum genome length (bp)')
    export_parser.add_argument('--limit', type=int, default=10000, help='Maximum genomes (default: 10000)')
    export_parser.add_argument(
        '--format',
        choices=['fasta', 'tsv', 'json'],
        default='fasta',
        help='Export format (default: fasta)'
    )
    export_parser.add_argument(
        '--output',
        required=True,
        help='Output file path'
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    # Create inspector
    inspector = DatabaseInspector(args.database)

    # Execute command
    if args.command == 'stats':
        inspector.print_stats(args.verbose)

    elif args.command == 'search':
        results = inspector.search(
            realm=args.realm,
            kingdom=args.kingdom,
            phylum=args.phylum,
            tax_class=args.tax_class,
            order=args.order,
            family=args.family,
            genus=args.genus,
            species=args.species,
            host=args.host,
            genome_type=args.genome_type,
            length_min=args.length_min,
            length_max=args.length_max,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            limit=args.limit
        )

        if not results:
            print("⚠️  No genomes match the search criteria.")
            sys.exit(0)

        print(f"Found {len(results)} genomes")
        print()

        if args.format == 'table':
            # Print as table
            print(f"{'Genome ID':<18} {'Species':<40} {'Type':<8} {'Length':>10} {'Family':<20}")
            print("-" * 110)
            for genome in results:
                species = genome['species_name'][:38] if genome['species_name'] else '-'
                family = genome.get('family', '-') or '-'
                family = family[:18] if family else '-'
                print(f"{genome['genome_id']:<18} {species:<40} {genome['genome_type']:<8} {genome['length']:>10,} {family:<20}")

        elif args.format == 'tsv':
            print("genome_id\tspecies_name\tgenome_type\tlength\tgc_content\trealm\tfamily\tgenus")
            for genome in results:
                print(f"{genome['genome_id']}\t{genome['species_name']}\t{genome['genome_type']}\t"
                      f"{genome['length']}\t{genome.get('gc_content', '')}\t{genome.get('realm', '')}\t"
                      f"{genome.get('family', '')}\t{genome.get('genus', '')}")

        elif args.format == 'json':
            print(json.dumps(results, indent=2))

    elif args.command == 'collection':
        info = inspector.collection_info(args.collection_id)

        if not info:
            print(f"❌ Error: Collection '{args.collection_id}' not found", file=sys.stderr)
            sys.exit(1)

        meta = info['metadata']

        print("=" * 70)
        print(f"Collection: {meta['name']}")
        print("=" * 70)
        print()
        print(f"ID:               {meta['collection_id']}")
        print(f"Body Site:        {meta['body_site']}")
        print(f"Sample Type:      {meta['sample_type']}")
        print(f"Health Status:    {meta['health_status']}")
        print(f"Genome Count:     {meta['genome_count']}")
        print(f"Phage Fraction:   {meta['phage_fraction']:.1%}")
        print(f"Curation Date:    {meta['curation_date']}")
        print()

        if meta['description']:
            print(f"Description:")
            print(f"  {meta['description']}")
            print()

        print(f"Top {args.top} Most Abundant Genomes:")
        print(f"{'Rank':<6} {'Abundance':>10} {'Genome ID':<18} {'Species':<40} {'Family':<20}")
        print("-" * 100)
        for i, genome in enumerate(info['genomes'][:args.top], 1):
            species = genome['species_name'][:38] if genome['species_name'] else '-'
            family = genome.get('family', '-') or '-'
            print(f"{i:<6} {genome['abundance']:>9.4%} {genome['genome_id']:<18} {species:<40} {family:<20}")
        print()

        if info['family_composition']:
            print("Taxonomic Composition (Top 10 Families):")
            for family, count in info['family_composition'].items():
                pct = count / meta['genome_count']
                print(f"  {family:40s} {count:>4} ({pct:>5.1%})")
            print()

    elif args.command == 'taxonomy':
        results = inspector.taxonomy_browse(args.rank, args.limit)

        if not results:
            print(f"⚠️  No taxonomy data at rank '{args.rank}'")
            sys.exit(0)

        print(f"Taxonomy Browser: {args.rank.capitalize()}")
        print("=" * 70)
        print()
        print(f"{'Rank':<4} {args.rank.capitalize():<45} {'Count':>10}")
        print("-" * 70)

        for i, (taxon, count) in enumerate(results, 1):
            taxon_display = taxon[:43] if taxon else '-'
            print(f"{i:<4} {taxon_display:<45} {count:>10,}")

        print()
        print(f"Showing top {len(results)} of {len(results)} {args.rank}s")

    elif args.command == 'export':
        print(f"Exporting genomes to {args.output}...")
        count = inspector.export_genomes(
            output_file=Path(args.output),
            format=args.format,
            realm=args.realm,
            family=args.family,
            genus=args.genus,
            species=args.species,
            host=args.host,
            genome_type=args.genome_type,
            length_min=args.length_min,
            length_max=args.length_max,
            limit=args.limit
        )

        if count > 0:
            print(f"✓ Exported {count} genomes to {args.output}")
        else:
            print("⚠️  No genomes exported")


if __name__ == '__main__':
    main()
