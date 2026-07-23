"""`viroforge benchmark` command: validate analysis pipelines against ground truth.

v1 implements Module 1 (QC): contamination removal and viral retention.
"""

from __future__ import annotations

import sys
from pathlib import Path

import json

from ..benchmarking import benchmark_qc, read_labels, read_names
from ..benchmarking.report import (
    assembly_to_markdown,
    taxonomy_to_markdown,
    to_markdown,
    write_reports,
)


def run_benchmark(args) -> int:
    if args.benchmark_command == "qc":
        return _run_qc(args)
    if args.benchmark_command == "assembly":
        return _run_assembly(args)
    if args.benchmark_command == "taxonomy":
        return _run_taxonomy(args)
    print("Usage: viroforge benchmark {qc,assembly,taxonomy} ...", file=sys.stderr)
    return 2


def _run_taxonomy(args) -> int:
    from ..benchmarking.taxonomy import (
        PARSERS,
        benchmark_taxonomy,
        benchmark_taxonomy_contigs,
        detect_format,
        parse_generic,
    )

    if not Path(args.ground_truth).exists():
        print(f"ERROR: file not found: {args.ground_truth}", file=sys.stderr)
        return 2
    if args.format not in PARSERS:
        print(f"ERROR: --format must be one of {sorted(PARSERS)}", file=sys.stderr)
        return 2

    # Auto-detect format if requested
    if args.format == "auto":
        detect_path = args.contig_taxonomy if args.mode == "contig-based" else args.pipeline_output
        if detect_path and Path(detect_path).exists():
            detected = detect_format(detect_path)
            if detected == "generic":
                # Unknown format — require user to specify columns
                print("ERROR: could not recognize the classification output format.\n\n"
                      "Your file must contain per-read classifications with numeric\n"
                      "NCBI taxids (not species names or abundance profiles).\n\n"
                      "Please re-run with --format generic and specify which columns\n"
                      "contain read IDs and taxids. For example:\n\n"
                      "  viroforge benchmark taxonomy --format generic \\\n"
                      "      --read-id-column 1 --taxid-column 3 ...\n",
                      file=sys.stderr)
                return 2
            args.format = detected
            print(f"Auto-detected format: {args.format}", file=sys.stderr)
        else:
            print("ERROR: --format auto requires a valid input file to inspect", file=sys.stderr)
            return 2
    tax_gt = json.loads(Path(args.ground_truth).read_text()).get("benchmarking", {}).get("taxonomy")
    if not tax_gt:
        print("ERROR: metadata has no benchmarking.taxonomy block. Regenerate the "
              "dataset with a current version to export per-genome taxonomy.",
              file=sys.stderr)
        return 2

    tree = None
    if args.taxdump_dir:
        from ..benchmarking.ncbi_tree import NcbiTree
        if not (Path(args.taxdump_dir) / "nodes.dmp").exists():
            print(f"ERROR: nodes.dmp not found in {args.taxdump_dir}", file=sys.stderr)
            return 2
        tree = NcbiTree.from_dir(args.taxdump_dir)

    if args.mode == "contig-based":
        missing = [f for f in ("contig_taxonomy", "contigs", "genomes")
                   if getattr(args, f) is None]
        if missing:
            print(f"ERROR: contig-based mode needs --{', --'.join(m.replace('_', '-') for m in missing)}",
                  file=sys.stderr)
            return 2
        for p in (args.contig_taxonomy, args.contigs, args.genomes):
            if not Path(p).exists():
                print(f"ERROR: file not found: {p}", file=sys.stderr)
                return 2
        if args.chimera_handling == "lca" and tree is None:
            print("ERROR: --chimera-handling lca requires --taxdump-dir", file=sys.stderr)
            return 2
        if args.format == "generic":
            assignments = parse_generic(args.contig_taxonomy, args.read_id_column, args.taxid_column)
        else:
            assignments = PARSERS[args.format](args.contig_taxonomy)
        metrics = benchmark_taxonomy_contigs(
            assignments, args.contigs, args.genomes, tax_gt,
            ncbi_tree=tree, chimera_handling=args.chimera_handling)
        kind = "taxonomy_contig"
    else:
        if not args.pipeline_output or not Path(args.pipeline_output).exists():
            print("ERROR: read-based mode needs --pipeline-output", file=sys.stderr)
            return 2
        if args.format == "generic":
            assignments = parse_generic(args.pipeline_output, args.read_id_column, args.taxid_column)
        else:
            assignments = PARSERS[args.format](args.pipeline_output)
        metrics = benchmark_taxonomy(assignments, tax_gt, ncbi_tree=tree)
        kind = "taxonomy"

    write_reports(metrics, json_path=args.output, md_path=args.markdown, kind=kind)
    from ..benchmarking.report import taxonomy_contig_to_markdown
    render = taxonomy_contig_to_markdown if kind == "taxonomy_contig" else taxonomy_to_markdown
    print(render(metrics))
    return 0 if metrics["reliable"] else 1


def _run_assembly(args) -> int:
    from ..benchmarking.assembly import benchmark_assembly

    for p in (args.contigs, args.genomes):
        if not Path(p).exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 2
    metadata = None
    if args.ground_truth:
        gt = Path(args.ground_truth)
        if not gt.exists():
            print(f"ERROR: file not found: {gt}", file=sys.stderr)
            return 2
        metadata = json.loads(gt.read_text())

    metrics = benchmark_assembly(args.contigs, args.genomes, metadata)
    write_reports(metrics, json_path=args.output, md_path=args.markdown, kind="assembly")
    print(assembly_to_markdown(metrics))
    return 0


def _run_qc(args) -> int:
    raw_paths = [Path(p) for p in args.raw_reads]
    cleaned_paths = [Path(p) for p in args.cleaned_reads]
    for p in raw_paths + cleaned_paths:
        if not p.exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 2

    keep_remove = None
    if args.keep_remove:
        keep_remove = {}
        for item in args.keep_remove:
            src, _, decision = item.partition(":")
            if decision not in ("keep", "remove"):
                print(f"ERROR: --keep-remove expects source:keep|remove, got {item}", file=sys.stderr)
                return 2
            keep_remove[src] = decision

    raw = read_labels(raw_paths)
    kept = read_names(cleaned_paths)
    metrics = benchmark_qc(raw, kept, keep_remove=keep_remove)

    write_reports(metrics, json_path=args.output, md_path=args.markdown)
    print(to_markdown(metrics))

    if not metrics["reliable"]:
        return 1  # unreliable match rate: signal via exit code
    return 0
