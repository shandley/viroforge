#!/usr/bin/env python3
"""Build the bacterial background reference set from RefSeq at build time.

ViroForge needs bacterial sequence to make bulk metagenomes (``--no-vlp``)
realistic: real bulk stool is 60-80% bacterial, where a virome-only community
leaves the viral fraction dominant no matter how far the contamination knobs are
turned up. See issue #37.

Two rules shape this script, both learned from PR #39, which committed 4.4M
lines of FASTA and a list of hand-written accessions containing at least one
wrong entry:

1. **No accessions in the repository.** Genomes are resolved by taxon NAME
   against RefSeq at build time. Nothing is recalled from memory, and there is
   no stored identifier list to go stale or be wrong. What the script writes is
   the accession it actually retrieved, recorded in each FASTA header for
   provenance.

2. **Fragments, not genomes.** Bacterial genomes are 3-6 Mbp. Only 10 kb slices
   are fetched, using Entrez ``seq_start``/``seq_stop``, so the whole genome is
   never downloaded or stored. This follows ``host_fragments.fasta``: 48
   fragments totalling 528 KB stands in for a 3 Gbp human genome.

The output is small enough to bundle, or to rebuild on demand. It is discovered
at runtime by ``viroforge/data/references/resolver.py``.

Usage:
    python scripts/curate_bacterial_background.py
    python scripts/curate_bacterial_background.py --communities gut oral
    python scripts/curate_bacterial_background.py --fragments-per-taxon 2 --dry-run
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).parent.parent / "viroforge" / "data" / "references"
OUTPUT_NAME = "bacterial_fragments.fasta"

# NCBI asks for no more than 3 requests/second without an API key.
REQUEST_DELAY_S = 0.4

FRAGMENT_LENGTH = 10000

# Assemblies contain N-padded gaps. A fragment above this much ambiguity is
# rejected and retried elsewhere in the genome.
MAX_N_FRACTION = 0.01
MAX_FRAGMENT_ATTEMPTS = 4

# Representative taxa per community. NAMES ONLY, deliberately: these are
# resolved against RefSeq at build time so no accession is ever committed. The
# genus must appear in the FASTA description because
# add_bacterial_background() matches fragments to a community by genus
# substring; keep these consistent with BACTERIAL_COMMUNITY_PROFILES in
# viroforge/core/contamination.py.
#
# Taxonomy moves. If a name stops resolving it has probably been reclassified
# (Bacteroides vulgatus -> Phocaeicola vulgatus, for instance); the script
# reports unresolved taxa rather than quietly writing a short file.
COMMUNITY_TAXA: dict[str, list[str]] = {
    "gut": [
        "Bacteroides fragilis",
        "Bacteroides thetaiotaomicron",
        "Faecalibacterium prausnitzii",
        "Escherichia coli",
        "Bifidobacterium longum",
        "Blautia wexlerae",
        "Roseburia intestinalis",
        "Akkermansia muciniphila",
        "Phocaeicola vulgatus",
    ],
    "oral": [
        "Streptococcus mutans",
        "Neisseria subflava",
        "Veillonella parvula",
        "Haemophilus parainfluenzae",
        "Rothia mucilaginosa",
        "Fusobacterium nucleatum",
    ],
    "skin": [
        "Cutibacterium acnes",
        "Staphylococcus epidermidis",
        "Corynebacterium tuberculostearicum",
    ],
    "respiratory": [
        "Streptococcus pneumoniae",
        "Haemophilus influenzae",
        "Moraxella catarrhalis",
        "Pseudomonas aeruginosa",
        "Staphylococcus aureus",
    ],
    "vaginal": [
        "Lactobacillus crispatus",
        "Lactobacillus iners",
        "Gardnerella vaginalis",
        "Prevotella bivia",
    ],
    "urinary": [
        "Escherichia coli",
        "Lactobacillus gasseri",
        "Streptococcus anginosus",
    ],
    "marine": [
        "Prochlorococcus marinus",
        "Synechococcus sp. WH 8102",
        "Alteromonas macleodii",
    ],
    "soil": [
        "Bradyrhizobium japonicum",
        "Streptomyces coelicolor",
        "Mycobacterium smegmatis",
        "Pseudomonas putida",
    ],
    "freshwater": [
        "Polynucleobacter necessarius",
        "Flavobacterium johnsoniae",
    ],
    "wastewater": [
        "Nitrosomonas europaea",
        "Acinetobacter baumannii",
        "Pseudomonas aeruginosa",
    ],
}


def resolve_refseq_genome(taxon: str) -> tuple[str, int] | None:
    """Find a RefSeq genomic sequence for a taxon name.

    Returns (accession, length_bp), or None if nothing suitable is found. The
    accession is discovered here rather than stored anywhere, which is the point.
    """
    from Bio import Entrez

    Entrez.email = "viroforge@example.com"

    # Restrict to RefSeq genomic DNA and prefer complete genomes, so we get a
    # single long contiguous record to slice rather than a draft assembly.
    query = (
        f'"{taxon}"[Organism] AND srcdb_refseq[PROP] AND biomol_genomic[PROP] '
        f'AND ("complete genome"[Title] OR "chromosome"[Title])'
    )
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=5,
                                sort="SLEN")
        result = Entrez.read(handle)
        handle.close()
    except Exception as e:
        logger.warning(f"  esearch failed for {taxon}: {e}")
        return None

    ids = result.get("IdList", [])
    if not ids:
        logger.warning(
            f"  no RefSeq complete genome found for '{taxon}' "
            "(likely reclassified, or absent from RefSeq)"
        )
        return None

    time.sleep(REQUEST_DELAY_S)
    try:
        handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
    except Exception as e:
        logger.warning(f"  esummary failed for {taxon}: {e}")
        return None

    # Longest record wins: most likely the main chromosome rather than a plasmid.
    best = None
    for s in summaries:
        length = int(s.get("Length", 0))
        acc = s.get("AccessionVersion") or s.get("Caption")
        if acc and length > FRAGMENT_LENGTH and (best is None or length > best[1]):
            best = (str(acc), length)

    if best is None:
        logger.warning(f"  no record longer than {FRAGMENT_LENGTH} bp for {taxon}")
    return best


def _fetch_slice(accession: str, start: int, end: int):
    """Fetch one region without downloading the surrounding genome."""
    from Bio import Entrez, SeqIO

    Entrez.email = "viroforge@example.com"
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text",
        seq_start=start + 1,  # Entrez is 1-based
        seq_stop=end,
    )
    try:
        return SeqIO.read(handle, "fasta")
    finally:
        handle.close()


def _n_fraction(record) -> float:
    seq = str(record.seq).upper()
    return seq.count("N") / len(seq) if seq else 1.0


def fetch_fragments(taxon: str, accession: str, genome_length: int,
                    community: str, n_fragments: int) -> list:
    """Fetch evenly spaced, low-ambiguity fragments from one genome.

    Assemblies contain N-padded gaps, and some records are long precisely
    because of them: the RefSeq E. coli record picked during development was
    7.6 Mbp with a run of Ns that swallowed a whole fragment. N-heavy reference
    sequence would generate reads full of Ns and read downstream as a QC bug,
    so each position is retried at an offset until a usable slice is found.
    """
    records = []

    # Spread fragments across the genome so a single biased region (a prophage,
    # an rRNA operon) cannot dominate the reference set.
    usable = genome_length - FRAGMENT_LENGTH
    if usable <= 0:
        return records
    step = usable // max(1, n_fragments)

    for i in range(n_fragments):
        label = f"{taxon.replace(' ', '_')}_{i:02d}"
        record = None

        for attempt in range(MAX_FRAGMENT_ATTEMPTS):
            # Walk forward by a fraction of the step on each retry, staying
            # inside this fragment's slot so fragments remain well spread.
            offset = (attempt * step) // MAX_FRAGMENT_ATTEMPTS
            start = min(i * step + offset, usable)
            end = start + FRAGMENT_LENGTH
            time.sleep(REQUEST_DELAY_S)
            try:
                candidate = _fetch_slice(accession, start, end)
            except Exception as e:
                logger.warning(f"  failed fragment {label}: {e}")
                continue

            n_frac = _n_fraction(candidate)
            if n_frac <= MAX_N_FRACTION:
                candidate.id = label
                # Genus first so the runtime community match finds it;
                # accession recorded because it was retrieved, not stored.
                candidate.description = (
                    f"{taxon} [{community}] [RefSeq {accession}:{start + 1}-{end}]"
                )
                record = candidate
                break

            logger.info(
                f"  {label} at {start + 1} is {n_frac:.0%} N, retrying further along"
            )

        if record is None:
            logger.warning(
                f"  no low-ambiguity fragment for {label} after "
                f"{MAX_FRAGMENT_ATTEMPTS} attempts; skipping"
            )
            continue
        records.append(record)

    return records


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--communities", nargs="+", default=sorted(COMMUNITY_TAXA),
                    choices=sorted(COMMUNITY_TAXA),
                    help="communities to build (default: all)")
    ap.add_argument("--fragments-per-taxon", type=int, default=3,
                    help="10 kb fragments per genome (default: 3)")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve taxa to accessions and stop, without fetching")
    args = ap.parse_args()

    # De-duplicate taxa shared between communities (E. coli is in gut and
    # urinary) while remembering every community each one serves.
    taxon_communities: dict[str, list[str]] = {}
    for community in args.communities:
        for taxon in COMMUNITY_TAXA[community]:
            taxon_communities.setdefault(taxon, []).append(community)

    logger.info(
        f"Resolving {len(taxon_communities)} taxa across "
        f"{len(args.communities)} communities against RefSeq"
    )

    all_records = []
    resolved, unresolved = [], []

    for taxon, communities in sorted(taxon_communities.items()):
        logger.info(f"{taxon} ({', '.join(communities)})")
        hit = resolve_refseq_genome(taxon)
        time.sleep(REQUEST_DELAY_S)
        if hit is None:
            unresolved.append(taxon)
            continue
        accession, length = hit
        logger.info(f"  -> {accession} ({length:,} bp)")
        resolved.append((taxon, accession, length))

        if args.dry_run:
            continue

        all_records.extend(fetch_fragments(
            taxon, accession, length, communities[0], args.fragments_per_taxon))

    print()
    print(f"Resolved {len(resolved)}/{len(taxon_communities)} taxa")
    if unresolved:
        print(f"UNRESOLVED ({len(unresolved)}): {', '.join(unresolved)}")
        print("  Check for reclassification before treating these as missing "
              "from RefSeq.")

    if args.dry_run:
        print("\nDry run: no sequences fetched.")
        return 0 if resolved else 1

    if not all_records:
        logger.error("No fragments fetched; refusing to write an empty reference")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / OUTPUT_NAME

    from Bio import SeqIO
    SeqIO.write(all_records, out, "fasta")

    total_bp = sum(len(r.seq) for r in all_records)
    n_bases = sum(str(r.seq).upper().count("N") for r in all_records)
    print(f"\nWrote {len(all_records)} fragments ({total_bp:,} bp, "
          f"{out.stat().st_size / 1024:.0f} KB) to {out}")
    print(f"Ambiguous bases: {n_bases:,} ({n_bases / total_bp:.3%})")
    print("Compare: host_fragments.fasta is 48 fragments / ~528 KB.")

    if n_bases / total_bp > MAX_N_FRACTION:
        logger.error(
            "Output exceeds the ambiguity threshold; reads generated from this "
            "reference would be N-heavy. Not usable as-is."
        )
        return 1

    if unresolved:
        print("\nSome taxa did not resolve, so coverage is incomplete. Rerun "
              "after correcting the names in COMMUNITY_TAXA.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
