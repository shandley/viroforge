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

2. **Fragments, not genomes.** Bacterial genomes are 3-6 Mbp. Only 50 kb slices
   are fetched, using Entrez ``seq_start``/``seq_stop``, so the whole genome is
   never downloaded or stored. 50 kb is long enough for realistic contig-length
   contamination and to naturally include prophage/CRISPR/ICE regions.

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
OUTPUT_NAME = "bacterial_fragments.fasta"  # overridden per domain

# NCBI asks for no more than 3 requests/second without an API key.
REQUEST_DELAY_S = 0.4

FRAGMENT_LENGTH = 50_000

# Assemblies contain N-padded gaps. A fragment above this much ambiguity is
# rejected and retried elsewhere in the genome.
MAX_N_FRACTION = 0.01
MAX_FRAGMENT_ATTEMPTS = 4

# Organelle and extrachromosomal records are not representative genomic
# background: they have their own base composition and gene content, and are
# tiny. The Trichoderma reesei search returned its 42 kb MITOCHONDRION as the
# top hit, which would have supplied only four non-overlapping fragments of
# atypical sequence.
EXCLUDED_TITLE_TERMS = (
    "mitochondrion", "mitochondrial", "chloroplast", "plastid",
    "apicoplast", "plasmid",
)

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
        # S. elongatus rather than a marine strain such as WH 8102: the
        # strain-qualified names do not resolve through an [Organism] search.
        # The genus and its GC profile are what matter for generating
        # background sequence, so the ecological mismatch is acceptable here.
        "Synechococcus elongatus",
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


# Mycobiome taxa. Fungal genomes are 12-40 Mbp and usually deposited as many
# per-chromosome records, which is fine: only 10 kb slices are fetched either way.
FUNGAL_TAXA: dict[str, list[str]] = {
    "gut": ["Candida albicans", "Saccharomyces cerevisiae"],
    "oral": ["Candida albicans"],
    "skin": ["Malassezia restricta"],
    "respiratory": ["Aspergillus fumigatus", "Candida albicans"],
    "vaginal": ["Candida albicans"],
    # Trichoderma reesei is deposited only as scaffolds plus a
    # mitochondrion, so Fusarium stands in as the second soil fungus.
    "soil": ["Aspergillus fumigatus", "Fusarium graminearum"],
    "wastewater": ["Candida albicans", "Aspergillus fumigatus"],
}

# Archaeome taxa. Methanogens in gut and wastewater; ammonia-oxidising
# Thaumarchaeota in soil and open ocean. Most body sites carry no appreciable
# archaeal population and are absent here rather than guessed at.
ARCHAEAL_TAXA: dict[str, list[str]] = {
    "gut": ["Methanobrevibacter smithii", "Methanosphaera stadtmanae"],
    "oral": ["Methanobrevibacter smithii"],
    "soil": ["Nitrososphaera viennensis"],
    "marine": ["Nitrosopumilus piranensis"],
    "wastewater": ["Methanosarcina barkeri", "Methanobrevibacter smithii"],
}

# Dietary background taxa. Food-derived plant viruses that survive digestion
# and appear in gut metagenomes. PMMoV is the most abundant RNA virus in human
# stool (Zhang et al. 2006). These genomes are small (5-10 kb), so whole
# genomes are fetched rather than fragments (min_length=1000 in DOMAINS).
# The community types here match DIETARY_COMMUNITY_PROFILES in
# contamination.py: all diets get Tobamovirus; plant-heavy diets get
# additional genera. CMV (Cucumovirus) was replaced with TBSV (Tombusvirus)
# because CMV is tripartite (3 RNA segments) and doesn't resolve as a single
# "complete genome" entry in RefSeq.
DIETARY_TAXA: dict[str, list[str]] = {
    "western": [
        "Pepper mild mottle virus",
        "Tobacco mosaic virus",
    ],
    "mediterranean": [
        "Pepper mild mottle virus",
        "Tobacco mosaic virus",
        "Potato virus Y",
    ],
    "plant_based": [
        "Pepper mild mottle virus",
        "Tobacco mosaic virus",
        "Potato virus Y",
        "Tomato bushy stunt virus",
    ],
    "high_fiber": [
        "Pepper mild mottle virus",
        "Potato virus Y",
    ],
}

# Curated "challenge regions" — bacterial genomic loci that cause false-positive
# viral calls in virome pipelines. Each entry is fetched as a single fragment
# (accession + start + length) and appended to the bacterial reference FASTA
# with a distinguishing prefix so benchmarking can track them.
#
# Three categories:
#   prophage_ : integrated phage genomes; assemblers produce chimeric contigs
#               crossing prophage-host junctions
#   crispr_   : CRISPR arrays containing phage-derived spacers; taxonomy
#               classifiers see viral k-mers in bacterial contigs
#   ice_      : integrative conjugative elements with phage-like genes
#               (integrases, structural proteins); misannotated as viral
#
# Coordinates are from RefSeq complete genomes. Each fragment is large enough
# (40-50 kb) to capture the element plus flanking bacterial sequence.
CHALLENGE_REGIONS: list[dict] = [
    # --- Prophage-containing regions (40 kb) ---
    # E. coli K-12 MG1655 prophages (U00096.3, 4,641,652 bp)
    # Coordinates from Wang et al. 2010 Nat Commun; 9 cryptic prophages
    {"accession": "U00096.3", "start": 242_000, "length": 40_000,
     "label": "prophage_ecoli_CP4-6", "community": "gut",
     "desc": "Escherichia coli K-12 CP4-6 prophage (34 kb lambdoid)"},
    {"accession": "U00096.3", "start": 544_000, "length": 40_000,
     "label": "prophage_ecoli_DLP12", "community": "gut",
     "desc": "Escherichia coli K-12 DLP12 defective prophage (21 kb)"},
    {"accession": "U00096.3", "start": 1_185_000, "length": 40_000,
     "label": "prophage_ecoli_e14", "community": "gut",
     "desc": "Escherichia coli K-12 e14 prophage (15 kb lambdoid)"},
    {"accession": "U00096.3", "start": 1_390_000, "length": 40_000,
     "label": "prophage_ecoli_Rac", "community": "gut",
     "desc": "Escherichia coli K-12 Rac prophage (23 kb, RecE/RecT)"},
    {"accession": "U00096.3", "start": 1_610_000, "length": 40_000,
     "label": "prophage_ecoli_Qin", "community": "gut",
     "desc": "Escherichia coli K-12 Qin prophage (20 kb, lysis/tail genes)"},
    {"accession": "U00096.3", "start": 2_736_000, "length": 40_000,
     "label": "prophage_ecoli_CP4-57", "community": "gut",
     "desc": "Escherichia coli K-12 CP4-57 prophage (22 kb)"},
    # S. aureus NCTC 8325 Sa3/phi13 prophage (NC_007795.1, 2,821,361 bp)
    {"accession": "NC_007795.1", "start": 460_000, "length": 40_000,
     "label": "prophage_saureus_phi13", "community": "respiratory",
     "desc": "Staphylococcus aureus NCTC 8325 Sa3/phi13 immune evasion prophage"},
    # P. aeruginosa PAO1 Pf4 filamentous prophage (NC_002516.2, 6,264,404 bp)
    # Pf4 at 785,311-797,747 per RefSeq annotation
    {"accession": "NC_002516.2", "start": 775_000, "length": 40_000,
     "label": "prophage_paerug_Pf4", "community": "respiratory",
     "desc": "Pseudomonas aeruginosa PAO1 Pf4 filamentous prophage"},

    # --- CRISPR-rich regions (50 kb) ---
    # E. coli K-12 Type I-E CRISPR-Cas (cas3 at 2,884,553)
    {"accession": "U00096.3", "start": 2_857_000, "length": 50_000,
     "label": "crispr_ecoli_IE", "community": "gut",
     "desc": "Escherichia coli K-12 Type I-E CRISPR-Cas + spacer array"},
    # S. thermophilus LMD-9 CRISPR3 (CP000419.1, 1,856,368 bp)
    {"accession": "CP000419.1", "start": 1_352_000, "length": 50_000,
     "label": "crispr_stherm_CR3", "community": "gut",
     "desc": "Streptococcus thermophilus LMD-9 CRISPR3 + cas genes"},
    # B. fragilis NCTC 9343 Type IB CRISPR (NC_003228.3, 5,205,140 bp)
    {"accession": "NC_003228.3", "start": 4_210_000, "length": 50_000,
     "label": "crispr_bfrag_IB", "community": "gut",
     "desc": "Bacteroides fragilis NCTC 9343 Type IB CRISPR-Cas"},
    # B. fragilis NCTC 9343 prophage region (also tests phage in gut anaerobe)
    {"accession": "NC_003228.3", "start": 2_710_000, "length": 50_000,
     "label": "prophage_bfrag_BF01", "community": "gut",
     "desc": "Bacteroides fragilis NCTC 9343 prophage region"},

    # --- ICE / defense island regions (50 kb) ---
    # B. thetaiotaomicron VPI-5482 CTnDOT-1 (NC_004663.1, 6,260,361 bp)
    # BT-CTn1 at 63-133 kb per comparative genomic studies
    {"accession": "NC_004663.1", "start": 63_000, "length": 50_000,
     "label": "ice_btheta_CTn1", "community": "gut",
     "desc": "Bacteroides thetaiotaomicron BT-CTn1 conjugative transposon"},
    # B. thetaiotaomicron VPI-5482 CTnDOT-2 at 2850-2930 kb
    {"accession": "NC_004663.1", "start": 2_850_000, "length": 50_000,
     "label": "ice_btheta_CTn2", "community": "gut",
     "desc": "Bacteroides thetaiotaomicron BT-CTn2 conjugative transposon"},
    # E. coli K-12 CPS-53/KpLE1 defense island (2,464-2,475 kb region)
    {"accession": "U00096.3", "start": 2_434_000, "length": 50_000,
     "label": "ice_ecoli_CPS53", "community": "gut",
     "desc": "Escherichia coli K-12 CPS-53/KpLE1 prophage remnant + RM cassette"},
]

DOMAINS: dict[str, dict] = {
    "bacterial": {"taxa": COMMUNITY_TAXA, "output": "bacterial_fragments.fasta"},
    "fungal": {"taxa": FUNGAL_TAXA, "output": "fungal_fragments.fasta"},
    "archaeal": {"taxa": ARCHAEAL_TAXA, "output": "archaeal_fragments.fasta"},
    "dietary": {"taxa": DIETARY_TAXA, "output": "dietary_fragments.fasta",
                "min_length": 1000},
}


def resolve_refseq_genome(taxon: str, min_length: int | None = None) -> tuple[str, int] | None:
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

    # Longest non-organelle record wins: most likely the main chromosome.
    # For small-genome domains (plant viruses ~6-10 kb), callers pass a lower
    # min_length so the resolver doesn't reject genomes shorter than FRAGMENT_LENGTH.
    threshold = min_length if min_length is not None else FRAGMENT_LENGTH
    best = None
    for s in summaries:
        title = (s.get("Title") or "").lower()
        if any(term in title for term in EXCLUDED_TITLE_TERMS):
            logger.info(f"  skipping {s.get('AccessionVersion')}: {title[:60]}")
            continue
        length = int(s.get("Length", 0))
        acc = s.get("AccessionVersion") or s.get("Caption")
        if acc and length > threshold and (best is None or length > best[1]):
            best = (str(acc), length)

    if best is None:
        logger.warning(f"  no record longer than {threshold} bp for {taxon}")
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


def fetch_challenge_regions() -> list:
    """Fetch curated challenge-region fragments from CHALLENGE_REGIONS.

    These are prophage, CRISPR, and ICE loci at specific coordinates in
    well-characterized RefSeq genomes. Unlike the generic fragments (evenly
    spaced), these target known problem regions for virome pipelines.
    """
    records = []
    for entry in CHALLENGE_REGIONS:
        label = entry["label"]
        accession = entry["accession"]
        start = entry["start"]
        length = entry["length"]
        end = start + length

        logger.info(f"  {label} from {accession}:{start + 1}-{end}")
        time.sleep(REQUEST_DELAY_S)
        try:
            record = _fetch_slice(accession, start, end)
        except Exception as e:
            logger.warning(f"  FAILED {label}: {e}")
            continue

        n_frac = _n_fraction(record)
        if n_frac > MAX_N_FRACTION:
            logger.warning(f"  {label} is {n_frac:.0%} N; skipping")
            continue

        record.id = label
        record.description = (
            f"{entry['desc']} [{entry['community']}] "
            f"[RefSeq {accession}:{start + 1}-{end}]"
        )
        records.append(record)

    return records


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--domain", default="bacterial", choices=sorted(DOMAINS),
                    help="which background to build (default: bacterial)")
    ap.add_argument("--communities", nargs="+", default=None,
                    help="communities to build (default: all for this domain)")
    ap.add_argument("--fragments-per-taxon", type=int, default=3,
                    help="50 kb fragments per genome (default: 3)")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve taxa to accessions and stop, without fetching")
    args = ap.parse_args()

    domain = DOMAINS[args.domain]
    taxa_map = domain["taxa"]
    if args.communities is None:
        args.communities = sorted(taxa_map)
    unknown = [c for c in args.communities if c not in taxa_map]
    if unknown:
        ap.error(f"{args.domain} has no taxa for {unknown}; "
                 f"available: {sorted(taxa_map)}")

    # De-duplicate taxa shared between communities (E. coli is in gut and
    # urinary) while remembering every community each one serves.
    taxon_communities: dict[str, list[str]] = {}
    for community in args.communities:
        for taxon in taxa_map[community]:
            taxon_communities.setdefault(taxon, []).append(community)

    logger.info(
        f"Resolving {len(taxon_communities)} {args.domain} taxa across "
        f"{len(args.communities)} communities against RefSeq"
    )

    all_records = []
    resolved, unresolved = [], []
    min_length = domain.get("min_length")

    for taxon, communities in sorted(taxon_communities.items()):
        logger.info(f"{taxon} ({', '.join(communities)})")
        hit = resolve_refseq_genome(taxon, min_length=min_length)
        time.sleep(REQUEST_DELAY_S)
        if hit is None:
            unresolved.append(taxon)
            continue
        accession, length = hit
        logger.info(f"  -> {accession} ({length:,} bp)")
        resolved.append((taxon, accession, length))

        if args.dry_run:
            continue

        if length <= FRAGMENT_LENGTH:
            # Small genomes (plant viruses, ~6-10 kb): fetch the whole sequence
            # rather than trying to cut fragments longer than the genome itself.
            time.sleep(REQUEST_DELAY_S)
            try:
                rec = _fetch_slice(accession, 0, length)
                rec.id = taxon.replace(" ", "_")
                rec.description = (
                    f"{taxon} [{communities[0]}] [RefSeq {accession}:1-{length}]"
                )
                if _n_fraction(rec) <= MAX_N_FRACTION:
                    all_records.append(rec)
                    logger.info(f"  fetched whole genome ({length:,} bp)")
                else:
                    logger.warning(f"  {accession} is {_n_fraction(rec):.0%} N; skipping")
            except Exception as e:
                logger.warning(f"  failed to fetch {accession}: {e}")
        else:
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

    # Append challenge regions (prophage, CRISPR, ICE) for bacterial domain
    n_challenge = 0
    if args.domain == "bacterial" and not args.dry_run:
        logger.info(f"\nFetching {len(CHALLENGE_REGIONS)} challenge-region fragments "
                     "(prophage, CRISPR, ICE)...")
        challenge = fetch_challenge_regions()
        all_records.extend(challenge)
        n_challenge = len(challenge)
        logger.info(f"  {n_challenge}/{len(CHALLENGE_REGIONS)} challenge regions fetched")

    if not all_records:
        logger.error("No fragments fetched; refusing to write an empty reference")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / domain["output"]

    from Bio import SeqIO
    SeqIO.write(all_records, out, "fasta")

    total_bp = sum(len(r.seq) for r in all_records)
    n_bases = sum(str(r.seq).upper().count("N") for r in all_records)
    n_generic = len(all_records) - n_challenge
    print(f"\nWrote {len(all_records)} fragments ({total_bp:,} bp, "
          f"{out.stat().st_size / 1024:.0f} KB) to {out}")
    if n_challenge:
        print(f"  {n_generic} generic (50 kb) + {n_challenge} challenge regions "
              f"(40-50 kb prophage/CRISPR/ICE)")
    print(f"Ambiguous bases: {n_bases:,} ({n_bases / total_bp:.3%})")
    print("Compare: host_fragments.fasta is 58 fragments / ~696 KB.")

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
