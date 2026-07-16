"""Taxonomy benchmarking metrics engine (Module 4, read-based v1).

Compares a read classifier's per-read NCBI taxid assignments against ViroForge
ground truth. The classifier's read IDs are the ViroForge read names, which embed
the source genome accession, so ground truth comes from the read ID plus the
per-genome taxonomy exported in the metadata (`benchmarking.taxonomy`) - no raw
reads or 500 MB database needed.

v1 scope: taxid-exact scoring (species/exact level; a higher-rank LCA assignment
counts as not-exact - rank-level metrics are a later addition) plus an abundance
profile comparison. Reads are stratified into known viruses (ICTV family assigned,
likely in reference databases) and dark/unknown-family viruses (viral dark matter),
because a classifier leaving dark matter unclassified is expected, not a failure.
"""

from __future__ import annotations

import re
from collections import Counter

_DUP = re.compile(r"_dup\d+$")
_TAXID_IN_NAME = re.compile(r"taxid\s*(\d+)")


def parse_genome_id(read_id: str) -> str:
    """Recover the ViroForge source genome id from a read name."""
    return _DUP.sub("", read_id).rsplit("_", 2)[0]


def _to_taxid(token: str):
    token = token.strip()
    if token.isdigit():
        t = int(token)
        return t if t != 0 else None
    m = _TAXID_IN_NAME.search(token)  # kraken2 --use-names: "Name (taxid 123)"
    return int(m.group(1)) if m else None


def parse_kraken2(path) -> dict[str, int | None]:
    """Kraken2 per-read output -> {read_id: taxid or None}.

    Columns: C/U, read_id, taxid (or 'name (taxid N)'), length, LCA map.
    """
    out: dict[str, int | None] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status, rid, taxid = parts[0], parts[1], parts[2]
            out[rid] = None if status == "U" else _to_taxid(taxid)
    return out


def parse_generic(path) -> dict[str, int | None]:
    """Generic TSV -> {read_id: taxid or None}. Columns: read_id, taxid[, ...]."""
    out: dict[str, int | None] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            out[parts[0]] = _to_taxid(parts[1])
    return out


PARSERS = {"kraken2": parse_kraken2, "generic": parse_generic}


def _bray_curtis(a: dict, b: dict) -> float:
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)
    den = sum(a.get(k, 0.0) + b.get(k, 0.0) for k in keys)
    return num / den if den else 0.0


def _pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / (vx * vy) ** 0.5 if vx and vy else None


def _spearman(xs, ys):
    if len(xs) < 3:
        return None
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


def _stratum_metrics(s: dict) -> dict:
    total = s["correct"] + s["unclassified"] + s["misclassified"]
    classified = s["correct"] + s["misclassified"]
    return {
        "n": total,
        "correct": s["correct"],
        "unclassified": s["unclassified"],
        "misclassified": s["misclassified"],
        "sensitivity": s["correct"] / total if total else None,       # correct / all
        "precision": s["correct"] / classified if classified else None,  # correct / classified
        "unclassified_rate": s["unclassified"] / total if total else None,
    }


RANKS = ("species", "genus", "family")


def _per_rank_metrics(items, tree, ranks) -> dict:
    """Per-rank precision/recall/F1 over the known-virus stratum.

    items are (assigned_taxid, true_taxid, is_known) tuples (one per read or
    contig). For each rank, a unit counts only if its true taxid has an ancestor
    at that rank. The assigned taxid is resolved to its ancestor at the rank;
    correct when it equals the true rank taxid, misclassified when it differs,
    unclassified when the classifier made no call or its call was too shallow to
    reach the rank. Dark matter is excluded, so correctly-unclassified novel
    content is not penalized.
    """
    stats = {r: {"correct": 0, "misclassified": 0, "unclassified": 0} for r in ranks}
    for assigned, true_taxid, is_known in items:
        if not is_known:
            continue
        for r in ranks:
            true_r = tree.rank_taxid(true_taxid, r)
            if true_r is None:
                continue  # ground truth has no taxid at this rank; can't score
            s = stats[r]
            if assigned is None:
                s["unclassified"] += 1
                continue
            asg_r = tree.rank_taxid(assigned, r)
            if asg_r is None:
                s["unclassified"] += 1  # call too shallow to reach this rank
            elif asg_r == true_r:
                s["correct"] += 1
            else:
                s["misclassified"] += 1
    out = {}
    for r, s in stats.items():
        total = s["correct"] + s["misclassified"] + s["unclassified"]
        classified = s["correct"] + s["misclassified"]
        precision = s["correct"] / classified if classified else None
        recall = s["correct"] / total if total else None
        f1 = (2 * precision * recall / (precision + recall)
              if precision and recall else None)
        out[r] = {**s, "n": total, "precision": precision, "recall": recall, "f1": f1}
    return out


def benchmark_taxonomy(
    assignments: dict[str, int | None],
    taxonomy_gt: dict,
    ncbi_tree=None,
    ranks=RANKS,
) -> dict:
    """Score per-read taxid assignments against the ground-truth taxonomy.

    Args:
        assignments: {read_id: assigned_taxid or None} from a classifier.
        taxonomy_gt: {genome_id: {ncbi_taxid, is_known, ...}} from the metadata.
        ncbi_tree: optional NcbiTree; enables per-rank (genus/family) metrics.
        ranks: taxonomic ranks to score when a tree is given.
    """
    items = []
    non_viral = 0
    for rid, assigned in assignments.items():
        gt = taxonomy_gt.get(parse_genome_id(rid))
        if gt is None:
            non_viral += 1
            continue
        items.append((assigned, gt.get("ncbi_taxid"), gt.get("is_known")))

    core = _score(items, ncbi_tree, ranks)
    return {
        "reliable": core["reliable"],
        "n_assignments": len(assignments),
        "n_viral_reads": core["n_viral"],
        "n_non_viral_reads": non_viral,
        "known_viruses": core["known_viruses"],
        "dark_matter": core["dark_matter"],
        "per_rank": core["per_rank"],
        "abundance_profile": core["abundance_profile"],
    }


def _score(items, ncbi_tree, ranks) -> dict:
    """Score (assigned_taxid, true_taxid, is_known) items into strata + profile.

    Shared by read-based and contig-based modes; the unit (read or contig) is
    whatever produced the items.
    """
    known = {"correct": 0, "unclassified": 0, "misclassified": 0}
    dark = {"correct": 0, "unclassified": 0, "misclassified": 0}
    true_counts: Counter = Counter()
    obs_counts: Counter = Counter()

    for assigned, true_taxid, is_known in items:
        stratum = known if is_known else dark
        true_counts[true_taxid] += 1
        if assigned is None:
            stratum["unclassified"] += 1
        else:
            obs_counts[assigned] += 1
            if assigned == true_taxid:
                stratum["correct"] += 1
            else:
                stratum["misclassified"] += 1

    n_viral = sum(true_counts.values())
    true_prof = {t: c / n_viral for t, c in true_counts.items()} if n_viral else {}
    obs_prof = {t: c / n_viral for t, c in obs_counts.items()} if n_viral else {}
    taxa = sorted(set(true_prof) | set(obs_prof), key=lambda t: (t is None, t))
    xs = [true_prof.get(t, 0.0) for t in taxa]
    ys = [obs_prof.get(t, 0.0) for t in taxa]
    return {
        "reliable": n_viral > 0,
        "n_viral": n_viral,
        "known_viruses": _stratum_metrics(known),
        "dark_matter": _stratum_metrics(dark),
        "per_rank": _per_rank_metrics(items, ncbi_tree, ranks) if ncbi_tree else None,
        "abundance_profile": {
            "n_taxa": len(taxa),
            "bray_curtis": _bray_curtis(true_prof, obs_prof) if n_viral else None,
            "pearson": _pearson(xs, ys),
            "spearman": _spearman(xs, ys),
            "mean_abs_error": (sum(abs(a - b) for a, b in zip(xs, ys)) / len(taxa)
                               if taxa else None),
        },
    }


def _lca(tree, a, b):
    """Lowest common ancestor taxid of two taxids (None if either is None)."""
    if a is None or b is None:
        return None
    b_anc = set(tree.path_to_root(b))
    for t in tree.path_to_root(a):
        if t in b_anc:
            return t
    return None


def benchmark_taxonomy_contigs(
    contig_assignments: dict[str, int | None],
    contigs_fasta,
    genome_fasta,
    taxonomy_gt: dict,
    ncbi_tree=None,
    chimera_handling: str = "exclude",
    ranks=RANKS,
) -> dict:
    """Contig-based taxonomy benchmark.

    Ground truth per contig is derived by aligning contigs to the true genomes
    (shared aligner): the contig's primary genome gives its true taxid. Contigs
    that align to no genome are novel; contigs aligning to a contaminant genome
    are excluded. Chimeric contigs are handled per `chimera_handling`:
    "exclude" (reported, not scored) or "lca" (true taxid = LCA of the two segment
    genomes, requires an NCBI tree).
    """
    from .align import align_contigs

    al = align_contigs(contigs_fasta, genome_fasta)
    items = []
    n_novel = n_contaminant = n_chimera = 0
    for c in al["contigs"]:
        assigned = contig_assignments.get(c["name"])  # None if classifier omitted it
        primary = c["primary_genome"]
        if primary is None:
            n_novel += 1
            continue
        gt = taxonomy_gt.get(primary)
        if gt is None:
            n_contaminant += 1
            continue
        if c["is_chimera"]:
            n_chimera += 1
            if chimera_handling == "lca" and ncbi_tree is not None:
                segs = c["chimera_segments"]
                t1 = taxonomy_gt.get(segs[0]["genome"], {}).get("ncbi_taxid")
                t2 = taxonomy_gt.get(segs[1]["genome"], {}).get("ncbi_taxid")
                items.append((assigned, _lca(ncbi_tree, t1, t2), gt.get("is_known")))
            # "exclude" (or lca without a tree): reported but not scored
            continue
        items.append((assigned, gt.get("ncbi_taxid"), gt.get("is_known")))

    core = _score(items, ncbi_tree, ranks)
    return {
        "reliable": core["reliable"],
        "n_contigs": len(al["contigs"]),
        "n_viral_contigs": core["n_viral"],
        "n_novel_contigs": n_novel,
        "n_contaminant_contigs": n_contaminant,
        "n_chimeric_contigs": n_chimera,
        "chimera_handling": chimera_handling,
        "known_viruses": core["known_viruses"],
        "dark_matter": core["dark_matter"],
        "per_rank": core["per_rank"],
        "abundance_profile": core["abundance_profile"],
    }
