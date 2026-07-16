"""Minimal NCBI taxonomy tree for rank-level taxonomy benchmarking.

Parses the NCBI taxdump (nodes.dmp, optionally names.dmp) to resolve any taxid to
its ancestor taxid at a named rank (species, genus, family, ...). This is what
higher-rank scoring needs: to know the family/genus of a classifier's assigned
taxid. No external dependency; the format is the standard NCBI taxdump.

Users already have the taxdump if they run Kraken2/Centrifuge. Download:
ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz (nodes.dmp, names.dmp).
"""

from __future__ import annotations

from pathlib import Path


class NcbiTree:
    def __init__(self, nodes_dmp, names_dmp=None):
        self.parent: dict[int, int] = {}
        self.rank: dict[int, str] = {}
        self._lineage_cache: dict[int, dict[str, int]] = {}
        with open(nodes_dmp) as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                if len(f) < 3:
                    continue
                tid, parent, rank = int(f[0]), int(f[1]), f[2]
                self.parent[tid] = parent
                self.rank[tid] = rank
        self.names: dict[int, str] = {}
        if names_dmp:
            with open(names_dmp) as fh:
                for line in fh:
                    f = [x.strip() for x in line.split("|")]
                    if len(f) >= 4 and f[3] == "scientific name":
                        self.names[int(f[0])] = f[1]

    @classmethod
    def from_dir(cls, taxdump_dir) -> "NcbiTree":
        d = Path(taxdump_dir)
        names = d / "names.dmp"
        return cls(d / "nodes.dmp", names if names.exists() else None)

    def lineage(self, taxid: int | None) -> dict[str, int]:
        """{rank: taxid} for the taxid's ancestors (nearest rank wins)."""
        if taxid is None:
            return {}
        if taxid in self._lineage_cache:
            return self._lineage_cache[taxid]
        out: dict[str, int] = {}
        t = taxid
        seen: set[int] = set()
        while t and t not in seen:
            seen.add(t)
            r = self.rank.get(t)
            if r and r != "no rank":
                out.setdefault(r, t)
            p = self.parent.get(t)
            if p is None or p == t:
                break
            t = p
        self._lineage_cache[taxid] = out
        return out

    def rank_taxid(self, taxid: int | None, rank: str) -> int | None:
        """The taxid's ancestor at `rank`, or None if it has none at that rank."""
        return self.lineage(taxid).get(rank)

    def path_to_root(self, taxid: int | None) -> list[int]:
        """Ordered ancestor taxids from `taxid` up to the root (inclusive)."""
        out: list[int] = []
        t = taxid
        seen: set[int] = set()
        while t and t not in seen:
            seen.add(t)
            out.append(t)
            p = self.parent.get(t)
            if p is None or p == t:
                break
            t = p
        return out
