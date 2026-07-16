#!/usr/bin/env python3
"""Populate missing higher taxonomy ranks from NCBI for genomes ViroForge's ICTV
name-match left completely unclassified.

Many RefSeq genomes (mostly modern bacteriophages) failed ICTV name-matching and
sit as Unknown at EVERY rank in the DB, even though NCBI classifies them to
realm/class/genus. NCBI (ICTV-synced) legitimately assigns NO family to most
Caudoviricetes genera - ICTV abolished the morphology-based phage families in
2021 - so this does NOT invent families; it fills the ranks NCBI actually provides
(realm..genus) and only sets family where NCBI has a real family (e.g. the handful
of eukaryotic viruses like Human bocavirus -> Parvoviridae that also failed the
name-match).

Two phases:
  --fetch   query NCBI E-utilities (batched efetch, taxonomy) for the target
            taxids and write taxid -> lineage to a tracked cache TSV.
  --apply   fill EMPTY taxonomy ranks in the DB from the cache. Never overwrites a
            non-empty rank and never touches `species`. Idempotent.

Usage:
  python3 scripts/enrich_taxonomy_from_ncbi.py --fetch
  python3 scripts/enrich_taxonomy_from_ncbi.py --apply           # audit
  python3 scripts/enrich_taxonomy_from_ncbi.py --apply --write    # write DB
"""
from __future__ import annotations

import argparse
import sqlite3
import time
import urllib.request
from pathlib import Path

RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus"]
CACHE = Path("data/reference_profiles/ncbi_lineage_cache.tsv")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def target_taxids(db: Path) -> list[str]:
    """Genomes that are Unknown at every rank (nothing to lose by filling)."""
    con = sqlite3.connect(str(db))
    rows = con.execute(
        """SELECT DISTINCT ncbi_taxid FROM taxonomy
           WHERE family='Unknown'
             AND (class IS NULL OR class IN ('Unknown',''))
             AND ncbi_taxid IS NOT NULL AND ncbi_taxid != ''"""
    ).fetchall()
    con.close()
    return [str(r[0]) for r in rows]


def _parse_efetch(xml: str) -> dict[str, dict]:
    """Map every requested taxid (incl. AkaTaxIds) -> {rank: name} from an efetch
    taxonomy response, using a real XML parser (nested <Taxon> in <LineageEx>)."""
    import xml.etree.ElementTree as ET

    out: dict[str, dict] = {}
    root = ET.fromstring(xml)
    for taxon in root.findall("Taxon"):  # top-level records only
        primary = taxon.findtext("TaxId")
        if not primary:
            continue
        aka_ids = [e.text for e in taxon.findall("AkaTaxIds/TaxId") if e.text]

        lineage: dict[str, str] = {}
        for anc in taxon.findall("LineageEx/Taxon"):
            rank = anc.findtext("Rank")
            name = anc.findtext("ScientificName")
            if rank in RANKS and name:
                lineage[rank] = name
        own_rank = taxon.findtext("Rank")
        own_name = taxon.findtext("ScientificName")
        if own_rank in RANKS and own_name:
            lineage.setdefault(own_rank, own_name)

        for k in [primary, *aka_ids]:
            out[k] = lineage
    return out


def fetch(taxids: list[str], batch: int = 180) -> dict[str, dict]:
    resolved: dict[str, dict] = {}
    for i in range(0, len(taxids), batch):
        chunk = taxids[i:i + batch]
        url = f"{EFETCH}?db=taxonomy&id={','.join(chunk)}&retmode=xml"
        with urllib.request.urlopen(url, timeout=60) as resp:
            xml = resp.read().decode("utf-8", "replace")
        resolved.update(_parse_efetch(xml))
        print(f"  fetched {min(i + batch, len(taxids))}/{len(taxids)} taxids "
              f"({sum(1 for t in chunk if t in resolved)}/{len(chunk)} resolved in batch)")
        time.sleep(0.4)  # be polite to E-utilities
    return resolved


def write_cache(resolved: dict[str, dict]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "w") as fh:
        fh.write("ncbi_taxid\t" + "\t".join(RANKS) + "\n")
        for tid, lin in sorted(resolved.items(), key=lambda kv: int(kv[0])):
            fh.write(tid + "\t" + "\t".join(lin.get(r, "") for r in RANKS) + "\n")
    print(f"cache -> {CACHE} ({len(resolved)} taxids)")


def load_cache() -> dict[str, dict]:
    import csv
    m: dict[str, dict] = {}
    with open(CACHE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m[row["ncbi_taxid"]] = {r: row[r] for r in RANKS if row.get(r)}
    return m


def cmd_fetch(args) -> None:
    taxids = target_taxids(args.db)
    print(f"{len(taxids)} unclassified taxids to resolve from NCBI")
    resolved = fetch(taxids)
    got_family = sum(1 for l in resolved.values() if l.get("family"))
    got_class = sum(1 for l in resolved.values() if l.get("class"))
    print(f"resolved: {len(resolved)} taxids | with class: {got_class} | with family: {got_family}")
    write_cache(resolved)


def cmd_apply(args) -> None:
    cache = load_cache()
    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    rows = con.execute(
        """SELECT genome_id, ncbi_taxid, realm, kingdom, phylum, class, order_name, family, genus
           FROM taxonomy WHERE family='Unknown'
             AND (class IS NULL OR class IN ('Unknown',''))
             AND ncbi_taxid IS NOT NULL AND ncbi_taxid != ''"""
    ).fetchall()

    col = {"order": "order_name"}  # DB column name for the 'order' rank
    updates = []
    stats = {"genomes": 0, "class_set": 0, "genus_set": 0, "family_set": 0, "no_class_after": 0}
    for r in rows:
        lin = cache.get(str(r["ncbi_taxid"]))
        if not lin:
            stats["no_class_after"] += 1
            continue
        sets = {}
        for rank in RANKS:
            dbcol = col.get(rank, rank)
            cur = r[dbcol]
            if (cur is None or cur in ("Unknown", "")) and lin.get(rank):
                sets[dbcol] = lin[rank]
        if not sets:
            continue
        stats["genomes"] += 1
        if "class" in sets:
            stats["class_set"] += 1
        if "genus" in sets:
            stats["genus_set"] += 1
        if "family" in sets:
            stats["family_set"] += 1
        updates.append((sets, r["genome_id"]))
    stats["no_class_after"] += sum(1 for r in rows if not cache.get(str(r["ncbi_taxid"]), {}).get("class"))

    print(f"target genomes: {len(rows)}")
    print(f"  would fill some rank for: {stats['genomes']}")
    print(f"  class set: {stats['class_set']} | genus set: {stats['genus_set']} | NEW family set: {stats['family_set']}")
    print(f"  still no class after enrichment (genuinely unclassified): {stats['no_class_after']}")

    if not args.write:
        print("\naudit only. re-run with --write to update the DB.")
        con.close()
        return

    for sets, gid in updates:
        assigns = ", ".join(f"{c}=?" for c in sets)
        con.execute(f"UPDATE taxonomy SET {assigns} WHERE genome_id=?", (*sets.values(), gid))
    con.commit()
    con.close()
    print(f"\napplied rank fills to {len(updates)} genomes.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    ap.add_argument("--fetch", action="store_true", help="query NCBI and write the cache TSV")
    ap.add_argument("--apply", action="store_true", help="apply cached lineage to the DB (audit unless --write)")
    ap.add_argument("--write", action="store_true", help="with --apply, actually write the DB")
    args = ap.parse_args()
    if args.fetch:
        cmd_fetch(args)
    if args.apply:
        cmd_apply(args)
    if not (args.fetch or args.apply):
        ap.error("pass --fetch and/or --apply")


if __name__ == "__main__":
    main()
