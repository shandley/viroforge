#!/usr/bin/env python3
"""Correct the viral composition of collections the biological-accuracy review
flagged as genuinely off (docs/BIOLOGICAL_ACCURACY_REVIEW.md).

Each correction only reweights genomes already in the collection or adds genomes
that already exist in the DB (real RefSeq accessions - no invented identifiers),
then the abundances are renormalized. Corrections are deterministic and idempotent:
re-running reproduces the same family/host shares and does not re-add genomes.

Two target modes per collection:
  family_shares  set listed families to explicit target abundance shares; unlisted
                 families split the remainder in proportion to their current weight.
  host_target    scale each host class (phage / eukaryotic) to a target share,
                 distributing within the class in proportion to current weight.

Targets are literature-informed (see the reference profile + the review report).

Audit-by-default; pass --apply to write. Backs up nothing itself - the caller
(setup-db) or user should back up the DB first.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from collections import defaultdict
from pathlib import Path

# --- corrections -----------------------------------------------------------
# Blood/lung: Anelloviridae should be the dominant family (Cebria-Mendoza 2021;
# Young 2015; Abbas 2017). Nasopharynx: eukaryote-dominant, Anellovirus-led, so
# add Anelloviridae (absent) and downweight phage (Wang 2016; Megremis 2023).
# Wastewater: bacteriophage-majority, enteric eukaryotic viruses a minority
# (Gulino 2020; Kuo 2023; Calusinska 2016).
CORRECTIONS = {
    17: {  # Blood/Plasma
        "family_shares": {
            "Anelloviridae": 0.60, "Orthoherpesviridae": 0.18, "Polyomaviridae": 0.08,
            "Parvoviridae": 0.04, "Papillomaviridae": 0.03, "Adenoviridae": 0.02,
        },
    },
    19: {  # Lower respiratory (lung)
        "family_shares": {
            "Anelloviridae": 0.45, "Orthoherpesviridae": 0.08, "Papillomaviridae": 0.06,
            "Coronaviridae": 0.04, "Orthomyxoviridae": 0.03, "Pneumoviridae": 0.03,
            "Adenoviridae": 0.03, "Picornaviridae": 0.02,
        },
    },
    4: {  # Nasopharynx
        "add": [("Anelloviridae", 6)],
        "family_shares": {
            "Anelloviridae": 0.42, "Orthoherpesviridae": 0.15, "Adenoviridae": 0.10,
            "Aliceevansviridae": 0.12, "Orthomyxoviridae": 0.06,
        },
    },
    9: {  # Wastewater
        "host_target": {"phage": 0.85, "eukaryotic": 0.15},
    },
}


def load_host_map(path: Path) -> dict[str, str]:
    m = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m[row["family"]] = row["host_type"]
    return m


def host_of(family: str | None, name: str, host_map: dict[str, str]) -> str:
    if family and family in host_map:
        return host_map[family]
    if "phage" in (name or "").lower():
        return "phage"
    return "unknown"


def _members(conn, cid):
    return conn.execute(
        """SELECT cg.genome_id, COALESCE(t.family,'Unknown') AS family,
                  g.genome_name AS name, cg.relative_abundance AS ab
           FROM collection_genomes cg JOIN genomes g ON cg.genome_id=g.genome_id
           LEFT JOIN taxonomy t ON cg.genome_id=t.genome_id
           WHERE cg.collection_id=?""", (cid,)).fetchall()


def apply_family_shares(rows, shares):
    """Return {genome_id: new_abundance} hitting the listed family shares; unlisted
    families split the remainder proportionally to current weight."""
    by_fam = defaultdict(list)
    for gid, fam, _name, ab in rows:
        by_fam[fam].append((gid, ab or 0.0))
    listed = {f: s for f, s in shares.items() if f in by_fam}
    remainder = max(0.0, 1.0 - sum(listed.values()))
    other_fams = {f: v for f, v in by_fam.items() if f not in listed}
    other_total = sum(ab for lst in other_fams.values() for _, ab in lst) or 1.0

    new = {}
    for fam, lst in by_fam.items():
        fam_tot = sum(ab for _, ab in lst) or 1.0
        if fam in listed:
            target = listed[fam]
        else:
            fam_weight = sum(ab for _, ab in lst)
            target = remainder * (fam_weight / other_total)
        for gid, ab in lst:
            new[gid] = target * (ab / fam_tot)
    return new


def apply_host_target(rows, targets, host_map):
    by_host = defaultdict(list)
    for gid, fam, name, ab in rows:
        by_host[host_of(fam, name, host_map)].append((gid, ab or 0.0))
    # unknown-host keeps a small residual proportional share of the leftover
    listed_total = sum(targets.values())
    leftover = max(0.0, 1.0 - listed_total)
    new = {}
    for host, lst in by_host.items():
        tot = sum(ab for _, ab in lst) or 1.0
        if host in targets:
            target = targets[host]
        else:
            target = leftover  # unknown-host bucket
        for gid, ab in lst:
            new[gid] = target * (ab / tot)
    return new


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    ap.add_argument("--properties", type=Path, default=Path("data/reference_profiles/family_properties.tsv"))
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    host_map = load_host_map(args.properties)
    conn = sqlite3.connect(str(args.db))

    for cid, spec in CORRECTIONS.items():
        name = conn.execute("SELECT collection_name FROM body_site_collections WHERE collection_id=?",
                            (cid,)).fetchone()
        name = name[0] if name else f"collection {cid}"
        print(f"\n[{cid}] {name}")

        # additions (real DB genomes not already present)
        for fam, k in spec.get("add", []):
            existing = conn.execute(
                "SELECT COUNT(*) FROM collection_genomes cg JOIN taxonomy t ON cg.genome_id=t.genome_id "
                "WHERE cg.collection_id=? AND t.family=?", (cid, fam)).fetchone()[0]
            if existing >= k:
                print(f"  add {fam}: already has {existing} (>= {k}), skip")
                continue
            need = k - existing
            cand = conn.execute(
                """SELECT g.genome_id FROM genomes g JOIN taxonomy t ON g.genome_id=t.genome_id
                   WHERE t.family=? AND g.genome_id NOT IN
                     (SELECT genome_id FROM collection_genomes WHERE collection_id=?)
                   ORDER BY g.genome_id LIMIT ?""", (fam, cid, need)).fetchall()
            print(f"  add {fam}: inserting {len(cand)} genomes (real DB accessions)")
            if args.apply:
                conn.executemany(
                    "INSERT INTO collection_genomes (collection_id, genome_id, relative_abundance, prevalence, abundance_rank) "
                    "VALUES (?,?,?,?,?)",
                    [(cid, gid, 0.01, 0.5, 999) for (gid,) in cand])

        rows = _members(conn, cid)
        if "family_shares" in spec:
            new = apply_family_shares(rows, spec["family_shares"])
        else:
            new = apply_host_target(rows, spec["host_target"], host_map)

        # report resulting family shares
        fam_share = defaultdict(float)
        for gid, fam, _n, _ab in rows:
            fam_share[fam] += new.get(gid, 0.0)
        top = sorted(fam_share.items(), key=lambda kv: -kv[1])[:4]
        print("  target family shares: " + ", ".join(f"{f} {s:.2f}" for f, s in top))

        if args.apply:
            conn.executemany("UPDATE collection_genomes SET relative_abundance=? WHERE collection_id=? AND genome_id=?",
                             [(ab, cid, gid) for gid, ab in new.items()])

    if args.apply:
        conn.commit()
        print("\napplied. (run renormalize_abundances.py next to normalize sums to 1.0)")
    else:
        print("\naudit only. re-run with --apply to write.")
    conn.close()


if __name__ == "__main__":
    main()
