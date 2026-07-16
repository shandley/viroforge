#!/usr/bin/env python3
"""Evaluate per-collection viral composition against literature expectations.

Two modes:

  observe   Compute the observed composition of every collection (phage vs
            eukaryote balance, nucleic-acid mix, dominant families) in three
            currencies (abundance-weighted, genome-count, prevalence-weighted),
            derived from the verified ICTV family property map rather than the
            unreliable genomes.genome_type column. Writes JSON; needs no profile.

  evaluate  Compare the observed composition against the tunable reference
            profile (data/reference_profiles/virome_composition.yaml) and emit a
            per-site scorecard with band status and severity. (Added once the
            profile exists.)

Usage:
  python3 scripts/evaluate_composition.py observe \
      --db viroforge/data/viral_genomes.db \
      --properties data/reference_profiles/family_properties.tsv \
      --out validation/observed_composition.json
"""
from __future__ import annotations

import argparse
import json
import sqlite3
from collections import defaultdict
from pathlib import Path


def load_property_map(path: Path) -> dict[str, dict]:
    """family -> {na_type, na_class, host_type, ...} from the ICTV-derived TSV."""
    import csv

    m: dict[str, dict] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m[row["family"]] = {
                "na_type": row["na_type"],
                "na_class": row["na_class"],
                "host_type": row["host_type"],
                "na_purity": float(row["na_purity"]),
                "host_purity": float(row["host_purity"]),
            }
    return m


def _pct(d: dict[str, float]) -> dict[str, float]:
    tot = sum(d.values())
    if tot <= 0:
        return {k: 0.0 for k in d}
    return {k: round(v / tot, 4) for k, v in d.items()}


def _looks_like_phage(name: str) -> bool:
    """Heuristic host call for genomes ICTV leaves unclassified. RefSeq phage
    names reliably contain 'phage' (or 'bacteriophage'/'prophage'); this recovers
    the large pool of bacteriophages carried as family='Unknown' so they are not
    miscounted as host-unknown."""
    n = (name or "").lower()
    return "phage" in n


def observe_collection(rows: list[dict], props: dict[str, dict]) -> dict:
    """rows: per-genome dicts with family, genome_name, relative_abundance,
    prevalence, class.

    Returns observed composition in three currencies. Host type comes from the
    ICTV family map; for genomes whose family is unmapped/Unknown, a name-based
    phage heuristic is applied before falling back to 'unknown'. Nucleic-acid
    class is only taken from the ICTV map (no name inference).
    """
    host_ab, host_n, host_prev = defaultdict(float), defaultdict(float), defaultdict(float)
    na_ab, na_n = defaultdict(float), defaultdict(float)
    fam_ab, fam_n = defaultdict(float), defaultdict(int)
    cls_ab = defaultdict(float)
    unknown_ab = 0.0
    unknown_n = 0
    inferred_phage_ab = 0.0

    for r in rows:
        fam = r["family"] or "Unknown"
        ab = r["relative_abundance"] or 0.0
        prev = r["prevalence"] if r["prevalence"] is not None else 1.0
        p = props.get(fam)
        na = p["na_class"] if p else "unknown"
        if p:
            host = p["host_type"]
        elif _looks_like_phage(r.get("genome_name", "")):
            host = "phage"
            inferred_phage_ab += ab
        else:
            host = "unknown"

        host_ab[host] += ab
        host_n[host] += 1
        host_prev[host] += prev
        na_ab[na] += ab
        na_n[na] += 1
        fam_ab[fam] += ab
        fam_n[fam] += 1
        cls_ab[r["class"] or "(none)"] += ab
        if fam == "Unknown" or p is None:
            unknown_ab += ab
            unknown_n += 1

    top_fam = sorted(fam_ab.items(), key=lambda kv: kv[1], reverse=True)
    top_cls = sorted(cls_ab.items(), key=lambda kv: kv[1], reverse=True)

    return {
        "n_genomes": len(rows),
        "unknown_family": {"abundance": round(unknown_ab, 4), "count": unknown_n},
        "inferred_phage": {"abundance": round(inferred_phage_ab, 4)},  # phage by name, family=Unknown
        "host_balance": {
            "abundance": _pct(host_ab),
            "count": _pct(host_n),
            "prevalence": _pct(host_prev),
        },
        "na_class_mix": {"abundance": _pct(na_ab), "count": _pct(na_n)},
        "top_families": [
            {"family": f, "abundance": round(a, 4), "count": fam_n[f]} for f, a in top_fam[:12]
        ],
        "top_classes": [{"class": c, "abundance": round(a, 4)} for c, a in top_cls[:8]],
    }


def _band_status(value: float, band: dict, strictness: float) -> dict:
    """Compare a value to a {low, expected, high} band scaled by strictness.

    strictness > 1 narrows the tolerated band (stricter); < 1 widens it. Returns
    status (within/below/above), the effective band, and a normalized deviation
    (0 inside the band, growing with distance outside it, in band-width units)."""
    lo, ex, hi = band["low"], band["expected"], band["high"]
    s = strictness if strictness > 0 else 1.0
    lo_eff = ex - (ex - lo) / s
    hi_eff = ex + (hi - ex) / s
    width = max(hi_eff - lo_eff, 1e-6)
    if value < lo_eff:
        return {"status": "below", "deviation": round((lo_eff - value) / width, 3),
                "effective_band": [round(lo_eff, 3), round(hi_eff, 3)]}
    if value > hi_eff:
        return {"status": "above", "deviation": round((value - hi_eff) / width, 3),
                "effective_band": [round(lo_eff, 3), round(hi_eff, 3)]}
    return {"status": "within", "deviation": 0.0,
            "effective_band": [round(lo_eff, 3), round(hi_eff, 3)]}


def _severity(deviation: float) -> str:
    if deviation == 0:
        return "ok"
    if deviation < 0.5:
        return "minor"
    if deviation < 1.5:
        return "moderate"
    return "major"


def evaluate_site(obs: dict, expect: dict, strictness: float) -> dict:
    """Score one site's observed composition against its expected bands.

    Balance and NA metrics are computed on the CLASSIFIED (known-host / known-NA)
    fraction so a high unknown-family fraction does not distort the biological
    judgment; unknown-family is instead flagged separately as a data-quality
    metric against its own ceiling."""
    findings = []
    hb = obs["host_balance"]["abundance"]
    na = obs["na_class_mix"]["abundance"]

    # phage fraction among classified (phage + eukaryotic)
    known_host = hb.get("phage", 0) + hb.get("eukaryotic", 0)
    if "phage_fraction" in expect and known_host > 0:
        v = hb.get("phage", 0) / known_host
        st = _band_status(v, expect["phage_fraction"], strictness)
        findings.append({"metric": "phage_fraction_of_classified", "observed": round(v, 3),
                         **st, "severity": _severity(st["deviation"]),
                         "cite": expect["phage_fraction"].get("cite", [])})

    # RNA fraction among classified NA
    known_na = na.get("DNA", 0) + na.get("RNA", 0) + na.get("RT", 0)
    if "rna_fraction" in expect and known_na > 0:
        v = na.get("RNA", 0) / known_na
        st = _band_status(v, expect["rna_fraction"], strictness)
        findings.append({"metric": "rna_fraction_of_classified", "observed": round(v, 3),
                         **st, "severity": _severity(st["deviation"]),
                         "cite": expect["rna_fraction"].get("cite", [])})

    # signature taxa presence / dominance
    present = {f["family"]: f for f in obs["top_families"]}
    all_fams_ab = {f["family"]: f["abundance"] for f in obs["top_families"]}
    for sig in expect.get("signature_taxa", []):
        fam = sig["family"]
        ab = all_fams_ab.get(fam, 0.0)
        role = sig.get("role", "present")
        rank = next((i + 1 for i, f in enumerate(obs["top_families"]) if f["family"] == fam), None)
        if ab <= 0:
            status, sev = "missing", ("major" if role == "dominant" else "moderate")
        elif role == "dominant" and rank != 1:
            # expected THE dominant family but it is not the most abundant one
            status, sev = "present_not_dominant", "moderate"
        else:
            status, sev = "present", "ok"
        findings.append({"metric": f"signature:{fam}", "role": role, "observed_abundance": round(ab, 3),
                         "observed_rank": rank, "status": status, "severity": sev,
                         "cite": sig.get("cite", [])})

    # unknown-family ceiling (data quality, not biology)
    if "max_unknown_family" in expect:
        v = obs["unknown_family"]["abundance"]
        ceil = expect["max_unknown_family"]["high"]
        ceil_eff = ceil / (strictness if strictness > 0 else 1.0)
        if v > ceil_eff:
            dev = round((v - ceil_eff) / max(ceil_eff, 1e-6), 3)
            findings.append({"metric": "unknown_family_fraction", "observed": round(v, 3),
                             "ceiling": round(ceil_eff, 3), "status": "above",
                             "deviation": dev, "severity": _severity(dev),
                             "kind": "data_quality", "cite": expect["max_unknown_family"].get("cite", [])})

    rank = ["ok", "minor", "moderate", "major"]
    bio = [f["severity"] for f in findings if f.get("kind") != "data_quality"]
    dq = [f["severity"] for f in findings if f.get("kind") == "data_quality"]
    bio_verdict = max(bio, key=rank.index, default="ok")
    dq_verdict = max(dq, key=rank.index, default="ok")
    # headline verdict reflects BIOLOGY; data-quality (taxonomy coverage) reported separately
    return {"site_verdict": bio_verdict, "data_quality_verdict": dq_verdict, "findings": findings}


def cmd_evaluate(args) -> None:
    import yaml

    props = load_property_map(args.properties)
    profile = yaml.safe_load(args.profile.read_text())
    strictness = float(profile.get("strictness", 1.0))
    sites = {str(k): v for k, v in (profile.get("sites") or {}).items()}

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    colls = con.execute("SELECT collection_id, collection_name FROM body_site_collections ORDER BY collection_id").fetchall()

    report = {"strictness": strictness, "sites": {}}
    for c in colls:
        cid = str(c["collection_id"])
        rows = [dict(r) for r in con.execute(
            """SELECT t.family AS family, t.class AS class, g.genome_name AS genome_name,
                      cg.relative_abundance AS relative_abundance, cg.prevalence AS prevalence
               FROM collection_genomes cg
               JOIN genomes g ON cg.genome_id = g.genome_id
               LEFT JOIN taxonomy t ON cg.genome_id = t.genome_id
               WHERE cg.collection_id = ?""", (c["collection_id"],)).fetchall()]
        obs = observe_collection(rows, props)
        obs["collection_name"] = c["collection_name"]
        site_expect = sites.get(cid, {}).get("expected")
        if site_expect:
            ev = evaluate_site(obs, site_expect, strictness)
        else:
            ev = {"site_verdict": "no_profile", "findings": []}
        report["sites"][cid] = {"name": c["collection_name"], "observed": obs, **ev}
    con.close()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2))
    print(f"evaluation (strictness={strictness}) for {len(report['sites'])} sites -> {args.out}\n")
    order = {"major": 0, "moderate": 1, "minor": 2, "ok": 3, "no_profile": 4}
    print(f"  {'id':>2} {'collection':36} {'biology':10} {'data-qual':10} flags")
    for cid, s in sorted(report["sites"].items(), key=lambda kv: order.get(kv[1]["site_verdict"], 9)):
        flags = [f["metric"] for f in s["findings"] if f["severity"] in ("major", "moderate")]
        dq = s.get("data_quality_verdict", "ok")
        print(f"  [{cid:>2}] {s['name'][:36]:36} {s['site_verdict']:10} {dq:10} {'; '.join(flags[:4])}")


def cmd_observe(args) -> None:
    props = load_property_map(args.properties)
    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    colls = con.execute(
        "SELECT collection_id, collection_name FROM body_site_collections ORDER BY collection_id"
    ).fetchall()

    out = {}
    for c in colls:
        rows = [
            dict(r)
            for r in con.execute(
                """SELECT t.family AS family, t.class AS class, g.genome_name AS genome_name,
                          cg.relative_abundance AS relative_abundance, cg.prevalence AS prevalence
                   FROM collection_genomes cg
                   JOIN genomes g ON cg.genome_id = g.genome_id
                   LEFT JOIN taxonomy t ON cg.genome_id = t.genome_id
                   WHERE cg.collection_id = ?""",
                (c["collection_id"],),
            ).fetchall()
        ]
        obs = observe_collection(rows, props)
        obs["collection_name"] = c["collection_name"]
        out[str(c["collection_id"])] = obs
    con.close()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"observed composition for {len(out)} collections -> {args.out}\n")

    # console summary
    print(f"{'id':>2} {'collection':34} {'phage%':>7} {'euk%':>6} {'unk%':>6} {'DNA%':>6} {'RNA%':>6}  top family")
    for cid, o in out.items():
        hb = o["host_balance"]["abundance"]
        na = o["na_class_mix"]["abundance"]
        tf = o["top_families"][0]["family"] if o["top_families"] else "-"
        print(
            f"{cid:>2} {o['collection_name'][:34]:34} "
            f"{hb.get('phage',0)*100:6.1f} {hb.get('eukaryotic',0)*100:5.1f} {hb.get('unknown',0)*100:5.1f} "
            f"{na.get('DNA',0)*100:5.1f} {na.get('RNA',0)*100:5.1f}  {tf}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    po = sub.add_parser("observe", help="dump observed composition (no profile needed)")
    po.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    po.add_argument("--properties", type=Path, default=Path("data/reference_profiles/family_properties.tsv"))
    po.add_argument("--out", type=Path, default=Path("validation/observed_composition.json"))
    po.set_defaults(func=cmd_observe)

    pe = sub.add_parser("evaluate", help="score observed composition vs the reference profile")
    pe.add_argument("--db", type=Path, default=Path("viroforge/data/viral_genomes.db"))
    pe.add_argument("--properties", type=Path, default=Path("data/reference_profiles/family_properties.tsv"))
    pe.add_argument("--profile", type=Path, default=Path("data/reference_profiles/virome_composition.yaml"))
    pe.add_argument("--out", type=Path, default=Path("validation/composition_scorecard.json"))
    pe.set_defaults(func=cmd_evaluate)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
