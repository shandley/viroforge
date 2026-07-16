#!/usr/bin/env python3
"""Assemble docs/composition_references.bib from the per-site literature dossiers.

Reads validation/dossiers/*.json, collects every citation, dedupes by PMID, and
writes a BibTeX file plus a pmid->citekey map (validation/citekeys.json) used to
wire citation keys into data/reference_profiles/virome_composition.yaml.
Citation keys are firstauthor+year, lowercased, ASCII, with a/b/... on collisions.
"""
from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path


def slug(author: str) -> str:
    a = unicodedata.normalize("NFKD", author).encode("ascii", "ignore").decode()
    a = a.split(",")[0].split()[0] if a.strip() else "anon"
    return re.sub(r"[^a-z]", "", a.lower()) or "anon"


def main() -> None:
    dossier_dir = Path("validation/dossiers")
    cites: dict[str, dict] = {}  # pmid -> citation
    for f in sorted(dossier_dir.glob("*.json")):
        data = json.loads(f.read_text())
        for site in data["sites"].values():
            for c in site.get("citations", []):
                pmid = str(c.get("pmid", "")).strip()
                if pmid and pmid not in cites:
                    cites[pmid] = c

    # assign keys
    used: dict[str, int] = {}
    keyed = []
    for pmid, c in sorted(cites.items(), key=lambda kv: (slug(kv[1]["first_author"]), kv[1]["year"])):
        base = f"{slug(c['first_author'])}{c['year']}"
        n = used.get(base, 0)
        key = base if n == 0 else f"{base}{chr(ord('a') + n)}"
        used[base] = n + 1
        keyed.append((key, pmid, c))

    bib_path = Path("docs/composition_references.bib")
    with open(bib_path, "w") as fh:
        fh.write("% ViroForge composition review bibliography\n")
        fh.write("% Auto-assembled by scripts/build_bibliography.py from the per-site dossiers.\n")
        fh.write("% Every entry to be confirmed against source via /verify-references.\n\n")
        for key, pmid, c in keyed:
            fh.write(f"@article{{{key},\n")
            fh.write(f"  title   = {{{c['title']}}},\n")
            fh.write(f"  author  = {{{c['first_author']} and others}},\n")
            fh.write(f"  journal = {{{c['journal']}}},\n")
            fh.write(f"  year    = {{{c['year']}}},\n")
            if c.get("doi"):
                fh.write(f"  doi     = {{{c['doi']}}},\n")
            fh.write(f"  pmid    = {{{pmid}}},\n")
            fh.write("}\n\n")

    keymap = {pmid: key for key, pmid, _ in keyed}
    Path("validation/citekeys.json").write_text(json.dumps(keymap, indent=2))
    print(f"{len(keyed)} unique references -> {bib_path}")
    print(f"pmid->citekey map -> validation/citekeys.json")


if __name__ == "__main__":
    main()
