# Bibliography verification (/verify-references)

50 unique references, 100 checks (each DOI via CrossRef + each PMID via PubMed).
Outcome: all 50 identifiers resolve to the correct papers.

- 94/100 checks VERIFIED directly.
- 6 flags were false positives from CrossRef/PubMed title markup (e.g. `<scp>DNA</scp>`,
  "double-stranded" vs "double stranded") plus 3 off-by-one publication years.
- 3 genuine year corrections applied at source (online-vs-print year):
  - Pride salivary virome (PMID 22158393): 2011 -> 2012
  - Young anellovirus blooms (PMID 25403800): 2014 -> 2015
  - Abbas perioperative lung virome (PMID 27731934): 2016 -> 2017
- After correction, the three re-verified VERIFIED.

No hallucinated or wrong-paper identifiers. Manifest: validation/references_manifest.yaml.
