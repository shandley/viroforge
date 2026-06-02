# Database Fix Scripts — Execution Order

## Overview

The ViroForge database (`viral_genomes.db`) is not tracked in git due to its size (~500 MB). Several fix scripts must be run **in order** after building or updating the database to correct genome curation issues.

These scripts are idempotent — running them multiple times is safe and produces the same result.

## Required Execution Order

**Run these in sequence. Each step depends on the previous one.**

### Step 1: Populate host_associations table

**PR #51** — `feat/populate-host-associations`

```bash
python scripts/populate_host_associations.py
```

- **What it does**: Parses host bacterium genus from 6,070 phage genome names and populates the `host_associations` table (previously empty — 0 rows)
- **Why it must be first**: Steps 2 and 3 rely on host association data to find correct replacement genomes. The `body_site_filter` in `curate_body_site_collections.py` also depends on this table.
- **Result**: 5,957 phage-host mappings with body site annotations (gut, oral, skin, respiratory, vaginal, marine, soil, etc.)
- **Verification**: `sqlite3 viroforge/data/viral_genomes.db "SELECT COUNT(*) FROM host_associations"` → should show 5957

### Step 2: Remove animal/plant viruses from human collections

**PR #41** — `fix/animal-virus-contamination`

```bash
python scripts/fix_animal_virus_contamination.py --apply
```

- **What it does**: Removes 108 animal viruses (bat, bovine, porcine, simian adenoviruses; penguin/koala herpesviruses; plant viruses) and PhiX174 from 14 human-associated collections. Replaces 82 with correct human viruses (human adenovirus, HSV, EBV, CMV, human TTV, influenza, SARS-CoV-2, etc.)
- **Why it must be second**: Uses genome name patterns to find human replacements. Does not depend on host_associations, but should run before Step 3 to avoid replacing animal viruses with wrong phages.
- **Result**: 0 animal/plant viruses in any human-associated collection
- **Verification**: `python scripts/fix_animal_virus_contamination.py --dry-run` → should show "Total animal/plant viruses found: 0"

### Step 3: Fix phage host specificity

**PR #50** — `fix/phage-host-specificity`

```bash
python scripts/fix_phage_host_specificity.py --apply
```

- **What it does**: Removes 84 non-site-appropriate phages (marine Vibrio phages in gut, plant Ralstonia phages in oral, E. coli phages in skin/vaginal, etc.) and replaces them with body-site-matched phages (Streptococcus phages for oral, Propionibacterium phages for skin, Bacteroides/Faecalibacterium phages for gut, etc.)
- **Why it must be third**: Depends on Step 2 having already removed animal viruses (to avoid double-counting). Uses host_associations data from Step 1 for some replacement logic.
- **Result**: 0 non-site-appropriate phages in any collection
- **Verification**: `python scripts/fix_phage_host_specificity.py --dry-run` → should show "Total: 0 wrong phages removed"

## Quick Copy-Paste

```bash
# Run all three fixes in order
python scripts/populate_host_associations.py
python scripts/fix_animal_virus_contamination.py --apply
python scripts/fix_phage_host_specificity.py --apply

# Verify all clean
python scripts/fix_animal_virus_contamination.py --dry-run
python scripts/fix_phage_host_specificity.py --dry-run
```

## What These Fix

| Issue | Affected Collections | Examples of Wrong Genomes |
|---|---|---|
| Animal viruses in human collections | 14 collections | Bat adenovirus, penguin herpesvirus, chicken anemia virus |
| PhiX174 lab spike-in as natural virus | 9 collections | GCF_000819615.1 in gut, oral, skin, etc. |
| Marine phages in human body sites | Gut, Oral, IBD, HIV+ | Vibrio phages (marine bacteria) |
| Plant phages in human body sites | Gut, Oral, Respiratory | Ralstonia, Xanthomonas phages (plant pathogens) |
| E. coli phages in non-gut sites | Oral, Skin, Respiratory, Vaginal | Coliphages, Enterobacteria phages |
| Insect phages in human body sites | Multiple | Spiroplasma phage (insect endosymbiont) |

## Root Cause Prevention

After running these fixes, the `curate_body_site_collections.py` script has been updated with `body_site_filter` parameters (PR #51) that use the populated `host_associations` table. Future collection rebuilds will automatically select site-appropriate phages.
