#!/usr/bin/env python3
"""
Build bacterial reference fragments for ViroForge.

Downloads representative genomes for dominant bacteria at each body site,
extracts random fragments, and saves them as FASTA files with abundance
metadata. Supports multiple body sites.

Usage:
    python scripts/build_bacterial_references.py                # Build all body sites
    python scripts/build_bacterial_references.py --site gut_human
    python scripts/build_bacterial_references.py --site vaginal_human
    python scripts/build_bacterial_references.py --list-sites

Literature references per body site:
- Gut: Qin et al. 2010 (MetaHIT); Lloyd-Price et al. 2019 (HMP2); Almeida et al. 2021 (UHGG)
- Vaginal: Ravel et al. 2011 (CST); Fettweis et al. 2019 (VIRGO/HMP)
- Oral: Dewhirst et al. 2010 (HOMD); Escapa et al. 2018 (HMP)
- Skin: Oh et al. 2014 (HMP skin); Byrd et al. 2018 (NISC)
- Respiratory: Dickson et al. 2016 (lung); Man et al. 2017 (nasopharynx)
- Blood: Potgieter et al. 2015; Castillo et al. 2019 (blood microbiome)
- Ocular: Ozkan et al. 2019; Doan et al. 2016 (ocular surface)
- Urinary: Hilt et al. 2014; Thomas-White et al. 2018 (female urinary)
- Lung: Dickson et al. 2016; Segal et al. 2016 (lower respiratory)
"""

import argparse
import os
import sys
import json
import random
import subprocess
from pathlib import Path

# ============================================================
# Body site bacterial species lists
# ============================================================
# Each entry: species name, RefSeq accession, relative abundance, prevalence
# Abundances are approximate and based on cited literature.

GUT_BACTERIA = [
    # Bacteroidetes (dominant in Western diet)
    {"species": "Bacteroides vulgatus", "accession": "GCF_000012825.1", "abundance": 0.065, "prevalence": 0.90},
    {"species": "Bacteroides uniformis", "accession": "GCF_000154205.1", "abundance": 0.055, "prevalence": 0.85},
    {"species": "Bacteroides thetaiotaomicron", "accession": "GCF_000011065.1", "abundance": 0.045, "prevalence": 0.80},
    {"species": "Bacteroides fragilis", "accession": "GCF_000009925.1", "abundance": 0.035, "prevalence": 0.75},
    {"species": "Bacteroides ovatus", "accession": "GCF_000154125.1", "abundance": 0.030, "prevalence": 0.70},
    {"species": "Bacteroides caccae", "accession": "GCF_000154845.1", "abundance": 0.020, "prevalence": 0.60},
    {"species": "Parabacteroides distasonis", "accession": "GCF_000012845.1", "abundance": 0.025, "prevalence": 0.70},
    {"species": "Parabacteroides merdae", "accession": "GCF_000154465.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Prevotella copri", "accession": "GCF_000157935.1", "abundance": 0.040, "prevalence": 0.35},
    {"species": "Alistipes putredinis", "accession": "GCF_000154525.1", "abundance": 0.020, "prevalence": 0.65},

    # Firmicutes - Lachnospiraceae
    {"species": "Roseburia intestinalis", "accession": "GCF_000169035.1", "abundance": 0.035, "prevalence": 0.75},
    {"species": "Eubacterium rectale", "accession": "GCF_000020605.1", "abundance": 0.040, "prevalence": 0.80},
    {"species": "Coprococcus comes", "accession": "GCF_000155815.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Dorea formicigenerans", "accession": "GCF_000156015.1", "abundance": 0.012, "prevalence": 0.60},
    {"species": "Dorea longicatena", "accession": "GCF_000155875.1", "abundance": 0.012, "prevalence": 0.55},
    {"species": "Blautia obeum", "accession": "GCF_000154065.1", "abundance": 0.018, "prevalence": 0.65},

    # Firmicutes - Ruminococcaceae
    {"species": "Faecalibacterium prausnitzii", "accession": "GCF_000162015.1", "abundance": 0.055, "prevalence": 0.90},
    {"species": "Ruminococcus bromii", "accession": "GCF_000156375.1", "abundance": 0.025, "prevalence": 0.70},
    {"species": "Subdoligranulum variabile", "accession": "GCF_000177015.2", "abundance": 0.015, "prevalence": 0.55},

    # Firmicutes - Clostridiaceae / other
    {"species": "Clostridium bolteae", "accession": "GCF_000155435.1", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Clostridium symbiosum", "accession": "GCF_000155205.1", "abundance": 0.008, "prevalence": 0.40},
    {"species": "Dialister invisus", "accession": "GCF_000157955.2", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Phascolarctobacterium succinatutens", "accession": "GCF_000424985.1", "abundance": 0.012, "prevalence": 0.50},

    # Firmicutes - Erysipelotrichaceae
    {"species": "Holdemanella biformis", "accession": "GCF_000156495.1", "abundance": 0.008, "prevalence": 0.35},

    # Firmicutes - Lactobacillaceae / Streptococcaceae
    {"species": "Streptococcus thermophilus", "accession": "GCF_000011825.1", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Lactobacillus ruminis", "accession": "GCF_000224985.2", "abundance": 0.003, "prevalence": 0.25},

    # Actinobacteria
    {"species": "Bifidobacterium longum", "accession": "GCF_000007525.1", "abundance": 0.030, "prevalence": 0.75},
    {"species": "Bifidobacterium adolescentis", "accession": "GCF_000010425.1", "abundance": 0.020, "prevalence": 0.60},
    {"species": "Bifidobacterium bifidum", "accession": "GCF_000167455.2", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Collinsella aerofaciens", "accession": "GCF_000169035.1", "abundance": 0.015, "prevalence": 0.55},
    {"species": "Eggerthella lenta", "accession": "GCF_000024265.1", "abundance": 0.008, "prevalence": 0.40},

    # Verrucomicrobia
    {"species": "Akkermansia muciniphila", "accession": "GCF_000020225.1", "abundance": 0.025, "prevalence": 0.65},

    # Proteobacteria (low abundance in healthy gut)
    {"species": "Escherichia coli", "accession": "GCF_000005845.2", "abundance": 0.010, "prevalence": 0.70},
    {"species": "Klebsiella pneumoniae", "accession": "GCF_000240185.1", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Sutterella wadsworthensis", "accession": "GCF_000154245.1", "abundance": 0.008, "prevalence": 0.40},
    {"species": "Bilophila wadsworthia", "accession": "GCF_000185705.2", "abundance": 0.006, "prevalence": 0.35},
    {"species": "Desulfovibrio piger", "accession": "GCF_000156375.1", "abundance": 0.004, "prevalence": 0.25},
    {"species": "Haemophilus parainfluenzae", "accession": "GCF_000191405.1", "abundance": 0.003, "prevalence": 0.20},

    # Fusobacteria
    {"species": "Fusobacterium nucleatum", "accession": "GCF_000007325.1", "abundance": 0.003, "prevalence": 0.20},

    # Archaea (methanogenic)
    {"species": "Methanobrevibacter smithii", "accession": "GCF_000016525.1", "abundance": 0.015, "prevalence": 0.50},

    # Additional common species
    {"species": "Odoribacter splanchnicus", "accession": "GCF_000190535.1", "abundance": 0.010, "prevalence": 0.45},
    {"species": "Barnesiella intestinihominis", "accession": "GCF_000311925.1", "abundance": 0.008, "prevalence": 0.35},
    {"species": "Megamonas hypermegale", "accession": "GCF_000428645.1", "abundance": 0.006, "prevalence": 0.25},
    {"species": "Mitsuokella multacida", "accession": "GCF_000178935.2", "abundance": 0.005, "prevalence": 0.30},
    {"species": "Veillonella parvula", "accession": "GCF_000024945.1", "abundance": 0.005, "prevalence": 0.35},
    {"species": "Megasphaera elsdenii", "accession": "GCF_000284795.1", "abundance": 0.004, "prevalence": 0.20},
    {"species": "Citrobacter freundii", "accession": "GCF_000648515.1", "abundance": 0.003, "prevalence": 0.15},
    {"species": "Enterococcus faecalis", "accession": "GCF_000007785.1", "abundance": 0.003, "prevalence": 0.25},
    {"species": "Ruminococcus gnavus", "accession": "GCF_000156515.1", "abundance": 0.018, "prevalence": 0.65},
    {"species": "Coprococcus eutactus", "accession": "GCF_000155435.1", "abundance": 0.012, "prevalence": 0.50},
    {"species": "Anaerostipes hadrus", "accession": "GCF_000332875.2", "abundance": 0.014, "prevalence": 0.55},
]

# Vaginal microbiome — Lactobacillus-dominated (CST I-III-like composite)
# Ravel et al. 2011 (PNAS 108:4680); Fettweis et al. 2019 (Nat Med 25:1012)
VAGINAL_BACTERIA = [
    # Lactobacillus spp. (dominant, 70-95% in healthy)
    {"species": "Lactobacillus crispatus", "accession": "GCF_000091685.1", "abundance": 0.35, "prevalence": 0.70},
    {"species": "Lactobacillus iners", "accession": "GCF_000160075.2", "abundance": 0.25, "prevalence": 0.85},
    {"species": "Lactobacillus jensenii", "accession": "GCF_000159215.2", "abundance": 0.10, "prevalence": 0.50},
    {"species": "Lactobacillus gasseri", "accession": "GCF_000014425.1", "abundance": 0.08, "prevalence": 0.40},

    # BV-associated (low abundance in healthy, elevated in dysbiosis)
    {"species": "Gardnerella vaginalis", "accession": "GCF_000025205.1", "abundance": 0.05, "prevalence": 0.60},
    {"species": "Atopobium vaginae", "accession": "GCF_000160015.2", "abundance": 0.03, "prevalence": 0.30},
    {"species": "Prevotella bivia", "accession": "GCF_000195635.1", "abundance": 0.025, "prevalence": 0.35},
    {"species": "Sneathia amnii", "accession": "GCF_000759755.1", "abundance": 0.015, "prevalence": 0.20},
    {"species": "Megasphaera genomosp type 1", "accession": "GCF_000284795.1", "abundance": 0.01, "prevalence": 0.15},
    {"species": "Mobiluncus curtisii", "accession": "GCF_000154405.1", "abundance": 0.01, "prevalence": 0.15},

    # Other common vaginal species
    {"species": "Streptococcus agalactiae", "accession": "GCF_000196055.1", "abundance": 0.02, "prevalence": 0.25},
    {"species": "Enterococcus faecalis", "accession": "GCF_000007785.1", "abundance": 0.01, "prevalence": 0.20},
    {"species": "Escherichia coli", "accession": "GCF_000005845.2", "abundance": 0.01, "prevalence": 0.25},
    {"species": "Ureaplasma parvum", "accession": "GCF_000006625.1", "abundance": 0.015, "prevalence": 0.30},
    {"species": "Mycoplasma hominis", "accession": "GCF_000085865.1", "abundance": 0.01, "prevalence": 0.15},
]

# Oral microbiome — Streptococcus-dominated saliva
# Dewhirst et al. 2010 (J Bacteriol 192:5002); Escapa et al. 2018 (Microbiome 6:173)
ORAL_BACTERIA = [
    # Streptococcus (dominant in saliva)
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.15, "prevalence": 0.95},
    {"species": "Streptococcus oralis", "accession": "GCF_000164555.2", "abundance": 0.10, "prevalence": 0.90},
    {"species": "Streptococcus salivarius", "accession": "GCF_000014205.1", "abundance": 0.08, "prevalence": 0.85},
    {"species": "Streptococcus sanguinis", "accession": "GCF_000014245.1", "abundance": 0.06, "prevalence": 0.80},
    {"species": "Streptococcus gordonii", "accession": "GCF_000017005.1", "abundance": 0.04, "prevalence": 0.70},
    {"species": "Streptococcus parasanguinis", "accession": "GCF_000227135.1", "abundance": 0.03, "prevalence": 0.65},

    # Haemophilus
    {"species": "Haemophilus parainfluenzae", "accession": "GCF_000191405.1", "abundance": 0.06, "prevalence": 0.80},

    # Neisseria
    {"species": "Neisseria subflava", "accession": "GCF_001027105.1", "abundance": 0.05, "prevalence": 0.75},
    {"species": "Neisseria mucosa", "accession": "GCF_000175855.2", "abundance": 0.03, "prevalence": 0.55},

    # Veillonella
    {"species": "Veillonella parvula", "accession": "GCF_000024945.1", "abundance": 0.05, "prevalence": 0.85},
    {"species": "Veillonella dispar", "accession": "GCF_000160055.1", "abundance": 0.03, "prevalence": 0.65},

    # Prevotella (common in saliva)
    {"species": "Prevotella melaninogenica", "accession": "GCF_000024805.1", "abundance": 0.04, "prevalence": 0.70},
    {"species": "Prevotella histicola", "accession": "GCF_000382765.1", "abundance": 0.02, "prevalence": 0.50},

    # Rothia
    {"species": "Rothia mucilaginosa", "accession": "GCF_000175615.1", "abundance": 0.04, "prevalence": 0.75},
    {"species": "Rothia dentocariosa", "accession": "GCF_000163955.2", "abundance": 0.02, "prevalence": 0.50},

    # Actinomyces
    {"species": "Actinomyces naeslundii", "accession": "GCF_000181595.1", "abundance": 0.03, "prevalence": 0.65},

    # Gemella
    {"species": "Gemella haemolysans", "accession": "GCF_000186725.1", "abundance": 0.02, "prevalence": 0.55},

    # Fusobacterium
    {"species": "Fusobacterium nucleatum", "accession": "GCF_000007325.1", "abundance": 0.03, "prevalence": 0.60},

    # Porphyromonas
    {"species": "Porphyromonas gingivalis", "accession": "GCF_000010505.1", "abundance": 0.015, "prevalence": 0.30},

    # Granulicatella
    {"species": "Granulicatella adiacens", "accession": "GCF_000186165.1", "abundance": 0.02, "prevalence": 0.50},
]

# Skin microbiome — Cutibacterium (P. acnes) dominated sebaceous sites
# Oh et al. 2014 (Genome Med 6:R152); Byrd et al. 2018 (Nat Genet 50:1341)
SKIN_BACTERIA = [
    # Cutibacterium (dominant at sebaceous sites)
    {"species": "Cutibacterium acnes", "accession": "GCF_000144405.1", "abundance": 0.35, "prevalence": 0.95},

    # Staphylococcus
    {"species": "Staphylococcus epidermidis", "accession": "GCF_000007645.1", "abundance": 0.15, "prevalence": 0.90},
    {"species": "Staphylococcus hominis", "accession": "GCF_000234795.1", "abundance": 0.05, "prevalence": 0.60},
    {"species": "Staphylococcus capitis", "accession": "GCF_000504445.1", "abundance": 0.03, "prevalence": 0.45},
    {"species": "Staphylococcus aureus", "accession": "GCF_000013425.1", "abundance": 0.02, "prevalence": 0.30},

    # Corynebacterium
    {"species": "Corynebacterium tuberculostearicum", "accession": "GCF_000215345.1", "abundance": 0.06, "prevalence": 0.55},
    {"species": "Corynebacterium jeikeium", "accession": "GCF_000006325.1", "abundance": 0.03, "prevalence": 0.35},
    {"species": "Corynebacterium kroppenstedtii", "accession": "GCF_000024185.1", "abundance": 0.02, "prevalence": 0.30},

    # Micrococcus
    {"species": "Micrococcus luteus", "accession": "GCF_000023205.1", "abundance": 0.04, "prevalence": 0.50},

    # Streptococcus (moist sites)
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.04, "prevalence": 0.45},

    # Brevibacterium
    {"species": "Brevibacterium linens", "accession": "GCF_000270005.1", "abundance": 0.02, "prevalence": 0.25},

    # Dermabacter
    {"species": "Dermabacter hominis", "accession": "GCF_000374645.1", "abundance": 0.02, "prevalence": 0.30},

    # Acinetobacter
    {"species": "Acinetobacter johnsonii", "accession": "GCF_000368385.1", "abundance": 0.03, "prevalence": 0.35},

    # Enhydrobacter
    {"species": "Enhydrobacter aerosaccus", "accession": "GCF_000368565.1", "abundance": 0.02, "prevalence": 0.25},

    # Pseudomonas (transient)
    {"species": "Pseudomonas aeruginosa", "accession": "GCF_000006765.1", "abundance": 0.01, "prevalence": 0.15},
]

# Respiratory (nasopharynx) — mixed flora
# Man et al. 2017 (Am J Respir Crit Care Med 196:1582); Dickson et al. 2016 (Ann Rev Physiol 78:481)
RESPIRATORY_BACTERIA = [
    # Streptococcus (dominant in nasopharynx)
    {"species": "Streptococcus pneumoniae", "accession": "GCF_000006885.1", "abundance": 0.10, "prevalence": 0.40},
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.08, "prevalence": 0.80},
    {"species": "Streptococcus oralis", "accession": "GCF_000164555.2", "abundance": 0.05, "prevalence": 0.65},
    {"species": "Streptococcus salivarius", "accession": "GCF_000014205.1", "abundance": 0.04, "prevalence": 0.55},

    # Moraxella (common in children)
    {"species": "Moraxella catarrhalis", "accession": "GCF_000092265.1", "abundance": 0.08, "prevalence": 0.45},

    # Haemophilus
    {"species": "Haemophilus influenzae", "accession": "GCF_000027305.1", "abundance": 0.08, "prevalence": 0.50},
    {"species": "Haemophilus parainfluenzae", "accession": "GCF_000191405.1", "abundance": 0.05, "prevalence": 0.60},

    # Staphylococcus
    {"species": "Staphylococcus aureus", "accession": "GCF_000013425.1", "abundance": 0.06, "prevalence": 0.35},
    {"species": "Staphylococcus epidermidis", "accession": "GCF_000007645.1", "abundance": 0.04, "prevalence": 0.50},

    # Corynebacterium
    {"species": "Corynebacterium pseudodiphtheriticum", "accession": "GCF_000235185.1", "abundance": 0.06, "prevalence": 0.55},
    {"species": "Corynebacterium propinquum", "accession": "GCF_000231325.2", "abundance": 0.03, "prevalence": 0.35},

    # Dolosigranulum
    {"species": "Dolosigranulum pigrum", "accession": "GCF_000382845.1", "abundance": 0.06, "prevalence": 0.50},

    # Neisseria
    {"species": "Neisseria subflava", "accession": "GCF_001027105.1", "abundance": 0.04, "prevalence": 0.45},

    # Prevotella
    {"species": "Prevotella melaninogenica", "accession": "GCF_000024805.1", "abundance": 0.04, "prevalence": 0.40},

    # Veillonella
    {"species": "Veillonella parvula", "accession": "GCF_000024945.1", "abundance": 0.04, "prevalence": 0.50},

    # Fusobacterium
    {"species": "Fusobacterium nucleatum", "accession": "GCF_000007325.1", "abundance": 0.02, "prevalence": 0.25},

    # Propionibacterium / Cutibacterium
    {"species": "Cutibacterium acnes", "accession": "GCF_000144405.1", "abundance": 0.05, "prevalence": 0.60},

    # Alloiococcus
    {"species": "Alloiococcus otitis", "accession": "GCF_000382865.1", "abundance": 0.03, "prevalence": 0.25},

    # Pseudomonas
    {"species": "Pseudomonas aeruginosa", "accession": "GCF_000006765.1", "abundance": 0.02, "prevalence": 0.15},
]

# Blood/Plasma — very low biomass, translocated species
# Potgieter et al. 2015 (Gut Microbes 6:321); Castillo et al. 2019 (Int J Mol Sci 20:5042)
BLOOD_BACTERIA = [
    # Proteobacteria (most common in blood)
    {"species": "Escherichia coli", "accession": "GCF_000005845.2", "abundance": 0.15, "prevalence": 0.30},
    {"species": "Pseudomonas aeruginosa", "accession": "GCF_000006765.1", "abundance": 0.08, "prevalence": 0.15},
    {"species": "Acinetobacter baumannii", "accession": "GCF_000746645.1", "abundance": 0.05, "prevalence": 0.10},

    # Firmicutes
    {"species": "Staphylococcus epidermidis", "accession": "GCF_000007645.1", "abundance": 0.12, "prevalence": 0.25},
    {"species": "Staphylococcus aureus", "accession": "GCF_000013425.1", "abundance": 0.08, "prevalence": 0.15},
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.08, "prevalence": 0.20},
    {"species": "Enterococcus faecalis", "accession": "GCF_000007785.1", "abundance": 0.05, "prevalence": 0.10},

    # Actinobacteria
    {"species": "Cutibacterium acnes", "accession": "GCF_000144405.1", "abundance": 0.12, "prevalence": 0.25},
    {"species": "Corynebacterium jeikeium", "accession": "GCF_000006325.1", "abundance": 0.05, "prevalence": 0.10},

    # Bacteroidetes (gut translocation)
    {"species": "Bacteroides fragilis", "accession": "GCF_000009925.1", "abundance": 0.06, "prevalence": 0.10},

    # Other
    {"species": "Klebsiella pneumoniae", "accession": "GCF_000240185.1", "abundance": 0.05, "prevalence": 0.10},
    {"species": "Serratia marcescens", "accession": "GCF_000513215.1", "abundance": 0.03, "prevalence": 0.05},
]

# Ocular surface — low biomass conjunctival flora
# Ozkan et al. 2019 (Exp Eye Res 187:107762); Doan et al. 2016 (Invest Ophthalmol Vis Sci 57:5116)
OCULAR_BACTERIA = [
    # Coagulase-negative Staphylococci (dominant)
    {"species": "Staphylococcus epidermidis", "accession": "GCF_000007645.1", "abundance": 0.20, "prevalence": 0.80},
    {"species": "Staphylococcus aureus", "accession": "GCF_000013425.1", "abundance": 0.05, "prevalence": 0.25},

    # Corynebacterium
    {"species": "Corynebacterium macginleyi", "accession": "GCF_000404265.1", "abundance": 0.15, "prevalence": 0.60},
    {"species": "Corynebacterium propinquum", "accession": "GCF_000231325.2", "abundance": 0.05, "prevalence": 0.30},

    # Cutibacterium
    {"species": "Cutibacterium acnes", "accession": "GCF_000144405.1", "abundance": 0.15, "prevalence": 0.65},

    # Streptococcus
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.08, "prevalence": 0.40},

    # Haemophilus
    {"species": "Haemophilus influenzae", "accession": "GCF_000027305.1", "abundance": 0.05, "prevalence": 0.20},

    # Pseudomonas
    {"species": "Pseudomonas aeruginosa", "accession": "GCF_000006765.1", "abundance": 0.04, "prevalence": 0.10},

    # Micrococcus
    {"species": "Micrococcus luteus", "accession": "GCF_000023205.1", "abundance": 0.06, "prevalence": 0.35},

    # Enhydrobacter
    {"species": "Enhydrobacter aerosaccus", "accession": "GCF_000368565.1", "abundance": 0.05, "prevalence": 0.25},

    # Acinetobacter
    {"species": "Acinetobacter johnsonii", "accession": "GCF_000368385.1", "abundance": 0.04, "prevalence": 0.20},

    # Sphingomonas
    {"species": "Sphingomonas paucimobilis", "accession": "GCF_000282435.1", "abundance": 0.03, "prevalence": 0.15},
]

# Urinary — female urinary microbiome
# Hilt et al. 2014 (J Clin Microbiol 52:871); Thomas-White et al. 2018 (J Urol 199:1497)
URINARY_BACTERIA = [
    # Lactobacillus (dominant in healthy female bladder)
    {"species": "Lactobacillus crispatus", "accession": "GCF_000091685.1", "abundance": 0.20, "prevalence": 0.60},
    {"species": "Lactobacillus iners", "accession": "GCF_000160075.2", "abundance": 0.12, "prevalence": 0.50},
    {"species": "Lactobacillus jensenii", "accession": "GCF_000159215.2", "abundance": 0.06, "prevalence": 0.30},
    {"species": "Lactobacillus gasseri", "accession": "GCF_000014425.1", "abundance": 0.05, "prevalence": 0.25},

    # Gardnerella
    {"species": "Gardnerella vaginalis", "accession": "GCF_000025205.1", "abundance": 0.08, "prevalence": 0.40},

    # Streptococcus
    {"species": "Streptococcus agalactiae", "accession": "GCF_000196055.1", "abundance": 0.06, "prevalence": 0.25},
    {"species": "Streptococcus anginosus", "accession": "GCF_000166295.2", "abundance": 0.04, "prevalence": 0.20},

    # Staphylococcus
    {"species": "Staphylococcus epidermidis", "accession": "GCF_000007645.1", "abundance": 0.05, "prevalence": 0.30},

    # Corynebacterium
    {"species": "Corynebacterium coyleae", "accession": "GCF_000382825.1", "abundance": 0.04, "prevalence": 0.20},

    # Enterococcus
    {"species": "Enterococcus faecalis", "accession": "GCF_000007785.1", "abundance": 0.05, "prevalence": 0.25},

    # Escherichia (UTI-associated)
    {"species": "Escherichia coli", "accession": "GCF_000005845.2", "abundance": 0.06, "prevalence": 0.30},

    # Prevotella
    {"species": "Prevotella bivia", "accession": "GCF_000195635.1", "abundance": 0.04, "prevalence": 0.20},

    # Atopobium
    {"species": "Atopobium vaginae", "accession": "GCF_000160015.2", "abundance": 0.03, "prevalence": 0.15},

    # Aerococcus
    {"species": "Aerococcus urinae", "accession": "GCF_000373925.1", "abundance": 0.04, "prevalence": 0.20},

    # Actinobaculum / Actinotignum
    {"species": "Actinotignum schaalii", "accession": "GCF_000382085.1", "abundance": 0.03, "prevalence": 0.15},
]

# Lung (lower respiratory) — low biomass, aspiration-derived
# Dickson et al. 2016 (Ann Rev Physiol 78:481); Segal et al. 2016 (mBio 7:e01285)
LUNG_BACTERIA = [
    # Oral-derived (dominant via microaspiration)
    {"species": "Streptococcus mitis", "accession": "GCF_000148585.2", "abundance": 0.12, "prevalence": 0.70},
    {"species": "Streptococcus pneumoniae", "accession": "GCF_000006885.1", "abundance": 0.08, "prevalence": 0.35},
    {"species": "Veillonella parvula", "accession": "GCF_000024945.1", "abundance": 0.08, "prevalence": 0.60},
    {"species": "Prevotella melaninogenica", "accession": "GCF_000024805.1", "abundance": 0.08, "prevalence": 0.55},
    {"species": "Haemophilus parainfluenzae", "accession": "GCF_000191405.1", "abundance": 0.06, "prevalence": 0.45},
    {"species": "Fusobacterium nucleatum", "accession": "GCF_000007325.1", "abundance": 0.04, "prevalence": 0.30},
    {"species": "Rothia mucilaginosa", "accession": "GCF_000175615.1", "abundance": 0.06, "prevalence": 0.50},
    {"species": "Neisseria subflava", "accession": "GCF_001027105.1", "abundance": 0.04, "prevalence": 0.35},
    {"species": "Gemella haemolysans", "accession": "GCF_000186725.1", "abundance": 0.03, "prevalence": 0.30},
    {"species": "Actinomyces naeslundii", "accession": "GCF_000181595.1", "abundance": 0.03, "prevalence": 0.25},

    # Respiratory-specific
    {"species": "Moraxella catarrhalis", "accession": "GCF_000092265.1", "abundance": 0.04, "prevalence": 0.20},
    {"species": "Haemophilus influenzae", "accession": "GCF_000027305.1", "abundance": 0.05, "prevalence": 0.30},
    {"species": "Pseudomonas aeruginosa", "accession": "GCF_000006765.1", "abundance": 0.03, "prevalence": 0.10},
    {"species": "Staphylococcus aureus", "accession": "GCF_000013425.1", "abundance": 0.03, "prevalence": 0.15},

    # Tropheryma (rare but lung-specific)
    {"species": "Tropheryma whipplei", "accession": "GCF_000016285.1", "abundance": 0.04, "prevalence": 0.20},

    # Cutibacterium
    {"species": "Cutibacterium acnes", "accession": "GCF_000144405.1", "abundance": 0.05, "prevalence": 0.40},

    # Corynebacterium
    {"species": "Corynebacterium pseudodiphtheriticum", "accession": "GCF_000235185.1", "abundance": 0.03, "prevalence": 0.25},
]

# ============================================================
# Body site registry
# ============================================================
BODY_SITE_BACTERIA = {
    "gut_human": {
        "species_list": GUT_BACTERIA,
        "description": "Top 50 human gut bacteria (healthy Western diet)",
        "source": "MetaHIT/HMP2/UHGG",
    },
    "vaginal_human": {
        "species_list": VAGINAL_BACTERIA,
        "description": "Vaginal microbiome (Lactobacillus-dominated, CST composite)",
        "source": "Ravel2011/Fettweis2019",
    },
    "oral_human": {
        "species_list": ORAL_BACTERIA,
        "description": "Oral microbiome (saliva, Streptococcus-dominated)",
        "source": "HOMD/HMP",
    },
    "skin_human": {
        "species_list": SKIN_BACTERIA,
        "description": "Skin microbiome (sebaceous sites, Cutibacterium-dominated)",
        "source": "Oh2014/Byrd2018",
    },
    "respiratory_human": {
        "species_list": RESPIRATORY_BACTERIA,
        "description": "Nasopharyngeal microbiome (mixed flora)",
        "source": "Man2017/Dickson2016",
    },
    "blood_human": {
        "species_list": BLOOD_BACTERIA,
        "description": "Blood microbiome (low biomass, translocated species)",
        "source": "Potgieter2015/Castillo2019",
    },
    "ocular_human": {
        "species_list": OCULAR_BACTERIA,
        "description": "Ocular surface microbiome (low biomass conjunctival flora)",
        "source": "Ozkan2019/Doan2016",
    },
    "urinary_human": {
        "species_list": URINARY_BACTERIA,
        "description": "Female urinary microbiome",
        "source": "Hilt2014/Thomas-White2018",
    },
    "lung_human": {
        "species_list": LUNG_BACTERIA,
        "description": "Lower respiratory microbiome (aspiration-derived)",
        "source": "Dickson2016/Segal2016",
    },
}

FRAGMENT_LENGTH = 50000
FRAGMENTS_PER_SPECIES = 35
RANDOM_SEED = 42


def download_genome(accession: str, output_dir: Path) -> Path:
    """Download genome from NCBI using datasets CLI or efetch."""
    fasta_path = output_dir / f"{accession}.fasta"
    if fasta_path.exists():
        return fasta_path

    # Try NCBI datasets CLI first
    try:
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(output_dir / f"{accession}.zip")
        ]
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)

        # Extract
        import zipfile
        zip_path = output_dir / f"{accession}.zip"
        with zipfile.ZipFile(zip_path) as zf:
            for name in zf.namelist():
                if name.endswith(".fna"):
                    with zf.open(name) as src, open(fasta_path, 'wb') as dst:
                        dst.write(src.read())
                    break
        zip_path.unlink()
        return fasta_path
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Fallback: efetch via Entrez
    try:
        from Bio import Entrez, SeqIO
        Entrez.email = "viroforge@example.com"

        # Get assembly info to find nucleotide accession
        handle = Entrez.esearch(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            assembly_id = record["IdList"][0]
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()

            # Get FTP path
            ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0].get('FtpPath_RefSeq', '')
            if ftp_path:
                import urllib.request
                fname = ftp_path.split('/')[-1]
                url = f"{ftp_path}/{fname}_genomic.fna.gz"
                gz_path = output_dir / f"{accession}.fna.gz"
                urllib.request.urlretrieve(url, gz_path)

                import gzip
                with gzip.open(gz_path, 'rt') as gz, open(fasta_path, 'w') as out:
                    out.write(gz.read())
                gz_path.unlink()
                return fasta_path
    except Exception as e:
        print(f"  Warning: Could not download {accession}: {e}")

    return None


def extract_fragments(fasta_path: Path, n_fragments: int, fragment_length: int,
                      species_name: str, rng: random.Random) -> list:
    """Extract random fragments from a genome FASTA."""
    from Bio import SeqIO

    # Read all sequences
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) >= fragment_length:
            sequences.append(record)

    if not sequences:
        print(f"  Warning: No sequences >= {fragment_length}bp in {fasta_path}")
        return []

    fragments = []
    species_tag = species_name.replace(" ", "_")

    for i in range(n_fragments):
        # Pick a random sequence weighted by length
        weights = [len(s.seq) - fragment_length for s in sequences]
        total_w = sum(weights)
        if total_w <= 0:
            continue

        seq_record = rng.choices(sequences, weights=weights, k=1)[0]
        max_start = len(seq_record.seq) - fragment_length
        start = rng.randint(0, max_start)
        end = start + fragment_length

        frag_seq = str(seq_record.seq[start:end])
        frag_id = f"bact_{species_tag}_{i:04d}"
        frag_desc = f"[{species_name}] {seq_record.id}:{start+1}-{end}"

        fragments.append((frag_id, frag_desc, frag_seq))

    return fragments


def build_site(body_site: str, output_base: Path, genome_cache: Path):
    """Build bacterial reference fragments for a single body site."""
    if body_site not in BODY_SITE_BACTERIA:
        print(f"ERROR: Unknown body site '{body_site}'")
        print(f"Available: {', '.join(BODY_SITE_BACTERIA.keys())}")
        return False

    site_info = BODY_SITE_BACTERIA[body_site]
    species_list = site_info["species_list"]

    output_dir = output_base
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_cache.mkdir(parents=True, exist_ok=True)

    rng = random.Random(RANDOM_SEED)

    print(f"\nBuilding bacterial reference fragments: {body_site}")
    print(f"  Description: {site_info['description']}")
    print(f"  Species: {len(species_list)}")
    print(f"  Fragments per species: {FRAGMENTS_PER_SPECIES}")
    print(f"  Fragment length: {FRAGMENT_LENGTH} bp")
    print()

    all_fragments = []
    metadata = []
    failed = []

    for i, bact in enumerate(species_list):
        species = bact["species"]
        accession = bact["accession"]
        print(f"  [{i+1}/{len(species_list)}] {species} ({accession})...", end=" ", flush=True)

        fasta_path = download_genome(accession, genome_cache)
        if fasta_path is None or not fasta_path.exists():
            print("FAILED (download)")
            failed.append(species)
            continue

        fragments = extract_fragments(fasta_path, FRAGMENTS_PER_SPECIES,
                                      FRAGMENT_LENGTH, species, rng)
        if not fragments:
            print(f"FAILED (no fragments)")
            failed.append(species)
            continue

        all_fragments.extend(fragments)
        metadata.append({
            "species": species,
            "accession": accession,
            "abundance": bact["abundance"],
            "prevalence": bact["prevalence"],
            "n_fragments": len(fragments),
            "fragment_length": FRAGMENT_LENGTH,
        })
        print(f"OK ({len(fragments)} fragments)")

    if not all_fragments:
        print(f"  ERROR: No fragments generated for {body_site}")
        return False

    # Write combined FASTA
    fasta_out = output_dir / f"{body_site}_fragments.fasta"
    with open(fasta_out, 'w') as f:
        for frag_id, frag_desc, frag_seq in all_fragments:
            f.write(f">{frag_id} {frag_desc}\n")
            for j in range(0, len(frag_seq), 80):
                f.write(frag_seq[j:j+80] + "\n")

    # Write metadata
    meta_out = output_dir / f"{body_site}_metadata.json"
    with open(meta_out, 'w') as f:
        json.dump({
            "body_site": body_site,
            "host_organism": "human",
            "description": site_info["description"],
            "source": site_info["source"],
            "n_species": len(metadata),
            "n_fragments_total": len(all_fragments),
            "fragment_length": FRAGMENT_LENGTH,
            "species": metadata,
        }, f, indent=2)

    # Write abundance table (TSV)
    abund_out = output_dir / f"{body_site}_abundances.tsv"
    with open(abund_out, 'w') as f:
        f.write("species\taccession\tabundance\tprevalence\n")
        for m in metadata:
            f.write(f"{m['species']}\t{m['accession']}\t{m['abundance']}\t{m['prevalence']}\n")

    print(f"\n  Output:")
    print(f"    FASTA: {fasta_out} ({len(all_fragments)} fragments)")
    print(f"    Metadata: {meta_out}")
    print(f"    Abundances: {abund_out}")
    print(f"    Successfully processed: {len(metadata)}/{len(species_list)} species")
    if failed:
        print(f"    Failed: {', '.join(failed)}")

    fasta_size = os.path.getsize(fasta_out)
    print(f"    FASTA size: {fasta_size / 1024 / 1024:.1f} MB")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Build bacterial reference fragments for ViroForge body sites"
    )
    parser.add_argument(
        "--site", type=str, default=None,
        help="Body site to build (e.g., gut_human, vaginal_human). "
             "Default: build all sites."
    )
    parser.add_argument(
        "--list-sites", action="store_true",
        help="List available body sites and exit"
    )
    args = parser.parse_args()

    if args.list_sites:
        print("Available body sites:")
        for site, info in BODY_SITE_BACTERIA.items():
            n_species = len(info["species_list"])
            print(f"  {site:25s} {n_species:3d} species  {info['description']}")
        return

    output_base = Path("viroforge/data/references/bacterial_fragments")
    genome_cache = output_base / "genome_cache"

    if args.site:
        build_site(args.site, output_base, genome_cache)
    else:
        # Build all sites
        print(f"Building bacterial references for {len(BODY_SITE_BACTERIA)} body sites")
        results = {}
        for site in BODY_SITE_BACTERIA:
            ok = build_site(site, output_base, genome_cache)
            results[site] = "OK" if ok else "FAILED"

        print("\n" + "=" * 60)
        print("Summary:")
        for site, status in results.items():
            print(f"  {site:25s} {status}")


if __name__ == "__main__":
    main()
