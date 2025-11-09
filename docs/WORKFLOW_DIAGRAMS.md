# ViroForge Workflow Diagrams

**Visual Guide to ViroForge Data Flow**

This document provides visual representations of how data flows through ViroForge, from genome database to final FASTQ files.

---

## Table of Contents

1. [Overview Workflow](#overview-workflow)
2. [FASTQ Generation Workflow](#fastq-generation-workflow)
3. [VLP Enrichment Workflow](#vlp-enrichment-workflow)
4. [Contamination Processing](#contamination-processing)
5. [Database-to-FASTQ Pipeline](#database-to-fastq-pipeline)

---

## Overview Workflow

### High-Level ViroForge Process

```
┌─────────────────────────────────────────────────────────────────┐
│                    VIROFORGE WORKFLOW                           │
└─────────────────────────────────────────────────────────────────┘

┌────────────────┐
│  RefSeq Viral  │  14,423 genomes
│    Genomes     │  + ICTV taxonomy
└────────┬───────┘
         │
         ▼
┌────────────────┐
│ SQLite Database│  Indexed, annotated
│  + Collections │  8 body site collections
└────────┬───────┘
         │
         ▼
┌────────────────┐
│ Generate FASTQ │  Command: scripts/generate_fastq_dataset.py
│   Script       │  Parameters: collection, coverage, protocol
└────────┬───────┘
         │
         ▼
┌────────────────────────────────────────────┐
│        VIRAL COMMUNITY COMPOSITION         │
│  - Load genomes from database              │
│  - Extract sequences + taxonomy            │
│  - Assign relative abundances              │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│        CONTAMINATION PROFILE               │
│  - Add host DNA (NCBI genomes)             │
│  - Add rRNA (ribosomal sequences)          │
│  - Add bacteria (gut microbiome)           │
│  - Add PhiX control                        │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│          MOCK COMPOSITION                  │
│  Viral fraction: 92.6% (realistic)         │
│  Contamination: 7.4%                       │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│        VLP ENRICHMENT (Optional)           │
│  - Estimate virion sizes                   │
│  - Apply size-based filtration             │
│  - Apply nuclease treatment                │
│  - Reduce contamination by type            │
│  Result: 99.3% viral (tangential flow)     │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│          PREPARE FASTA FOR ISS             │
│  - Write genomes with abundances           │
│  - Format: >genome_id|abundance            │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│        INSILICOSEQ FASTQ GENERATION        │
│  - Generate paired-end reads               │
│  - Platform-specific error models          │
│  - Coverage-based read count               │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│             OUTPUT FILES                   │
│  ├── fasta/*.fasta (reference)             │
│  ├── fastq/*_R1.fastq (forward reads)      │
│  ├── fastq/*_R2.fastq (reverse reads)      │
│  └── metadata/*.json (ground truth)        │
└────────────────────────────────────────────┘
```

---

## FASTQ Generation Workflow

### From Collection to Reads

```
INPUT: Collection ID (9 = Gut Virome, 134 genomes)

┌───────────────────────────────────────────────────────────┐
│  STEP 1: Load Collection from Database                   │
│                                                           │
│  SELECT g.genome_id, g.sequence, g.length, g.gc_content, │
│         t.family, cg.relative_abundance                  │
│  FROM genomes g                                          │
│  JOIN collection_genomes cg ON g.genome_id = cg.genome_id│
│  JOIN taxonomy t ON g.genome_id = t.genome_id           │
│  WHERE cg.collection_id = 9                              │
└─────────────┬─────────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 2: Add Contamination Profile                     │
│                                                         │
│  Based on --contamination-level:                       │
│  ┌──────────────┬─────────────┬────────────────────┐  │
│  │ Level        │ Initial %   │ Components         │  │
│  ├──────────────┼─────────────┼────────────────────┤  │
│  │ clean        │ 0.7%        │ 0.3% host DNA      │  │
│  │              │             │ 0.2% bacteria      │  │
│  │              │             │ 0.1% rRNA          │  │
│  │              │             │ 0.1% PhiX          │  │
│  ├──────────────┼─────────────┼────────────────────┤  │
│  │ realistic    │ 7.4%        │ 3.0% host DNA      │  │
│  │              │             │ 2.5% bacteria      │  │
│  │              │             │ 1.5% rRNA          │  │
│  │              │             │ 0.4% PhiX          │  │
│  ├──────────────┼─────────────┼────────────────────┤  │
│  │ heavy        │ 27.1%       │ 12% host DNA       │  │
│  │              │             │ 10% bacteria       │  │
│  │              │             │ 4% rRNA            │  │
│  │              │             │ 1.1% PhiX          │  │
│  └──────────────┴─────────────┴────────────────────┘  │
│                                                         │
│  Result: Mixed composition (viral + contaminants)      │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 3: Apply VLP Enrichment (if --vlp-protocol set)  │
│                                                         │
│  For each sequence:                                     │
│    1. Estimate virion size from genome properties      │
│    2. Calculate retention probability (size-based)     │
│    3. Apply nuclease treatment (if not encapsidated)   │
│    4. Reduce abundance based on protocol efficiency    │
│                                                         │
│  Tangential Flow Example:                              │
│  ┌─────────────────┬──────────┬──────────┬──────────┐ │
│  │ Sequence Type   │ Before   │ After    │ Change   │ │
│  ├─────────────────┼──────────┼──────────┼──────────┤ │
│  │ Viral genomes   │ 92.6%    │ 99.3%    │ +7.2%    │ │
│  │ Host DNA        │ 3.0%     │ 0.11%    │ -96.4%   │ │
│  │ rRNA            │ 1.5%     │ 0.16%    │ -89.3%   │ │
│  │ Bacteria        │ 2.5%     │ 0.03%    │ -98.8%   │ │
│  │ PhiX            │ 0.4%     │ 0.24%    │ -40.0%   │ │
│  └─────────────────┴──────────┴──────────┴──────────┘ │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 4: Calculate Read Count                          │
│                                                         │
│  total_genome_length = sum(genome.length * abundance)  │
│  n_reads = (total_length * coverage) / (2 * read_len)  │
│                                                         │
│  Example (Gut, 10x, 150bp reads):                      │
│    total_length ≈ 3,500,000 bp                         │
│    n_reads = (3.5M * 10) / (2 * 150) = 116,667 pairs  │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 5: Write Abundance File for ISS                  │
│                                                         │
│  Format: >genome_id abundance                          │
│                                                         │
│  Example:                                               │
│  >GCF_000001405.40 0.0234                              │
│  ATCGATCGATCG...                                       │
│  >GCF_000002435.3 0.0156                               │
│  GCTAGCTAGCTA...                                       │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 6: Run InSilicoSeq                               │
│                                                         │
│  Command:                                               │
│  iss generate                                           │
│    --genomes reference.fasta                            │
│    --abundance abundances.txt                           │
│    --model {novaseq, miseq, hiseq}                      │
│    --n_reads 116667                                     │
│    --output output_prefix                               │
│                                                         │
│  ISS generates:                                         │
│    - output_prefix_R1.fastq (forward reads)             │
│    - output_prefix_R2.fastq (reverse reads)             │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 7: Export Ground Truth Metadata                  │
│                                                         │
│  Creates 3 files:                                       │
│  1. metadata.json - Complete information               │
│     {                                                   │
│       "generation_info": {...},                         │
│       "collection": {...},                              │
│       "configuration": {...},                           │
│       "enrichment_statistics": {...},                   │
│       "sequences": [...]                                │
│     }                                                   │
│                                                         │
│  2. composition.tsv - Abundance table                   │
│     genome_id | name | abundance | family | ...        │
│                                                         │
│  3. abundances.txt - ISS format                         │
│     genome_id abundance                                 │
└─────────────────────────────────────────────────────────┘

OUTPUT: Complete FASTQ dataset with ground truth
```

---

## VLP Enrichment Workflow

### Phase 5 Enhanced VLP Modeling

```
INPUT: Mock Composition (92.6% viral, 7.4% contamination)

┌─────────────────────────────────────────────────────────┐
│  FOR EACH SEQUENCE in Composition                      │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 1: Classify Sequence Type                        │
│                                                         │
│  If viral genome:                                       │
│    → Estimate virion size                               │
│    → Apply size-based filtration                        │
│    → Protected from nuclease (encapsidated)             │
│                                                         │
│  If contaminant:                                        │
│    → Determine contaminant type                         │
│    → Apply type-specific reduction                      │
│    → May be nuclease-susceptible                        │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼ (viral branch)
┌─────────────────────────────────────────────────────────┐
│  VIRAL SEQUENCE PROCESSING                              │
│                                                         │
│  1. Estimate Virion Diameter                            │
│    ┌────────────────────────────────────────────────┐  │
│    │ Genome Type │ Formula                          │  │
│    ├────────────────────────────────────────────────┤  │
│    │ dsDNA       │ 2.34 + 0.48*log10(length)       │  │
│    │ ssDNA       │ 2.34 + 0.48*log10(length) - 0.2 │  │
│    │ ssRNA       │ 2.34 + 0.48*log10(length) - 0.3 │  │
│    │ dsRNA       │ 2.34 + 0.48*log10(length)       │  │
│    └────────────────────────────────────────────────┘  │
│                                                         │
│    Example: 40kb dsDNA phage                            │
│    diameter = 2.34 + 0.48*log10(40000) = 4.57 (72nm)   │
│                                                         │
│  2. Apply Filtration Retention Curve                    │
│                                                         │
│    Tangential Flow (0.2 μm pore, sigmoid curve):       │
│    retention = 1 / (1 + exp(-k*(diameter - 3.5)))      │
│    where k = 0.008 (steepness)                         │
│                                                         │
│    ┌─────────────────────────────────────┐            │
│    │ Diameter (log nm) │ Retention       │            │
│    ├─────────────────────────────────────┤            │
│    │ 1.5 (27 nm, PhiX) │ 60%             │            │
│    │ 2.0 (100 nm)      │ 90%             │            │
│    │ 2.5 (316 nm)      │ 98%             │            │
│    │ 3.0 (1 μm)        │ 99.5%           │            │
│    └─────────────────────────────────────┘            │
│                                                         │
│  3. Adjust Abundance                                    │
│    new_abundance = old_abundance * retention           │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼ (contamination branch)
┌─────────────────────────────────────────────────────────┐
│  CONTAMINATION REDUCTION                                │
│                                                         │
│  Type-Specific Reduction Factors:                       │
│  ┌──────────────┬──────────────────────────────────┐  │
│  │ Contam Type  │ Tangential Flow Efficiency       │  │
│  ├──────────────┼──────────────────────────────────┤  │
│  │ Host DNA     │ 96.4% removal                    │  │
│  │ (free)       │ Mechanism: DNase digestion       │  │
│  │              │ Efficiency: 98% nuclease         │  │
│  │              │ + 50% from filtration            │  │
│  ├──────────────┼──────────────────────────────────┤  │
│  │ rRNA         │ 89.3% removal                    │  │
│  │ (debris)     │ Mechanism: RNase + filtration    │  │
│  │              │ Efficiency: 90% nuclease         │  │
│  │              │ Partial membrane protection      │  │
│  ├──────────────┼──────────────────────────────────┤  │
│  │ Bacteria     │ 98.8% removal                    │  │
│  │ (cells)      │ Mechanism: Size filtration       │  │
│  │              │ Bacteria 1-5 μm, pore 0.2 μm    │  │
│  │              │ 30% lysed → nuclease-susceptible │  │
│  ├──────────────┼──────────────────────────────────┤  │
│  │ PhiX         │ 40% loss                         │  │
│  │ (27 nm)      │ Mechanism: Small virus           │  │
│  │              │ Treated like viral genome        │  │
│  │              │ Size-based retention only        │  │
│  └──────────────┴──────────────────────────────────┘  │
│                                                         │
│  Stochastic Variation:                                  │
│    reduction_factor = protocol_factor * N(1.0, 0.1)    │
│    (10% coefficient of variation)                       │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  RENORMALIZE ABUNDANCES                                 │
│                                                         │
│  total = sum(all abundances after reduction)           │
│  for each sequence:                                     │
│    final_abundance = abundance / total                  │
│                                                         │
│  Result: Abundances sum to 1.0                          │
└─────────────┬───────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────────────────────┐
│  CALCULATE ENRICHMENT STATISTICS                        │
│                                                         │
│  viral_fraction_before = sum(viral abundances before)  │
│  viral_fraction_after = sum(viral abundances after)    │
│  viral_enrichment = after / before                      │
│                                                         │
│  For each contaminant type:                             │
│    reduction_pct = 1 - (after / before)                │
│                                                         │
│  Example Output:                                        │
│  {                                                      │
│    "viral_fraction_before": 0.926,                      │
│    "viral_fraction_after": 0.993,                       │
│    "contamination_reduction": {                         │
│      "overall_reduction": 0.912,                        │
│      "host_dna_removal": 0.964,                         │
│      "rrna_removal": 0.893,                             │
│      "bacteria_removal": 0.988,                         │
│      "phix_retention": 0.600                            │
│    }                                                    │
│  }                                                      │
└─────────────────────────────────────────────────────────┘

OUTPUT: Enriched Composition (99.3% viral, 0.7% contamination)
```

---

## Contamination Processing

### How Contamination is Added and Reduced

```
┌─────────────────────────────────────────────────────────┐
│  CONTAMINATION LIFECYCLE                                │
└─────────────────────────────────────────────────────────┘

STAGE 1: Viral Collection (Pure)
┌────────────────────────────────┐
│  Gut Virome Collection         │
│  134 viral genomes             │
│  100% viral                    │
│  Total: ~3.5 Mb                │
└─────────────┬──────────────────┘
              │
              ▼
STAGE 2: Add Contamination (Mock Sample)
┌────────────────────────────────────────────────────────┐
│  Contamination Sources:                                │
│                                                        │
│  1. Host DNA (3.0% of total)                           │
│     Source: NCBI human genome (chr1 random segments)  │
│     Size: ~120 kb total                                │
│     GC: 40-42%                                         │
│                                                        │
│  2. Reagent Bacteria (2.5%)                            │
│     Species: E. coli, Pseudomonas, Bacillus           │
│     Source: RefSeq bacterial genomes                   │
│     Size: ~100 kb total                                │
│                                                        │
│  3. rRNA (1.5%)                                        │
│     Types: 16S, 18S, 23S, 28S                          │
│     Source: SILVA/RDP databases                        │
│     Size: ~60 kb total                                 │
│                                                        │
│  4. PhiX Control (0.4%)                                │
│     Genome: PhiX174 complete (5,386 bp)                │
│     Source: NCBI NC_001422.1                           │
│     Encapsidated (27 nm virion)                        │
└─────────────┬──────────────────────────────────────────┘
              │
              ▼
STAGE 3: Mock Composition (Pre-VLP)
┌────────────────────────────────┐
│  Total composition:            │
│  ├─ Viral: 92.6%               │
│  ├─ Host DNA: 3.0%             │
│  ├─ Bacteria: 2.5%             │
│  ├─ rRNA: 1.5%                 │
│  └─ PhiX: 0.4%                 │
│                                │
│  Represents "input" to VLP     │
│  extraction in real sample     │
└─────────────┬──────────────────┘
              │
              ▼ (if VLP enrichment applied)
STAGE 4: VLP Processing
┌────────────────────────────────────────────────────────┐
│  Protocol: Tangential Flow Filtration                 │
│                                                        │
│  Step 1: Filtration (0.2 μm pore)                     │
│  ┌──────────────────────────────────────────────────┐ │
│  │ Component    │ Before   │ Filtration  │ After   │ │
│  ├──────────────────────────────────────────────────┤ │
│  │ Viral        │ 92.6%    │ 85% retained│ 94.8%   │ │
│  │ Host DNA     │ 3.0%     │ 50% removed │ 1.8%    │ │
│  │ Bacteria     │ 2.5%     │ 95% removed │ 0.15%   │ │
│  │ rRNA         │ 1.5%     │ 60% removed │ 0.72%   │ │
│  │ PhiX         │ 0.4%     │ 40% removed │ 0.29%   │ │
│  └──────────────────────────────────────────────────┘ │
│                                                        │
│  Step 2: Nuclease Treatment (DNase/RNase)             │
│  ┌──────────────────────────────────────────────────┐ │
│  │ Component    │ Before   │ Nuclease    │ After   │ │
│  ├──────────────────────────────────────────────────┤ │
│  │ Viral        │ 94.8%    │ Protected   │ 99.3%   │ │
│  │ Host DNA     │ 1.8%     │ 98% removed │ 0.04%   │ │
│  │ Bacteria     │ 0.15%    │ 30% lysed   │ 0.03%   │ │
│  │ rRNA         │ 0.72%    │ 90% removed │ 0.16%   │ │
│  │ PhiX         │ 0.29%    │ Protected   │ 0.29%   │ │
│  └──────────────────────────────────────────────────┘ │
└─────────────┬──────────────────────────────────────────┘
              │
              ▼
STAGE 5: VLP-Enriched Composition
┌────────────────────────────────┐
│  Final composition:            │
│  ├─ Viral: 99.35%              │
│  ├─ Host DNA: 0.11%            │
│  ├─ Bacteria: 0.03%            │
│  ├─ rRNA: 0.16%                │
│  └─ PhiX: 0.24%                │
│                                │
│  Total contamination: 0.54%    │
│  Reduction: 92.7%              │
│                                │
│  ViromeQC score: 0.993 (PASS)  │
└────────────────────────────────┘
```

---

## Database-to-FASTQ Pipeline

### Complete End-to-End Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                  DATABASE LAYER                                 │
└─────────────────────────────────────────────────────────────────┘

    ┌──────────────────────────────────────────────────────────┐
    │  viroforge/data/viral_genomes.db (SQLite)                │
    │                                                          │
    │  Tables:                                                 │
    │  ┌────────────────┬──────────┬────────────────────────┐ │
    │  │ genomes        │ 14,423   │ RefSeq viral genomes   │ │
    │  │ taxonomy       │ 14,423   │ ICTV taxonomy (53.9%)  │ │
    │  │ collections    │ 8        │ Body site collections  │ │
    │  │ coll_genomes   │ 1,198    │ Collection membership  │ │
    │  └────────────────┴──────────┴────────────────────────┘ │
    └────────────────────────┬─────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  COLLECTION LOADER                              │
│  (scripts/generate_fastq_dataset.py::CollectionLoader)         │
└─────────────────────────────────────────────────────────────────┘

    Load Collection
    ───────────────────────────────────────────────────────────
    SELECT g.genome_id, g.name, g.sequence, g.length,
           g.gc_content, t.family, t.genus, t.species,
           cg.relative_abundance
    FROM genomes g
    JOIN collection_genomes cg ON g.genome_id = cg.genome_id
    JOIN taxonomy t ON g.genome_id = t.genome_id
    WHERE cg.collection_id = ?
    ───────────────────────────────────────────────────────────
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  COMPOSITION BUILDER                            │
│  (scripts/generate_fastq_dataset.py::FASTQGenerator)           │
└─────────────────────────────────────────────────────────────────┘

    Create Viral Community
    ──────────────────────────────────────────────────────
    viral_genomes = []
    for row in db_results:
        genome = {
            'genome_id': row['genome_id'],
            'sequence': row['sequence'],
            'abundance': row['relative_abundance'],
            'taxonomy': {...}
        }
        viral_genomes.append(genome)
    ──────────────────────────────────────────────────────
                             │
                             ▼
    Add Contamination Profile
    ──────────────────────────────────────────────────────
    contamination = create_contamination_profile(
        level='realistic',  # 7.4% total
        body_site='gut'
    )

    # Combine viral + contamination
    all_sequences = viral_genomes + contamination
    ──────────────────────────────────────────────────────
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  VLP ENRICHMENT MODULE                          │
│  (viroforge/enrichment/vlp.py::VLPEnrichment)                   │
└─────────────────────────────────────────────────────────────────┘

    Apply VLP Protocol
    ──────────────────────────────────────────────────────
    vlp = VLPEnrichment(
        protocol=VLPProtocol.tangential_flow_standard()
    )

    enriched_sequences, stats = vlp.apply_enrichment(
        all_sequences
    )
    ──────────────────────────────────────────────────────
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  FASTA WRITER                                   │
│  (scripts/generate_fastq_dataset.py::write_fasta)              │
└─────────────────────────────────────────────────────────────────┘

    Write Reference FASTA
    ──────────────────────────────────────────────────────
    with open('collection.fasta', 'w') as f:
        for seq in enriched_sequences:
            f.write(f">{seq['genome_id']}\n")
            f.write(f"{seq['sequence']}\n")
    ──────────────────────────────────────────────────────
                             │
    Write ISS Abundance File
    ──────────────────────────────────────────────────────
    with open('abundances.txt', 'w') as f:
        for seq in enriched_sequences:
            f.write(f"{seq['genome_id']} {seq['abundance']}\n")
    ──────────────────────────────────────────────────────
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  INSILICOSEQ WRAPPER                            │
│  (viroforge/simulators/illumina.py)                             │
└─────────────────────────────────────────────────────────────────┘

    Call ISS
    ──────────────────────────────────────────────────────
    cmd = [
        'iss', 'generate',
        '--genomes', 'collection.fasta',
        '--abundance_file', 'abundances.txt',
        '--model', 'novaseq',
        '--n_reads', '116667',
        '--output', 'output_prefix'
    ]
    subprocess.run(cmd)
    ──────────────────────────────────────────────────────
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  OUTPUT FILES                                   │
└─────────────────────────────────────────────────────────────────┘

    Directory Structure:
    ──────────────────────────────────────────────────────
    output_dir/
    ├── fasta/
    │   └── collection_9.fasta                 (3.5 MB)
    ├── fastq/
    │   ├── collection_9_R1.fastq             (35 MB)
    │   └── collection_9_R2.fastq             (35 MB)
    └── metadata/
        ├── collection_9_metadata.json        (450 KB)
        ├── collection_9_composition.tsv      (25 KB)
        └── collection_9_abundances.txt       (5 KB)
    ──────────────────────────────────────────────────────

COMPLETE: Ready for pipeline benchmarking!
```

---

## Summary

These diagrams illustrate:

1. **Overview**: High-level data flow from database to FASTQ
2. **FASTQ Generation**: Step-by-step process with examples
3. **VLP Enrichment**: Detailed biological modeling
4. **Contamination**: How contamination is added and reduced
5. **Database Pipeline**: Complete technical implementation

**Key Insight**: ViroForge provides complete traceability from RefSeq genomes through VLP enrichment to final reads, with every step validated against literature.

---

**For More Details**:
- [FASTQ Generation Guide](PHASE4_FASTQ_GENERATION.md)
- [VLP Protocol Comparison](VLP_PROTOCOL_COMPARISON_TUTORIAL.md)
- [VLP Integration Guide](VLP_CONTAMINATION_INTEGRATION.md)
