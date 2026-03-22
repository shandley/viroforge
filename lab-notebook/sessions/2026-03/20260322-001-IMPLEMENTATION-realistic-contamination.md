# Realistic Contamination Reference Sequences

**Date**: 2026-03-22
**Session Type**: IMPLEMENTATION
**Phase**: Post-13A Enhancement
**Status**: Complete

---

## Objective

Replace synthetic random-GC contamination sequences with real reference sequences detectable by standard QC tools (SortMeRNA, fastp, BBDuk, Kraken2).

## Background

Integration testing of ViroForge with virome-qc revealed that ViroForge's contamination model generates random DNA sequences with correct GC content but no biological k-mer structure. This makes contamination invisible to real k-mer-based QC tools, limiting ViroForge's usefulness for QC pipeline validation.

Specific gaps identified:
1. rRNA contamination uses synthetic sequences, not detectable by rRNA screeners
2. Host DNA uses random sequences, not detectable by host removal tools
3. No adapter read-through modeling for testing adapter trimming
4. PhiX sequence was synthetic instead of the real NC_001422.1

## Implementation Summary

### Bundled Reference Data (`viroforge/data/references/`)

| File | Content | Size |
|------|---------|------|
| `phix174.fasta` | NC_001422.1 genome | 5.4 KB |
| `rrna_representatives.fasta` | 23 rRNA sequences (E. coli, human, gut bacteria, archaea, mouse) | 41 KB |
| `host_fragments.fasta` | 48 human genome fragments from T2T-CHM13v2.0 | 481 KB |
| `adapters.fasta` | 11 Illumina TruSeq/Nextera adapter sequences | <1 KB |

Total bundled data: ~528 KB.

### Reference Resolver (`viroforge/data/references/resolver.py`)

Auto-discovers reference files with priority chain:
1. User-supplied path (explicit argument)
2. Environment variable (VIROFORGE_PHIX_GENOME, VIROFORGE_RRNA_DB, etc.)
3. Bundled references (shipped with package)
4. Synthetic fallback (old behavior, with log warning)

### Contamination Module Changes (`viroforge/core/contamination.py`)

- `add_phix_control()`: Now auto-resolves bundled PhiX174
- `add_rrna_contamination()`: Now auto-resolves bundled rRNA, samples with replacement when bundled set < requested count
- `add_host_contamination()`: Now auto-resolves bundled host fragments or full genome
- `create_contamination_profile()`: New `use_real_references` parameter
- `create_rna_contamination_profile()`: New `use_real_references` parameter

### Adapter Read-Through (`viroforge/simulators/adapters.py`)

Post-processes FASTQ files to inject adapter contamination at 3' ends:
- Supports TruSeq and Nextera adapter types
- Configurable rate (fraction of reads affected)
- Random insert size determines adapter length per read
- High quality scores for adapter bases (Q30-37)

### CLI Args (scripts/generate_fastq_dataset.py)

- `--adapter-rate FLOAT`: Fraction of reads with adapter read-through (default: 0.0)
- `--adapter-type {truseq,nextera}`: Adapter type (default: truseq)
- `--host-genome PATH`: Override bundled host fragments
- `--rrna-database PATH`: Override bundled rRNA
- `--no-real-contaminants`: Force synthetic sequences (old behavior)

### Curation Script (`scripts/curate_reference_sequences.py`)

Downloads and bundles reference FASTA files from NCBI. Serves as provenance record with all accession numbers documented. Idempotent and can regenerate all bundled references.

## Testing

23 new tests in `tests/test_realistic_contamination.py`:
- 8 resolver tests (priority chain, env vars, fallbacks)
- 2 PhiX real vs synthetic tests
- 3 rRNA real vs synthetic tests
- 2 host DNA real vs synthetic tests
- 2 create_contamination_profile tests
- 6 adapter read-through tests (injection, detection, rate, types, lengths)

All 60 contamination-related tests passing (23 new + 16 VLP + 21 RNA).

## Key Design Decisions

1. **Bundled over downloaded**: Reference files ship with the package (~528 KB) rather than requiring a download step. Small enough to commit to git.
2. **Backward compatible**: Default behavior changed (real references used automatically) but `--no-real-contaminants` preserves old behavior for reproducibility.
3. **Adapter rate defaults to 0**: Adapter contamination is opt-in to avoid breaking existing workflows.
4. **Lazy imports**: Resolver is imported inside functions to avoid circular imports and keep the module lightweight.
