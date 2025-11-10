# Phase 10: Long-Read Sequencing - Integration Architecture

**Date**: 2025-11-10
**Status**: Design Phase
**Implementation**: viroforge/simulators/longread.py

---

## Architecture Overview

The long-read simulator integration follows ViroForge's established patterns:

```
ViroForge Architecture (Existing)
├── simulators/
│   ├── illumina.py          # InSilicoSeq wrapper (short reads)
│   └── longread.py          # PBSIM3 wrapper (NEW)
├── enrichment/
│   └── vlp.py               # VLP enrichment (updated for long reads)
├── workflows/
│   └── rna_virome.py        # RNA workflow (compatible with PBSIM3)
└── core/
    ├── community.py         # Viral community
    └── contamination.py     # Contamination profiles
```

### Integration Points

1. **New Module**: `viroforge/simulators/longread.py`
   - PacBio HiFi simulation
   - Nanopore simulation
   - Follows `illumina.py` design patterns

2. **Updated Modules**:
   - `viroforge/enrichment/vlp.py`: Reduced size bias for long reads
   - `scripts/generate_fastq_dataset.py`: Add `--platform {pacbio-hifi,nanopore}` flags

3. **Dependencies** (add to environment):
   - `pbsim3` (bioconda)
   - `pbccs` (PacBio consensus calling)
   - `samtools` (BAM handling)

---

## Class Design: LongReadSimulator

### Design Pattern

Following `illumina.py` structure:
- Top-level function: `generate_long_reads()`
- Helper functions: `_write_genome_fasta()`, `_run_pbsim3()`, etc.
- Check function: `check_pbsim3_installed()`
- Validation: Ground truth creation
- Cleanup: Temporary file management

### Configuration Classes (Dataclass Pattern)

Following `rna_virome.py` pattern for platform-specific parameters:

```python
from dataclasses import dataclass
from enum import Enum
from typing import Optional

class LongReadPlatform(Enum):
    """Long-read sequencing platforms."""
    PACBIO_HIFI = "pacbio-hifi"
    NANOPORE = "nanopore"

@dataclass
class PacBioHiFiConfig:
    """
    Configuration for PacBio HiFi simulation.

    PacBio HiFi uses circular consensus sequencing (CCS) to generate
    high-accuracy reads (QV20+, >99.9% accuracy) from multiple passes
    of the same molecule.
    """
    passes: int = 10              # Number of passes (default: 10)
    min_passes: int = 3           # Minimum passes for CCS
    accuracy: str = "QSHMM-RSII"  # Error model (QSHMM-RSII, QSHMM-SEQUEL)
    read_length_mean: int = 15000 # Mean read length (bp)
    read_length_sd: int = 5000    # Read length std dev
    clr_error_rate: float = 0.15  # CLR error rate before CCS (15%)

@dataclass
class NanoporeConfig:
    """
    Configuration for Nanopore simulation.

    Oxford Nanopore generates ultra-long reads with characteristic
    homopolymer errors and quality-length relationships.
    """
    chemistry: str = "R10.4"      # Nanopore chemistry version
    read_length_mean: int = 20000 # Mean read length (bp)
    read_length_sd: int = 10000   # Read length std dev
    error_rate: float = 0.05      # Base error rate (5%)
    hp_del_bias: int = 6          # Homopolymer deletion bias
    quality_mean: int = 10        # Mean quality score
```

### Main Function Signature

```python
def generate_long_reads(
    composition,                     # MockViromeComposition object
    output_prefix: str,             # Output file prefix
    platform: LongReadPlatform,     # Platform: PacBio HiFi or Nanopore
    depth: float = 10.0,            # Sequencing depth (coverage)
    platform_config: Optional[Union[PacBioHiFiConfig, NanoporeConfig]] = None,
    validate_output: bool = True,   # Validate FASTQ after generation
    random_seed: Optional[int] = None,
    keep_temp_files: bool = False
) -> Dict[str, Path]:
    """
    Generate long reads from mock virome composition.

    Two-step process for PacBio HiFi:
    1. PBSIM3 generates CLR (continuous long reads) with multi-pass
    2. PacBio ccs generates HiFi consensus reads

    Single-step for Nanopore:
    1. PBSIM3 generates ONT reads with homopolymer errors

    Returns:
        Dict with 'reads', 'ground_truth', optionally 'temp_*' paths
    """
```

---

## Implementation Details

### 1. PacBio HiFi Workflow

**Two-Step Process** (following PBSIM3 design):

```python
def _simulate_pacbio_hifi(
    genomes_fasta: Path,
    abundance_file: Path,
    output_prefix: str,
    depth: float,
    config: PacBioHiFiConfig,
    seed: Optional[int]
) -> Path:
    """
    Simulate PacBio HiFi reads (two-step process).

    Step 1: PBSIM3 generates CLR with multiple passes
    Step 2: PacBio ccs generates HiFi consensus
    """

    # Step 1: Generate CLR with PBSIM3
    clr_bam = _run_pbsim3_clr(
        genomes_fasta=genomes_fasta,
        abundance_file=abundance_file,
        output_prefix=output_prefix,
        depth=depth,
        passes=config.passes,
        error_model=config.accuracy,
        read_length_mean=config.read_length_mean,
        read_length_sd=config.read_length_sd,
        seed=seed
    )

    # Step 2: Generate HiFi consensus with ccs
    hifi_fastq = _run_pacbio_ccs(
        clr_bam=clr_bam,
        output_prefix=output_prefix,
        min_passes=config.min_passes
    )

    return hifi_fastq
```

**PBSIM3 Command (Step 1)**:
```bash
pbsim --strategy wgs \
      --method qshmm \
      --qshmm QSHMM-RSII.model \
      --depth 10 \
      --genome genomes.fasta \
      --pass-num 10 \
      --length-mean 15000 \
      --length-sd 5000 \
      --accuracy-mean 0.85 \
      --prefix output \
      --seed 42
```

**PacBio CCS Command (Step 2)**:
```bash
ccs clr.bam hifi.fastq.gz \
    --min-passes 3 \
    --min-rq 0.99 \
    --log-level INFO
```

### 2. Nanopore Workflow

**Single-Step Process**:

```python
def _simulate_nanopore(
    genomes_fasta: Path,
    abundance_file: Path,
    output_prefix: str,
    depth: float,
    config: NanoporeConfig,
    seed: Optional[int]
) -> Path:
    """
    Simulate Nanopore reads (single-step process).
    """

    # PBSIM3 with ONT error model
    fastq = _run_pbsim3_nanopore(
        genomes_fasta=genomes_fasta,
        abundance_file=abundance_file,
        output_prefix=output_prefix,
        depth=depth,
        chemistry=config.chemistry,
        error_rate=config.error_rate,
        hp_del_bias=config.hp_del_bias,
        read_length_mean=config.read_length_mean,
        read_length_sd=config.read_length_sd,
        seed=seed
    )

    return fastq
```

**PBSIM3 Command (Nanopore)**:
```bash
pbsim --strategy wgs \
      --method errhmm \
      --errhmm ERRHMM-ONT.model \
      --depth 10 \
      --genome genomes.fasta \
      --length-mean 20000 \
      --length-sd 10000 \
      --hp-del-bias 6 \
      --prefix output \
      --seed 42
```

### 3. Helper Functions

Following `illumina.py` patterns:

```python
def check_pbsim3_installed() -> bool:
    """Check if PBSIM3 is installed."""
    try:
        result = subprocess.run(
            ['pbsim', '--version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False

def check_pbccs_installed() -> bool:
    """Check if PacBio ccs is installed (required for HiFi)."""
    try:
        result = subprocess.run(
            ['ccs', '--version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False

def _write_genome_fasta(composition, output_path: Path) -> Dict[str, str]:
    """
    Write genomes to multi-FASTA (reuse from illumina.py).
    Same format as short-read workflow.
    """
    # REUSE EXISTING FUNCTION FROM illumina.py

def _write_abundance_file(composition, output_path: Path) -> None:
    """
    Write abundance file for PBSIM3.

    PBSIM3 uses different format than InSilicoSeq:
    >genome_id abundance
    ATCG...

    NOT the TSV format! We'll need to modify genomes FASTA headers.
    """
    # IMPLEMENTATION NEEDED (different from illumina.py)

def _create_ground_truth_mapping(
    composition,
    output_dir: Path,
    genome_types: Dict[str, str],
    platform: LongReadPlatform
) -> Path:
    """
    Create ground truth file (reuse from illumina.py with platform field).
    """
    # EXTEND EXISTING FUNCTION FROM illumina.py
```

### 4. VLP Modeling Updates

**Key Change**: Long reads less affected by size bias

```python
# viroforge/enrichment/vlp.py

class VLPEnrichment:

    def apply_enrichment(
        self,
        genomes: List[Genome],
        protocol: VLPProtocol,
        read_type: str = "short"  # NEW PARAMETER: "short" or "long"
    ) -> List[Genome]:
        """
        Apply VLP enrichment with read-type-specific bias.

        For short reads (150-300bp):
            - Strong size bias during sequencing
            - Apply full size-based filtration

        For long reads (10-30kb):
            - Size bias only from VLP enrichment step
            - Minimal read-length-dependent bias
            - Larger genomes slightly favored in fragmentation
        """

        if read_type == "short":
            # Existing logic (full size bias)
            return self._apply_short_read_enrichment(genomes, protocol)
        else:  # long reads
            # Reduced size bias (primarily from VLP filtration)
            return self._apply_long_read_enrichment(genomes, protocol)

    def _apply_long_read_enrichment(
        self,
        genomes: List[Genome],
        protocol: VLPProtocol
    ) -> List[Genome]:
        """
        Apply VLP enrichment for long reads.

        Size bias sources:
        1. VLP capsid size filtration (20-200nm typically captured)
        2. Shearing/fragmentation (10-30kb fragments favor larger genomes slightly)
        3. Minimal sequencing bias (unlike short reads)
        """
        # Implementation: reduce size bias by 50-70% compared to short reads
```

---

## File Structure

```
viroforge/
├── simulators/
│   ├── __init__.py          # Export generate_long_reads
│   ├── illumina.py          # Existing short-read simulator
│   └── longread.py          # NEW: Long-read simulator
│       ├── check_pbsim3_installed()
│       ├── check_pbccs_installed()
│       ├── generate_long_reads()        # Main entry point
│       ├── _simulate_pacbio_hifi()      # PacBio HiFi workflow
│       ├── _simulate_nanopore()         # Nanopore workflow
│       ├── _run_pbsim3_clr()           # PBSIM3 CLR generation
│       ├── _run_pbsim3_nanopore()      # PBSIM3 Nanopore
│       ├── _run_pacbio_ccs()           # PacBio ccs consensus
│       ├── _write_genome_fasta()       # Reuse from illumina.py
│       ├── _write_abundance_file()     # PBSIM3-specific format
│       └── _create_ground_truth()      # Extended from illumina.py
└── enrichment/
    └── vlp.py               # Updated: apply_enrichment(read_type="long")
```

---

## Integration with generate_fastq_dataset.py

**Add Platform Selection**:

```python
# scripts/generate_fastq_dataset.py

def main():
    parser = argparse.ArgumentParser(...)

    # EXISTING
    parser.add_argument(
        '--platform',
        choices=['novaseq', 'miseq', 'hiseq', 'nextseq',
                 'pacbio-hifi', 'nanopore'],  # NEW OPTIONS
        default='novaseq',
        help='Sequencing platform'
    )

    # NEW: Long-read specific parameters
    lr_group = parser.add_argument_group('long-read options')
    lr_group.add_argument(
        '--depth',
        type=float,
        default=10.0,
        help='Sequencing depth (coverage) for long reads (default: 10x)'
    )
    lr_group.add_argument(
        '--pacbio-passes',
        type=int,
        default=10,
        help='Number of passes for PacBio HiFi (default: 10)'
    )
    lr_group.add_argument(
        '--ont-chemistry',
        choices=['R9.4', 'R10.4'],
        default='R10.4',
        help='Nanopore chemistry version (default: R10.4)'
    )

    # ... existing arguments ...

    args = parser.parse_args()

    # Route to appropriate simulator
    if args.platform in ['pacbio-hifi', 'nanopore']:
        # Use long-read simulator
        from viroforge.simulators.longread import generate_long_reads

        # Determine platform
        platform = LongReadPlatform.PACBIO_HIFI if args.platform == 'pacbio-hifi' \
                   else LongReadPlatform.NANOPORE

        # Generate reads
        output = generate_long_reads(
            composition=composition,
            output_prefix=args.output,
            platform=platform,
            depth=args.depth,
            random_seed=args.seed
        )
    else:
        # Use existing short-read simulator
        from viroforge.simulators.illumina import generate_reads
        # ... existing logic ...
```

---

## Ground Truth Format

Extended from short-read format to include platform-specific metadata:

```tsv
genome_id       genome_type     taxonomy        length  gc_content      abundance       source          platform        read_type
NC_123456       viral           Norovirus       7500    0.48            0.15            viral_community pacbio-hifi     long
NC_789012       bacterial       E. coli K12     4600000 0.51            0.05            contamination   pacbio-hifi     long
```

---

## Error Handling

Following established patterns:

```python
def generate_long_reads(...):
    """Generate long reads with comprehensive error handling."""

    # 1. Check dependencies
    if not check_pbsim3_installed():
        raise RuntimeError(
            "PBSIM3 not installed. Install with: conda install bioconda::pbsim3"
        )

    if platform == LongReadPlatform.PACBIO_HIFI and not check_pbccs_installed():
        raise RuntimeError(
            "PacBio ccs not installed. Install with: conda install bioconda::pbccs"
        )

    # 2. Validate inputs
    if depth <= 0:
        raise ValueError(f"Depth must be positive, got {depth}")

    # 3. Run with try/finally for cleanup
    temp_dir = tempfile.mkdtemp(prefix='viroforge_longread_')
    try:
        # ... simulation logic ...
    finally:
        if not keep_temp_files:
            shutil.rmtree(temp_dir)
```

---

## Testing Strategy

### Unit Tests

```python
# tests/test_longread_simulator.py

def test_check_pbsim3_installed():
    """Test PBSIM3 installation check."""

def test_pacbio_hifi_config_defaults():
    """Test PacBio HiFi config defaults."""

def test_nanopore_config_defaults():
    """Test Nanopore config defaults."""

def test_write_abundance_file():
    """Test PBSIM3-specific abundance file format."""
```

### Integration Tests

```python
def test_generate_pacbio_hifi_reads():
    """Test complete PacBio HiFi workflow (dry run if tools not installed)."""

def test_generate_nanopore_reads():
    """Test complete Nanopore workflow (dry run if tools not installed)."""
```

---

## Documentation Updates

1. **README.md**: Add long-read examples
2. **Tutorial**: `docs/LONGREAD_TUTORIAL.md`
3. **API docs**: Document new functions
4. **Dependencies**: Update environment.yml

---

## Timeline

### Week 1-2: PacBio HiFi Implementation
- Day 1-2: Create `longread.py` with PacBio HiFi support
- Day 3-4: Integrate with `generate_fastq_dataset.py`
- Day 5-7: Testing, debugging, validation

### Week 2-3: Nanopore Implementation
- Day 8-9: Add Nanopore simulation to `longread.py`
- Day 10-11: Update VLP modeling for long reads
- Day 12-14: Testing, integration tests

### Week 3-4: Documentation & Polish
- Day 15-17: Write tutorials and documentation
- Day 18-19: Create examples and benchmarks
- Day 20-21: Final testing, bug fixes

---

## Success Criteria

- ✅ PBSIM3 wrapper functional for both platforms
- ✅ PacBio HiFi generates >99.9% accuracy reads
- ✅ Nanopore generates reads with characteristic errors
- ✅ VLP modeling updated for long reads
- ✅ Integration tests pass
- ✅ Documentation complete
- ✅ Tutorial examples work end-to-end

---

## Next Steps

1. **Implement** `viroforge/simulators/longread.py`
2. **Update** `viroforge/enrichment/vlp.py` for long-read support
3. **Integrate** with `scripts/generate_fastq_dataset.py`
4. **Test** with real collections
5. **Document** usage and examples
