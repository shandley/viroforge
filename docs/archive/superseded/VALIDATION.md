# ViroForge Validation Framework

**Quality Assurance and Data Integrity**

This document describes the validation framework built into ViroForge to ensure data integrity and prevent common errors in sequence generation and file output.

---

## Table of Contents

1. [Overview](#overview)
2. [Why Validation Matters](#why-validation-matters)
3. [Validation Functions](#validation-functions)
4. [Using Validation](#using-validation)
5. [Best Practices](#best-practices)
6. [FASTQ Generation Guidelines](#fastq-generation-guidelines)

---

## Overview

The ViroForge validation framework provides comprehensive quality checks for:

- **Sequence validation** - Ensure only valid DNA characters (ATCGN)
- **Length validation** - Detect mismatches between metadata and actual sequence length
- **Abundance validation** - Verify abundances sum correctly and are in valid range
- **FASTQ validation** - **CRITICAL** for preventing sequence/quality length mismatches
- **File integrity** - Verify files are complete and not truncated
- **GC content** - Check calculated vs expected GC content

All validation is **optional** by default but **strongly recommended** for production use.

---

## Why Validation Matters

### Common FASTQ Errors (That Validation Prevents)

1. **Sequence length ≠ Quality length**
   ```
   @read1
   ATCG        # 4 bases
   +
   III         # 3 quality scores ← ERROR!
   ```
   **Prevention:** `validate_fastq_record()` catches this before writing

2. **File truncation**
   - Process killed mid-write
   - Disk full
   - Buffer not flushed

   **Prevention:** `verify_fastq_file()` confirms file integrity after writing

3. **Invalid characters**
   - Non-DNA characters in sequences
   - Invalid quality scores (ASCII outside 33-126)

   **Prevention:** `validate_sequence()` and `validate_fastq_record()` check characters

---

## Validation Functions

### Sequence Validation

#### `validate_sequence(seq, allow_n=True)`

Validates DNA sequence contains only valid characters.

```python
from viroforge.utils import validate_sequence

# Valid sequences
validate_sequence("ATCG")      # ✅ OK
validate_sequence("ATCGN")     # ✅ OK (N allowed by default)
validate_sequence("atcg")      # ✅ OK (case insensitive)

# Invalid sequences
validate_sequence("ATCGX")     # ❌ ValueError: Invalid DNA characters: {'X'}
validate_sequence("ATCGN", allow_n=False)  # ❌ ValueError
```

#### `validate_sequence_length(genome)`

Validates genome.length matches len(genome.sequence).

```python
from viroforge.utils import validate_sequence_length
from viroforge.core import ViralGenome

genome = ViralGenome("test", "ATCG", "Test")
validate_sequence_length(genome)  # ✅ OK

# Manually break it
genome.length = 100
validate_sequence_length(genome)  # ❌ ValueError: Length mismatch
```

### Abundance Validation

#### `validate_abundances(genomes, tolerance=1e-6)`

Validates abundances sum to approximately 1.0.

```python
from viroforge.utils import validate_abundances

class MockGenome:
    def __init__(self, abundance):
        self.abundance = abundance

genomes = [MockGenome(0.5), MockGenome(0.5)]
validate_abundances(genomes)  # ✅ OK (sum = 1.0)

genomes = [MockGenome(0.5), MockGenome(0.4)]
validate_abundances(genomes, warn_only=True)  # ⚠️ Warning logged
validate_abundances(genomes, warn_only=False) # ❌ ValueError
```

#### `validate_abundance_range(abundance, min_val=0.0, max_val=1.0)`

Validates abundance is within valid range.

```python
from viroforge.utils import validate_abundance_range

validate_abundance_range(0.5)    # ✅ OK
validate_abundance_range(-0.1)   # ❌ ValueError
validate_abundance_range(1.5)    # ❌ ValueError
```

### FASTQ Validation (CRITICAL)

#### `validate_fastq_record(read_id, sequence, quality)`

**MUST be called before writing every FASTQ record!**

```python
from viroforge.utils import validate_fastq_record

# Valid record
validate_fastq_record("read1", "ATCG", "IIII")  # ✅ OK

# Length mismatch - MOST COMMON ERROR
validate_fastq_record("read1", "ATCG", "III")
# ❌ ValueError: FASTQ LENGTH MISMATCH for read1:
#    Sequence: 4 bp
#    Quality:  3 chars
#    These MUST be equal!

# Invalid sequence characters
validate_fastq_record("read1", "ATCGX", "IIIII")
# ❌ ValueError: Invalid DNA characters

# Invalid quality scores
validate_fastq_record("read1", "ATCG", "III\x00")
# ❌ ValueError: Invalid quality score (null character)
```

#### `verify_fastq_file(file_path)`

Verifies FASTQ file integrity after writing.

```python
from viroforge.utils import verify_fastq_file

# After writing FASTQ file
n_reads = verify_fastq_file("output.fastq")
print(f"Verified {n_reads} reads")  # ✅ OK if file is valid

# Detects:
# - Length mismatches
# - Truncated files (incomplete records)
# - Malformed headers (missing @)
# - Malformed separators (missing +)
```

### GC Content Validation

#### `validate_gc_content(genome, tolerance=5.0)`

Validates actual GC content matches expected.

```python
from viroforge.utils import validate_gc_content

genome = ViralGenome("test", "ATCG", "Test")  # 50% GC
validate_gc_content(genome)  # ✅ OK (actual matches calculated)

# Manually set wrong value
genome.gc_content = 10.0
validate_gc_content(genome, warn_only=True)   # ⚠️ Warning
validate_gc_content(genome, warn_only=False)  # ❌ ValueError
```

### File Validation

#### `validate_file_not_empty(file_path)`

Validates file exists and has content.

```python
from viroforge.utils import validate_file_not_empty

validate_file_not_empty("output.fasta")  # ✅ OK
validate_file_not_empty("missing.fasta") # ❌ FileNotFoundError
validate_file_not_empty("empty.fasta")   # ❌ ValueError: File is empty
```

#### `validate_output_directory(dir_path, create=False)`

Validates output directory exists and is writable.

```python
from viroforge.utils import validate_output_directory

validate_output_directory("output", create=True)
# ✅ Creates directory if needed
```

### Batch Validation

#### `validate_genome_collection(genomes, ...)`

Validates entire genome collection at once.

```python
from viroforge.utils import validate_genome_collection

validate_genome_collection(
    genomes,
    check_sequences=True,   # Check all sequences valid
    check_lengths=True,     # Check length metadata
    check_abundances=True,  # Check abundances sum to 1.0
    check_gc=False          # Optionally check GC content
)
```

---

## Using Validation

### Built-in Validation Methods

Both `ViralGenome` and `ContaminantGenome` have `.validate()` methods:

```python
from viroforge.core import ViralGenome

# Create genome
genome = ViralGenome("test", "ATCG", "Test", abundance=0.5)

# Validate it (optional, but recommended)
genome.validate()  # ✅ Runs all default checks

# Selective validation
genome.validate(
    check_sequence=True,   # Check DNA characters
    check_length=True,     # Check length consistency
    check_abundance=True,  # Check abundance range
    check_gc=False         # Skip GC check (faster)
)
```

### Validation in Workflows

#### Development/Testing

Use validation **liberally** during development:

```python
# Create genomes
for i in range(100):
    genome = create_genome(i)
    genome.validate()  # Catch errors immediately

# Validate collections
validate_genome_collection(all_genomes)
```

#### Production Pipelines

**Recommended approach:**

```python
# 1. Validate inputs
validate_file_not_empty(input_fasta)

# 2. Create data
community = create_body_site_profile('gut', n_genomes=100)

# 3. Validate composition
validate_genome_collection(community.genomes)
validate_abundances(community.genomes)

# 4. Generate reads (when implemented)
# CRITICAL: validate EVERY FASTQ record
for read in generated_reads:
    validate_fastq_record(read.id, read.seq, read.qual)  # ← MUST DO THIS
    write_fastq_record(read)

# 5. Verify output
n_reads = verify_fastq_file(output_path)
logger.info(f"Successfully generated and verified {n_reads} reads")
```

---

## Best Practices

### 1. Always Validate FASTQ Records

**CRITICAL:** Never write a FASTQ record without validation:

```python
# ❌ BAD - No validation
for read in reads:
    write_fastq_record(read)

# ✅ GOOD - Validate before writing
for read in reads:
    validate_fastq_record(read.id, read.seq, read.qual)
    write_fastq_record(read)
```

### 2. Verify Files After Writing

```python
# Write file
with open(output_path, 'w') as f:
    for read in reads:
        validate_fastq_record(read.id, read.seq, read.qual)
        f.write(f"@{read.id}\n{read.seq}\n+\n{read.qual}\n")
    f.flush()  # Ensure buffered data is written

# Verify file integrity
n_reads = verify_fastq_file(output_path)
```

### 3. Use Defensive File Writing

```python
import tempfile
import shutil

# Write to temporary file first
with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
    tmp_path = tmp.name
    try:
        for read in reads:
            validate_fastq_record(read.id, read.seq, read.qual)
            tmp.write(f"@{read.id}\n{read.seq}\n+\n{read.qual}\n")
        tmp.flush()

        # Verify temp file
        verify_fastq_file(tmp_path)

        # Atomic rename (only if successful)
        shutil.move(tmp_path, final_path)

    except Exception as e:
        Path(tmp_path).unlink()  # Clean up on error
        raise
```

### 4. Handle Large Datasets Efficiently

```python
# Don't load all reads in memory
BATCH_SIZE = 10000

with open(output_path, 'w') as f:
    for batch_num in range(0, total_reads, BATCH_SIZE):
        reads = generate_read_batch(batch_num, BATCH_SIZE)

        for read in reads:
            validate_fastq_record(read.id, read.seq, read.qual)
            f.write(f"@{read.id}\n{read.seq}\n+\n{read.qual}\n")

        reads.clear()  # Free memory

    f.flush()

# Verify file
verify_fastq_file(output_path)
```

### 5. Use Appropriate Tolerance

```python
# Abundances may not sum to exactly 1.0 due to floating point
validate_abundances(genomes, tolerance=1e-6)  # ✅ Reasonable

# GC content rounding is normal
validate_gc_content(genome, tolerance=5.0, warn_only=True)  # ✅ Warn, don't error
```

---

## FASTQ Generation Guidelines

When implementing FASTQ generation (future feature), follow these guidelines:

### 1. Generate Quality Scores Correctly

```python
def generate_quality_scores(length, mean_quality=30):
    """Generate quality scores matching sequence length."""
    # Quality scores MUST be same length as sequence
    quality = some_distribution(length, mean_quality)

    # Validate
    assert len(quality) == length, "Quality length mismatch!"

    return quality
```

### 2. Handle Read Trimming Carefully

```python
# If trimming reads, trim BOTH sequence AND quality
if trim_read:
    trimmed_seq = seq[:trim_length]
    trimmed_qual = qual[:trim_length]  # ← DON'T FORGET THIS

    # Validate after trimming
    validate_fastq_record(read_id, trimmed_seq, trimmed_qual)
```

### 3. Validate At Multiple Stages

```python
# Stage 1: After generation
read = generate_read(genome)
validate_fastq_record(read.id, read.seq, read.qual)  # ← Here

# Stage 2: After any modifications
if add_errors:
    read = introduce_errors(read)
    validate_fastq_record(read.id, read.seq, read.qual)  # ← Here too

# Stage 3: Before writing
validate_fastq_record(read.id, read.seq, read.qual)  # ← And here
write_fastq_record(read)

# Stage 4: After writing file
verify_fastq_file(output_path)  # ← Final check
```

### 4. Test With Edge Cases

```python
# Test validation with edge cases
test_cases = [
    ("short", "AT", "II"),                    # 2 bp read
    ("long", "A" * 300, "I" * 300),          # 300 bp read
    ("with_n", "ATCGN", "IIIII"),            # N in sequence
    ("low_qual", "ATCG", "!!!!"),            # Phred+33 = 0
    ("high_qual", "ATCG", "~~~~"),           # Phred+33 = 93
]

for read_id, seq, qual in test_cases:
    validate_fastq_record(read_id, seq, qual)
```

---

## Summary

### Required for Production

- ✅ **Always** validate FASTQ records before writing
- ✅ **Always** verify FASTQ files after writing
- ✅ Use defensive file I/O (temp files, atomic rename)
- ✅ Validate abundances sum correctly
- ✅ Check sequences contain only valid characters

### Recommended for Quality

- ✅ Validate genome collections during testing
- ✅ Check length consistency
- ✅ Verify GC content (with warnings)
- ✅ Test with edge cases

### Performance Considerations

- ✅ Validation adds minimal overhead (<1% for most operations)
- ✅ For large datasets, validate in batches
- ✅ Skip GC validation if not needed (saves computation)
- ✅ Use `warn_only=True` for non-critical checks

---

## Need Help?

- **Questions:** Open an issue on GitHub
- **Bug reports:** Include validation error messages
- **Feature requests:** Suggest additional validators

**Remember:** The few milliseconds spent on validation can save hours of debugging corrupted files!

---

*Last updated: 2025-01-30*
*ViroForge validation framework version: 0.1.0*
