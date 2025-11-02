# ViroForge Validation Test Suite

**Purpose**: Pipeline-independent validation of ViroForge-generated datasets

**Audience**: Virome research community - any user generating synthetic datasets

---

## Overview

The validation test suite provides comprehensive quality control for ViroForge-generated FASTQ datasets, independent of any specific analysis pipeline (Hecatomb, Viromeome, etc.). It verifies data integrity, format correctness, and compositional accuracy.

### What It Tests

**9 Independent Validation Tests:**
1. Directory structure
2. Required files presence
3. Metadata integrity
4. FASTA reference validation
5. FASTQ format compliance
6. Abundance calculations
7. Coverage and read counts
8. Quality score distributions
9. Compositional characteristics

---

## Quick Start

### Basic Usage

```bash
# Validate a dataset
python scripts/validate_fastq_dataset.py \
    --dataset data/benchmark_datasets/collection_16_cov10x_novaseq_vlp

# Verbose output with detailed statistics
python scripts/validate_fastq_dataset.py \
    --dataset my_dataset \
    --verbose

# Save validation report to file
python scripts/validate_fastq_dataset.py \
    --dataset my_dataset \
    --report validation_report.txt
```

### Example Output

```
================================================================================
VIROFORGE DATASET VALIDATION REPORT
================================================================================
Dataset: data/benchmark_datasets/collection_16_cov10x_novaseq_vlp

STATUS: ✓ PASSED - All validation tests passed

================================================================================
```

---

## Validation Tests Explained

### Test 1: Directory Structure

**Checks:**
- Dataset directory exists
- Required subdirectories present: `fasta/`, `fastq/`, `metadata/`

**Why It Matters:** Ensures complete dataset with all components

**Failure Example:**
```
✗ Missing directories: fastq, metadata
```

---

### Test 2: Required Files

**Checks:**
- FASTA reference file exists
- R1 and R2 FASTQ files exist
- Metadata JSON file exists
- Composition TSV and abundances files (warnings if missing)

**Why It Matters:** Ensures all necessary files for analysis and validation

**Failure Example:**
```
✗ No R2 FASTQ file found in fastq/
```

---

### Test 3: Metadata Integrity

**Checks:**
- Valid JSON structure
- Required fields present ('genomes' at minimum)
- Genomes list is not empty
- Each genome has required fields: genome_id, length, abundance

**Why It Matters:** Ground truth metadata is essential for benchmarking

**Success Output:**
```
✓ Metadata valid (134 genomes)
```

**Failure Example:**
```
✗ Invalid JSON: Expecting property name enclosed in double quotes
```

---

### Test 4: FASTA Reference Validation

**Checks:**
- FASTA file is parseable
- Contains sequences
- No empty sequences
- Genome IDs extracted correctly

**Why It Matters:** Reference genomes are needed for read mapping validation

**Success Output:**
```
✓ FASTA valid (134 sequences)
```

**Failure Example:**
```
✗ Empty sequences found: GCF_123456789.1, GCF_987654321.1
```

---

### Test 5: FASTQ Format Compliance

**Checks:**
- Valid 4-line record structure (@header, sequence, +, quality)
- Sequence and quality lengths match **for every read**
- No empty sequences
- Read pair counts match (R1 == R2)
- Headers start with '@'
- Plus lines start with '+'

**Why It Matters:** Invalid FASTQ causes pipeline failures and data corruption

**Success Output:**
```
✓ FASTQ format valid (138,688 read pairs)
  Read length: 150.0 bp (range: 150-150)
```

**Failure Examples:**
```
✗ gut_virome_R1.fastq line 1234: Length mismatch (seq=150, qual=149)
✗ Read count mismatch: R1=100,000, R2=99,999
```

**This is the most critical test** - it prevents the exact FASTQ issues previously encountered (sequence/quality length mismatches, truncated files).

---

### Test 6: Abundance Validation

**Checks:**
- All genomes have abundance values
- No negative abundances
- No abundances > 1.0
- Abundances sum to ~1.0 (within 0.99-1.01)
- Genomes in metadata match genomes in FASTA

**Why It Matters:** Ensures ground truth is mathematically valid

**Success Output:**
```
✓ Abundances valid (sum=1.000000)
```

**Failure Examples:**
```
✗ Negative abundance for GCF_123456789.1: -0.05
✗ Abundances sum to 1.234567 (expected ~1.0)
✗ Genomes in FASTA but not metadata: 5
```

---

### Test 7: Coverage and Read Count Validation

**Checks:**
- Read count matches expected from coverage calculation
- Actual coverage calculated from total bases
- Read length matches specification

**Why It Matters:** Verifies generation parameters were applied correctly

**Success Output:**
```
✓ Read count matches coverage (expected=138,688, actual=138,688)
✓ Actual coverage: 10.00x (target: 10.0x)
```

**Warning Example:**
```
⚠ Read count deviation: expected=277,377, actual=138,688 (ratio=0.500)
```

**Note:** Warnings are issued if reads deviate >5% from expected. This may occur due to InSilicoSeq's read counting methodology but doesn't indicate a critical error.

---

### Test 8: Quality Score Validation

**Checks:**
- Quality scores in valid Phred33 range (ASCII 33-126)
- No out-of-range values
- Mean quality calculated

**Why It Matters:** Invalid quality scores break many tools

**Success Output:**
```
✓ Quality scores valid (range: 33-73, mean: 62.5)
```

**Failure Examples:**
```
✗ Quality scores below Phred33 minimum (33): 25
✗ Quality scores above Phred33 maximum (126): 127
```

---

### Test 9: Compositional Validation

**Checks:**
- Shannon diversity calculated
- Diversity reasonable for collection type (warnings if unexpected)
- Family and genus counts tabulated
- Taxonomic representation assessed

**Why It Matters:** Ensures compositional realism

**Success Output:**
```
✓ Shannon diversity: 3.456
  Families: 15
  Genera: 42
```

**Warning Examples:**
```
⚠ High diversity for skin collection: 4.2
⚠ Low diversity for marine collection: 1.8
```

---

## Use Cases

### 1. Post-Generation Quality Control

After generating a dataset, immediately validate it:

```bash
# Generate dataset
python scripts/generate_fastq_dataset.py \
    --collection-id 9 \
    --output gut_virome \
    --coverage 10

# Validate
python scripts/validate_fastq_dataset.py \
    --dataset gut_virome \
    --verbose
```

### 2. Batch Dataset Validation

Validate all generated datasets:

```bash
for dataset in data/benchmark_datasets/*/; do
    echo "Validating: $dataset"
    python scripts/validate_fastq_dataset.py --dataset "$dataset"
done
```

### 3. Pre-Analysis Verification

Before running through your pipeline, verify dataset integrity:

```bash
# Validate
python scripts/validate_fastq_dataset.py \
    --dataset received_dataset \
    --report validation_report.txt

# Check exit code
if [ $? -eq 0 ]; then
    echo "Dataset valid, proceeding with analysis"
    your_pipeline received_dataset
else
    echo "Dataset validation failed, check validation_report.txt"
    exit 1
fi
```

### 4. Continuous Integration Testing

Add to CI/CD pipeline:

```yaml
# .github/workflows/validate.yml
- name: Validate generated datasets
  run: |
    python scripts/validate_fastq_dataset.py \
        --dataset test_dataset \
        --report validation_report.txt
```

---

## What Gets Validated vs. Not Validated

### ✓ Validated

- **File format correctness** - FASTQ, FASTA, JSON structures
- **Data integrity** - No corruption, complete records
- **Mathematical validity** - Abundances, lengths, counts
- **Metadata completeness** - Ground truth fields present
- **Quality score validity** - Proper Phred33 encoding

### ✗ Not Validated

- **Biological realism** - Whether composition reflects nature
- **Pipeline compatibility** - Specific tool requirements
- **Read mapping accuracy** - That reads map to correct genomes
- **Taxonomic correctness** - ICTV taxonomy accuracy
- **Assembly quality** - How well reads assemble

The validation suite focuses on **data quality**, not biological accuracy. It ensures datasets are technically correct and usable, not that they perfectly represent natural viromes.

---

## Interpreting Results

### Exit Codes

- `0` - All tests passed (may have warnings)
- `1` - One or more tests failed (errors present)

### Status Indicators

- `✓ PASSED` - All tests passed, no errors
- `⚠ WARNING` - Tests passed but warnings issued
- `✗ FAILED` - One or more tests failed

### When to Be Concerned

**Critical (Must Fix):**
- Invalid FASTQ format
- Sequence/quality length mismatches
- Missing files
- Corrupted metadata
- Invalid abundance values

**Warnings (Review, May Be Acceptable):**
- Read count deviation (may be ISS behavior)
- Unexpected diversity (may be collection-specific)
- Missing optional files (composition TSV, abundances TXT)

---

## Common Issues and Solutions

### Issue: "Length mismatch (seq=150, qual=149)"

**Cause:** FASTQ record has different sequence and quality lengths

**Fix:** Regenerate dataset - this indicates corruption during generation

**Prevention:** This validation prevents this issue from propagating

---

### Issue: "Genomes in FASTA but not metadata"

**Cause:** Metadata and FASTA got out of sync

**Fix:** Regenerate dataset with correct metadata export

**Prevention:** Atomic generation process ensures sync

---

### Issue: "Read count deviation: ratio=0.500"

**Cause:** InSilicoSeq may count reads differently than expected

**Status:** Warning only - dataset is still usable

**Action:** Note actual coverage, adjust if needed for next generation

---

### Issue: "Abundances sum to 1.234567"

**Cause:** Abundance normalization error

**Fix:** Regenerate with corrected abundance calculation

**Prevention:** Validation catches this before use

---

## Advanced Usage

### Custom Validation Thresholds

Modify `scripts/validate_fastq_dataset.py` to adjust tolerances:

```python
# Line 387: Abundance sum tolerance
if not (0.99 <= total <= 1.01):  # Change to 0.98-1.02 for more tolerance

# Line 439: Read count deviation tolerance
if not (0.95 <= ratio <= 1.05):  # Change to 0.90-1.10 for more tolerance
```

### Integration with Other Tools

**Example: Validate then run Hecatomb**

```bash
#!/bin/bash
if python scripts/validate_fastq_dataset.py --dataset $1; then
    hecatomb run --reads $1/fastq/*_R{1,2}.fastq --outdir results
else
    echo "Validation failed, not running Hecatomb"
    exit 1
fi
```

**Example: Python wrapper**

```python
import subprocess
import sys

def validate_and_process(dataset_path):
    # Validate
    result = subprocess.run(
        ['python', 'scripts/validate_fastq_dataset.py', '--dataset', dataset_path],
        capture_output=True
    )

    if result.returncode == 0:
        print("Dataset valid, processing...")
        # Your processing code here
    else:
        print("Validation failed:")
        print(result.stdout.decode())
        sys.exit(1)
```

---

## For Pipeline Developers

If you're developing a virome analysis pipeline, you can use this validation suite to:

1. **Test input data** - Verify datasets before analysis
2. **Benchmark your pipeline** - Use validated datasets for testing
3. **Develop your own tests** - Extend this suite for pipeline-specific checks

### Extending the Validator

Add custom tests by creating a new method in `DatasetValidator` class:

```python
def validate_custom_check(self) -> bool:
    """Your custom validation logic"""
    # Your code here
    if everything_ok:
        self.log_info("✓ Custom check passed")
        return True
    else:
        self.log_error("✗ Custom check failed")
        return False
```

Then add to `validate_all()`:

```python
def validate_all(self) -> bool:
    # ... existing tests ...

    # Test 10: Custom check
    logger.info("Test 10: Custom validation...")
    if not self.validate_custom_check():
        return False
```

---

## Technical Details

### Performance

- **Mouse Gut (22 genomes, 8K reads)**: <1 second
- **Gut Virome (134 genomes, 139K reads)**: <2 seconds
- **Marine (448 genomes, 620K reads)**: ~5 seconds

Validation is fast enough to run routinely on all generated datasets.

### Sampling Strategy

For performance, some tests use sampling:
- Quality score validation: First 10,000 reads
- Read length statistics: First 1,000 reads

This provides sufficient coverage while maintaining speed.

### File Format Support

- **FASTQ**: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
- **FASTA**: `.fasta`, `.fa`
- **Metadata**: `.json`, `.tsv`, `.txt`

---

## Support

**Questions**: scott.handley@wustl.edu

**Issues**: https://github.com/shandley/viroforge/issues

**Documentation**: See `docs/` directory

---

## Citation

If you use ViroForge validation in your research, please cite:

```
[ViroForge citation will go here upon publication]
```

---

**Last Updated**: 2025-11-01

**Version**: 0.3.0
