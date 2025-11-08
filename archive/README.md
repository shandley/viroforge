# ViroForge Archive

This directory contains archived files from ViroForge development that are no longer actively used but retained for historical reference.

## Archive Structure

### old-tests/
Historical test scripts that have been superseded by the comprehensive test suite in tests/:
- test_read_generation.py - Early end-to-end test script (Oct 2025)
  - Superseded by: tests/test_fastq_integration.py and tests/test_integration_workflow.py

## Archive Policy

Files are archived when:
1. They are superseded by newer implementations
2. They represent early development iterations
3. They may contain useful reference code but are no longer maintained
4. They are kept for audit trail and historical understanding

## Current Testing

Active test suite is located in tests/ directory:
- tests/test_amplification.py - Amplification bias tests
- tests/test_artifacts.py - Platform artifact tests
- tests/test_enrichment.py - VLP enrichment tests
- tests/test_fastq_integration.py - FASTQ workflow integration tests
- tests/test_genome_database.py - Database tests
- tests/test_integration_workflow.py - End-to-end workflow tests
- tests/test_validation.py - Validation framework tests
- tests/test_vlp_contamination.py - VLP/contamination integration tests
- tests/test_vlp_integration.py - VLP integration tests

Run tests with: `pytest tests/ -v`
