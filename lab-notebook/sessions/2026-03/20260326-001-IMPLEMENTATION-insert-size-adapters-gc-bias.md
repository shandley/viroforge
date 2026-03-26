# Insert-Size-Driven Adapters and GC-Biased Duplicates

**Date**: 2026-03-26
**Session Type**: IMPLEMENTATION
**Status**: Complete

---

## Objective

Upgrade adapter and duplicate injection for virome-qc expected range calibration.

## Changes

### Insert-size-driven adapter contamination
- New function `add_insert_size_adapters()` in adapters.py
- Adapter rate is emergent from insert size distribution vs read length
- Chimeric adapter support (internal adapter from ligation events)
- CLI: --mean-insert-size, --insert-size-sd, --chimera-rate

### GC/length-biased PCR duplicates
- Template selection weighted by GC content and fragment length
- Moderate GC and shorter fragments preferentially duplicated
- gc_length_bias parameter (default True)

### Tests
- 7 new tests (54 total), all passing
- Insert size: short inserts -> high adapter rate, long inserts -> low rate
- Chimera injection and header tagging
- GC bias affects template selection
