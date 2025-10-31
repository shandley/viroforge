# Topic: [Topic Name]

**Status**: [Not Started|In Development|Complete|On Hold]
**Phase**: [1|2|3|4]
**Owner**: Scott Handley + Claude
**Last Updated**: YYYY-MM-DD

---

## Overview

[Brief description of the feature/concept and its importance for ViroForge]

[Why this matters for realistic virome simulation]

---

## Design Evolution

### YYYY-MM-DD: [Session Title] (Session NNN)

**Decision**: [What was decided]

**Rationale**: [Why]

**Literature Support**: [Citations]

### YYYY-MM-DD: [Next Session] (Session NNN)

[Progression of design thinking]

---

## Current Implementation

**Status**: [Development stage]

**Files**:
- `path/to/file.py` (XXX lines)
- `tests/test_file.py` (XXX lines)

**Key Classes/Functions**:
- `ClassName`: Purpose
- `function_name()`: Purpose

**API Example**:
```python
# Example usage
from viroforge.module import Class

obj = Class(
    parameter1=value,
    parameter2=value
)
result = obj.method()
```

**Configuration Options**:
```yaml
# Example YAML configuration
feature:
  parameter1: value
  parameter2: value
```

---

## Literature Support

**Foundational Papers**:

1. **Author et al. (Year) Journal** - Title
   - DOI: XX.XXXX/citation
   - Key finding 1
   - Key finding 2
   - Impact on our design: [explanation]

2. **Author et al. (Year) Journal** - Title
   - DOI: XX.XXXX/citation
   - [Key findings]

**Supporting Evidence**:

- Author (Year): [Finding relevant to implementation]
- Author (Year): [Finding relevant to validation]

**Validation Ranges from Literature**:

| Parameter | Literature Range | ViroForge Implementation | Source |
|-----------|-----------------|--------------------------|--------|
| [Param 1] | [Range] | [Our value] | Author Year |
| [Param 2] | [Range] | [Our value] | Author Year |

---

## Biological Validation

**Expected Behavior**:
- [Biological process 1] should result in [outcome]
- [Biological process 2] should result in [outcome]

**Literature-Based Validation**:
- [Metric]: Expected [range] (Source: [Citation])
- [Metric]: Expected [range] (Source: [Citation])

**Testing Results**:
- [Metric]: Achieved [value] ✅ Within expected range
- [Metric]: Achieved [value] ⚠️ [Explanation if outside range]

---

## Implementation Details

### Key Algorithms

**Algorithm 1: [Name]**
```python
# Pseudocode or key code snippet
def algorithm():
    # Implementation
    pass
```

**Biological Justification**: [Why this algorithm models the real biology]

**Literature Support**: [Citation and finding]

### Key Parameters

| Parameter | Type | Default | Range | Literature Support |
|-----------|------|---------|-------|-------------------|
| param1 | float | 0.95 | 0.8-0.99 | Author Year |
| param2 | int | 40 | 10-50 | Author Year |

---

## Testing and Validation

**Unit Tests**:
- `test_feature_basic()`: Tests basic functionality
- `test_feature_edge_cases()`: Tests boundary conditions
- `test_feature_literature_ranges()`: Validates against literature

**Integration Tests**:
- [Test description]

**Validation Status**:
- ✅ Unit tests passing (XX/XX)
- ✅ Literature ranges validated
- ⏳ Integration tests (in progress)

---

## Open Questions

- [ ] Question 1: [What we still need to figure out]
- [ ] Question 2: [Research needed]
- [x] Question 3: [Resolved - see Session NNN]

---

## Related Sessions

**Design and Planning**:
- Session NNN (YYYYMMDD-NNN-DESIGN-*): Initial design
- Session NNN (YYYYMMDD-NNN-LITERATURE-*): Literature review

**Implementation**:
- Session NNN (YYYYMMDD-NNN-IMPLEMENTATION-*): Core implementation
- Session NNN (YYYYMMDD-NNN-TESTING-*): Validation

**Decisions**:
- Session NNN (YYYYMMDD-NNN-DECISION-*): [Key decision]

---

## Related Topics

- `topic-name-1.md` (interacts with this feature)
- `topic-name-2.md` (depends on this feature)

---

## Publication Notes

**Methods Section Content**:
- [What should go in methods for publication]
- [Algorithm description]
- [Parameter justification]

**Results to Report**:
- [Validation metrics]
- [Comparison to literature]

**Figures Needed**:
- [ ] Figure 1: [Description]
- [ ] Figure 2: [Description]

---

**Version History**:
- v0.1 (YYYY-MM-DD): Initial design (Session NNN)
- v0.2 (YYYY-MM-DD): Implementation complete (Session NNN)
- v1.0 (YYYY-MM-DD): Tested and validated (Session NNN)
