---
entry_id: YYYYMMDD-NNN-TYPE-short-description
date: YYYY-MM-DD
type: [DESIGN|IMPLEMENTATION|TESTING|LITERATURE|DECISION|INTEGRATION|STRATEGIC|PUBLICATION|BUGFIX|REVIEW]
status: [in-progress|complete|deferred]
phase: [1|2|3|4]
week: [1-12 for Phase 2]

author: Scott Handley + Claude

references:
  literature:
    - doi:XX.XXXX/journal.citation  # Author Year - Title
  prior_sessions:
    - YYYYMMDD-NNN
  related_topics:
    - topic-name
  code_files:
    - path/to/file.py

tags:
  - tag1
  - tag2

key_decisions:
  - Brief description of critical decision 1
  - Brief description of critical decision 2

commits:
  - abc1234  # Commit message summary

raw_data: raw-data/YYYYMMDD-NNN/
---

# [Descriptive Title]

**Date**: [Month DD, YYYY]
**Phase**: [N], Week [N]
**Status**: [Status]

## Goals

[What we intended to accomplish this session]

## Background

[Context from previous work, literature, or discussions]

## Approach

[How we approached the problem]

## Key Decisions

### Decision 1: [Decision Name]

**Rationale**: [Why this decision was made]

**Approach**:
```python
# Example code or pseudocode
```

**Alternatives Considered**:
- Option A (rejected - reason)
- Option B (rejected - reason)

**Literature Support**:
- [Citation] reports [finding]
- Supports our approach because [reason]

### Decision 2: [...]

## Implementation Notes

[Key implementation details, files created/modified]

## Testing

[Tests written, validation performed, results]

## Literature References

**Key Papers**:
1. Author et al. (Year) Journal - Title
   - Key finding 1
   - Key finding 2

2. Author et al. (Year) Journal - Title
   - Key finding 1

## Results

[What was accomplished, what worked, what didn't]

## Next Steps

- [ ] Immediate follow-up action 1
- [ ] Open question to address
- [ ] Future work item

## Related Documents

- `path/to/doc.md` lines XX-YY (description)
- `lab-notebook/topics/topic-name.md` (updated with this work)

---

**Status**: [Status]
**Impact**: [Brief statement of impact on project]
