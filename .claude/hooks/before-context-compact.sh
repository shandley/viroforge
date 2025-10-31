#!/bin/bash
# Hook: Runs before context window compaction
# Purpose: Ensure critical information is documented before memory loss

cat << 'EOF'

⚠️  CONTEXT COMPACTION IMMINENT ⚠️

BEFORE THIS SESSION ENDS:

📝 Document Important Decisions:
   • Design decisions → lab-notebook/sessions/YYYY-MM/YYYYMMDD-NNN-DESIGN-*.md
   • Literature findings → lab-notebook/literature/*.md
   • Implementation notes → Update relevant topic doc in lab-notebook/topics/

✅ Check Progress:
   • Are design decisions documented? → Session entry with type DESIGN/DECISION
   • Are tests documented? → Session entry with type TESTING
   • Is ground truth validated? → Document validation results
   • Are topic docs updated? → lab-notebook/topics/*.md

🚨 Critical Items:
   • Update lab-notebook/INDEX.md with today's entries
   • Commit lab notebook changes separately
   • Update .claude.md if major progress made
   • Update relevant topic docs in lab-notebook/topics/

💡 Publication Preparation:
   • Document methods-quality details for publication
   • Note any literature support for design decisions
   • Record validation results for Results section

📚 Literature:
   • Document any papers reviewed in lab-notebook/literature/
   • Record parameter ranges and biological justifications
   • Note any conflicts or gaps in literature

EOF

# Check if we have uncommitted lab notebook changes
if git diff --name-only | grep -q "lab-notebook/"; then
    cat << 'WARN'

⚠️  UNCOMMITTED LAB NOTEBOOK CHANGES DETECTED
   Consider committing lab notebook updates before session ends:
   git add lab-notebook/
   git commit -m "docs: update lab notebook [session topic]"

WARN
fi

# Check if INDEX.md is up to date
today=$(date +%Y%m%d)
today_entries=$(find lab-notebook/sessions -name "${today}-*.md" 2>/dev/null | wc -l | xargs)
if [ $today_entries -gt 0 ]; then
    if ! git diff --name-only | grep -q "lab-notebook/INDEX.md"; then
        if ! git diff --cached --name-only | grep -q "lab-notebook/INDEX.md"; then
            cat << 'INDEX'

💡 INDEX.md REMINDER
   Today's entries exist but INDEX.md not modified
   Consider updating lab-notebook/INDEX.md with entry summaries

INDEX
        fi
    fi
fi
