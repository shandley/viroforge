#!/bin/bash
# Hook: Runs after user submits prompt, before Claude responds
# Purpose: Reinforce ViroForge philosophy and documentation practices

cat << 'EOF'

🧬 VIROFORGE MISSION 🧬

CORE PHILOSOPHY - Biological Accuracy First:
• Literature-validated - Every simulation parameter must match published virome biology
• Community-focused - Lab-agnostic design that serves the entire field
• Publication-ready - Systematic documentation for Nature Biotechnology-level rigor
• Ground truth - Complete, accurate metadata for benchmarking validation

CRITICAL QUESTIONS:
❓ Is this biologically accurate? (Check literature)
❓ Is this lab-agnostic? (Works for any protocol)
❓ Are we documenting this properly? (Publication preparation)
❓ Is the ground truth accurate? (Validation depends on it)

FOR EVERY FEATURE:
1. ✓ Literature review (what do real viromes look like?)
2. ✓ Design (flexible, composable, lab-agnostic)
3. ✓ Implementation (clean, tested, documented)
4. ✓ Validation (matches literature ranges)
5. ✓ Documentation (methods-quality writeup)
6. ✓ Ground truth (complete metadata)

📖 See docs/IMPLEMENTATION_PLAN.md for Phase 2 details

EOF

# Lab notebook suggestion
USER_MESSAGE="$1"

if echo "$USER_MESSAGE" | grep -qiE "implement|design|test|integrate|validate|literature|review"; then
    today=$(date +%Y%m%d)
    recent_entries=$(find lab-notebook/sessions -name "${today}-*.md" 2>/dev/null | wc -l | xargs)

    if [ $recent_entries -eq 0 ]; then
        cat << 'LABEOF'

💡 LAB NOTEBOOK REMINDER
   Consider creating a lab notebook entry for this work:
   Format: sessions/YYYY-MM/YYYYMMDD-NNN-TYPE-description.md
   Types: DESIGN, IMPLEMENTATION, TESTING, LITERATURE, DECISION, INTEGRATION, STRATEGIC, PUBLICATION

LABEOF
    fi
fi
