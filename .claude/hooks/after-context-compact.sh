#!/bin/bash
# Hook: Runs after context window compaction
# Purpose: Re-ground in mission and restore context quickly

cat << 'EOF'

🔄 CONTEXT RESTORED - VIROFORGE SESSION CONTINUES

📋 PROJECT STATUS:
   • Phase 1: ~80% complete (FASTQ generation production-ready)
   • Phase 2: In progress (12-week timeline - see docs/IMPLEMENTATION_PLAN.md)
   • Current focus: Check lab-notebook/INDEX.md for active work

🎯 MISSION:
   • Generate realistic synthetic virome data with complete ground truth
   • Lab-agnostic, community-focused design
   • Publication-quality implementation and documentation
   • Literature-validated biological accuracy

📚 KEY DOCUMENTS:
   • .claude.md - Concise session context
   • lab-notebook/INDEX.md - Current status and recent work
   • docs/IMPLEMENTATION_PLAN.md - Phase 2 detailed plan (12 weeks)
   • lab-notebook/sessions/YYYY-MM/ - Recent session notes
   • lab-notebook/topics/ - Topic evolution and design decisions

🔬 CONTEXT RECOVERY:
   1. Read lab-notebook/INDEX.md "Active Status" section
   2. Check "In-Progress" sessions
   3. Review relevant topic docs in lab-notebook/topics/
   4. Continue from where we left off

CORE PRINCIPLES (DO NOT LOSE):
• Every feature must be lab-agnostic (flexible frameworks, not hard-coded protocols)
• Every parameter must have literature support
• Complete ground truth metadata is non-negotiable
• Documentation is for publication preparation (Nature Biotechnology level)

EOF

# Show recent activity if available
if [ -d "lab-notebook/sessions" ]; then
    recent=$(find lab-notebook/sessions -name "*.md" -type f -mtime -7 2>/dev/null | wc -l | xargs)
    if [ $recent -gt 0 ]; then
        echo ""
        echo "📅 RECENT ACTIVITY (last 7 days): $recent entries"
        find lab-notebook/sessions -name "*.md" -type f -mtime -7 2>/dev/null | sort -r | head -5 | while read file; do
            echo "     • $(basename "$file")"
        done
    fi
fi

echo ""
