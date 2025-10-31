#!/bin/bash
# Hook: Runs at the start of each Claude Code session
# Purpose: Display lab notebook status and current phase

cat << 'EOF'

📔 VIROFORGE LAB NOTEBOOK
EOF

if [ -d "lab-notebook" ]; then
    # Count entries
    total_entries=$(find lab-notebook/sessions -name "*.md" 2>/dev/null | wc -l | xargs)
    echo "   Total entries: $total_entries"

    # Today's entries
    today=$(date +%Y%m%d)
    today_entries=$(find lab-notebook/sessions -name "${today}-*.md" 2>/dev/null | wc -l | xargs)
    if [ $today_entries -gt 0 ]; then
        echo "   Today's entries: $today_entries"
        find lab-notebook/sessions -name "${today}-*.md" 2>/dev/null | while read file; do
            echo "     • $(basename "$file")"
        done
    else
        echo "   Today's entries: 0"
    fi

    # Current phase - extract from INDEX.md
    if [ -f "lab-notebook/INDEX.md" ]; then
        current_phase=$(grep -E "^\*\*Phase [0-9]" lab-notebook/INDEX.md | head -1 2>/dev/null)
        if [ -n "$current_phase" ]; then
            echo ""
            echo "🚧 $current_phase"
        fi
    fi

    # In-progress sessions
    in_progress=$(grep -l "^status: in-progress" lab-notebook/sessions/**/*.md 2>/dev/null | wc -l | xargs)
    if [ $in_progress -gt 0 ]; then
        echo ""
        echo "🔬 IN-PROGRESS: $in_progress"
        grep -l "^status: in-progress" lab-notebook/sessions/**/*.md 2>/dev/null | while read file; do
            echo "     • $(basename "$file" .md)"
        done
    fi

    # Check for INDEX.md
    if [ ! -f "lab-notebook/INDEX.md" ]; then
        echo ""
        echo "⚠️  INDEX.md not found - should be created"
    fi
else
    echo "   ⚠️  Lab notebook directory not found"
fi

echo ""
