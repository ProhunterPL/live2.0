#!/bin/bash
# Quick script to check Phase 2B results on AWS

echo "Phase 2B Results Check"
echo "====================="
echo ""

# Check main reports
echo "📄 Main Reports:"
echo "----------------"
ls -lh results/phase2b_additional/*.md 2>/dev/null || echo "  No reports found"

echo ""
echo "📊 Formamide Debug:"
echo "-------------------"
ls -lh results/phase2b_additional/formamide_debug/formamide_debug_report.md 2>/dev/null || echo "  No debug report"

echo ""
echo "🧪 Simulation Results:"
echo "---------------------"
echo "Miller-Urey:"
ls -ld results/phase2b_additional/miller_urey_extended/ 2>/dev/null | wc -l
ls results/phase2b_additional/miller_urey_extended/ 2>/dev/null | wc -l

echo "Hydrothermal:"
ls -ld results/phase2b_additional/hydrothermal_extended/ 2>/dev/null | wc -l
ls results/phase2b_additional/hydrothermal_extended/ 2>/dev/null | wc -l

echo "Formamide:"
ls -ld results/phase2b_additional/formamide_extended/ 2>/dev/null | wc -l
ls results/phase2b_additional/formamide_extended/ 2>/dev/null | wc -l

echo ""
echo "📈 Quick Summary:"
echo "------------------"
if [ -f results/phase2b_additional/phase2b_analysis_report.md ]; then
    echo "First 50 lines of analysis:"
    head -50 results/phase2b_additional/phase2b_analysis_report.md
fi

