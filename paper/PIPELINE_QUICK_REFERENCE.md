# 🚀 Analysis Pipeline - Quick Reference

**Status**: ✅ Complete  
**Ready for**: AWS data

---

## 📋 ONE COMMAND Solution

### **Master Script** (Recommended)

```bash
# Process everything at once!
python scripts/process_phase2b_for_paper.py \
    --input results/phase2b_additional

# Expected time: 10-30 minutes
# Output: All data, figures, and tables!
```

**This runs**:
1. ✅ Autocatalysis detection
2. ✅ Complexity metrics
3. ✅ Figure generation (4 figures)
4. ✅ Table generation (3 tables)

---

## 🔧 Individual Scripts (Advanced)

### **1. Analysis**

```bash
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data

# Time: ~10-20 min
# Output: JSON files with all statistics
```

**Generates**:
- `miller_urey_extended_analysis.json`
- `hydrothermal_extended_analysis.json`
- `formamide_extended_analysis.json`
- `scenario_comparison.json`
- `summary_table.csv`

---

### **2. Figures**

```bash
python scripts/generate_all_figures.py \
    --data paper/results_data \
    --output paper/figures

# Time: ~2-5 min
# Output: 4 publication-ready figures (300 DPI)
```

**Generates**:
- `figure3_molecular_diversity.png` (4-panel)
- `figure4_reaction_networks.png` (4-panel)
- `figure5_autocatalytic_cycles.png` (4-panel)
- `figure6_novel_molecules.png` (4-panel)

---

### **3. Tables**

```bash
python scripts/generate_all_tables.py \
    --data paper/results_data \
    --output paper/tables

# Time: ~1-2 min
# Output: Tables in CSV, LaTeX, Markdown
```

**Generates**:
- `table5_hub_molecules` (.csv, .tex, .md)
- `table6_novel_molecules` (.csv, .tex, .md)
- `tableS2_network_metrics` (.csv, .tex, summary)

---

## 📊 Output Structure

```
paper/
├── results_data/
│   ├── miller_urey_extended_analysis.json
│   ├── hydrothermal_extended_analysis.json
│   ├── formamide_extended_analysis.json
│   ├── scenario_comparison.json
│   ├── summary_table.csv
│   ├── figure_data.json
│   └── latex_snippets.txt
│
├── figures/
│   ├── figure3_molecular_diversity.png
│   ├── figure4_reaction_networks.png
│   ├── figure5_autocatalytic_cycles.png
│   └── figure6_novel_molecules.png
│
└── tables/
    ├── table5_hub_molecules.csv
    ├── table5_hub_molecules.tex
    ├── table5_hub_molecules.md
    ├── table6_novel_molecules.csv
    ├── table6_novel_molecules.tex
    ├── table6_novel_molecules.md
    ├── tableS2_network_metrics.csv
    ├── tableS2_network_metrics.tex
    └── tableS2_summary.csv
```

---

## ⚡ Quick Commands

### **Test Pipeline** (before AWS completes):

```bash
# Test with dummy data
python backend/sim/analysis/autocatalysis_detector.py
python backend/sim/core/complexity_metrics.py

# Generate mockup figures
python scripts/generate_all_figures.py \
    --data paper/results_data \
    --output paper/figures_mockup
```

### **Skip Steps** (if already done):

```bash
# Skip analysis (if already run)
python scripts/process_phase2b_for_paper.py \
    --input results/phase2b_additional \
    --skip-analysis

# Skip figures
python scripts/process_phase2b_for_paper.py \
    --input results/phase2b_additional \
    --skip-figures

# Skip tables
python scripts/process_phase2b_for_paper.py \
    --input results/phase2b_additional \
    --skip-tables
```

---

## 📝 What Data You Get

### **Autocatalysis Statistics**:
```json
{
  "total_cycles": 127,
  "cycles_per_run_mean": 4.2,
  "cycle_types": {
    "direct": 18,
    "indirect": 95,
    "hypercycle": 14
  },
  "amplification_stats": {
    "mean": 3.7,
    "median": 3.2,
    "max": 12.8
  }
}
```

### **Complexity Metrics**:
```json
{
  "shannon_entropy": 3.42,
  "species_richness": 142,
  "self_organization_index": 0.67,
  "phases": {
    "exploration": {"growth_rate": 0.000012},
    "diversification": {"growth_rate": 0.000034},
    "consolidation": {"growth_rate": 0.000008}
  }
}
```

---

## 🔍 Filling Paper Sections

### **Results 3.3 - Autocatalytic Cycles**:

1. Open: `paper/results_data/latex_snippets.txt`
2. Copy autocatalysis snippet
3. Paste into `manuscript_draft.tex` Section 3.3
4. Replace [XX] with real numbers from JSON

### **Discussion 4.1 - Emergent Complexity**:

1. Open: `paper/results_data/miller_urey_extended_analysis.json`
2. Get `complexity.phases` data
3. Fill three-phase pattern description
4. Add quantitative growth rates

### **Discussion 4.3 - Autocatalysis**:

1. Open: `paper/results_data/scenario_comparison.json`
2. Get `autocatalysis_comparison` data
3. Fill cycle frequencies
4. Add amplification factors

---

## ⏱️ Timeline

| Step | Time | Output |
|------|------|--------|
| Download AWS data | 5 min | Local directory |
| Run master script | 20 min | All outputs |
| Review outputs | 30 min | Verify quality |
| Fill Results 3.3 | 30 min | Section complete |
| Fill Discussion 4.1 | 30 min | Section complete |
| Fill Discussion 4.3 | 30 min | Section complete |
| **TOTAL** | **~3 hours** | **Paper sections complete!** |

---

## 🚨 Troubleshooting

### **"Module not found"**:
```bash
# Activate environment
.\activate_live_env.ps1  # Windows
# OR
source venv/bin/activate  # Linux/Mac
```

### **"No data files"**:
```bash
# Check AWS results downloaded
ls results/phase2b_additional/

# Should see:
# - miller_urey_extended/
# - hydrothermal_extended/
# - formamide_extended/
```

### **"Figure generation failed"**:
```bash
# Install matplotlib-venn if missing
pip install matplotlib-venn

# Rerun just figures
python scripts/generate_all_figures.py \
    --data paper/results_data \
    --output paper/figures
```

---

## ✅ Success Checklist

After running pipeline, verify:

- [ ] `paper/results_data/` has 5+ JSON files
- [ ] `paper/figures/` has 4 PNG files (300 DPI)
- [ ] `paper/tables/` has 6+ files (CSV + LaTeX)
- [ ] `latex_snippets.txt` has copy-pasteable text
- [ ] Figures look good (open and inspect)
- [ ] Tables have reasonable numbers

If all ✅ → Ready to fill paper sections!

---

## 🎯 Final Step: Filling Paper

```bash
# 1. Open manuscript
code paper/manuscript_draft.tex

# 2. Find Section 3.3
# 3. Copy from latex_snippets.txt
# 4. Replace [XX] with real numbers

# 5. Find Discussion 4.1
# 6. Add phase data from complexity metrics

# 7. Find Discussion 4.3
# 8. Add cycle data from autocatalysis analysis

# 9. Save and compile LaTeX
pdflatex paper/manuscript_draft.tex

# 10. Review PDF
# 11. Celebrate! 🎉
```

---

**Status**: ✅ Pipeline ready  
**Next**: Wait for AWS → Run master script → Fill paper  
**ETA**: 3 hours from AWS completion to paper sections complete

**Questions?** See `TIER1_IMPLEMENTATION_GUIDE.md` for details.

