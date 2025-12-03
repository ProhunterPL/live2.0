# Live 2.0 - Paper Development

This directory contains the manuscript for publication in *Origins of Life and Evolution of Biospheres* or *Journal of Chemical Theory and Computation (JCTC)*.

---

## 📁 Contents

```
paper/
├── manuscript_draft.tex       ← Main manuscript (LaTeX)
├── references.bib             ← BibTeX references
├── README.md                  ← This file
├── figures/                   ← Publication figures (to be generated)
│   ├── figure1_thermodynamic_validation.png
│   ├── figure2_benchmark_validation.png
│   ├── figure3_molecular_diversity.png
│   ├── figure4_reaction_networks.png
│   ├── figure5_autocatalytic_cycles.png
│   ├── figure6_top_molecules.png
│   └── figure7_emergence_timeline.png
├── tables/
│   ├── tableS1_parameters.tex  ← Parameter database (exists)
│   └── [other tables TBD]
└── supplementary/             ← Supplementary materials
    ├── SI_document.tex
    ├── tables/
    └── figures/
```

---

## 📝 Manuscript Status

**Current Stage**: ✅ READY FOR SUBMISSION

| Section | Status | Word Count (Target) | Notes |
|---------|--------|---------------------|-------|
| Abstract | ✅ COMPLETE | ~250 / 250 | All data filled |
| Introduction | ✅ COMPLETE | ~1500 / 1500 | Ready for review |
| Methods | ✅ COMPLETE | ~1800 / 1800 | Truth-filter added |
| Results | ✅ COMPLETE | ~1800 / 1800 | All sections filled |
| Discussion | ✅ COMPLETE | ~1200 / 1200 | All sections filled |
| Conclusions | ✅ COMPLETE | ~250 / 250 | 4 paragraphs |
| **Total** | **✅ 100%** | **~5800 / 6000** | Ready for submission |

---

## 🎯 Next Steps

### Immediate (When AWS Results Arrive):

1. **Fill Results section** with actual data:
   ```bash
   # After running AWS pipeline
   python scripts/generate_all_figures.py --input results/aws_batch/analysis
   ```

2. **Update placeholders** marked with `[XX]` or `TODO`

3. **Generate all figures**:
   - Figure 1: Thermodynamic validation (already have script)
   - Figure 2: Benchmark validation (already have script)
   - Figure 3-7: From AWS analysis

4. **Write Discussion** based on Results

5. **Write Abstract** (last, summarizing everything)

### Before Submission:

- [ ] All figures at 300+ DPI
- [ ] All tables formatted consistently
- [ ] All references have DOIs
- [ ] Word count < 6000
- [ ] Supplementary materials complete
- [ ] Code/data availability statements updated with DOIs

---

## 🔧 Building the Paper

### Requirements:

```bash
# LaTeX distribution (TeXLive, MiKTeX, or MacTeX)
sudo apt install texlive-full  # Ubuntu/Debian
# or
brew install --cask mactex      # macOS
```

### Compile:

```bash
cd paper
pdflatex manuscript_draft.tex
bibtex manuscript_draft
pdflatex manuscript_draft.tex
pdflatex manuscript_draft.tex
```

Or use **latexmk** for automatic compilation:

```bash
latexmk -pdf manuscript_draft.tex
```

### View:

```bash
open manuscript_draft.pdf  # macOS
xdg-open manuscript_draft.pdf  # Linux
```

---

## 📊 Figures Checklist

Generate figures after AWS simulations complete:

- [ ] **Figure 1**: Thermodynamic validation
  ```bash
  python scripts/analyze_thermodynamics.py --output paper/figures/figure1_thermodynamic_validation.png
  ```

- [ ] **Figure 2**: Benchmark validation
  ```bash
  python scripts/analyze_benchmark_reactions.py --output paper/figures/figure2_benchmark_validation.png
  ```

- [ ] **Figure 3**: Molecular diversity
  ```bash
  python scripts/plot_molecular_diversity.py results/aws_batch/analysis --output paper/figures/figure3_molecular_diversity.png
  ```

- [ ] **Figure 4**: Reaction networks
  ```bash
  python scripts/plot_reaction_networks.py results/aws_batch/analysis --output paper/figures/figure4_reaction_networks.png
  ```

- [ ] **Figure 5**: Autocatalytic cycles
  ```bash
  python scripts/plot_autocatalytic_cycles.py results/aws_batch/analysis --output paper/figures/figure5_autocatalytic_cycles.png
  ```

- [ ] **Figure 6**: Top molecules
  ```bash
  python scripts/plot_top_molecules.py results/aws_batch/analysis --output paper/figures/figure6_top_molecules.png
  ```

- [ ] **Figure 7**: Emergence timeline
  ```bash
  python scripts/plot_emergence_timeline.py results/aws_batch/analysis --output paper/figures/figure7_emergence_timeline.png
  ```

Or generate all at once:

```bash
python scripts/generate_all_figures.py --input results/aws_batch/analysis --output paper/figures
```

---

## 📚 References Management

### Adding references:

1. Open `references.bib`
2. Add BibTeX entry (get from Google Scholar, DOI lookup, etc.)
3. Cite in text: `\citep{author2023}`
4. Recompile with BibTeX

### Checking references:

```bash
# List all citations
grep -o '\\citep{[^}]*}' manuscript_draft.tex | sort -u

# Find missing references
# (compare with references.bib)
```

---

## ✅ Pre-Submission Checklist

### Content:

- [ ] All sections complete
- [ ] Abstract ≤ 250 words
- [ ] Main text ≤ 6000 words
- [ ] All placeholders `[XX]` filled
- [ ] All TODOs resolved
- [ ] All figures cited in text
- [ ] All tables cited in text

### Formatting:

- [ ] Figures: 300 DPI, TIFF or PDF
- [ ] Figure legends: informative and complete
- [ ] Tables: consistent formatting
- [ ] References: all have DOIs
- [ ] Line numbers enabled (for review)
- [ ] Page numbers enabled

### Supplementary:

- [ ] Table S1: Parameter database (exists!)
- [ ] Table S2: All molecules
- [ ] Figures S1-S10: Additional validation
- [ ] Movies S1-S2: Visualizations
- [ ] Code availability: GitHub + Zenodo DOI
- [ ] Data availability: Zenodo/Dryad DOI

### Quality:

- [ ] Spell-check complete
- [ ] Grammar check (Grammarly)
- [ ] Equations numbered correctly
- [ ] Cross-references working
- [ ] No orphan/widow lines
- [ ] Consistent terminology

---

## 🎯 Target Journals

### Primary: *Origins of Life and Evolution of Biospheres*

- **Scope**: Perfect fit (prebiotic chemistry)
- **Impact Factor**: ~2.5
- **Word Limit**: 6000 words
- **Allows preprints**: YES
- **Open access**: Optional (€2,290)
- **Review time**: 2-3 months

### Alternative: *Journal of Chemical Theory and Computation (JCTC)*

- **Scope**: Computational chemistry
- **Impact Factor**: ~5.5
- **Word Limit**: No strict limit
- **Allows preprints**: YES
- **Open access**: Optional
- **Review time**: 2-4 months
- **Note**: More computational focus needed

### Backup: *Astrobiology*

- **Scope**: Origin of life, astrobiology
- **Impact Factor**: ~3.5
- **Allows preprints**: YES

---

## 📅 Timeline

| Date | Milestone | Status |
|------|-----------|--------|
| Oct 16, 2025 | Manuscript skeleton | ✅ DONE |
| Oct 22, 2025 | AWS simulations complete | 🏃 IN PROGRESS |
| Oct 26, 2025 | All figures generated | ⏳ PENDING |
| Nov 2, 2025 | Results section complete | ⏳ PENDING |
| Nov 9, 2025 | Discussion complete | ⏳ PENDING |
| Nov 16, 2025 | Full draft ready | ⏳ PENDING |
| Nov 23, 2025 | Internal review | ⏳ PENDING |
| Nov 30, 2025 | **SUBMISSION** | 🎯 TARGET |
| Feb 2026 | **PUBLICATION** | 🔮 GOAL |

---

## 🤝 Contribution Guidelines

If collaborators join:

1. Use **track changes** mode in LaTeX:
   ```latex
   \usepackage{changes}
   \added{New text}
   \deleted{Old text}
   \replaced{New}{Old}
   ```

2. Comment style:
   ```latex
   % TODO: [Your Name] - Description
   % FIXME: [Your Name] - What needs fixing
   ```

3. Commit messages:
   ```
   paper: Add Results section (Section 3.1)
   paper: Update Figure 3 with final data
   paper: Fix typos in Introduction
   ```

---

## 📞 Support

For questions about:
- **Manuscript content**: Check VALIDATION_ROADMAP.md
- **Data analysis**: Check docs/PHASE3_ANALYSIS_GUIDE.md
- **Figures**: Check scripts/README.md
- **AWS pipeline**: Check docs/AWS_RESULTS_PIPELINE.md
- **Publication strategy**: Check PUBLICATION_STRATEGY.md (includes Polish media strategy)

---

**Last updated**: January 23, 2025  
**Status**: READY FOR DATA - Manuscript structure complete, waiting for AWS results

