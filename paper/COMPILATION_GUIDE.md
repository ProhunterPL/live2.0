# LaTeX Compilation Guide

**Date**: 2025-01-23  
**Manuscript**: `manuscript_draft.tex`

---

## 📋 Prerequisites

### Option 1: Local LaTeX Installation

**Windows:**
- Install [MiKTeX](https://miktex.org/download) (recommended)
- Or [TeX Live](https://www.tug.org/texlive/windows.html)

**macOS:**
```bash
brew install --cask mactex
```

**Linux:**
```bash
sudo apt-get install texlive-full  # Ubuntu/Debian
```

### Option 2: Online Compilation

Use online LaTeX compilers:
- [Overleaf](https://www.overleaf.com) (recommended - free, collaborative)
- [ShareLaTeX](https://www.sharelatex.com)
- [LaTeX Base](https://latexbase.com)

---

## 🔧 Compilation Steps

### Standard Compilation (3-pass)

```bash
cd paper

# Pass 1: Initial compilation
pdflatex manuscript_draft.tex

# Pass 2: Process bibliography
bibtex manuscript_draft

# Pass 3: Final compilation (2x for cross-references)
pdflatex manuscript_draft.tex
pdflatex manuscript_draft.tex
```

### Using latexmk (Automatic)

```bash
cd paper
latexmk -pdf manuscript_draft.tex
```

This automatically runs all necessary passes.

---

## ✅ Verification Checklist

After compilation, verify:

1. **No LaTeX Errors**
   - Check console output for errors
   - Look for "Error" messages (warnings are usually OK)

2. **Cross-References Resolved**
   - Open PDF and check:
     - All `\ref{fig:...}` show figure numbers (not "??")
     - All `\ref{tab:...}` show table numbers (not "??")
     - All `\ref{sec:...}` show section numbers (not "??")

3. **Figures Present**
   - All 6 figures should appear in PDF
   - Check figure quality (300 DPI)

4. **Tables Present**
   - Table 5 (Hub Molecules) should appear
   - Table 6 (Novel Molecules) should appear

5. **Bibliography**
   - All citations should show as numbers
   - References section should list all cited papers

---

## 🐛 Common Issues & Solutions

### Issue 1: Missing Figures

**Error**: `! LaTeX Error: File 'figureX.png' not found.`

**Solution**:
- Verify all figure files exist in `paper/figures/`
- Check paths in `\includegraphics{}` commands
- Use relative paths: `figures/figure1_thermodynamic_validation.png`

### Issue 2: Missing Bibliography

**Error**: Citations show as `[?]` or `[author?]`

**Solution**:
- Run `bibtex manuscript_draft` after first pdflatex pass
- Verify `references.bib` exists and contains all citations
- Run pdflatex twice more after bibtex

### Issue 3: Cross-References Not Resolved

**Error**: References show as `??`

**Solution**:
- Run pdflatex multiple times (usually 2-3 passes needed)
- Check that all `\label{}` commands are present
- Verify label names match reference names exactly

### Issue 4: Missing Packages

**Error**: `! LaTeX Error: File 'package.sty' not found.`

**Solution**:
- Install missing packages via package manager
- For MiKTeX: Package Manager will auto-install on first use
- For TeX Live: `tlmgr install <package>`

---

## 📦 Required LaTeX Packages

The manuscript uses these packages (should be auto-installed):

```latex
\usepackage[utf8]{inputenc}    % UTF-8 encoding
\usepackage{amsmath}            % Math environments
\usepackage{graphicx}           % Graphics inclusion
\usepackage{hyperref}           % Hyperlinks
\usepackage{natbib}             % Bibliography
\usepackage{lineno}             % Line numbers (draft mode)
\usepackage[margin=1in]{geometry} % Page margins
```

Plus table packages (booktabs, etc.) - should be in standard distributions.

---

## 🎯 Quick Compilation Test

To quickly test if compilation works:

```bash
cd paper

# Test compilation (will show errors if any)
pdflatex -interaction=nonstopmode manuscript_draft.tex 2>&1 | grep -i error

# If no errors, continue with full compilation
bibtex manuscript_draft
pdflatex manuscript_draft.tex
pdflatex manuscript_draft.tex
```

---

## 📄 Output Files

After successful compilation, you should have:

- `manuscript_draft.pdf` - Final PDF (for submission)
- `manuscript_draft.aux` - Auxiliary file (can delete)
- `manuscript_draft.log` - Compilation log (can delete)
- `manuscript_draft.bbl` - Bibliography (can delete)
- `manuscript_draft.blg` - BibTeX log (can delete)
- `manuscript_draft.out` - Hyperref output (can delete)

**Keep only**: `manuscript_draft.pdf` and source files (`.tex`, `.bib`)

---

## 🌐 Online Compilation (Overleaf)

### Steps:

1. **Create Account**: Sign up at https://www.overleaf.com (free)

2. **Create New Project**:
   - Click "New Project" → "Upload Project"
   - Upload entire `paper/` directory

3. **Compile**:
   - Click "Recompile" button
   - Overleaf will automatically run all passes

4. **Download PDF**:
   - Click "Download" → "PDF"

### Advantages:
- ✅ No local installation needed
- ✅ Automatic compilation
- ✅ Real-time error detection
- ✅ Collaborative editing
- ✅ Version control

---

## ✅ Pre-Submission PDF Check

Before submitting, verify PDF contains:

- [ ] Title page with author information
- [ ] Abstract (250 words)
- [ ] All 5 main sections
- [ ] All 6 figures (numbered correctly)
- [ ] All 2 tables (numbered correctly)
- [ ] All citations as numbers
- [ ] References section at end
- [ ] Page numbers
- [ ] Line numbers (draft mode)

---

## 📝 Notes

- **Draft Mode**: Manuscript uses `\linenumbers` for review
- **Final Version**: Remove `\linenumbers` before final submission if journal doesn't require it
- **File Size**: PDF should be ~2-5 MB (with figures)
- **Quality**: All figures are 300 DPI (publication quality)

---

**Last Updated**: 2025-01-23  
**Status**: Ready for compilation

