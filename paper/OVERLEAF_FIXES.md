# Overleaf Compilation Fixes

**Date**: 2025-01-23  
**Status**: ✅ Fixes Applied

---

## 🔧 Problems Fixed

### 1. Missing Packages ✅
**Problem**: Tables use `\toprule`, `\midrule`, `\bottomrule` but `booktabs` package not loaded

**Fix Applied**: Added `\usepackage{booktabs}` to manuscript

### 2. Missing Math Fonts ✅
**Problem**: `\mathbb{R}` requires `amsfonts` package

**Fix Applied**: Added `\usepackage{amsfonts}` and `\usepackage{amssymb}`

### 3. Bibliography Style Mismatch ✅
**Problem**: Using `naturemag` style (author-year) with `\citep{}` (numerical citations)

**Fix Applied**: 
- Changed to `plainnat` style (numerical citations)
- Added `[numbers]` option to natbib: `\usepackage[numbers]{natbib}`

### 4. Wrong File Being Compiled ⚠️
**Problem**: Overleaf may be compiling `manuscript_draft_filled.tex` instead of `manuscript_draft.tex`

**Action Required**: 
- In Overleaf, make sure main file is set to `manuscript_draft.tex`
- Go to: Menu → Settings → Main document → Select `manuscript_draft.tex`

---

## 📋 Steps to Fix in Overleaf

### Step 1: Set Correct Main File
1. In Overleaf, click **Menu** (top left)
2. Go to **Settings**
3. Under **Main document**, select `manuscript_draft.tex`
4. Click **Save**

### Step 2: Recompile
1. Click **Recompile** button
2. Check if errors are resolved

### Step 3: Verify Tables
1. Scroll to end of PDF
2. Look for "Tables" section
3. Tables 5 and 6 should be visible

---

## ✅ Changes Made to manuscript_draft.tex

### Added Packages:
```latex
\usepackage{amsfonts}  % For \mathbb{R}
\usepackage{amssymb}  % For additional math symbols
\usepackage{booktabs}  % For \toprule, \midrule, \bottomrule in tables
```

### Changed natbib:
```latex
% Before:
\usepackage{natbib}

% After:
\usepackage[numbers]{natbib}  % Use numerical citations
```

### Changed Bibliography Style:
```latex
% Before:
\bibliographystyle{naturemag}

% After:
\bibliographystyle{plainnat}  % Numerical citation style
```

---

## 🐛 If Tables Still Don't Appear

### Check 1: Table Files Exist
- Verify `tables/table5_hub_molecules.tex` exists
- Verify `tables/table6_novel_molecules.tex` exists

### Check 2: Table Location
Tables are in section "Tables" at the end of document (after References, after Figures).

### Check 3: Compilation Order
Make sure to run:
1. `pdflatex` (first pass)
2. `bibtex` (process bibliography)
3. `pdflatex` (second pass - resolves citations)
4. `pdflatex` (third pass - resolves cross-references)

Overleaf should do this automatically, but you can force it by clicking "Recompile" multiple times.

### Check 4: Unicode Characters
Tables contain Unicode subscripts (CH₂O, H₂CO₃, etc.). If these cause issues:
- Overleaf should handle UTF-8 correctly
- If not, may need to use LaTeX math mode: `CH$_2$O`

---

## 📝 Verification Checklist

After fixes, verify:
- [ ] No natbib errors
- [ ] No "undefined control sequence" errors
- [ ] No "missing $ inserted" errors
- [ ] Tables visible in PDF
- [ ] All citations show as numbers [1], [2], etc.
- [ ] All cross-references resolved (no "??")

---

**Status**: ✅ Fixes applied to `manuscript_draft.tex`  
**Next**: Set main file in Overleaf and recompile

