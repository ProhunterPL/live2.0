# Overleaf Compilation Instructions

**Date**: 2025-01-23  
**Status**: ✅ All fixes applied to `manuscript_draft.tex`

---

## ⚠️ IMPORTANT: Set Main File

**CRITICAL**: Overleaf may be compiling the wrong file!

### Fix This First:
1. In Overleaf, click **Menu** (☰ icon, top left)
2. Go to **Settings**
3. Under **Main document**, select **`manuscript_draft.tex`**
4. Click **Save**

**DO NOT** use `manuscript_draft_filled.tex` - that's an old file with errors!

---

## ✅ Fixes Already Applied

I've fixed all issues in `manuscript_draft.tex`:

### 1. Added Missing Packages
```latex
\usepackage{amsfonts}  % For \mathbb{R}
\usepackage{amssymb}   % For math symbols
\usepackage{booktabs}  % For table formatting (\toprule, etc.)
```

### 2. Fixed Bibliography Style
```latex
% Changed from:
\usepackage{natbib}
\bibliographystyle{naturemag}

% To:
\usepackage[numbers]{natbib}  % Numerical citations
\bibliographystyle{plainnat}  % Compatible style
```

---

## 🔄 Compilation Steps in Overleaf

1. **Set Main File** (see above) ⚠️ CRITICAL
2. Click **Recompile** button
3. Wait for compilation to complete
4. Check for errors in log

### Expected Result:
- ✅ No natbib errors
- ✅ No "undefined control sequence" errors
- ✅ Tables visible in PDF (at end, after Figures)
- ✅ All citations as numbers [1], [2], etc.

---

## 📊 Where to Find Tables

Tables are located at the **end of the document**:
1. Main text (Introduction → Methods → Results → Discussion → Conclusions)
2. References section
3. **Figures section** (6 figures)
4. **Tables section** ← Tables 5 and 6 are here
5. Supplementary Information section

**Scroll to the very end of the PDF** to see tables!

---

## 🐛 If Tables Still Don't Appear

### Check 1: File Paths
In Overleaf, verify these files exist:
- `tables/table5_hub_molecules.tex`
- `tables/table6_novel_molecules.tex`

### Check 2: Compilation Log
Look for errors like:
- `File 'table5_hub_molecules.tex' not found`
- `Undefined control sequence \toprule`

If you see these, the fixes may not have been uploaded yet.

### Check 3: Force Full Recompile
1. Click **Menu** → **Compiler** → Select **pdfLaTeX**
2. Click **Recompile** 3-4 times
3. This ensures all passes complete (pdflatex → bibtex → pdflatex → pdflatex)

---

## 📝 Quick Checklist

- [ ] Main file set to `manuscript_draft.tex` (NOT `manuscript_draft_filled.tex`)
- [ ] All files uploaded to Overleaf
- [ ] Recompiled after setting main file
- [ ] No errors in compilation log
- [ ] Scrolled to end of PDF to see tables
- [ ] All citations show as numbers

---

## 💡 Tip: Verify Main File

To check which file is main:
- Look at the file list in Overleaf
- The main file has a **star icon** ⭐ next to it
- If `manuscript_draft_filled.tex` has the star, click on `manuscript_draft.tex` and select "Set as Main Document"

---

**Status**: ✅ Ready to compile  
**Action**: Set main file and recompile in Overleaf

