# Overleaf Compilation Fixes - Complete

**Date**: 2025-01-23  
**Status**: ✅ All fixes applied

---

## ✅ Fixes Applied

### 1. Unicode Subscripts in Tables ✅
**Problem**: LaTeX doesn't support Unicode subscripts (₂, ₃, etc.) directly

**Fix**: Converted all Unicode subscripts to LaTeX math mode:
- `CH₂O` → `CH$_2$O`
- `C₈H₁₂N₂O₃` → `C$_8$H$_{12}$N$_2$O$_3$`
- etc.

**Files Fixed**:
- `tables/table5_hub_molecules.tex` ✅
- `tables/table6_novel_molecules.tex` ✅

### 2. Unicode Sigma (Σ) ✅
**Problem**: Unicode character Σ in Shannon entropy formula

**Fix**: Changed to LaTeX math mode:
- `H = -Σp_i log(p_i)` → `H = $-\sum p_i \log(p_i)$`

**File Fixed**: `manuscript_draft.tex` line 426 ✅

### 3. Unicode Plus-Minus (±) ✅
**Problem**: Unicode ± characters throughout text

**Fix**: Converted all to LaTeX math mode:
- `56.2 ± 8.6` → `56.2 $\pm$ 8.6`

**Locations Fixed**:
- Abstract ✅
- Results Section 3.1 ✅
- Results Section 3.3 ✅
- Discussion Section 4.2 ✅
- Discussion Section 4.3 ✅

### 4. Bibliography Style ✅
**Problem**: natbib error - bibliography not compatible

**Fix**: 
- Changed `\bibliographystyle{naturemag}` → `\bibliographystyle{plainnat}`
- Added `[numbers]` option: `\usepackage[numbers]{natbib}`

**File Fixed**: `manuscript_draft.tex` ✅

### 5. BibTeX Entry Error ✅
**Problem**: `hu2019taichi` has both `volume` and `number` fields

**Fix**: Removed `number` field (kept `volume`)

**File Fixed**: `references.bib` ✅

### 6. Missing Packages ✅
**Problem**: Missing packages for math symbols and table formatting

**Fix**: Added packages:
```latex
\usepackage{amsfonts}  % For \mathbb{R}
\usepackage{amssymb}   % For math symbols
\usepackage{booktabs}  % For \toprule, \midrule, \bottomrule
```

**File Fixed**: `manuscript_draft.tex` ✅

---

## 📋 Files Modified

1. ✅ `manuscript_draft.tex`
   - Added packages (amsfonts, amssymb, booktabs)
   - Fixed natbib configuration
   - Fixed Unicode ± characters (5 locations)
   - Fixed Unicode Σ character

2. ✅ `tables/table5_hub_molecules.tex`
   - Converted all Unicode subscripts to LaTeX math mode

3. ✅ `tables/table6_novel_molecules.tex`
   - Converted all Unicode subscripts to LaTeX math mode

4. ✅ `references.bib`
   - Fixed hu2019taichi entry (removed number field)

---

## 🔄 Next Steps in Overleaf

### Step 1: Upload Updated Files
Make sure these updated files are in Overleaf:
- `manuscript_draft.tex` ✅
- `tables/table5_hub_molecules.tex` ✅
- `tables/table6_novel_molecules.tex` ✅
- `references.bib` ✅

### Step 2: Set Main File
- Menu → Settings → Main document → `manuscript_draft.tex`

### Step 3: Clean and Recompile
1. Click **Menu** → **Clear Cache**
2. Click **Recompile**
3. If errors persist, try:
   - Menu → **Compiler** → Select **pdfLaTeX**
   - Click **Recompile** 3-4 times (to ensure all passes)

### Step 4: Verify
- [ ] No natbib errors
- [ ] No "undefined control sequence" errors
- [ ] No Unicode character errors
- [ ] Tables visible in PDF (at end)
- [ ] All citations show as numbers [1], [2], etc.

---

## 🐛 If Errors Persist

### Error: "Undefined control sequence \mathbb{R}"
**Solution**: Make sure `\usepackage{amsfonts}` is in the preamble (it is now)

### Error: "Is \usepackage{booktabs} missing?"
**Solution**: Make sure `\usepackage{booktabs}` is in the preamble (it is now)

### Error: Tables still not visible
**Solution**: 
1. Scroll to very end of PDF
2. Tables are after Figures section
3. Check compilation log for table file errors

### Error: natbib still complaining
**Solution**:
1. Delete `.aux` file in Overleaf
2. Recompile from scratch
3. Make sure `plainnat` style is used

---

## ✅ Expected Result

After fixes, compilation should show:
- ✅ No critical errors
- ✅ Only minor warnings (overfull hbox, etc. - these are OK)
- ✅ Tables visible at end of PDF
- ✅ All citations as numbers
- ✅ All cross-references resolved

---

**Status**: ✅ All fixes applied  
**Action**: Upload updated files to Overleaf and recompile

