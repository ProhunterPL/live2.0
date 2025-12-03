# Overleaf Upload Instructions - Fix Unicode Errors

**Date**: 2025-01-23  
**Problem**: Overleaf still sees Unicode characters in tables  
**Solution**: Files are fixed locally, need to upload to Overleaf

---

## ⚠️ CRITICAL: Upload Updated Files

The files on your local machine are **already fixed**, but Overleaf is using **old versions**.

### Files That MUST Be Uploaded:

1. **`tables/table5_hub_molecules.tex`** ✅ Fixed locally
2. **`tables/table6_novel_molecules.tex`** ✅ Fixed locally
3. **`manuscript_draft.tex`** ✅ Fixed locally
4. **`references.bib`** ✅ Fixed locally

---

## 📤 How to Upload to Overleaf

### Method 1: Replace Files (Recommended)

1. **In Overleaf**, go to the file list
2. **Delete** the old versions:
   - Right-click `tables/table5_hub_molecules.tex` → Delete
   - Right-click `tables/table6_novel_molecules.tex` → Delete
3. **Upload** new versions:
   - Click **Upload** button (top of file list)
   - Select `table5_hub_molecules.tex` from your local `paper/tables/` folder
   - Select `table6_novel_molecules.tex` from your local `paper/tables/` folder
4. **Verify** file contents:
   - Open `table5_hub_molecules.tex` in Overleaf
   - Check line 9: Should show `CH$_2$O` (NOT `CH₂O`)
   - If you see `CH₂O` with Unicode, the wrong file was uploaded

### Method 2: Copy-Paste Content

1. **Open** local file: `paper/tables/table5_hub_molecules.tex`
2. **Select All** (Ctrl+A) and **Copy** (Ctrl+C)
3. **In Overleaf**, open `tables/table5_hub_molecules.tex`
4. **Select All** and **Paste** (replace entire content)
5. **Save** (Ctrl+S)
6. **Repeat** for `table6_novel_molecules.tex`

---

## ✅ Verification Checklist

After uploading, verify in Overleaf:

### Check table5_hub_molecules.tex:
- [ ] Line 9: `CH$_2$O` (with `$_2$`, NOT Unicode ₂)
- [ ] Line 11: `NH$_3$` (with `$_3$`, NOT Unicode ₃)
- [ ] Line 12: `H$_2$CO$_3$` (with `$_2$` and `$_3$`, NOT Unicode)

### Check table6_novel_molecules.tex:
- [ ] Line 9: `C$_8$H$_{12}$N$_2$O$_3$` (with `$_8$`, `$_{12}$`, etc., NOT Unicode)
- [ ] No Unicode subscripts visible anywhere

### Quick Test:
- In Overleaf, open the table file
- Search for `₂` (Unicode subscript 2)
- If found → Wrong file uploaded
- If not found → Correct file ✅

---

## 🔄 After Uploading

1. **Clear Cache** in Overleaf:
   - Menu → **Clear Cache**
   
2. **Recompile**:
   - Click **Recompile** button
   - Wait for compilation

3. **Check Logs**:
   - Should see **NO** Unicode character errors
   - Should see **NO** "undefined control sequence" errors for tables

---

## 📋 What Should Be in Files

### table5_hub_molecules.tex (Line 9):
```latex
CH$_2$O & Formaldehyde & 28 & 0.420000 & All & Central building block \\
```

**NOT**:
```latex
CH₂O & Formaldehyde & 28 & 0.420000 & All & Central building block \\
```

### table6_novel_molecules.tex (Line 9):
```latex
1 & C$_8$H$_{12}$N$_2$O$_3$ & 184 & 7.800000 & Formamide & 342000 \\
```

**NOT**:
```latex
1 & C₈H₁₂N₂O₃ & 184 & 7.800000 & Formamide & 342000 \\
```

---

## 🐛 If Still Seeing Errors

### Problem: Still seeing Unicode errors after upload

**Solution 1**: Check file encoding
- Make sure files are saved as **UTF-8** (should be default)
- In Overleaf, check if file shows correct content

**Solution 2**: Force refresh
- Delete file in Overleaf
- Upload again from local machine
- Clear cache
- Recompile

**Solution 3**: Manual edit in Overleaf
- Open table file in Overleaf
- Find any Unicode subscript (₂, ₃, etc.)
- Replace manually:
  - `₂` → `$_2$`
  - `₃` → `$_3$`
  - `₄` → `$_4$`
  - etc.

---

## ✅ Expected Result

After correct upload and recompilation:
- ✅ No Unicode character errors
- ✅ Tables compile successfully
- ✅ Tables visible in PDF
- ✅ All subscripts render correctly

---

**Status**: Files fixed locally ✅  
**Action**: Upload fixed files to Overleaf

