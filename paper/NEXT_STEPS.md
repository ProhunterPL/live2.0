# Next Steps - Manuscript Submission

**Date**: 2025-01-23  
**Status**: ✅ Manuscript ready, awaiting compilation and final checks

---

## ✅ Completed Steps

1. ✅ **Final Review** - Complete
2. ✅ **Data Consistency** - Fixed
3. ✅ **All Figures Generated** - 6/6 figures present
4. ✅ **All Tables Present** - 2/2 tables present
5. ✅ **Technical Validation** - 0 errors

---

## 📋 Remaining Steps

### Step 1: Compile PDF (5-10 minutes)

**Option A: Local Compilation** (if LaTeX installed)
```bash
cd paper
pdflatex manuscript_draft.tex
bibtex manuscript_draft
pdflatex manuscript_draft.tex
pdflatex manuscript_draft.tex
```

**Option B: Online Compilation** (recommended if no LaTeX)
1. Go to https://www.overleaf.com
2. Create free account
3. Upload `paper/` directory
4. Click "Recompile"
5. Download PDF

**See**: `COMPILATION_GUIDE.md` for detailed instructions

---

### Step 2: Verify PDF (5 minutes)

Check that PDF contains:
- [ ] Title page with correct author info
- [ ] Abstract (250 words)
- [ ] All 5 main sections
- [ ] All 6 figures (numbered 1-6)
- [ ] All 2 tables (numbered 5-6)
- [ ] All citations as numbers (not "??")
- [ ] References section at end
- [ ] Page numbers visible
- [ ] Line numbers (draft mode)

---

### Step 3: Final Grammar Check (30 minutes)

**Recommended Tools:**
- [Grammarly](https://www.grammarly.com) (free version sufficient)
- [LanguageTool](https://languagetool.org) (open source)
- Microsoft Word (if available - paste text)

**Check for:**
- [ ] Typos and spelling errors
- [ ] Grammar mistakes
- [ ] Consistency in terminology
- [ ] Number formatting (commas, decimals)
- [ ] Scientific notation consistency

---

### Step 4: Journal-Specific Formatting (10 minutes)

**For Origins of Life and Evolution of Biospheres:**
- [ ] Check word limit (6000-8000 words typical)
- [ ] Verify abstract ≤ 250 words ✅
- [ ] Check figure format (300 DPI) ✅
- [ ] Verify reference style (natbib) ✅
- [ ] Check if line numbers required (currently enabled)

**For JCTC (Alternative):**
- [ ] Check word limit (8000-10000 words typical)
- [ ] Prepare Supporting Information
- [ ] Verify code availability statement

---

### Step 5: Prepare Submission Package (10 minutes)

**Required Files:**
- [ ] `manuscript_draft.pdf` - Main manuscript
- [ ] All figure files (6 PNG files, 300 DPI)
- [ ] All table files (if journal requires separate)
- [ ] Cover letter (if required)
- [ ] Author information form (if required)

**Optional Files:**
- [ ] Supplementary materials
- [ ] Data availability statement
- [ ] Conflict of interest statement

---

### Step 6: Upload to Journal (15 minutes)

1. **Access Submission System**
   - Origins of Life: https://www.springer.com/journal/11084
   - JCTC: https://pubs.acs.org/journal/jctcce

2. **Fill Submission Form**
   - Author information
   - Corresponding author
   - Keywords
   - Suggested reviewers (if optional)
   - Cover letter (if required)

3. **Upload Files**
   - Main manuscript PDF
   - Figures (if separate upload required)
   - Tables (if separate upload required)
   - Supplementary materials (if any)

4. **Review and Submit**
   - Check all information
   - Verify file uploads
   - Submit!

---

## 📝 Post-Submission Tasks

### Immediate (After Submission)

1. **Save Submission Confirmation**
   - Download confirmation email
   - Note submission ID/reference number
   - Save submission date

2. **Update Documentation**
   - Update `SUBMISSION_CHECKLIST.md` with submission date
   - Create submission log entry

### Within 1 Week

1. **Upload Data to Zenodo**
   - Follow `DATA_REPOSITORY_TEMPLATE.md`
   - Get DOI
   - Update manuscript Data Availability section (for revision if needed)

2. **Prepare for Review**
   - Anticipate 2-3 month review time
   - Prepare response template for potential revisions

---

## 🎯 Quick Checklist

**Before Submission:**
- [ ] PDF compiled successfully
- [ ] All figures visible in PDF
- [ ] All tables visible in PDF
- [ ] All cross-references resolved (no "??")
- [ ] Grammar check completed
- [ ] Author information correct
- [ ] Acknowledgments complete
- [ ] Data availability statement ready (with placeholders)

**During Submission:**
- [ ] All required fields filled
- [ ] All files uploaded
- [ ] Submission confirmation received

**After Submission:**
- [ ] Confirmation saved
- [ ] Submission date logged
- [ ] Plan for data repository upload

---

## ⏱️ Time Estimate

| Task | Time |
|------|------|
| Compile PDF | 5-10 min |
| Verify PDF | 5 min |
| Grammar check | 30 min |
| Journal formatting | 10 min |
| Prepare package | 10 min |
| Upload & submit | 15 min |
| **Total** | **~1.5 hours** |

---

## 📞 Support

**If Issues Arise:**

1. **LaTeX Compilation Errors**
   - Check `COMPILATION_GUIDE.md`
   - Use Overleaf for online compilation
   - Check LaTeX log file for details

2. **Missing Files**
   - Run `python paper/check_manuscript.py`
   - Verify all files in `paper/` directory

3. **Journal Requirements**
   - Check journal website for author guidelines
   - Contact journal editorial office if unclear

---

**Status**: Ready for Step 1 (Compilation)  
**Next Action**: Compile PDF using local LaTeX or Overleaf

