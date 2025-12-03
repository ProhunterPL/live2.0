# Grammar Check Guide - Step 2

**Date**: 2025-01-23  
**Status**: PDF Compiled Successfully ✅  
**Next**: Final Grammar Check

---

## 📝 About the NIST Link Issue

**Link**: `https://cccbdb.nist.gov/` (Reference [15] in manuscript)

**Status**: Error 503 (likely temporary server issue)

**Action**: ✅ **No action needed**
- The link is correctly formatted in `references.bib`
- Journal reviewers understand that external links may be temporarily unavailable
- NIST databases are well-known and the URL is correct
- This is a standard reference format for databases

**Note**: If the error persists, you can add a note in the reference:
```bibtex
note={National Institute of Standards and Technology. URL accessed: [date]}
```

But this is **not required** for submission - the URL format is correct.

---

## 🔍 Grammar Check - Step 2

### Option 1: Grammarly (Recommended)

**Steps**:
1. Go to https://www.grammarly.com
2. Create free account (if needed)
3. Click "New" → "Upload" or paste text
4. Upload the PDF or copy text from Overleaf

**What to check**:
- [ ] Spelling errors
- [ ] Grammar mistakes
- [ ] Punctuation
- [ ] Clarity and readability
- [ ] Scientific terminology consistency

**Note**: Grammarly may flag some scientific terms as errors (e.g., "prebiotic", "autocatalytic") - these are **correct** in scientific context.

---

### Option 2: LanguageTool (Free, Open Source)

**Steps**:
1. Go to https://languagetool.org
2. Paste text or upload document
3. Select language: English (US) or English (UK)
4. Review suggestions

**Advantages**:
- Free and open source
- Good for scientific writing
- Can be installed locally

---

### Option 3: Microsoft Word

**Steps**:
1. Copy text from Overleaf PDF
2. Paste into Word
3. Run spell-check (F7)
4. Review grammar suggestions

**Note**: Word may not handle LaTeX formatting well - use PDF text or plain text.

---

### Option 4: Manual Review (Recommended for Scientific Papers)

**Checklist**:

#### Spelling & Typos
- [ ] All scientific terms spelled correctly
- [ ] Names (Miller-Urey, Kauffman, etc.) correct
- [ ] Chemical formulas correct (CH₂O, HCN, etc.)
- [ ] Numbers and units consistent

#### Grammar
- [ ] Subject-verb agreement
- [ ] Tense consistency (past tense for methods/results)
- [ ] Article usage (a/an/the)
- [ ] Plural/singular consistency

#### Scientific Writing Style
- [ ] Passive voice used appropriately (Methods section)
- [ ] Active voice used appropriately (Discussion section)
- [ ] Technical terms defined on first use
- [ ] Abbreviations defined on first use

#### Consistency
- [ ] Terminology consistent throughout:
  - "prebiotic chemistry" vs "prebiotic Chemistry"
  - "autocatalytic cycle" vs "autocatalytic Cycle"
  - "Miller-Urey" vs "Miller–Urey" (hyphen vs en-dash)
- [ ] Number formatting consistent:
  - "500,000" vs "500 000" vs "5×10⁵"
  - Decimals: "0.1" vs ".1"
- [ ] Units consistent:
  - "amu" vs "Da" (atomic mass units)
  - "a.u." vs "au" (atomic units)

#### Common Issues to Watch For

1. **Comma usage**:
   - ✅ "We conducted 30 simulations, each running for 500,000 steps."
   - ❌ "We conducted 30 simulations each running for 500,000 steps."

2. **Hyphenation**:
   - ✅ "prebiotic chemistry" (no hyphen)
   - ✅ "far-from-equilibrium" (hyphenated)
   - ✅ "self-sustaining" (hyphenated)

3. **Scientific notation**:
   - ✅ "10⁶" or "10^6" or "1×10⁶"
   - ✅ Consistent throughout

4. **Citations**:
   - ✅ "Miller-Urey experiment [12]"
   - ✅ "as shown in [12]"
   - ✅ Consistent format

---

## 📋 Specific Sections to Review

### Abstract
- [ ] Word count: ~250 words ✅
- [ ] All key findings mentioned
- [ ] No undefined abbreviations

### Introduction
- [ ] All citations present
- [ ] Historical context clear
- [ ] Research gap identified

### Methods
- [ ] Technical details clear
- [ ] All parameters defined
- [ ] Reproducibility information present

### Results
- [ ] All numbers accurate
- [ ] Statistical tests have p-values
- [ ] Figures referenced correctly

### Discussion
- [ ] Interpretations supported by data
- [ ] Limitations acknowledged
- [ ] Future work mentioned

### Conclusions
- [ ] Key findings summarized
- [ ] Significance stated
- [ ] No new information introduced

---

## ⚠️ Common Scientific Writing Errors

### 1. Tense Issues
- **Methods**: Past tense ✅ "We conducted simulations..."
- **Results**: Past tense ✅ "We detected 2,315 species..."
- **Discussion**: Present tense ✅ "These results suggest..."

### 2. Article Usage
- ✅ "the Miller-Urey experiment" (specific)
- ✅ "a prebiotic scenario" (general)
- ✅ "the origin of life" (specific concept)

### 3. Number Formatting
- ✅ "2,315 species" (comma for thousands)
- ✅ "500,000 steps" (comma for thousands)
- ✅ "0.1%" (decimal point)

### 4. Chemical Formulas
- ✅ "CH₂O" (subscripts in math mode: `CH$_2$O`)
- ✅ "H₂CO₃" (subscripts: `H$_2$CO$_3$`)
- ✅ "C₈H₁₂N₂O₃" (subscripts: `C$_8$H$_{12}$N$_2$O$_3$`)

---

## 🎯 Quick Grammar Check Workflow

1. **Extract text from PDF** (5 min)
   - Copy text from Overleaf PDF
   - Or use PDF text extraction tool

2. **Run automated check** (10 min)
   - Grammarly or LanguageTool
   - Review flagged issues
   - Accept/reject suggestions

3. **Manual scientific review** (15 min)
   - Check terminology consistency
   - Verify number formatting
   - Review chemical formulas
   - Check citation format

4. **Final read-through** (10 min)
   - Read entire manuscript
   - Check flow and clarity
   - Verify all sections complete

**Total time**: ~40 minutes

---

## ✅ After Grammar Check

Once grammar check is complete:

1. **Document any changes made**
   - Note major corrections
   - Keep track of terminology decisions

2. **Update Overleaf** (if changes made)
   - Make corrections in LaTeX source
   - Recompile PDF
   - Verify changes appear correctly

3. **Proceed to Step 3**: Journal Selection & Submission

---

## 📝 Notes

- **Scientific terminology**: Some words may be flagged by grammar checkers but are correct in scientific context (e.g., "prebiotic", "autocatalytic", "hypercycle")
- **LaTeX formatting**: Grammar checkers may not understand LaTeX commands - focus on the text content
- **Citations**: Format is handled by BibTeX - grammar checkers may flag citation format, but this is correct
- **Figures/Tables**: Captions should be checked separately

---

**Status**: Ready for Grammar Check  
**Estimated Time**: 30-40 minutes  
**Next Step**: After grammar check → Journal Submission

