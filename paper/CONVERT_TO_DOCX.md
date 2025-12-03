# Convert LaTeX to DOCX for Submission

**Date**: 2025-01-23  
**Issue**: Journal prefers editable .docx files over PDF

---

## 🎯 Quick Solution

### Option 1: Overleaf Export (Easiest) ⭐ RECOMMENDED

1. **Go to Overleaf project**
   - Open your manuscript in Overleaf
   
2. **Export to DOCX**
   - Click "Menu" (top left)
   - Select "Download" → "Source"
   - Or use "Download" → "DOCX" (if available)
   
3. **Alternative in Overleaf**:
   - Click "Menu" → "Compiler" → "LaTeX"
   - Then "Download" → "DOCX" (if option appears)

**Note**: Overleaf may not have direct DOCX export. If not available, use Option 2.

---

### Option 2: Pandoc Conversion (Best Quality)

#### Step 1: Install Pandoc

**Windows:**
```powershell
# Using Chocolatey
choco install pandoc

# Or download from: https://pandoc.org/installing.html
```

**macOS:**
```bash
brew install pandoc
```

**Linux:**
```bash
sudo apt-get install pandoc  # Ubuntu/Debian
```

#### Step 2: Convert LaTeX to DOCX

```bash
cd paper

# Basic conversion
pandoc manuscript_draft.tex -o manuscript_draft.docx --bibliography=references.bib --citeproc

# Better conversion with math support
pandoc manuscript_draft.tex -o manuscript_draft.docx \
  --bibliography=references.bib \
  --citeproc \
  --mathml \
  --standalone
```

#### Step 3: Post-Processing

After conversion, open in Word and:
- [ ] Check all figures are embedded
- [ ] Verify tables format correctly
- [ ] Check math equations render properly
- [ ] Verify citations format correctly
- [ ] Check cross-references

---

### Option 3: Online Converter

**Tools:**
- [CloudConvert](https://cloudconvert.com/latex-to-docx) - Upload .tex file
- [Zamzar](https://www.zamzar.com/convert/latex-to-docx/) - Free conversion
- [Online-Convert](https://www.online-convert.com/) - Multiple format support

**Limitations:**
- May not handle complex LaTeX well
- Bibliography may need manual fixing
- Math equations may not convert perfectly

---

### Option 4: Manual Conversion (If Others Fail)

1. **Compile PDF in Overleaf**
2. **Open PDF in Word** (Word 2013+ can import PDF)
   - File → Open → Select PDF
   - Word will convert PDF to editable format
3. **Clean up formatting**
   - Fix tables
   - Re-insert figures if needed
   - Fix citations format

---

## ⚠️ Important Notes

### What May Not Convert Well:

1. **Complex LaTeX Commands**
   - Custom commands may not work
   - Some packages may not be supported

2. **Bibliography**
   - May need to use `--citeproc` flag
   - Or manually format citations in Word

3. **Math Equations**
   - Use `--mathml` for better math support
   - Or convert to images

4. **Figures**
   - May need to re-insert manually
   - Check all figures are visible

5. **Tables**
   - LaTeX tables may need reformatting
   - Check table alignment

---

## ✅ Recommended Workflow

### Best Approach:

1. **Try Overleaf Export First** (if available)
   - Easiest, maintains formatting

2. **If Not Available, Use Pandoc**
   ```bash
   pandoc manuscript_draft.tex -o manuscript_draft.docx \
     --bibliography=references.bib \
     --citeproc \
     --mathml \
     --standalone
   ```

3. **Open in Word and Fix**
   - Check formatting
   - Fix any issues
   - Save as .docx

4. **Verify Before Submission**
   - All figures present
   - All tables formatted
   - All citations working
   - Math equations readable

---

## 🔧 Troubleshooting

### Issue: Pandoc Not Found

**Solution**: Install pandoc (see Option 2, Step 1)

### Issue: Bibliography Not Converting

**Solution**: 
```bash
# Use citeproc
pandoc manuscript_draft.tex -o manuscript_draft.docx \
  --bibliography=references.bib \
  --citeproc
```

### Issue: Math Equations Broken

**Solution**:
```bash
# Use mathml for better math support
pandoc manuscript_draft.tex -o manuscript_draft.docx --mathml
```

### Issue: Figures Missing

**Solution**: 
- Open DOCX in Word
- Manually insert figures from `paper/figures/` directory
- Or use PDF import method (Option 4)

---

## 📋 Final Checklist

Before submitting DOCX:

- [ ] All text converted correctly
- [ ] All figures embedded and visible
- [ ] All tables formatted correctly
- [ ] All citations formatted (numbers or author-year, depending on journal)
- [ ] Math equations readable
- [ ] Cross-references working (or removed if broken)
- [ ] Page numbers present
- [ ] Line numbers removed (if journal doesn't want them)
- [ ] File size reasonable (<10 MB)

---

## 🎯 Quick Command (Copy-Paste)

If you have pandoc installed:

```bash
cd paper
pandoc manuscript_draft.tex -o manuscript_draft.docx --bibliography=references.bib --citeproc --mathml --standalone
```

Then open `manuscript_draft.docx` in Word, check formatting, and save.

---

**Status**: Ready for conversion  
**Recommended**: Try Overleaf export first, then pandoc if needed

