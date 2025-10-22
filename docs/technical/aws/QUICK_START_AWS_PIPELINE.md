# AWS Pipeline - Quick Start Guide

**When AWS tests complete, run this ONE command:**

```bash
bash scripts/aws_pipeline.sh <your-aws-ip> ~/.ssh/aws_key.pem
```

**Example:**
```bash
bash scripts/aws_pipeline.sh 54.123.45.67 ~/.ssh/aws_key.pem
```

This will automatically:
1. ✅ Download all results from AWS (~5-15 min)
2. ✅ Extract molecules from each simulation (~10-20 min)
3. ✅ Build reaction networks (~5 min)
4. ✅ Detect autocatalytic cycles (~3 min)
5. ✅ Generate all figures (~5 min)
6. ✅ Create comprehensive report

**Total time: ~30-45 minutes for 24 simulations**

---

## Then Update Paper

After pipeline completes:

```bash
# 1. Copy figures to paper
cp results/aws_batch/analysis/figures/*.png paper/figures/

# 2. Open analysis report
cat results/aws_batch/analysis/batch_analysis_report.json | python -m json.tool

# 3. Fill Results section in manuscript
cd paper
# Edit manuscript_draft.tex - fill [XX] placeholders with data from report

# 4. Compile paper
pdflatex manuscript_draft.tex
bibtex manuscript_draft
pdflatex manuscript_draft.tex
pdflatex manuscript_draft.tex

# 5. Review PDF
open manuscript_draft.pdf
```

---

## Next Steps

1. Fill Results section (~2-3 days)
2. Write Discussion (~1 week)
3. Write Abstract (~1 day)
4. Submit! (~Nov 30, 2025)

---

**For detailed instructions, see: `docs/AWS_RESULTS_PIPELINE.md`**  
**For paper guidelines, see: `paper/README.md`**

