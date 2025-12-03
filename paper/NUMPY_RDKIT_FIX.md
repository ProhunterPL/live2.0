# Fix NumPy 2.x / RDKit Compatibility Issue

**Date**: 2025-01-23  
**Problem**: RDKit compiled with NumPy 1.x cannot run with NumPy 2.3.3

---

## 🔴 Problem

```
AttributeError: _ARRAY_API not found
A module that was compiled using NumPy 1.x cannot be run in NumPy 2.3.3
```

**Cause**: RDKit was compiled against NumPy 1.x, but environment has NumPy 2.3.3

---

## ✅ Solutions

### Solution 1: Downgrade NumPy (Quickest) ⭐ RECOMMENDED

```bash
pip install 'numpy<2'
```

**Then verify**:
```bash
python -c "import numpy; print(numpy.__version__)"  # Should be 1.x
python -c "from rdkit import Chem; print('RDKit OK')"  # Should work
```

### Solution 2: Use Simple Script (No RDKit)

**For molecular structures**, use the simple version that doesn't require RDKit:

```bash
python scripts/generate_molecular_structures_simple.py \
    --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures_panel.png
```

This script:
- ✅ Doesn't require RDKit
- ✅ Works with NumPy 2.x
- ✅ Queries PubChem via REST API
- ✅ Generates text-based molecular structures panel

### Solution 3: Skip Matcher (For Now)

If you just need to submit quickly:

1. **Skip molecular structures** for now (can add later)
2. **Generate other figures** (thermodynamic, benchmark, network)
3. **Add structures in revision** if journal requests

---

## 🎯 Recommended Workflow

### For Quick Submission:

1. **Fix NumPy** (if you want structures):
   ```bash
   pip install 'numpy<2'
   ```

2. **Or use simple script** (no fix needed):
   ```bash
   python scripts/generate_molecular_structures_simple.py \
       --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
       --output paper/figures/molecular_structures_panel.png
   ```

3. **Generate other figures** (these don't need RDKit):
   ```bash
   python scripts/generate_paper_figures_from_real_data.py \
       --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
       --output-dir paper/figures \
       --skip-structures  # Skip if NumPy issue persists
   ```

---

## 📋 Updated Script Status

### `generate_paper_figures_from_real_data.py`
- ✅ **Fixed**: Lazy imports, handles NumPy/RDKit issues gracefully
- ✅ **Thermodynamic**: Works (no RDKit needed)
- ✅ **Benchmark**: Works (no RDKit needed)
- ⚠️ **Structures**: Requires NumPy <2 or use simple script
- ✅ **Network**: Works (no RDKit needed)

### `generate_molecular_structures_simple.py` (NEW)
- ✅ **No RDKit dependency**
- ✅ **Works with NumPy 2.x**
- ✅ **Queries PubChem via REST API**
- ✅ **Generates publication-ready panel**

---

## ✅ Action Plan

### Option A: Fix Environment (Best for Full Features)

```bash
# 1. Downgrade NumPy
pip install 'numpy<2'

# 2. Verify
python -c "from rdkit import Chem; print('OK')"

# 3. Run full script
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures
```

### Option B: Use Simple Script (No Fix Needed)

```bash
# 1. Generate structures with simple script
python scripts/generate_molecular_structures_simple.py \
    --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures_panel.png

# 2. Generate other figures (skip structures)
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures \
    --skip-structures
```

### Option C: Skip Structures for Now

```bash
# Generate everything except structures
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures \
    --skip-structures
```

---

## 🎯 Recommendation

**For fastest submission**:
- Use **Option B** (simple script) - no environment changes needed
- Or use **Option C** (skip structures) - add later if needed

**For best results**:
- Use **Option A** (fix NumPy) - full features available

---

**Status**: ✅ **FIXED** (script handles gracefully)  
**Alternative**: ✅ **Simple script available** (no RDKit needed)

