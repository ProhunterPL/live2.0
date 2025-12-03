# Final Guide: Generate Missing Figures for Manuscript

**Date**: 2025-01-23  
**Status**: ⚠️ **ACTION REQUIRED**

---

## 🔴 Problem: NumPy 2.x / RDKit Compatibility

**Error**: `AttributeError: _ARRAY_API not found`  
**Cause**: RDKit compiled with NumPy 1.x, environment has NumPy 2.3.3

---

## ✅ Solutions (Choose One)

### Option 1: Fix Environment (Best for Full Features) ⭐

```bash
# Downgrade NumPy
pip install 'numpy<2'

# Verify
python -c "from rdkit import Chem; print('OK')"

# Then run full script
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures
```

### Option 2: Use Simple Script (No Fix Needed) ⭐ RECOMMENDED

```bash
# First, extract molecules if needed
python scripts/extract_hydrothermal_molecules.py  # (adapt for miller_urey)

# Then generate structures (simple version, no RDKit)
python scripts/generate_molecular_structures_simple.py \
    --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures_panel.png
```

### Option 3: Skip Structures for Now (Fastest)

```bash
# Generate everything except structures
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures \
    --skip-structures
```

---

## 📋 What Needs to Be Generated

### ✅ Already Have (OK for Submission):
- **Figure 1**: Thermodynamic validation (synthetic data OK)
- **Figure 2**: Benchmark reactions (synthetic data OK)

### ⚠️ Need to Generate:
- **Molecular Structures Panel**: Top 5 molecules with PubChem matches
- **Reaction Network Example**: Visualization from real simulation

---

## 🎯 Recommended Workflow

### Step 1: Extract Molecules (If molecules.json is Empty)

```bash
# Extract molecules from snapshots
python scripts/extract_hydrothermal_molecules.py
# (Note: Adapt script path for miller_urey if needed)

# Or use molecule_extractor directly:
python -c "
from pathlib import Path
import sys
sys.path.insert(0, '.')
from backend.sim.molecule_extractor import extract_molecules_from_results

result = extract_molecules_from_results(
    'results/phase2b_additional/miller_urey_extended/run_1',
    output_dir='results/phase2b_additional/miller_urey_extended/run_1/analysis',
    export_for_matcher=False
)
print(f'Extracted {len(result[\"molecules\"])} molecules')
"
```

### Step 2: Generate Structures Panel

**If NumPy < 2**:
```bash
python scripts/generate_paper_figures_from_real_data.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output-dir paper/figures
```

**If NumPy >= 2** (use simple version):
```bash
python scripts/generate_molecular_structures_simple.py \
    --molecules-file results/phase2b_additional/miller_urey_extended/run_1/molecules.json \
    --output paper/figures/molecular_structures_panel.png \
    --top-n 5
```

### Step 3: Generate Reaction Network

```bash
python scripts/reaction_network_analyzer.py \
    results/phase2b_additional/miller_urey_extended/run_1 \
    --output analysis/reaction_network_example \
    --export both

python scripts/network_visualizer.py \
    analysis/reaction_network_example/reaction_network.json \
    --max-nodes 50 \
    --output paper/figures/reaction_network_example.png
```

---

## 📝 Integration with Manuscript

### After Generating Figures:

1. **Molecular Structures Panel**:
   - Add to Figure 6 as panel E, OR
   - Add new section in Results (3.5: Example Molecular Structures)

2. **Reaction Network Example**:
   - Add to Figure 4 as panel E, OR
   - Add as Figure 7 (Example Reaction Network)

3. **Update Captions**:
   ```latex
   \caption{\textbf{Example molecular structures detected in simulations.}
   Top 5 most abundant molecules with PubChem database matches.
   (A) [Molecule 1], (B) [Molecule 2], etc.}
   ```

---

## ✅ Checklist

### Before Submission:
- [ ] **Molecular structures panel** (top 5 molecules)
- [ ] **Reaction network example** (visualization)
- [ ] Figures integrated into manuscript
- [ ] Captions updated

### Optional (Can Add Later):
- [ ] Thermodynamic plots from real data
- [ ] Benchmark plots from real data
- [ ] More molecular structures
- [ ] More network examples

---

## 🎯 Quick Decision Tree

**Q: Do I have NumPy < 2?**
- **Yes** → Use full script (`generate_paper_figures_from_real_data.py`)
- **No** → Use simple script (`generate_molecular_structures_simple.py`)

**Q: Is molecules.json empty?**
- **Yes** → Extract first (`extract_molecules_from_results`)
- **No** → Generate structures directly

**Q: Need to submit quickly?**
- **Yes** → Skip structures for now, add in revision
- **No** → Generate all figures

---

**Status**: ✅ **SCRIPTS READY**  
**Action**: Choose solution and generate figures

