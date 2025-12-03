# Zenodo Deposit Guide - Live 2.0 Data Repository

**Date**: 2025-12-03  
**Purpose**: Create Zenodo deposit for simulation data and obtain DOI for manuscript

---

## 📋 Pre-Deposit Checklist

### Files to Include

1. **Simulation Results** (from `results/phase2b_additional/`)
   - All `run_*/` directories with:
     - `results.json` (simulation metadata)
     - `molecules.json` (detected molecules)
     - `snapshots/` (molecular snapshots)
     - `checkpoints/` (simulation checkpoints)

2. **Analysis Scripts** (from `scripts/`)
   - `run_phase2_full.py` (main simulation runner)
   - `analyze_results.py` (results analysis)
   - `reaction_network_analyzer.py` (network analysis)
   - `molecule_extractor.py` (molecule extraction)

3. **Configuration Files** (from `aws_test/configs/`)
   - `phase2_miller_urey_extended_SUPER_FAST.yaml`
   - `phase2_hydrothermal_extended_SUPER_FAST.yaml`
   - `phase2_formamide_extended_SUPER_FAST.yaml`

4. **Documentation** (from `docs/`)
   - `INDEX.md` (main documentation)
   - `NAVIGATION_GUIDE.md` (project structure)
   - `VALIDATION_ROADMAP.md` (validation methods)

5. **README for Zenodo**
   - Create `ZENODO_README.md` (see template below)

---

## 📝 Zenodo Deposit Metadata

### Title
```
Simulation Data for "Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach"
```

### Authors
- **Name**: Michał Klawikowski
- **Affiliation**: Live 2.0
- **ORCID**: (if available)

### Description
```
This dataset contains simulation results, analysis scripts, and configuration files for the prebiotic chemistry simulations described in Klawikowski (2025). The dataset includes:

1. **Simulation Results**: 30 complete simulation runs (18 Miller-Urey, 17 Hydrothermal, 8 Formamide) with 500,000 steps each, including:
   - Molecular snapshots at 50,000-step intervals
   - Detected molecules with formulas and abundances
   - Reaction network data
   - Autocatalytic cycle detections

2. **Analysis Scripts**: Python scripts for:
   - Running simulations (run_phase2_full.py)
   - Extracting molecules from snapshots (molecule_extractor.py)
   - Analyzing reaction networks (reaction_network_analyzer.py)
   - Statistical analysis (analyze_results.py)

3. **Configuration Files**: YAML configuration files for all three prebiotic scenarios (Miller-Urey, Hydrothermal, Formamide)

4. **Documentation**: Complete project documentation including validation methods, data structure, and usage instructions

**Data Format**: JSON (results, molecules, snapshots), YAML (configurations), Python (scripts), Markdown (documentation)

**Total Size**: ~XXX MB (estimate)

**License**: MIT License (code) / CC-BY-4.0 (data)

**Related Publication**: Klawikowski, M. (2025). "Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach" (submitted)

**Repository**: https://github.com/ProhunterPL/live2.0
```

### Keywords
```
prebiotic chemistry, origin of life, molecular dynamics, autocatalysis, reaction networks, Miller-Urey, hydrothermal vents, formamide, computational chemistry, Taichi
```

### License
- **Code**: MIT License
- **Data**: CC-BY-4.0

### Version
```
1.0.0
```

### Publication Date
```
2025-12-03 (or date of deposit)
```

### Related Identifiers
- **Type**: IsSupplementTo
- **Identifier**: (GitHub repository URL: https://github.com/ProhunterPL/live2.0)

---

## 📄 ZENODO_README.md Template

Create this file in the root of your deposit:

```markdown
# Live 2.0 Simulation Data Repository

**DOI**: [To be assigned by Zenodo]

**Related Publication**: Klawikowski, M. (2025). "Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach" (submitted)

**Repository**: https://github.com/ProhunterPL/live2.0

---

## Dataset Contents

### 1. Simulation Results (`results/`)

Contains 30 complete simulation runs organized by scenario:

- `phase2b_additional/miller_urey_extended/run_1/` through `run_18/`
- `phase2b_additional/hydrothermal_extended/run_1/` through `run_17/`
- `phase2b_additional/formamide_extended/run_1/` through `run_8/`

Each run directory contains:
- `results.json`: Simulation metadata, statistics, and summary
- `molecules.json`: Detected molecules with formulas and abundances
- `snapshots/`: Molecular snapshots at 50,000-step intervals (step_*.json)
- `checkpoints/`: Simulation checkpoints at 100,000-step intervals

### 2. Analysis Scripts (`scripts/`)

- `run_phase2_full.py`: Main simulation runner
- `analyze_results.py`: Results analysis and statistics
- `reaction_network_analyzer.py`: Reaction network analysis
- `molecule_extractor.py`: Molecule extraction from snapshots

### 3. Configuration Files (`aws_test/configs/`)

YAML configuration files for all three scenarios:
- `phase2_miller_urey_extended_SUPER_FAST.yaml`
- `phase2_hydrothermal_extended_SUPER_FAST.yaml`
- `phase2_formamide_extended_SUPER_FAST.yaml`

### 4. Documentation (`docs/`)

- `INDEX.md`: Main documentation index
- `NAVIGATION_GUIDE.md`: Project structure guide
- `VALIDATION_ROADMAP.md`: Validation methods

---

## Data Format

### results.json
```json
{
  "simulation": {
    "n_particles": 1000,
    "max_steps": 500000,
    "dt": 0.001,
    ...
  },
  "statistics": {
    "total_molecules": 1234,
    "unique_species": 567,
    ...
  }
}
```

### molecules.json
```json
[
  {
    "formula": "CH2O",
    "abundance": 1234,
    "first_seen": 50000,
    ...
  },
  ...
]
```

### Snapshots (step_*.json)
```json
{
  "step": 50000,
  "bonds": [[0, 1], [1, 2], ...],
  "attributes": [[mass, charge_x, charge_y, charge_z], ...],
  ...
}
```

---

## Usage

### Running Simulations
```bash
python scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \
    --output results/test_run \
    --seed 100 \
    --steps 500000
```

### Extracting Molecules
```bash
python scripts/molecule_extractor.py \
    --results-dir results/phase2b_additional/miller_urey_extended/run_1 \
    --output molecules_extracted.json
```

### Analyzing Networks
```bash
python scripts/reaction_network_analyzer.py \
    --results-dirs results/phase2b_additional/miller_urey_extended/run_* \
    --output network_analysis/
```

---

## Citation

If you use this dataset, please cite:

```
Klawikowski, M. (2025). Simulation Data for "Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach". Zenodo. https://doi.org/[DOI]
```

And the related publication:

```
Klawikowski, M. (2025). Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach. [Journal], [Volume], [Pages]. https://doi.org/[Article DOI]
```

---

## License

- **Code**: MIT License (see LICENSE file)
- **Data**: CC-BY-4.0 (Creative Commons Attribution 4.0 International)

---

## Contact

- **Author**: Michał Klawikowski
- **Email**: klawikowski@klawikowski.pl
- **Repository**: https://github.com/ProhunterPL/live2.0
```

---

## 🚀 Step-by-Step Deposit Instructions

### 1. Prepare Files

```bash
# Create a clean directory for Zenodo deposit
mkdir zenodo_deposit
cd zenodo_deposit

# Copy simulation results (selective - may need to compress large files)
cp -r ../results/phase2b_additional/ .

# Copy scripts
cp -r ../scripts/ .

# Copy configs
cp -r ../aws_test/configs/ .

# Copy documentation
cp -r ../docs/ .

# Copy README
cp ../ZENODO_README.md .
```

### 2. Create Zenodo Account

1. Go to https://zenodo.org
2. Sign up / Log in
3. Click "Upload" → "New upload"

### 3. Fill Metadata

- **Title**: Use template above
- **Description**: Use template above
- **Authors**: Add Michał Klawikowski
- **Keywords**: Use template above
- **License**: Select MIT (code) and CC-BY-4.0 (data)
- **Version**: 1.0.0
- **Publication Date**: Today's date

### 4. Upload Files

- Upload all files from `zenodo_deposit/` directory
- Zenodo will create a ZIP automatically
- Check file sizes (Zenodo limit: 50 GB per deposit)

### 5. Reserve DOI

- Click "Reserve DOI" to get DOI before publication
- Copy DOI for manuscript

### 6. Publish

- Review all metadata
- Click "Publish"
- Copy final DOI

---

## 📝 Update Manuscript

After obtaining DOI, update `manuscript_draft.tex`:

```latex
All simulation data, analysis code, and visualization scripts are publicly available at \url{https://github.com/ProhunterPL/live2.0} (DOI: [GitHub DOI - to be assigned]). Raw simulation outputs are deposited at Zenodo (DOI: 10.5281/zenodo.XXXXXXX).
```

---

## ✅ Checklist

- [ ] All files prepared
- [ ] ZENODO_README.md created
- [ ] Zenodo account created
- [ ] Metadata filled
- [ ] Files uploaded
- [ ] DOI reserved
- [ ] Deposit published
- [ ] DOI added to manuscript
- [ ] Manuscript updated

---

**Status**: ✅ **GUIDE READY**  
**Next Step**: Prepare files and create Zenodo deposit

