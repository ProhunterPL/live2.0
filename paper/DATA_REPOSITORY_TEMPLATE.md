# Data Repository Template - Zenodo/Dryad

**Purpose**: Template for preparing simulation data for public repository (Zenodo/Dryad)

**Status**: Template ready, awaiting data package preparation

---

## 📦 Repository Structure

```
live2.0-phase2b-data/
├── README.md                    # This file (repository description)
├── LICENSE                      # License file (MIT/CC-BY-4.0)
├── CITATION.cff                 # Citation metadata
├── metadata.json                # Repository metadata
│
├── raw_data/                    # Raw simulation outputs
│   ├── miller_urey_extended/
│   │   ├── run_1/
│   │   │   ├── results.json
│   │   │   ├── molecules.json
│   │   │   ├── snapshots/       # Every 50K steps
│   │   │   └── simulation.log
│   │   ├── run_2/
│   │   └── ... (runs 1-18)
│   ├── hydrothermal_extended/
│   │   └── ... (runs 1-17)
│   └── formamide_extended/
│       └── ... (runs 1-8)
│
├── processed_data/              # Processed analysis results
│   ├── summary_statistics.json
│   ├── scenario_comparison.json
│   ├── autocatalytic_cycles.json
│   ├── network_topology.json
│   └── novel_molecules.json
│
├── analysis_scripts/            # Code for reproducing analysis
│   ├── analyze_phase2b_complete.py
│   ├── extract_molecules.py
│   ├── detect_autocatalysis.py
│   └── generate_figures.py
│
├── documentation/
│   ├── DATA_DESCRIPTION.md      # Detailed data description
│   ├── FILE_FORMATS.md          # File format specifications
│   └── REPRODUCTION_GUIDE.md    # How to reproduce results
│
└── supplementary/
    ├── tables/                  # All tables from paper
    ├── figures/                 # All figures from paper
    └── validation/               # Validation results
```

---

## 📋 Repository Metadata (Zenodo)

### Basic Information

**Title**: 
```
Live 2.0 Phase 2B: Prebiotic Chemistry Simulation Data
```

**Description**:
```
This dataset contains raw and processed results from 30 independent prebiotic chemistry simulations (Phase 2B) across three scenarios: Miller-Urey reducing atmosphere, alkaline hydrothermal vents, and formamide-rich environments. Each simulation ran for 500,000 steps (~140 hours simulated time) and generated molecular diversity, reaction networks, and autocatalytic cycle data. The dataset includes complete particle configurations, bond graphs, molecular inventories, and statistical analysis results used in the publication "Emergent Molecular Complexity in Prebiotic Chemistry Simulations: A Physics-Based Approach".
```

**Keywords**:
```
prebiotic chemistry, origin of life, molecular dynamics, autocatalysis, emergent complexity, Miller-Urey, hydrothermal vents, formamide, reaction networks, chemical evolution
```

**Creators**:
- Name: Michał Klawikowski
- Affiliation: Live 2.0
- ORCID: [to be added if available]

**Contributors**: (if applicable)

**Version**: 1.0.0

**Publication Date**: [to be set after submission]

**License**: 
- Code: MIT License
- Data: CC-BY-4.0 (Creative Commons Attribution 4.0 International)

**Related Publications**:
- Manuscript DOI: [to be added after publication]
- GitHub Repository: https://github.com/ProhunterPL/live2.0

---

## 📄 README.md Template

```markdown
# Live 2.0 Phase 2B: Prebiotic Chemistry Simulation Data

## Overview

This repository contains simulation data from Phase 2B of the Live 2.0 project, 
consisting of 30 independent prebiotic chemistry simulations across three scenarios.

## Dataset Description

### Scenarios
1. **Miller-Urey** (18 runs): Reducing atmosphere conditions
2. **Hydrothermal** (17 runs): Alkaline hydrothermal vent conditions  
3. **Formamide** (8 runs): Formamide-rich environment conditions

### Data Contents
- Raw simulation outputs (results.json, molecules.json, snapshots)
- Processed analysis results (statistics, network topology, autocatalytic cycles)
- Analysis scripts for reproducing results
- Documentation and file format specifications

## File Formats

See `documentation/FILE_FORMATS.md` for detailed format specifications.

## Citation

If you use this dataset, please cite:

```
Klawikowski, M. (2025). Live 2.0 Phase 2B: Prebiotic Chemistry Simulation Data. 
Zenodo. https://doi.org/[DOI]
```

## License

- Code: MIT License
- Data: CC-BY-4.0

## Contact

Michał Klawikowski  
Email: klawikowski@klawikowski.pl  
GitHub: https://github.com/ProhunterPL/live2.0
```

---

## 📊 Data Package Checklist

### Before Upload

- [ ] All raw simulation outputs collected (43 runs total)
- [ ] Processed analysis results generated
- [ ] Analysis scripts included and tested
- [ ] Documentation complete (README, file formats, reproduction guide)
- [ ] License files added (MIT for code, CC-BY-4.0 for data)
- [ ] Citation metadata (CITATION.cff) created
- [ ] File sizes checked (Zenodo limit: 50 GB per dataset)
- [ ] Sensitive data removed (if any)
- [ ] Data validated (no corrupted files)

### Zenodo Upload Steps

1. **Create Zenodo Account** (if not exists)
   - Sign up at https://zenodo.org
   - Link ORCID if available

2. **Create New Upload**
   - Click "New Upload"
   - Select "Dataset" as upload type

3. **Fill Metadata**
   - Use metadata from `metadata.json` template above
   - Add all creators, keywords, description
   - Set license (CC-BY-4.0)
   - Set publication date (after manuscript acceptance)

4. **Upload Files**
   - Upload entire `live2.0-phase2b-data/` directory
   - Or create ZIP archive first (if large)
   - Wait for upload completion

5. **Reserve DOI**
   - Click "Reserve DOI" to get permanent DOI
   - Save DOI for manuscript update

6. **Publish**
   - Review all metadata
   - Click "Publish" (makes dataset public)
   - Or "Save" for draft (can publish later)

---

## 🔗 Integration with Manuscript

After obtaining Zenodo DOI:

1. Update `paper/manuscript_draft.tex`:
   ```latex
   All simulation data, analysis code, and visualization scripts are 
   publicly available at \url{https://github.com/ProhunterPL/live2.0} 
   (DOI: 10.5281/zenodo.XXXXXXX). Raw simulation outputs are deposited 
   at Zenodo (DOI: 10.5281/zenodo.XXXXXXX).
   ```

2. Update `paper/SUBMISSION_CHECKLIST.md`:
   - Mark data availability as complete
   - Add DOI to checklist

3. Add to manuscript acknowledgments (if required by journal)

---

## 📝 Notes

- **File Size**: Total dataset size ~5-10 GB (estimate)
- **Compression**: Consider ZIP compression for large files
- **Embargo**: No embargo needed (data can be public immediately)
- **Versioning**: Zenodo supports versioning (v1.0.0, v1.1.0, etc.)
- **Updates**: Can update dataset after publication if needed

---

**Created**: 2025-01-23  
**Status**: Template ready  
**Next Step**: Prepare actual data package and upload to Zenodo

