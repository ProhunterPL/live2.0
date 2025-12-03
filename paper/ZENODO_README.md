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

