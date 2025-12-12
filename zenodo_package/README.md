# Simulation Data for: Emergent molecular complexity in physics based simulations of prebiotic chemistry

**Authors**: Michał Klawikowski  
**Affiliation**: Independent Researcher, Pruszcz Gdański, Poland  
**Publication**: Submitted to Discover Life (Springer Nature)  
**Date**: 2025-12-12

## Dataset Description

This dataset contains raw simulation outputs from 43 independent prebiotic chemistry simulations conducted as part of Phase 2B extended runs. The simulations explore three distinct prebiotic scenarios: Miller-Urey reducing atmosphere, alkaline hydrothermal vents, and formamide-rich environments.

## Dataset Structure

```
zenodo_package/
├── README.md (this file)
├── configs/
│   ├── phase2_miller_urey_extended_SUPER_FAST.yaml
│   ├── phase2_hydrothermal_extended_SUPER_FAST.yaml
│   └── phase2_formamide_extended_SUPER_FAST.yaml
├── results/
│   ├── miller_urey_extended/
│   │   ├── run_1/
│   │   │   ├── results.json
│   │   │   ├── molecules.json
│   │   │   └── snapshots/ (representative samples)
│   │   └── ... (runs 2-18)
│   ├── hydrothermal_extended/
│   │   └── ... (runs 1-17)
│   └── formamide_extended/
│       └── ... (runs 1-8)
```

## Simulation Parameters

- **Total runs**: 43 (18 Miller-Urey + 17 Hydrothermal + 8 Formamide)
- **Simulation steps**: 500,000 per run
- **Simulated time**: ~140 hours per run
- **Particles**: 360-400 atoms per simulation
- **Box size**: 100 × 100 Å
- **Temperature**: 298 K (Miller-Urey, Formamide) or 373 K (Hydrothermal)

## File Formats

### results.json
Contains simulation metadata, statistics, and summary:
- Simulation parameters (n_particles, max_steps, dt, etc.)
- Final statistics (total molecules, unique species, etc.)
- Energy conservation metrics
- Completion status

### molecules.json
Contains all detected molecular species:
- Molecular formulas
- Abundances over time
- First detection step
- SMILES representations (where available)
- PubChem matches (where available)

### snapshots/
Periodic state snapshots saved every 50,000 steps:
- Particle positions and velocities
- Bond connectivity
- Molecular structures
- System state at specific time points

## Configuration Files

YAML configuration files defining simulation parameters for each scenario:
- `phase2_miller_urey_extended_SUPER_FAST.yaml`: Miller-Urey reducing atmosphere
- `phase2_hydrothermal_extended_SUPER_FAST.yaml`: Alkaline hydrothermal vents
- `phase2_formamide_extended_SUPER_FAST.yaml`: Formamide-rich environment

## Data Usage

This dataset supports the findings reported in the manuscript:
- Molecular diversity analysis (2,315 unique species detected)
- Autocatalytic cycle detection (769,315 cycles identified)
- Reaction network topology analysis
- Scenario comparison (Miller-Urey vs. Hydrothermal vs. Formamide)

## Code Availability

Simulation code and analysis scripts are available at:
https://github.com/ProhunterPL/live2.0

## Citation

If you use this dataset, please cite:
Klawikowski, M. (2025). Emergent molecular complexity in physics based simulations of prebiotic chemistry. *Discover Life* (submitted).

## License

This dataset is provided under [specify license - e.g., CC BY 4.0]

## Contact

Michał Klawikowski  
Email: klawikowski@klawikowski.pl

