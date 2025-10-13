# Phase 2: Open-Ended Experiments

**Status**: 📋 Ready to Execute  
**Duration**: 2 weeks (estimated)  
**Goal**: Generate novel results for publication

---

## Overview

Phase 2 runs long-term simulations under three prebiotic scenarios to:
1. **Generate novel molecules** (target: 100+)
2. **Identify compounds** using PubChem Matcher v2
3. **Analyze reaction networks** and autocatalytic cycles
4. **Validate findings** with quantum chemistry (DFT)

---

## Scenarios

### 1. Miller-Urey (1953 Classic)

**Conditions**:
- **Atmosphere**: CH₄, NH₃, H₂O, H₂ (reducing)
- **Temperature**: 298K (room temperature)
- **Energy source**: Electrical discharge (pulses every 1000 steps)
- **Duration**: 10⁷ steps × 10 runs

**Expected Products**:
- Glycine (NH₂CH₂COOH)
- Alanine (NH₂CH(CH₃)COOH)
- Formaldehyde (CH₂O)
- Hydrogen cyanide (HCN)
- Formic acid (HCOOH)

**Configuration**: `configs/phase2_miller_urey.yaml`

### 2. Hydrothermal Vent

**Conditions**:
- **Fluids**: H₂, H₂S, CO₂, H₂O, NH₃
- **Temperature**: 373K (100°C), variable 323-423K (50-150°C)
- **pH**: 10.0 (alkaline)
- **Catalysts**: FeS, FeS₂ (iron-sulfur minerals)
- **Duration**: 10⁷ steps × 10 runs

**Expected Products**:
- Formic acid (HCOOH)
- Acetic acid (CH₃COOH)
- Pyruvic acid (CH₃COCOOH)
- Iron-sulfur clusters ([FeS]ₙ)
- Thiols (R-SH)

**Configuration**: `configs/phase2_hydrothermal.yaml`

### 3. Formamide-Rich

**Conditions**:
- **Solvent**: HCONH₂ (formamide, 60%), H₂O (20%), NH₃ + HCN (20%)
- **Temperature**: 323K (50°C)
- **Energy source**: UV radiation (pulses every 500 steps)
- **Catalysts**: TiO₂ (photocatalyst), ZnS, montmorillonite clay
- **Duration**: 10⁷ steps × 10 runs

**Expected Products**:
- **Nucleobases**: Adenine (C₅H₅N₅), Guanine, Cytosine, Uracil
- **Amino acids**: Glycine
- **Sugars**: Glycolaldehyde (C₂H₄O₂)
- Carbodiimide (HN=C=NH)

**Configuration**: `configs/phase2_formamide.yaml`

---

## Workflow

### Step 1: Run Batch Simulations

```bash
# Run all scenarios (30 simulations total)
python scripts/run_phase2_batch.py --all --runs 10

# Or run individual scenarios
python scripts/run_phase2_batch.py --scenario miller_urey --runs 10
python scripts/run_phase2_batch.py --scenario hydrothermal --runs 10
python scripts/run_phase2_batch.py --scenario formamide --runs 10
```

**Expected Duration**:
- ~2 hours per run
- 10 runs per scenario
- 3 scenarios
- **Total**: ~60 hours (2.5 days of compute time)

**Parallelization**: Can run multiple scenarios simultaneously on different machines/GPUs.

### Step 2: Analyze Results

```bash
# Analyze all scenarios
python scripts/analyze_phase2_results.py --all

# Or analyze individual scenarios
python scripts/analyze_phase2_results.py --scenario miller_urey
```

**Analysis includes**:
- Extract unique molecules from all runs
- Match to PubChem using MatcherV2 (ML + multi-metric)
- Build reaction networks
- Identify autocatalytic cycles
- Generate statistics and comparisons

### Step 3: Generate Figures

```bash
# Generate publication-ready figures
python scripts/plot_phase2_figures.py
```

**Figures to generate**:
- **Figure 5**: Molecule diversity across scenarios
- **Figure 6**: Reaction network visualization
- **Figure 7**: Autocatalytic cycles identified
- **Figure S2**: Size distribution of molecules
- **Figure S3**: Temporal evolution of complexity

### Step 4: DFT Validation (Top 5 Molecules)

For the 5 most interesting novel molecules:
1. Export structures to .xyz format
2. Run DFT calculations (B3LYP/6-31G* or similar)
3. Validate:
   - Geometry optimization
   - Vibrational frequencies (real, not imaginary)
   - Energy comparison with simulation
   - Stability assessment

**Tools**: Gaussian, ORCA, or Psi4

---

## Expected Deliverables

### Data
- ✅ 30 simulation trajectories (3 scenarios × 10 runs)
- ✅ 100+ unique molecules catalog
- ✅ PubChem matches for all molecules
- ✅ Reaction networks (JSON format)
- ✅ Autocatalytic cycles list

### Analysis
- ✅ Statistical comparison of scenarios
- ✅ Novelty assessment (high-confidence matches vs unknowns)
- ✅ Complexity metrics (molecule size, bond diversity)
- ✅ Reaction kinetics analysis

### Validation
- ✅ Top 5 molecules validated with DFT
- ✅ Thermodynamic feasibility confirmed
- ✅ Chemical plausibility verified

### Figures (Publication-Ready)
- ✅ Figure 5: Molecule diversity
- ✅ Figure 6: Reaction networks
- ✅ Figure 7: Autocatalytic cycles
- ✅ Supplementary figures (S2, S3)

---

## Success Criteria

### Minimum Requirements
- [ ] 30/30 simulations completed successfully
- [ ] 100+ unique molecules identified
- [ ] 50+ molecules matched to PubChem (confidence > 0.5)
- [ ] 10+ autocatalytic cycles found
- [ ] 5/5 top molecules validated with DFT

### Excellence Criteria
- [ ] 150+ unique molecules
- [ ] 20+ high-confidence PubChem matches
- [ ] 5+ novel molecules not in PubChem
- [ ] 20+ autocatalytic cycles
- [ ] Evidence of emergent complexity

---

## Timeline

| Week | Tasks | Deliverables |
|------|-------|--------------|
| **Week 1** | Run simulations (30 total) | Raw data, trajectories |
| **Week 2** | Analysis + validation | Molecule catalog, figures, DFT validation |

**Parallel work**:
- Simulations can run overnight/weekend
- Analysis can start as soon as first runs complete
- DFT validation can be done for top candidates as identified

---

## Technical Details

### Computational Requirements

**Per Simulation**:
- Steps: 10⁷ (10 million)
- Particles: 2000
- Estimated time: ~2 hours (GPU)
- Memory: ~2 GB RAM
- Disk: ~500 MB per run

**Total**:
- 30 runs
- ~60 hours compute time
- ~60 GB RAM (if parallel)
- ~15 GB disk space

**Recommended**: 
- NVIDIA GPU (RTX 3070 or better)
- 16 GB RAM
- SSD storage

### Output Format

Each simulation produces:
```
results/phase2/{scenario}/{run_id}/
├── final_state.json       # Final positions, velocities
├── molecules.json         # All molecules detected
├── reactions.json         # All reactions recorded
├── metrics.csv            # Time series metrics
├── snapshots/             # Periodic state snapshots
│   ├── step_0000000.json
│   ├── step_0050000.json
│   └── ...
└── validation_log.json    # Thermodynamic validation
```

### Analysis Output

```
results/phase2/
├── miller_urey_analysis.json
├── hydrothermal_analysis.json
├── formamide_analysis.json
├── phase2_comparison.json
└── figures/
    ├── figure5_diversity.png
    ├── figure6_networks.png
    └── figure7_cycles.png
```

---

## Integration with Phase 1

Phase 2 builds directly on Phase 1 infrastructure:

| Phase 1 Component | Phase 2 Usage |
|-------------------|---------------|
| Thermodynamic validation | Continuous monitoring during runs |
| Physics database | Literature parameters used |
| Benchmark reactions | Expected products tracked |
| MatcherV2 | Molecule identification |

**Synergy**: Phase 1 validation ensures Phase 2 results are scientifically rigorous!

---

## Troubleshooting

### Common Issues

1. **Simulation crashes**
   - Check thermodynamic validation alerts
   - Verify GPU memory not exhausted
   - Reduce `n_particles` if needed

2. **No molecules forming**
   - Increase `pulse_energy` for energy injection
   - Check `bond_formation_threshold`
   - Verify initial composition loaded correctly

3. **MatcherV2 fails**
   - Check internet connection (PubChem API)
   - Verify RDKit installation
   - Fall back to similarity-only matching

4. **Long runtime**
   - Use GPU acceleration (Taichi CUDA)
   - Reduce `max_steps` for testing
   - Parallelize scenarios across machines

---

## Next Steps After Phase 2

Once Phase 2 is complete:
1. Review results with Phase 2 completion checklist
2. Begin Phase 3: Paper writing
3. Generate all paper figures (7+ main + supplementary)
4. Write Methods, Results, Discussion sections

---

## Files

### Configurations
- `configs/phase2_miller_urey.yaml`
- `configs/phase2_hydrothermal.yaml`
- `configs/phase2_formamide.yaml`

### Scripts
- `scripts/run_phase2_batch.py` - Batch simulation runner
- `scripts/analyze_phase2_results.py` - Results analyzer (uses MatcherV2)
- `scripts/plot_phase2_figures.py` - Figure generator (TODO)

### Documentation
- This file: `docs/PHASE2_EXPERIMENTS.md`
- Phase 1 summary: `docs/PHASE1_COMPLETION_SUMMARY.md`
- Overall roadmap: `docs/VALIDATION_ROADMAP.md`

---

**Status**: 📋 **Ready to Execute**  
**Prerequisite**: Phase 1 complete ✅  
**Next**: Run simulations and analyze results!

*Last updated: October 13, 2025*

