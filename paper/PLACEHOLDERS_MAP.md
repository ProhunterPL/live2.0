# 📋 Placeholders Map - Mapping Analysis Results to Paper

**Data**: 2025-11-28  
**Status**: Ready for data filling  
**Purpose**: Complete mapping of all placeholders in manuscript to analysis outputs

---

## 🎯 Overview

This document maps every placeholder `[XX]`, `[YY]`, etc. in `manuscript_draft.tex` to specific data outputs from `analyze_phase2b_complete.py`.

**Total Placeholders**: ~50+  
**Data Source**: `paper/results_data/` (after running analysis)  
**Update Strategy**: Fill placeholders → regenerate manuscript

---

## 📊 Abstract Placeholders

| Placeholder | Location | Data Source | Description |
|-------------|----------|-------------|-------------|
| `[XX]` (unique species) | Abstract, Results line 37 | `summary_table.csv` → `total_unique_species` | Total unique molecular species across all 30 runs |
| `[XX]` (autocatalytic cycles) | Abstract, Results line 37 | `scenario_comparison.json` → `autocatalysis.total_cycles` | Total autocatalytic cycles detected |
| `[XX]` (hub molecules) | Abstract, Results line 37 | `table5_hub_molecules.csv` → count rows | Number of hub molecules (degree > threshold) |
| `[XX]%` (benchmark accuracy) | Abstract, Results line 37 | `benchmark_validation.json` → `overall_accuracy` | Benchmark reaction validation accuracy |

**Action**: After analysis, extract these 4 values and fill Abstract.

---

## 📈 Results Section 3.1: Molecular Diversity

### Paragraph 1: Global Statistics

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (total unique species) | `summary_table.csv` → `total_unique_species` | Total across all 30 runs |
| `[YY]` (max heavy atoms) | `summary_table.csv` → `max_heavy_atoms` | Largest molecule detected |
| `[XXX,XXX]` (steepest accumulation) | `figure3_data.json` → `steepest_phase_end_step` | Step where accumulation rate peaks |

### Paragraph 2: Scenario Comparison

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX ± YY]` (Miller-Urey species) | `miller_urey_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD across 10 runs |
| `[AA ± BB]` (Hydrothermal species) | `hydrothermal_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD across 10 runs |
| `[CC ± DD]` (Formamide species) | `formamide_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD across 8 runs |
| `[P]` (Kruskal-Wallis p-value) | `scenario_comparison.json` → `diversity.kruskal_wallis_p` | Statistical test p-value |

### Paragraph 3: Size Distributions

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X]` (Miller-Urey median) | `miller_urey_extended_analysis.json` → `size_distribution.median` | Median atoms per molecule |
| `[Y]-[Z]` (Miller-Urey IQR) | `miller_urey_extended_analysis.json` → `size_distribution.q25`-`q75` | Interquartile range |
| `[A]` (Formamide median) | `formamide_extended_analysis.json` → `size_distribution.median` | Median atoms per molecule |
| `[MAX]` (largest molecule atoms) | `summary_table.csv` → `max_heavy_atoms` | Maximum heavy atoms |
| `[scenario]` (where largest found) | `summary_table.csv` → `max_molecule_scenario` | Scenario containing largest molecule |

### Paragraph 4: Shannon Entropy

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X.X ± Y.Y]` (Miller-Urey entropy) | `miller_urey_extended_analysis.json` → `entropy.final_mean ± final_std` | Shannon entropy at 500K steps |
| `[A.A ± B.B]` (Hydrothermal entropy) | `hydrothermal_extended_analysis.json` → `entropy.final_mean ± final_std` | Shannon entropy at 500K steps |
| `[C.C ± D.D]` (Formamide entropy) | `formamide_extended_analysis.json` → `entropy.final_mean ± final_std` | Shannon entropy at 500K steps |

### Paragraph 5: Overlap Analysis

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]%` (scenario-specific) | `scenario_comparison.json` → `overlap.scenario_specific_percent` | Percentage unique to one scenario |
| `[YY]%` (shared across all) | `scenario_comparison.json` → `overlap.shared_all_percent` | Percentage in all 3 scenarios |
| `[list examples]` (core molecules) | `scenario_comparison.json` → `overlap.shared_molecules` | List of molecules in all scenarios |

---

## 🔗 Results Section 3.2: Reaction Network Topology

### Paragraph 1: Network Construction

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X.X ± Y.Y]` (avg path length) | `scenario_comparison.json` → `network.avg_path_length_mean ± std` | Average shortest path length |
| `[A.A ± B.B]` (clustering coeff) | `scenario_comparison.json` → `network.clustering_mean ± std` | Average clustering coefficient |

### Paragraph 2: Hub Molecules

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (formaldehyde degree) | `table5_hub_molecules.csv` → row where molecule="CH2O" → degree | Formaldehyde connectivity |
| `[YY]` (HCN degree) | `table5_hub_molecules.csv` → row where molecule="HCN" → degree | HCN connectivity |
| `[ZZ]` (ammonia degree) | `table5_hub_molecules.csv` → row where molecule="NH3" → degree | Ammonia connectivity |

### Paragraph 3: Topology Comparison

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X.X ± Y.Y]` (Miller-Urey γ) | `miller_urey_extended_analysis.json` → `network.power_law_exponent_mean ± std` | Power-law exponent |
| `[A.A ± B.B]` (Hydrothermal γ) | `hydrothermal_extended_analysis.json` → `network.power_law_exponent_mean ± std` | Power-law exponent |
| `[C.C ± D.D]` (Formamide γ) | `formamide_extended_analysis.json` → `network.power_law_exponent_mean ± std` | Power-law exponent |

### Paragraph 4: Network Metrics

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X.X ± Y.Y]` (Formamide avg degree) | `formamide_extended_analysis.json` → `network.avg_degree_mean ± std` | Average node degree |
| `[A.A ± B.B]` (Miller-Urey avg degree) | `miller_urey_extended_analysis.json` → `network.avg_degree_mean ± std` | Average node degree |
| `[L_f]` (Formamide path length) | `formamide_extended_analysis.json` → `network.avg_path_length_mean` | Average path length |
| `[L_m]` (Miller-Urey path length) | `miller_urey_extended_analysis.json` → `network.avg_path_length_mean` | Average path length |

---

## 🔄 Results Section 3.3: Autocatalytic Cycles

### Paragraph 1: Detection Method

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (total unique cycles) | `scenario_comparison.json` → `autocatalysis.total_unique_cycles` | Total distinct cycles across all runs |

### Paragraph 2: Frequency by Scenario

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[X ± Y]` (Miller-Urey cycles/run) | `miller_urey_extended_analysis.json` → `autocatalysis.cycles_per_run_mean ± std` | Mean cycles per run |
| `[A ± B]` (Hydrothermal cycles/run) | `hydrothermal_extended_analysis.json` → `autocatalysis.cycles_per_run_mean ± std` | Mean cycles per run |
| `[C ± D]` (Formamide cycles/run) | `formamide_extended_analysis.json` → `autocatalysis.cycles_per_run_mean ± std` | Mean cycles per run |
| `[P]` (p-value) | `scenario_comparison.json` → `autocatalysis.fisher_exact_p` | Statistical test p-value |

### Paragraph 3: Cycle Types

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]%` (simple loops) | `scenario_comparison.json` → `autocatalysis.simple_loops_percent` | 2-3 node cycles |
| `[YY]%` (medium loops) | `scenario_comparison.json` → `autocatalysis.medium_loops_percent` | 4-6 node cycles |
| `[ZZ]%` (complex networks) | `scenario_comparison.json` → `autocatalysis.complex_loops_percent` | >6 node cycles |
| `[N]` (direct autocatalysis) | `scenario_comparison.json` → `autocatalysis.direct_count` | A + B → 2A instances |

### Paragraph 4: Amplification Factors

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[min]` (min amplification) | `scenario_comparison.json` → `autocatalysis.amplification_min` | Minimum amplification factor |
| `[max]` (max amplification) | `scenario_comparison.json` → `autocatalysis.amplification_max` | Maximum amplification factor |
| `[X.X]` (median amplification) | `scenario_comparison.json` → `autocatalysis.amplification_median` | Median amplification |
| `[Y.Y]-[Z.Z]` (IQR) | `scenario_comparison.json` → `autocatalysis.amplification_q25`-`q75` | Interquartile range |
| `[molecule names]` (strongest) | `scenario_comparison.json` → `autocatalysis.top_amplifiers` | List of top amplifying molecules |
| `[scenario]` (where strongest) | `scenario_comparison.json` → `autocatalysis.top_amplifier_scenario` | Scenario with strongest amplification |
| `[XX]` (fold amplification) | `scenario_comparison.json` → `autocatalysis.max_amplification_factor` | Maximum fold increase |
| `[time]` (steps) | `scenario_comparison.json` → `autocatalysis.max_amplification_time` | Time to reach max amplification |

### Paragraph 5: Formose-like Cycles

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[N]` (formamide runs with formose) | `formamide_extended_analysis.json` → `autocatalysis.formose_like_count` | Number of runs with formose-like cycles |
| `[XX]` (glycolaldehyde amplification) | `formamide_extended_analysis.json` → `autocatalysis.glycolaldehyde_max_amplification` | Maximum glycolaldehyde fold increase |
| `[time]` (steps) | `formamide_extended_analysis.json` → `autocatalysis.glycolaldehyde_amplification_time` | Time to reach amplification |

---

## 🆕 Results Section 3.4: Novel Molecules

### Paragraph 1: Novel Molecule Definition

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (potentially novel) | `scenario_comparison.json` → `novelty.total_potentially_novel` | Total novel species detected |

### Paragraph 2: Novelty Distribution

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[YY]%` (novel percentage) | `scenario_comparison.json` → `novelty.novel_percent` | Percentage of total species that are novel |
| `[M]` (median mass) | `scenario_comparison.json` → `novelty.median_mass_amu` | Median molecular mass |
| `[min]-[max]` (mass range) | `scenario_comparison.json` → `novelty.mass_range` | Min-max mass range |
| `[XXX,XXX]` (novel detection time) | `scenario_comparison.json` → `novelty.median_detection_step` | Median step when novel molecules first appear |
| `[YYY,YYY]` (known detection time) | `scenario_comparison.json` → `novelty.known_median_detection_step` | Median step for known molecules |

### Paragraph 3: Top Novel Molecules

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[SMILES]` (most complex) | `table6_novel_molecules.csv` → row 1 → SMILES | SMILES string of most complex |
| `[mass]` (most complex mass) | `table6_novel_molecules.csv` → row 1 → mass_amu | Mass in amu |
| `[scenario]` (where found) | `table6_novel_molecules.csv` → row 1 → scenario | Scenario name |
| `[N]` (run number) | `table6_novel_molecules.csv` → row 1 → run_id | Run ID |
| `[XXX,XXX]` (detection step) | `table6_novel_molecules.csv` → row 1 → step_detected | Step when detected |
| `[X]` (heavy atoms) | `table6_novel_molecules.csv` → row 1 → heavy_atoms | Number of heavy atoms |
| `[feature description]` | `table6_novel_molecules.csv` → row 1 → features | Structural features |

### Paragraph 4: Formation Pathways

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[N]` (longest pathway steps) | `scenario_comparison.json` → `novelty.longest_pathway_steps` | Number of reaction steps |
| `[M]` (intermediate molecules) | `scenario_comparison.json` → `novelty.longest_pathway_intermediates` | Number of intermediates |
| `[molecules]` (common intermediates) | `scenario_comparison.json` → `novelty.common_intermediates` | List of frequent intermediates |

### Paragraph 5: Scenario Specificity

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]%` (formamide unique) | `formamide_extended_analysis.json` → `novelty.scenario_unique_percent` | Percentage unique to formamide |
| `[YY]%` (Miller-Urey unique) | `miller_urey_extended_analysis.json` → `novelty.scenario_unique_percent` | Percentage unique to Miller-Urey |
| `[ZZ]%` (Hydrothermal unique) | `hydrothermal_extended_analysis.json` → `novelty.scenario_unique_percent` | Percentage unique to hydrothermal |

---

## 💬 Discussion Section Placeholders

### Section 4.1: Emergent Complexity

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (species from <10 types) | `summary_table.csv` → `total_unique_species` | Total species detected |
| `[YY]` (fold increase) | `summary_table.csv` → `diversity_fold_increase` | Fold increase from starting molecules |

### Section 4.2: Scenario-Specific Chemistry

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX ± YY]` (Formamide species) | `formamide_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD |
| `[AA ± BB]` (Miller-Urey species) | `miller_urey_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD |
| `[CC ± DD]` (Hydrothermal species) | `hydrothermal_extended_analysis.json` → `diversity.mean_species ± std_species` | Mean ± SD |
| `[ZZ]%` (scenario-specific) | `scenario_comparison.json` → `overlap.scenario_specific_percent` | Percentage unique to one scenario |

### Section 4.3: Autocatalysis

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (cycles in all scenarios) | `scenario_comparison.json` → `autocatalysis.total_unique_cycles` | Total cycles detected |
| `[X ± Y]` (cycles per run) | `scenario_comparison.json` → `autocatalysis.cycles_per_run_mean ± std` | Mean cycles per run |

---

## 📊 Figure Captions Placeholders

### Figure 1: Thermodynamic Validation

| Placeholder | Data Source | Description |
|-------------|-------------|-------------|
| `[XX]` (χ² p-value) | `validation_results.json` → `maxwell_boltzmann.chi2_pvalue` | Chi-squared test p-value |
| `[XX]%` (entropy compliance) | `validation_results.json` → `entropy.compliance_percent` | Percentage of steps with ΔS ≥ 0 |

### Figure 2: Benchmark Validation

All data from `benchmark_validation.json` → individual reaction results.

---

## 🔧 How to Fill Placeholders

### Step 1: Run Analysis
```bash
python scripts/analyze_phase2b_complete.py \
    --input results/phase2b_additional \
    --output paper/results_data
```

### Step 2: Extract Values
Use this script to extract all placeholder values:
```bash
python scripts/extract_placeholder_values.py \
    --data paper/results_data \
    --output paper/placeholder_values.json
```

### Step 3: Fill Manuscript
```bash
python scripts/fill_manuscript_placeholders.py \
    --template paper/manuscript_draft.tex \
    --values paper/placeholder_values.json \
    --output paper/manuscript_filled.tex
```

---

## ✅ Checklist

After analysis completes:

- [ ] Run `analyze_phase2b_complete.py`
- [ ] Verify all JSON files in `paper/results_data/`
- [ ] Extract placeholder values
- [ ] Fill Abstract (4 placeholders)
- [ ] Fill Results 3.1 (15+ placeholders)
- [ ] Fill Results 3.2 (10+ placeholders)
- [ ] Fill Results 3.3 (15+ placeholders)
- [ ] Fill Results 3.4 (10+ placeholders)
- [ ] Fill Discussion (5+ placeholders)
- [ ] Fill Figure captions (2+ placeholders)
- [ ] Verify all `[XX]` replaced
- [ ] Check consistency across sections
- [ ] Generate final manuscript

---

**Status**: Ready for data  
**Next**: Run analysis → Extract values → Fill placeholders  
**ETA**: ~2 hours after analysis completes

