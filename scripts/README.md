# Scripts Directory

**Complete toolset for Live 2.0 simulation and analysis**

---

## 📁 Directory Structure

### Phase 2: Simulation Execution

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_phase2_full.py` | Main Phase 2 simulation runner | `python scripts/run_phase2_full.py --config configs/phase2_miller_urey.yaml` |
| `phase2_master.py` | Orchestrates multiple simulations | `python scripts/phase2_master.py --mode full --scenarios all` |
| `start_overnight_test.ps1` | Launch overnight Miller-Urey test | `.\scripts\start_overnight_test.ps1` |
| `watch_and_analyze.ps1` | Monitor simulation + auto-analyze | `.\scripts\watch_and_analyze.ps1 -ResultDir results/overnight_test` |

---

### Phase 3: Analysis Tools

| Script | Purpose | Output |
|--------|---------|--------|
| `quick_analyze.py` | Fast molecule extraction | `analysis/summary.txt` |
| `compare_scenarios.py` | Compare multiple scenarios | `analysis/scenario_comparison.txt` |
| `reaction_network_analyzer.py` | Build reaction networks | `analysis/reaction_network/` |
| `autocatalytic_detector.py` | Find autocatalytic cycles | `analysis/autocatalytic_cycles/` |
| `network_visualizer.py` | Create visualizations | `analysis/visualizations/` |

---

### Utilities

| Script | Purpose |
|--------|---------|
| `validate_parameters.py` | Validate physics parameters database |
| `analyze_diagnostics.py` | Analyze diagnostic outputs |

---

## 🚀 Quick Start Workflows

### 1. Run Single Simulation

```bash
# Miller-Urey test (100K steps, ~2-4 hours)
.\scripts\start_overnight_test.ps1
```

### 2. Monitor Running Simulation

```bash
# Auto-analyze when complete
.\scripts\watch_and_analyze.ps1 -ResultDir "results/overnight_test_2025-10-13_18-17-09"
```

### 3. Analyze Completed Simulation

```bash
# Quick analysis
python scripts/quick_analyze.py results/overnight_test --full

# Full pipeline
python scripts/quick_analyze.py results/overnight_test && \
python scripts/reaction_network_analyzer.py results/overnight_test --export both && \
python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json && \
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json \
    --interactive
```

### 4. Run Full Phase 2 Pipeline (30 simulations)

```bash
# All 3 scenarios × 10 runs each
python scripts/phase2_master.py --mode full --scenarios all
```

### 5. Compare Scenarios

```bash
# After multiple runs complete
python scripts/compare_scenarios.py \
    results/phase2/miller_urey \
    results/phase2/hydrothermal \
    results/phase2/formamide
```

---

## 📊 Typical Analysis Workflow

```
1. Run simulation
   └─> scripts/run_phase2_full.py
   └─> results/simulation_xxx/

2. Quick check
   └─> scripts/quick_analyze.py
   └─> analysis/summary.txt

3. Build network
   └─> scripts/reaction_network_analyzer.py
   └─> analysis/reaction_network/*.{json,graphml}

4. Find cycles
   └─> scripts/autocatalytic_detector.py
   └─> analysis/autocatalytic_cycles/*.json

5. Visualize
   └─> scripts/network_visualizer.py
   └─> analysis/visualizations/*.{png,html}

6. Compare scenarios (if multiple runs)
   └─> scripts/compare_scenarios.py
   └─> analysis/scenario_comparison/*.txt
```

---

## 🎯 Command Reference

### Quick Analysis
```bash
# Basic
python scripts/quick_analyze.py results/run_name

# With PubChem matching
python scripts/quick_analyze.py results/run_name --full

# Batch mode
python scripts/quick_analyze.py results/ --recursive
```

### Network Analysis
```bash
# Single run
python scripts/reaction_network_analyzer.py results/run_name

# Merge multiple runs
python scripts/reaction_network_analyzer.py results/run1 results/run2 --merge

# Export GraphML for Gephi
python scripts/reaction_network_analyzer.py results/run_name --export graphml
```

### Cycle Detection
```bash
# Standard
python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json

# Detailed search
python scripts/autocatalytic_detector.py network.json --max-cycle-length 15 --detailed
```

### Visualization
```bash
# All visualizations
python scripts/network_visualizer.py analysis/reaction_network/reaction_network.json

# With cycles
python scripts/network_visualizer.py network.json --cycles cycles.json

# Interactive HTML
python scripts/network_visualizer.py network.json --interactive

# Limit network size
python scripts/network_visualizer.py network.json --max-nodes 50
```

### Scenario Comparison
```bash
# Compare 2+ scenarios
python scripts/compare_scenarios.py scenario1/ scenario2/ scenario3/

# Custom output
python scripts/compare_scenarios.py scenario1/ scenario2/ --output analysis/comparison
```

---

## 📚 Documentation

- **Quick Reference**: [ANALYSIS_QUICK_REF.md](ANALYSIS_QUICK_REF.md)
- **Full Guide**: [../docs/PHASE3_ANALYSIS_GUIDE.md](../docs/PHASE3_ANALYSIS_GUIDE.md)
- **Phase 2 Usage**: [../docs/PHASE2_USAGE.md](../docs/PHASE2_USAGE.md)
- **Roadmap**: [../docs/VALIDATION_ROADMAP.md](../docs/VALIDATION_ROADMAP.md)

---

## 🔧 Requirements

### Python Dependencies
```bash
pip install -r requirements.txt
```

### Additional (for visualization)
```bash
pip install matplotlib networkx
```

---

## 💡 Tips

1. **Always run network analysis before cycle detection**
   - Cycle detector needs network JSON as input

2. **Use `--recursive` for batch processing**
   - Processes all subdirectories automatically

3. **Export to GraphML for custom layouts**
   - Open in Gephi/Cytoscape for publication figures

4. **Use interactive HTML for exploration**
   - Share with collaborators via browser

5. **Monitor long simulations with watch_and_analyze.ps1**
   - Automatically runs analysis when complete

---

## 🐛 Troubleshooting

### Script not found
```bash
# Make sure you're in project root
cd C:\Users\user\Desktop\live2.0
```

### Module import errors
```bash
# Reinstall dependencies
pip install -r requirements.txt
```

### matplotlib/networkx not found
```bash
# Install visualization dependencies
pip install matplotlib networkx
```

### No results found
```bash
# Check simulation completed
ls results/run_name/
Get-Content results/run_name/stderr.log -Tail 20
```

---

## 🎓 Examples

### Example 1: Quick overnight test analysis
```powershell
# After test completes
python scripts/quick_analyze.py results/overnight_test_2025-10-13_18-17-09 --full
```

### Example 2: Full network analysis
```powershell
python scripts/reaction_network_analyzer.py results/overnight_test --export both
python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json
python scripts/network_visualizer.py `
    analysis/reaction_network/reaction_network.json `
    --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json `
    --interactive
```

### Example 3: Batch scenario analysis
```powershell
# Analyze all Miller-Urey runs
python scripts/reaction_network_analyzer.py results/phase2/miller_urey/* --merge --output analysis/miller_urey_merged
```

---

**Questions?** See full documentation in `docs/` or open an issue on GitHub.

