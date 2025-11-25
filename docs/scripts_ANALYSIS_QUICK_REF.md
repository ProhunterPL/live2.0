# Phase 2 & 3 Analysis - Quick Reference

**One-page reference for all analysis tools**

---

## 🎯 After Simulation Completes

### 1️⃣ Quick Check (1 minute)
```bash
python scripts/quick_analyze.py results/overnight_test
```
→ `analysis/summary.txt`

### 2️⃣ Build Network (5 minutes)
```bash
python scripts/reaction_network_analyzer.py results/overnight_test --export both
```
→ `analysis/reaction_network/*.{json,graphml,txt}`

### 3️⃣ Find Cycles (10 minutes)
```bash
python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json
```
→ `analysis/autocatalytic_cycles/*.{json,txt}`

### 4️⃣ Visualize (5 minutes)
```bash
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json \
    --interactive
```
→ `analysis/visualizations/*.{png,html}`

---

## 📊 Compare Multiple Scenarios

```bash
python scripts/compare_scenarios.py \
    results/phase2/miller_urey \
    results/phase2/hydrothermal \
    results/phase2/formamide
```
→ `analysis/scenario_comparison/*.{txt,json}`

---

## 🔄 Complete Pipeline (One Command)

### Single Simulation
```bash
# Analyze everything
python scripts/quick_analyze.py results/overnight_test --full && \
python scripts/reaction_network_analyzer.py results/overnight_test --export both && \
python scripts/autocatalytic_detector.py analysis/reaction_network/reaction_network.json && \
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json \
    --interactive
```

### Batch Analysis (All Scenarios)
```bash
# PowerShell
foreach ($scenario in @('miller_urey','hydrothermal','formamide')) {
    Write-Host "Analyzing $scenario..."
    python scripts/reaction_network_analyzer.py results/phase2/$scenario/* --merge --output analysis/$scenario/network
    python scripts/autocatalytic_detector.py analysis/$scenario/network/reaction_network.json --output analysis/$scenario/cycles
    python scripts/network_visualizer.py analysis/$scenario/network/reaction_network.json --cycles analysis/$scenario/cycles/autocatalytic_cycles.json --output analysis/$scenario/viz --interactive
}
python scripts/compare_scenarios.py results/phase2/* --output analysis/comparison
```

---

## 📝 Key Output Files

| File | Description | Tool |
|------|-------------|------|
| `molecules.txt` | All molecules found | quick_analyze |
| `pubchem_matches.txt` | PubChem matches | quick_analyze --full |
| `reaction_network.json` | Network structure | reaction_network_analyzer |
| `reaction_network.graphml` | For Gephi/Cytoscape | reaction_network_analyzer |
| `network_analysis.txt` | Network statistics | reaction_network_analyzer |
| `autocatalytic_cycles.json` | All cycles | autocatalytic_detector |
| `autocatalytic_report.txt` | Cycle summary | autocatalytic_detector |
| `degree_distribution.png` | Network metrics | network_visualizer |
| `network_topology.png` | Network graph | network_visualizer |
| `autocatalytic_cycles.png` | Cycle diagrams | network_visualizer |
| `network_interactive.html` | Interactive view | network_visualizer --interactive |
| `scenario_comparison.txt` | Scenario stats | compare_scenarios |

---

## 🎨 For Publication Figures

### High-Quality Network Graph
```bash
python scripts/network_visualizer.py network.json --max-nodes 50
```
→ 300 DPI PNG files

### Custom Layout in Gephi
```bash
python scripts/reaction_network_analyzer.py results/overnight_test --export graphml
# Then open reaction_network.graphml in Gephi
```

### Cycle Analysis Summary
```bash
python scripts/autocatalytic_detector.py network.json --detailed
```
→ `autocatalytic_report.txt` (formatted for paper)

---

## 🐛 Common Issues

**"No molecules found"**
```bash
# Check simulation output
ls results/overnight_test/snapshots/
Get-Content results/overnight_test/stderr.log -Tail 20
```

**"Network file not found"**
```bash
# Run network analysis first
python scripts/reaction_network_analyzer.py results/overnight_test
```

**"matplotlib not installed"**
```bash
pip install matplotlib networkx
```

---

## 📚 Full Documentation

- **Detailed Guide**: `docs/PHASE3_ANALYSIS_GUIDE.md`
- **Phase 2 Usage**: `docs/PHASE2_USAGE.md`
- **Roadmap**: `docs/VALIDATION_ROADMAP.md`

---

## ⚡ Performance Tips

1. **Use --recursive for batch**
   ```bash
   python scripts/quick_analyze.py results/ --recursive
   ```

2. **Limit network size**
   ```bash
   python scripts/network_visualizer.py network.json --max-nodes 100
   ```

3. **Run in order**
   - Always: quick_analyze → network_analyzer → autocatalytic_detector → visualizer

---

**Need help?** Check `docs/PHASE3_ANALYSIS_GUIDE.md` for detailed usage.

