# Phase 3 Analysis Tools - Usage Guide

**Version**: 1.0  
**Date**: October 13, 2025  
**Status**: Production Ready

---

## 📋 Overview

This guide covers the Phase 3 analysis infrastructure for processing Phase 2 simulation results, including:

1. **Quick Analysis** - Fast molecule extraction and summarization
2. **Scenario Comparison** - Compare multiple experimental conditions
3. **Reaction Network Analysis** - Build and analyze chemical reaction networks
4. **Autocatalytic Cycle Detection** - Identify self-amplifying reaction loops
5. **Network Visualization** - Create publication-quality figures

---

## 🚀 Quick Start

### 1. Analyze Single Simulation Run

```bash
# Extract molecules and generate summary
python scripts/quick_analyze.py results/overnight_test_2025-10-13_18-17-09

# With PubChem matching (slower, more detailed)
python scripts/quick_analyze.py results/overnight_test_2025-10-13_18-17-09 --full
```

**Output**:
- `analysis/molecules.txt` - List of all molecules found
- `analysis/pubchem_matches.txt` - PubChem database matches (if --full)
- `analysis/summary.txt` - Quick summary statistics

---

### 2. Compare Multiple Scenarios

```bash
# Compare all 3 scenarios
python scripts/compare_scenarios.py \
    results/phase2/miller_urey \
    results/phase2/hydrothermal \
    results/phase2/formamide

# Output to specific directory
python scripts/compare_scenarios.py \
    results/phase2/* \
    --output analysis/scenario_comparison
```

**Output**:
- `scenario_comparison.txt` - Human-readable comparison
- `comparison_data.json` - Structured data for further analysis

**Metrics Compared**:
- Unique molecules produced
- Average molecule size
- Total reactions
- Reaction rates per step
- Diversity indices

---

### 3. Reaction Network Analysis

```bash
# Analyze single simulation
python scripts/reaction_network_analyzer.py results/overnight_test

# Merge multiple runs into one network
python scripts/reaction_network_analyzer.py \
    results/phase2/miller_urey/* \
    --merge \
    --output analysis/miller_urey_network

# Export formats
python scripts/reaction_network_analyzer.py \
    results/overnight_test \
    --export both  # json, graphml, or both
```

**Output**:
- `reaction_network.json` - Network data structure
- `reaction_network.graphml` - GraphML for Gephi/Cytoscape
- `network_analysis.txt` - Statistics and key molecules

**Network Metrics**:
- Number of molecules and reactions
- Degree distributions (in/out/total)
- Source molecules (no incoming reactions)
- Sink molecules (no outgoing reactions)
- Hub molecules (high connectivity)

---

### 4. Autocatalytic Cycle Detection

```bash
# Detect autocatalytic cycles
python scripts/autocatalytic_detector.py \
    analysis/reaction_network/reaction_network.json

# Adjust search parameters
python scripts/autocatalytic_detector.py \
    analysis/reaction_network/reaction_network.json \
    --max-cycle-length 15 \
    --detailed

# Output to specific directory
python scripts/autocatalytic_detector.py \
    analysis/reaction_network/reaction_network.json \
    --output analysis/autocatalytic_cycles
```

**Output**:
- `autocatalytic_cycles.json` - All detected cycles
- `autocatalytic_report.txt` - Human-readable summary

**Cycle Types Detected**:

1. **Direct Autocatalysis**
   - A + B → 2A (molecule catalyzes its own production)
   - Example: A + B → A + A

2. **Indirect Autocatalysis**
   - A → B → C → A (cycle with amplification)
   - Example: A + X → B, B + Y → C + C, C → A + A

3. **Hypercycles**
   - A catalyzes B, B catalyzes C, C catalyzes A
   - Mutual catalytic support

4. **RAF Sets** (Reflexively Autocatalytic Food-generated)
   - Complex autocatalytic networks
   - Self-sustaining reaction sets

---

### 5. Network Visualization

```bash
# Create all visualizations
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json

# Include autocatalytic cycles
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --cycles analysis/autocatalytic_cycles/autocatalytic_cycles.json

# Generate interactive HTML
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --interactive

# Control graph size
python scripts/network_visualizer.py \
    analysis/reaction_network/reaction_network.json \
    --max-nodes 50
```

**Output**:
- `degree_distribution.png` - In/out/total degree histograms
- `network_topology.png` - Network graph layout
- `autocatalytic_cycles.png` - Cycle diagrams (if cycles provided)
- `statistics_summary.png` - Key metrics summary
- `network_interactive.html` - Interactive browser visualization

---

## 📊 Complete Analysis Pipeline

### Full Workflow (After Simulation Completes)

```bash
# 1. Quick analysis - extract molecules
python scripts/quick_analyze.py results/overnight_test --full

# 2. Build reaction network
python scripts/reaction_network_analyzer.py \
    results/overnight_test \
    --output analysis/network \
    --export both

# 3. Detect autocatalytic cycles
python scripts/autocatalytic_detector.py \
    analysis/network/reaction_network.json \
    --output analysis/cycles

# 4. Create visualizations
python scripts/network_visualizer.py \
    analysis/network/reaction_network.json \
    --cycles analysis/cycles/autocatalytic_cycles.json \
    --interactive
```

### Batch Analysis (Multiple Scenarios)

```bash
# Analyze all scenarios
for scenario in miller_urey hydrothermal formamide; do
    echo "Analyzing $scenario..."
    
    # Quick analysis
    python scripts/quick_analyze.py \
        results/phase2/$scenario \
        --recursive
    
    # Network analysis
    python scripts/reaction_network_analyzer.py \
        results/phase2/$scenario/* \
        --merge \
        --output analysis/$scenario/network
    
    # Cycle detection
    python scripts/autocatalytic_detector.py \
        analysis/$scenario/network/reaction_network.json \
        --output analysis/$scenario/cycles
    
    # Visualization
    python scripts/network_visualizer.py \
        analysis/$scenario/network/reaction_network.json \
        --cycles analysis/$scenario/cycles/autocatalytic_cycles.json \
        --output analysis/$scenario/visualizations \
        --interactive
done

# Compare all scenarios
python scripts/compare_scenarios.py \
    results/phase2/miller_urey \
    results/phase2/hydrothermal \
    results/phase2/formamide \
    --output analysis/comparison
```

---

## 📈 Output Files Reference

### Quick Analysis
```
analysis/
├── molecules.txt          # All molecules found
├── pubchem_matches.txt    # PubChem matches (with --full)
└── summary.txt            # Summary statistics
```

### Reaction Network
```
analysis/reaction_network/
├── reaction_network.json     # Network structure
├── reaction_network.graphml  # GraphML format
└── network_analysis.txt      # Statistics report
```

### Autocatalytic Cycles
```
analysis/autocatalytic_cycles/
├── autocatalytic_cycles.json  # All cycles detected
└── autocatalytic_report.txt   # Human-readable summary
```

### Visualizations
```
analysis/visualizations/
├── degree_distribution.png    # Degree histograms
├── network_topology.png       # Network graph
├── autocatalytic_cycles.png   # Cycle diagrams
├── statistics_summary.png     # Key metrics
└── network_interactive.html   # Interactive HTML
```

### Scenario Comparison
```
analysis/scenario_comparison/
├── scenario_comparison.txt    # Comparison report
└── comparison_data.json       # Structured data
```

---

## 🔬 Advanced Usage

### Custom Network Analysis

```python
from scripts.reaction_network_analyzer import ReactionNetworkAnalyzer
from pathlib import Path

# Load and analyze
analyzer = ReactionNetworkAnalyzer(
    [Path("results/run1"), Path("results/run2")],
    Path("analysis/custom")
)

analyzer.load_results()
stats = analyzer.analyze()

# Access network
network = analyzer.network
sources = network.get_sources()
hubs = network.get_hubs(threshold=10)

# Custom export
analyzer.export_json()
```

### Custom Cycle Detection

```python
from scripts.autocatalytic_detector import AutocatalyticDetector
from pathlib import Path

# Load network
detector = AutocatalyticDetector(
    Path("analysis/network/reaction_network.json"),
    Path("analysis/custom_cycles")
)

detector.load_network()

# Detect specific types
direct = detector.detect_direct_autocatalysis()
indirect = detector.detect_indirect_autocatalysis(max_length=15)
hypercycles = detector.detect_hypercycles()

# Filter by criteria
large_cycles = [c for c in indirect if c.size() >= 5]
amplifying = [c for c in direct if c.amplification_factor > 1.5]
```

---

## 🎯 Tips & Best Practices

### Performance Optimization

1. **Use `--recursive` for batch analysis**
   ```bash
   python scripts/quick_analyze.py results/ --recursive
   ```

2. **Limit network size for visualization**
   ```bash
   python scripts/network_visualizer.py network.json --max-nodes 50
   ```

3. **Run network analysis before cycle detection**
   - Cycle detection requires network JSON as input
   - Always run network analysis first

### For Publication Figures

1. **High-resolution exports**
   - All PNG files are 300 DPI (publication quality)
   - Use GraphML for Gephi/Cytoscape custom layouts

2. **Interactive HTML for exploration**
   - Share with collaborators via browser
   - No software installation required

3. **Combine multiple visualizations**
   - Use `--interactive` for overview
   - Create custom layouts in Gephi for paper figures

---

## 🐛 Troubleshooting

### "No molecules found"
- Check that simulation completed successfully
- Verify `results.json` or snapshot files exist
- Run with `--recursive` if analyzing directory

### "Network file not found"
- Run `reaction_network_analyzer.py` before `autocatalytic_detector.py`
- Check output paths match between tools

### "matplotlib not installed"
- Install: `pip install matplotlib networkx`
- Required for visualization tools

### Large networks (>1000 molecules)
- Use `--max-nodes` to limit visualization
- Consider splitting into sub-networks
- Use GraphML export for external tools

---

## 📚 See Also

- [VALIDATION_ROADMAP.md](VALIDATION_ROADMAP.md) - Overall project status
- [PHASE2_USAGE.md](PHASE2_USAGE.md) - Running simulations
- [PHASE2_TECHNICAL.md](PHASE2_TECHNICAL.md) - Technical details

---

**Questions?** Open an issue or see documentation in `docs/`

