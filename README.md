# LIVE 2.0 - Open-Ended Prebiotic Chemistry Simulator 🧬

<div align="center">

**GPU-accelerated emergent chemistry simulation with real-time molecular discovery**

[![CI Tests](https://github.com/klawi/live2.0/workflows/CI%20Tests/badge.svg)](https://github.com/klawi/live2.0/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![React](https://img.shields.io/badge/react-18+-61DAFB.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/fastapi-0.104+-009688.svg)](https://fastapi.tiangolo.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[Features](#-features) • [Quick Start](#-quick-start) • [Documentation](#-documentation) • [Architecture](#-architecture)

</div>

---

## 🌟 Features

### Core Simulation
- 🔬 **Emergent Chemistry** - Novel molecular structures emerge from particle interactions
- ⚡ **GPU Acceleration** - Taichi-powered simulation engine for real-time performance
- 🎯 **Novelty Detection** - Automatic discovery and cataloging of unique substances
- 🔄 **Open Chemistry Mode** - Continuous particle attributes with emergent binding
- 📊 **Real-time Metrics** - Live monitoring of novelty rate, particle count, and system health

### New: PubChem Matcher 🧪
- 🔍 **One-Click Matching** - Compare discovered clusters with real molecules from PubChem
- 🧬 **SMILES Conversion** - Automatic topology-to-SMILES translation
- 📁 **Multi-Format Export** - Generate .mol and .xyz files with 3D coordinates
- 🖼️ **Visual Comparison** - Side-by-side panels of LIVE vs PubChem structures
- ⚗️ **Chemical Integration** - Use exported files in PyMOL, Avogadro, ChemDraw, and more

### Web Interface
- 🎨 **Modern UI** - React + TypeScript with real-time WebSocket streaming
- 📈 **Live Visualization** - Particle rendering, energy fields, and bond networks
- 💾 **Snapshot System** - Save and restore interesting simulation states
- 🎛️ **Interactive Controls** - Adjust parameters and inject energy in real-time

## 🚀 Quick Start

### Prerequisites
```bash
Python 3.9+
Node.js 18+
CUDA-capable GPU (recommended)
```

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/live2.0.git
cd live2.0
```

2. **Install dependencies**
```bash
# Backend dependencies
pip install -r requirements.txt

# Frontend dependencies
cd frontend
npm install
cd ..
```

### Running

**Start backend:**
```bash
cd backend
python -m api.server
```

**Start frontend:**
```bash
cd frontend
npm run dev
```

**Access:**
- 🌐 Frontend: http://localhost:3000
- 🔌 API: http://localhost:8001
- 📚 API Docs: http://localhost:8001/docs

## ⚡ Performance Optimization

### Check Your Backend

Want to know if your simulation is running on GPU or CPU?

```bash
python scripts/check_current_backend.py
```

This will show:
- Available backends (CUDA, Vulkan, CPU)
- Current active backend
- Performance recommendations

### GPU vs CPU Benchmark

To test which backend is fastest for your system:

**Windows:**
```powershell
.\run_benchmark.ps1
```

**Linux/Mac:**
```bash
python tests/benchmark_gpu_vs_cpu.py
```

The benchmark will:
- ✅ Test GPU (CUDA) if available
- ✅ Test CPU with different thread counts
- ✅ Measure both simulation and visualization performance
- ✅ Provide clear recommendations

**Expected results:**
- **GPU (CUDA)**: Best for visualization and large particle counts (10-50x faster)
- **CPU (many threads)**: Can be faster for smaller simulations with complex chemistry
- Results depend on: particle count, GPU model, CPU core count

**Note from Phase 2B testing:** On AWS instances with 96 vCPUs, CPU was faster than GPU for our specific workload (chemistry-heavy, not visualization-heavy). Your results may vary!

### 🚀 NEW: Hybrid GPU+CPU Mode

Get the best of both worlds - GPU for physics, CPU for chemistry!

```powershell
.\run_hybrid_test.ps1
```

**Hybrid mode uses:**
- 🎮 **GPU (CUDA)**: Fast particle physics in main thread
- 🧮 **CPU (Python)**: Complex chemistry analysis in background thread
- 🔄 **Async communication**: No blocking, smooth real-time experience

**Benefits:**
- ✅ GPU runs at full speed (chemistry doesn't block)
- ✅ CPU handles complex logic (graphs, branching) - its strength!
- ✅ 5-10x faster than Pure GPU for chemistry-heavy simulations
- ✅ Smooth 30+ FPS visualization with full analysis

**When to use:**
- You have NVIDIA GPU
- Visualization shows slow "Bonds/Clusters" timing (>200ms)
- You want real-time UI with full chemistry analysis

See [HYBRID_GPU_CPU_GUIDE.md](docs/HYBRID_GPU_CPU_GUIDE.md) for usage instructions.

## 📖 Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- **[README.md](docs/README.md)** - Detailed project overview
- **[README_MATCHER.md](docs/README_MATCHER.md)** - PubChem matcher guide
- **[QUICK_START.md](docs/QUICK_START.md)** - Getting started guide
- **[INSTALLATION.md](docs/INSTALLATION.md)** - Installation instructions
- **[protocol.md](docs/protocol.md)** - API protocol specification
- **[SCIENTIFIC_OVERVIEW.md](docs/SCIENTIFIC_OVERVIEW.md)** - Scientific background

## 🏗️ Architecture

```
live2.0/
├── requirements.txt     # All Python dependencies
├── backend/             # Python simulation engine
│   ├── sim/            # Core simulation modules
│   │   ├── core/       # Particle system, binding, graphs
│   │   └── io/         # Snapshot management
│   ├── api/            # FastAPI server + WebSocket
│   └── tests/          # Backend tests
│
├── logs/                # Application logs
│
├── matcher/             # PubChem matching system
│   ├── chem.py          # RDKit + PubChem integration
│   ├── compose.py       # Image composition
│   ├── matcher.py       # CLI tool
│   └── watcher.py       # File system watcher
│
├── frontend/            # React TypeScript UI
│   ├── src/
│   │   ├── components/  # React components
│   │   └── lib/         # API client
│   └── package.json
│
├── tests/               # Test suite
├── docs/                # Documentation
├── matches/             # PubChem match results
└── docker/              # Docker configuration
```

## 🧪 Usage Example

### 1. Start Simulation
```typescript
// Frontend automatically creates simulation on load
// Or use API directly:
POST /simulation/create
```

### 2. Monitor Discoveries
```typescript
// View novel substances in Novelty Panel
GET /simulation/{id}/novel-substances?count=20
```

### 3. Match to PubChem
```typescript
// Click download icon on any substance
POST /simulation/{id}/substance/{substance_id}/match

// Returns SMILES, PubChem match, and file paths
{
  "smiles": "CCCCCC",
  "pubchem_match": {
    "cid": 8058,
    "name": "hexane",
    "formula": "C6H14"
  },
  "files": {
    "mol": "matches/cluster_XXX.mol",
    "xyz": "matches/cluster_XXX.xyz",
    "panel": "matches/cluster_XXX_match.png"
  }
}
```

### 4. Use Exported Files
```bash
# Visualize in PyMOL
pymol matches/cluster_XXX.xyz

# Convert to other formats
obabel matches/cluster_XXX.mol -O output.pdb

# Load in RDKit
from rdkit import Chem
mol = Chem.MolFromMolFile("matches/cluster_XXX.mol")
```

## 🎯 Key Concepts

### Open Chemistry Mode
- Particles with continuous mass and charge attributes
- Emergent bond formation based on energy and distance
- Novel substance discovery through graph isomorphism detection
- Energy-driven mutations and structural evolution

### Novelty Detection
- Automatic cataloging of unique molecular graphs
- Real-time complexity and density metrics
- Deduplication using Weisfeiler-Lehman graph hashing
- Discovery timeline with timestamps

### PubChem Integration
- Automatic SMILES generation from cluster topology
- Similarity search (80% threshold) in PubChem database
- Chemical file format export (MOL, XYZ)
- 3D coordinate generation with ETKDG + UFF optimization

## 🛠️ Development

### Running Tests
```bash
# Backend tests
cd backend
pytest tests/

# Frontend tests
cd frontend
npm test
```

### Code Quality
```bash
# Backend formatting
black backend/sim backend/api
isort backend/sim backend/api

# Frontend linting
cd frontend
npm run lint
```

## 📊 Performance

- **GPU Acceleration**: Taichi-based particle physics (~10,000 particles)
- **Real-time Streaming**: msgpack binary protocol over WebSocket
- **Spatial Hashing**: O(1) neighbor finding for collision detection
- **Adaptive Rendering**: Dynamic frame rate adjustment
- **Optimized Novelty Detection**: Caching and incremental updates

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- **[Taichi](https://www.taichi-lang.org/)** - GPU-accelerated computation framework
- **[FastAPI](https://fastapi.tiangolo.com/)** - Modern web framework for APIs
- **[RDKit](https://www.rdkit.org/)** - Cheminformatics and machine learning toolkit
- **[PubChem](https://pubchem.ncbi.nlm.nih.gov/)** - Chemical substance database
- **[React](https://reactjs.org/)** - UI framework
- **[NetworkX](https://networkx.org/)** - Graph analysis library

## 📚 Citation

If you use LIVE 2.0 in your research, please cite:

```bibtex
@software{live2_simulator,
  title={LIVE 2.0: Open-Ended Prebiotic Chemistry Simulator},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/live2.0}
}
```

## 🔬 Scientific Background

LIVE 2.0 explores emergent chemistry through:
- **Bottom-up self-organization** of molecular structures
- **Novelty-driven evolution** of chemical complexity
- **Open-ended discovery** of reaction pathways
- **Prebiotic chemistry** simulation scenarios

For detailed scientific context, see [SCIENTIFIC_OVERVIEW.md](docs/SCIENTIFIC_OVERVIEW.md).

---

<div align="center">

**[Report Bug](https://github.com/yourusername/live2.0/issues)** • **[Request Feature](https://github.com/yourusername/live2.0/issues)** • **[Discussions](https://github.com/yourusername/live2.0/discussions)**

Made with ❤️ by the LIVE 2.0 Team

</div>

