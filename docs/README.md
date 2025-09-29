# Live 2.0 - Open-Ended Prebiotic Simulator

A GPU-accelerated 2D simulation system capable of generating new particles, substances, and reactions through emergent chemistry.

## Features

- **Open Chemistry Mode**: Emergent particle interactions with continuous attributes
- **Preset Prebiotic Mode**: Validated chemical reactions for testing
- **Real-time Visualization**: WebSocket streaming with heatmaps and particle rendering
- **Novelty Detection**: Automatic discovery and cataloging of new molecular structures
- **Snapshot System**: Save and restore simulation states
- **GPU Acceleration**: Taichi-based simulation engine
- **Modern UI**: React frontend with TypeScript

## Architecture

```
live2/
├── backend/           # Python simulation engine
│   ├── sim/          # Core simulation modules
│   │   ├── core/     # Grid, particles, potentials, binding
│   │   └── io/       # Snapshot management
│   └── api/          # FastAPI server with WebSocket
├── frontend/         # React TypeScript frontend
│   ├── src/
│   │   ├── components/  # UI components
│   │   └── lib/         # API client and utilities
└── docker/           # Containerization
```

## Quick Start

### Prerequisites

- Python 3.9+
- Node.js 18+
- CUDA-capable GPU (recommended)

### Backend Setup

```bash
cd backend
pip install -r requirements.txt
python -m api.server
```

### Frontend Setup

```bash
cd frontend
npm install
npm run dev
```

### Access

- Frontend: http://localhost:3000
- API: http://localhost:8000
- API Docs: http://localhost:8000/docs

## Usage

1. **Start Simulation**: The frontend automatically creates and starts a simulation
2. **Monitor Metrics**: View particle count, novelty rate, and health score
3. **Explore Discoveries**: Browse novel substances and molecular structures
4. **Save Snapshots**: Capture interesting simulation states
5. **Adjust Parameters**: Modify simulation settings in real-time

## Simulation Modes

### Open Chemistry (Default)
- Continuous particle attributes (mass, charge)
- Emergent binding and unbinding
- Novel substance discovery
- Energy-driven mutations

### Preset Prebiotic
- Predefined chemical species
- Known reaction pathways
- Validation and testing
- Concentration field visualization

## API Endpoints

- `POST /simulation/create` - Create new simulation
- `GET /simulation/{id}/status` - Get simulation status
- `POST /simulation/{id}/start` - Start simulation
- `POST /simulation/{id}/pause` - Pause simulation
- `POST /simulation/{id}/stop` - Stop simulation
- `GET /simulation/{id}/novel-substances` - Get discoveries
- `POST /simulation/{id}/snapshot/save` - Save snapshot
- `POST /simulation/{id}/snapshot/load` - Load snapshot

## WebSocket Streaming

Real-time data stream at `/simulation/{id}/stream`:
- Particle positions and attributes
- Energy field data
- Bond information
- Cluster data
- Metrics updates

## Configuration

### Simulation Parameters

```python
{
  "grid_height": 256,
  "grid_width": 256,
  "max_particles": 10000,
  "dt": 0.01,
  "energy_decay": 0.95,
  "binding_threshold": 0.8,
  "novelty_window": 100
}
```

### Energy Sources

- Spatial energy impulses
- Continuous energy fields
- Configurable intensity and radius
- Temporal decay

## Development

### Backend Development

```bash
cd backend
python -m pytest tests/
python -m black sim/ api/
python -m isort sim/ api/
```

### Frontend Development

```bash
cd frontend
npm run lint
npm run build
npm run preview
```

### Testing

```bash
# Backend tests
cd backend
pytest tests/

# Frontend tests
cd frontend
npm test
```

## Performance

- **GPU Acceleration**: Taichi-based particle simulation
- **Binary Streaming**: msgpack encoding for efficient data transfer
- **Spatial Hashing**: O(1) neighbor finding
- **Adaptive Timestep**: Dynamic simulation speed
- **Memory Management**: Efficient particle lifecycle

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## License

MIT License - see LICENSE file for details

## Citation

If you use Live 2.0 in your research, please cite:

```bibtex
@software{live2_simulator,
  title={Live 2.0: Open-Ended Prebiotic Simulator},
  author={Your Name},
  year={2024},
  url={https://github.com/ProhunterPL/live2.0}
}
```

## Acknowledgments

- Taichi framework for GPU computing
- FastAPI for web API
- React for frontend
- NetworkX for graph analysis
