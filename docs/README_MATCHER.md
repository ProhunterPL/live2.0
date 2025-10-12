# LIVE 2.0 → PubChem Matcher

Automatic cluster-to-molecule matching system for LIVE 2.0 simulation.

## Overview

This system automatically compares LIVE 2.0 simulation clusters with real molecules from the PubChem database. 

**Fully integrated with backend** - no separate processes or manual file transfers needed!

When you click the download icon on a cluster, it:

1. **Converts** the cluster topology to a SMILES molecular representation
2. **Searches** PubChem for the most similar real molecule
3. **Generates** a side-by-side comparison panel
4. **Exports** to standard chemical formats (.mol, .xyz)
5. **Provides** detailed chemical information about the match

### Key Features

✅ **One-click operation** - just click the icon in the UI  
✅ **No file management** - everything handled by backend  
✅ **Automatic format conversion** - MOL and XYZ files generated  
✅ **3D coordinates** - ETKDG + UFF optimization  
✅ **PubChem integration** - automatic similarity search  
✅ **Visual comparison** - side-by-side panels  
✅ **Optional standalone** - CLI and watcher still available for batch processing

## Quick Start

### 1. Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

The requirements include all necessary dependencies for the matcher:
- `rdkit-pypi>=2022.9.5` - Chemistry toolkit
- `requests>=2.32.0` - HTTP client for PubChem API
- `Pillow>=10.4.0` - Image processing
- `watchdog>=4.0.0` - File system monitoring (optional, for standalone watcher)

### 2. Start the Backend

Make sure your LIVE 2.0 backend is running:

```bash
python backend/api/server.py
```

The matcher functionality is **automatically integrated** into the backend API - no separate process needed!

### 3. Match Clusters from Frontend

1. Start your LIVE 2.0 simulation
2. Navigate to the "Novelty Detection" panel
3. Click the download icon (📥) on any novel substance
4. Wait a few seconds - the backend will:
   - Convert the cluster to SMILES
   - Query PubChem for similar molecules
   - Generate comparison panel and chemical format files
   - Display the match results

### 4. View Results

Check the `matches/` folder for:
- `cluster_TIMESTAMP_match.png` - Side-by-side comparison panel
- `cluster_TIMESTAMP_match.json` - Detailed match information
- `cluster_TIMESTAMP.mol` - LIVE cluster in MDL Molfile format
- `cluster_TIMESTAMP.xyz` - LIVE cluster in XYZ format with 3D coordinates
- `cluster_TIMESTAMP_pubchem.mol` - PubChem match in MOL format (if found)
- `cluster_TIMESTAMP_pubchem.xyz` - PubChem match in XYZ format (if found)

## File Structure

```
live2.0/
├── requirements.txt       # All Python dependencies (includes matcher)
├── backend/               # Simulation engine
├── matches/               # Output folder (generated panels)
│   ├── cluster_XXX_match.png
│   ├── cluster_XXX_match.json
│   ├── cluster_XXX.mol
│   ├── cluster_XXX.xyz
│   ├── cluster_XXX_pubchem.mol
│   └── cluster_XXX_pubchem.xyz
├── matcher/               # Matcher system code
│   ├── __init__.py
│   ├── chem.py           # RDKit & PubChem integration
│   ├── compose.py        # Panel composition
│   ├── matcher.py        # CLI tool (optional)
│   └── watcher.py        # File watcher (optional)
└── tests/                 # All test files
```

## Data Format

### Input JSON (cluster_XXX.json)

The cluster JSON file must contain topology information:

```json
{
  "id": "cluster_2025-10-11T12-34-56",
  "nodes": [
    {
      "id": 0,
      "label": "A",
      "mass": 1.0,
      "charge": [0.0, 0.0, 0.0],
      "energy": 59.4
    },
    {
      "id": 1,
      "label": "A",
      "mass": 1.0,
      "charge": [0.0, 0.0, 0.0],
      "energy": 56.2
    }
  ],
  "bonds": [
    {"a": 0, "b": 1, "order": 1},
    {"a": 1, "b": 2, "order": 1}
  ],
  "metadata": {
    "size": 8,
    "bonds": 10,
    "density": 0.357,
    "avg_mass": 1.0,
    "complexity": 12.5
  }
}
```

### Output JSON (cluster_XXX_match.json)

The result contains SMILES and PubChem match information:

```json
{
  "cluster_id": "cluster_2025-10-11T12-34-56",
  "smiles_live": "CCCCCC",
  "pubchem_top": {
    "cid": 8058,
    "name": "hexane",
    "formula": "C6H14",
    "mw": 86.18,
    "smiles": "CCCCCC",
    "inchikey": "VLKZOEOYAKHREP-UHFFFAOYSA-N"
  },
  "panel_path": "matches/cluster_2025-10-11T12-34-56_match.png",
  "cluster_metadata": {
    "size": 6,
    "bonds": 5,
    "density": 0.333
  }
}
```

## API Usage

### Backend Endpoint

The matching functionality is exposed via REST API:

```http
POST /simulation/{simulation_id}/substance/{substance_id}/match
```

**Response:**

```json
{
  "success": true,
  "cluster_id": "SUB_abc123",
  "smiles": "CCCCCC",
  "timestamp": "2025-10-11_12-34-56",
  "files": {
    "panel": "matches/cluster_2025-10-11_12-34-56_match.png",
    "live_png": "matches/cluster_2025-10-11_12-34-56_live.png",
    "mol": "matches/cluster_2025-10-11_12-34-56.mol",
    "xyz": "matches/cluster_2025-10-11_12-34-56.xyz",
    "pubchem_png": "matches/cluster_2025-10-11_12-34-56_pubchem.png",
    "pubchem_mol": "matches/cluster_2025-10-11_12-34-56_pubchem.mol",
    "pubchem_xyz": "matches/cluster_2025-10-11_12-34-56_pubchem.xyz"
  },
  "pubchem_match": {
    "cid": 8058,
    "name": "hexane",
    "formula": "C6H14",
    "mw": 86.18,
    "smiles": "CCCCCC",
    "inchikey": "VLKZOEOYAKHREP-UHFFFAOYSA-N"
  }
}
```

### Frontend Usage

In TypeScript/JavaScript:

```typescript
import { SimulationAPI } from './lib/api'

const api = new SimulationAPI()

// Match a substance to PubChem
const result = await api.matchSubstanceToPubchem(simulationId, substanceId)

console.log(`SMILES: ${result.smiles}`)
if (result.pubchem_match) {
  console.log(`Match: ${result.pubchem_match.name} (CID ${result.pubchem_match.cid})`)
}
```

## Manual/Standalone Usage (Optional)

If you want to process clusters independently of the backend:

### Using the Watcher

```bash
python matcher/watcher.py
```

Place PNG + JSON files in `clustershots/` folder and the watcher will process them automatically.

### Using the CLI

```bash
python matcher/matcher.py clustershots/cluster_XXX.png
```

This will process the single cluster and generate results in `matches/`.

## How It Works

### Complete Flow

1. **User clicks download icon** in Novelty Panel
2. **Frontend calls** `POST /simulation/{id}/substance/{substance_id}/match`
3. **Backend retrieves** substance from simulation catalog
4. **Backend processes:**
   - Converts topology to RDKit Mol
   - Generates SMILES representation
   - Queries PubChem for similar molecules
   - Exports to MOL and XYZ formats
   - Renders comparison images
   - Composes side-by-side panel
5. **Backend returns** match results to frontend
6. **User sees** success message with match details

### 1. Topology to SMILES Conversion

The matcher converts cluster topology to molecular representation:

- **Nodes** → **Atoms**: Maps particles to chemical elements
- **Bonds** → **Bonds**: Preserves connection topology
- **Heuristics**: Uses node degree to guess element types:
  - Degree 1 → Hydrogen (H)
  - Degree 2 → Oxygen (O)
  - Degree ≥3 → Carbon (C)

### 2. Chemical File Format Export

After conversion, the matcher automatically exports to standard formats:

- **MOL format (.mol)**: MDL Molfile, widely supported in chemistry software
- **XYZ format (.xyz)**: Simple Cartesian coordinates with 3D geometry
  - Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) for 3D structure generation
  - UFF (Universal Force Field) optimization for realistic geometry
  - Compatible with visualization tools like VMD, PyMOL, Avogadro, Jmol

### 3. PubChem Similarity Search

Once converted to SMILES, the matcher:

1. Queries PubChem's similarity API (threshold=80%)
2. Retrieves the top matching compound
3. Fetches detailed properties (name, formula, weight, etc.)
4. Exports PubChem molecule to the same formats

### 4. Panel Composition

Creates a side-by-side comparison:

- **Left**: LIVE 2.0 cluster (rendered from topology)
- **Right**: PubChem molecule (retrieved structure)
- **Metadata**: Size, bonds, density, chemical name

## Advanced Configuration

### Custom Element Mapping

To improve SMILES accuracy, you can enhance the element mapping in `matcher/chem.py`:

```python
def choose_symbol(deg: int, given_label: str | None = None, energy: float = 0.0) -> str:
    # Add energy-based heuristics
    if energy > 100.0:
        return "O"  # High energy → oxygen
    if energy > 50.0:
        return "N"  # Medium energy → nitrogen
    return "C"  # Default → carbon
```

### Adjusting Similarity Threshold

In `matcher/matcher.py`, change the PubChem threshold:

```python
top = pubchem_similar_top(smiles, threshold=70)  # Lower = more permissive
```

## API Endpoints

### Match Substance to PubChem

Main endpoint for cluster matching (see "API Usage" section above):

```http
POST /simulation/{simulation_id}/substance/{substance_id}/match
```

### Get Substance Details

Retrieve detailed topology for a specific substance (used internally by the match endpoint):

```http
GET /simulation/{simulation_id}/substance/{substance_id}/details
```

**Response:**

```json
{
  "id": "SUB_abc123",
  "nodes": [...],
  "bonds": [...],
  "metadata": {...}
}
```

## Troubleshooting

### "No PubChem match found"

- The cluster structure is too unusual or not chemically valid
- Try lowering the similarity threshold (default: 80%)
- Check that the cluster has at least 2 nodes and 1 bond

### "Missing metadata JSON"

- Ensure both `.png` and `.json` files are saved with the same name
- The watcher waits up to 4 seconds for the JSON file
- Check file permissions in the `clustershots/` folder

### Import Errors

If you get `ModuleNotFoundError: No module named 'rdkit'`:

```bash
# Try conda instead
conda install -c conda-forge rdkit

# Or verify pip installation
pip show rdkit-pypi
```

### PubChem API Timeout

- PubChem API can be slow for complex queries
- The timeout is set to 25 seconds per request
- Check your internet connection

## Using Exported Chemical Files

The matcher exports structures in standard chemical file formats that can be opened in various chemistry software:

### MOL Files (.mol)

MOL files can be imported into:
- **ChemDraw** - Chemical drawing software
- **Marvin Sketch** - Structure editor
- **RDKit** - Cheminformatics toolkit
- **OpenBabel** - Chemical toolbox for conversion

Example Python usage:
```python
from rdkit import Chem

# Load MOL file
mol = Chem.MolFromMolFile("matches/cluster_XXX.mol")
print(f"Atoms: {mol.GetNumAtoms()}")
print(f"SMILES: {Chem.MolToSmiles(mol)}")
```

### XYZ Files (.xyz)

XYZ files contain 3D coordinates and can be visualized in:
- **VMD** (Visual Molecular Dynamics) - Molecular visualization
- **PyMOL** - 3D molecular graphics
- **Avogadro** - Molecular editor and visualizer
- **Jmol** - Java-based molecular viewer

Example visualization in PyMOL:
```bash
pymol matches/cluster_XXX.xyz
```

Example loading in Python:
```python
# Read XYZ file
with open("matches/cluster_XXX.xyz") as f:
    lines = f.readlines()
    num_atoms = int(lines[0])
    comment = lines[1]
    for i in range(2, 2 + num_atoms):
        symbol, x, y, z = lines[i].split()
        print(f"{symbol}: ({x}, {y}, {z})")
```

### Conversion to Other Formats

Use OpenBabel to convert to additional formats:
```bash
# Convert to PDB format
obabel matches/cluster_XXX.mol -O cluster.pdb

# Convert to SDF format
obabel matches/cluster_XXX.mol -O cluster.sdf

# Convert to SMILES
obabel matches/cluster_XXX.mol -osmi
```

## Examples

### Example 1: Simple Linear Chain

**Cluster**: 6 atoms in a chain

```
A-A-A-A-A-A
```

**SMILES**: `CCCCCC`

**PubChem Match**: Hexane (CID 8058)

### Example 2: Branched Structure

**Cluster**: 5 atoms with branching

```
  A
  |
A-A-A-A
```

**SMILES**: `CC(C)CC`

**PubChem Match**: Isopentane (CID 6556)

### Example 3: Cyclic Structure

**Cluster**: 6 atoms in a ring

```
 A---A
 |   |
 A   A
 |   |
 A---A
```

**SMILES**: `C1CCCCC1`

**PubChem Match**: Cyclohexane (CID 8078)

## Limitations

1. **Element Mapping**: The heuristic mapping is approximate. Real chemistry uses specific elements, while LIVE 2.0 uses generic particles.

2. **Bond Orders**: The system defaults to single bonds. Real molecules have single, double, and triple bonds.

3. **3D Structure**: LIVE 2.0 is 2D, but real molecules are 3D. The matcher ignores stereochemistry.

4. **Unusual Structures**: Highly unusual or physically impossible structures may not match anything in PubChem.

## Future Enhancements

- [ ] Machine learning-based element prediction
- [ ] Bond order inference from particle properties
- [ ] Local Tanimoto similarity scoring
- [ ] Materials Project integration for periodic structures
- [ ] Direct backend upload (no manual file placement)
- [ ] Real-time matching preview in UI

## Contributing

To add features or fix bugs:

1. Modify the relevant module in `matcher/`
2. Test with sample cluster files
3. Update this README with new functionality

## References

- [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [SMILES Notation](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)

## License

Part of the LIVE 2.0 project. See main project LICENSE for details.

## Contact

For questions or issues, please open an issue in the LIVE 2.0 repository.

---

**Happy Matching! 🧪🔬**

