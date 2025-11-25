# Physics Parameters Database

**Date**: October 12, 2025  
**Status**: ✅ Day 1 Complete - Infrastructure Ready

## Overview

The Physics Database provides **literature-cited parameters** for all physical interactions in Live 2.0 simulations. This ensures scientific rigor and traceability to peer-reviewed sources.

## Architecture

### Core Components

1. **`backend/sim/core/physics_db.py`**
   - `Citation` - Dataclass for scientific citations
   - `BondParameters` - Morse potential parameters with metadata
   - `VanDerWaalsParameters` - Lennard-Jones parameters
   - `PhysicsDatabase` - Central database manager

2. **`data/physics_parameters.json`**
   - JSON database of all parameters
   - Follows schema in `data/physics_parameters_schema.json`
   - Version controlled

3. **`scripts/validate_parameters.py`**
   - Validation script for database integrity
   - Checks citations, physical plausibility, completeness
   - Generates validation reports

## Data Structures

### Bond Parameters

Morse potential: `V(r) = D_e * (1 - exp(-a*(r - r_e)))^2`

```python
@dataclass
class BondParameters:
    atom_pair: Tuple[str, str]  # ('C', 'C')
    bond_order: int              # 1, 2, or 3
    D_e: float                   # Dissociation energy (kJ/mol)
    r_e: float                   # Equilibrium length (Angstrom)
    a: float                     # Width parameter (1/Angstrom)
    k_spring: float              # Force constant (kJ/mol/A²)
    source: Citation             # Literature citation
    confidence: str              # 'high'/'medium'/'low'
    method: str                  # 'experimental'/'DFT'/etc
```

### Van der Waals Parameters

Lennard-Jones potential: `V(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]`

```python
@dataclass
class VanDerWaalsParameters:
    atom_type: str      # 'C', 'N', 'O', etc
    epsilon: float      # Well depth (kJ/mol)
    sigma: float        # Zero-crossing (Angstrom)
    source: Citation    # Literature citation
    method: str         # 'UFF'/'OPLS'/etc
```

## Usage

### Loading Database

```python
from backend.sim.core.physics_db import PhysicsDatabase

db = PhysicsDatabase('data/physics_parameters.json')

# Get bond parameters
params = db.get_bond_parameters('C', 'C', order=1)
print(f"C-C bond: D_e = {params.D_e} kJ/mol")
print(f"Citation: {params.source.format_apa()}")

# Get VDW parameters (with combination rules)
epsilon, sigma = db.get_vdw_parameters('C', 'N')
```

### Adding Parameters

```python
from backend.sim.core.physics_db import BondParameters, Citation

# Create citation
source = Citation(
    doi='10.1201/9781420007282',
    authors=['Luo, Y.-R.'],
    title='Comprehensive Handbook of Chemical Bond Energies',
    journal='CRC Press',
    year=2007
)

# Create bond parameters
bond = BondParameters(
    atom_pair=('C', 'O'),
    bond_order=2,  # C=O double bond
    D_e=745.0,     # kJ/mol (strong!)
    r_e=1.21,      # Angstrom
    a=2.3,         # 1/Angstrom
    source=source,
    confidence='high',
    method='experimental'
)

db.add_bond_parameters(bond)
db.save()
```

### Validation

```bash
# Validate database
python scripts/validate_parameters.py --db data/physics_parameters.json

# Strict mode (warnings = errors)
python scripts/validate_parameters.py --db data/physics_parameters.json --strict
```

## Validation Checks

The validation script checks:

1. **Citations**
   - All parameters have DOI
   - Complete author/year/title
   
2. **Physical Plausibility**
   - D_e in [50, 1000] kJ/mol
   - r_e in [0.5, 3.0] Angstrom
   - Higher bond orders → higher D_e
   
3. **Completeness**
   - Important bonds (C-C, C-H, C-N, C-O, etc.)
   - Important atoms (H, C, N, O, S, P)
   
4. **Consistency**
   - k_spring matches Morse parameters
   - No duplicates

## Data Sources

### Recommended Sources

**Bond Parameters:**
- NIST Chemistry WebBook: https://webbook.nist.gov/
- CCCBDB: https://cccbdb.nist.gov/
- Luo (2007): doi:10.1201/9781420007282

**Van der Waals:**
- UFF (Rappé et al. 1992): doi:10.1021/ja00051a040
- OPLS force fields
- AMBER force fields

### Citation Format

Use APA style:
```python
citation = Citation(
    doi='10.1021/ja00051a040',
    authors=['Rappé, A. K.', 'Casewit, C. J.', 'Colwell, K. S.', 
             'Goddard III, W. A.', 'Skiff, W. M.'],
    title='UFF, a full periodic table force field',
    journal='Journal of the American Chemical Society',
    year=1992
)
```

## Export for Paper

### LaTeX Table

```python
db.export_table_for_paper('paper/tables/tableS1_parameters.tex', format='latex')
```

Generates:
```latex
\begin{table}[h]
\centering
\caption{Physical Parameters from Literature}
\begin{tabular}{llccccl}
\hline
Atoms & Order & $D_e$ (kJ/mol) & $r_e$ (\AA) & ... \\
\hline
C--C & 1 & 348.0 & 1.540 & 1.800 & experimental & \cite{...} \\
\hline
\end{tabular}
\end{table}
```

### CSV for Supplementary

```python
db.export_table_for_paper('data/supplementary/parameters.csv', format='csv')
```

## Integration with Simulation

### Current Status
🟡 **In Progress** - Day 1 complete, integration pending

### Next Steps (Day 4)

1. Modify `PotentialSystem` to use database:
```python
# backend/sim/core/potentials.py
class PotentialSystem:
    def __init__(self, config, physics_db=None):
        self.db = physics_db
        self.use_db = config.use_physics_db
```

2. Add config flag:
```yaml
# config.yaml
use_physics_db: true
physics_db_path: "data/physics_parameters.json"
```

3. Fallback to config values if parameter not found

## Statistics

### Current Database (Example)
- Bonds: 1
- VDW: 1
- Citations: 2
- Methods: experimental (100%)

### Target (End of Week 2)
- Bonds: 50+
- VDW: 10 (H, C, N, O, S, P, F, Cl, Br, I)
- Citations: 30+
- Methods: experimental (70%), DFT (20%), fitted (10%)

## Roadmap

### ✅ Day 1: Infrastructure (COMPLETE)
- [x] PhysicsDatabase class
- [x] JSON schema
- [x] Validation script
- [x] Example data
- [x] Documentation

### 📋 Day 2-3: Data Collection (TODO)
- [ ] Collect 50+ bond parameters from literature
- [ ] Collect 10 VDW parameters
- [ ] Verify all citations
- [ ] Validate physical plausibility

### 📋 Day 4: Integration (TODO)
- [ ] Modify PotentialSystem
- [ ] Add config flags
- [ ] Test simulations with DB vs config
- [ ] Migration guide

## References

1. Luo, Y.-R. (2007). *Comprehensive Handbook of Chemical Bond Energies*. CRC Press. doi:10.1201/9781420007282

2. Rappé, A. K., et al. (1992). UFF, a full periodic table force field. *Journal of the American Chemical Society*, 114(25), 10024-10035. doi:10.1021/ja00051a040

3. NIST Chemistry WebBook. https://webbook.nist.gov/

4. CCCBDB: Computational Chemistry Comparison and Benchmark Database. https://cccbdb.nist.gov/

## Notes

- All parameters must have literature citations
- Prefer experimental > high-level theory (CCSD(T)) > DFT > fitted
- Document any assumptions or approximations
- Version control the database JSON file
- Update `last_updated` timestamp on changes

---

**Last Updated**: October 12, 2025  
**Status**: ✅ Day 1 Complete  
**Next**: Day 2 - Literature data collection

