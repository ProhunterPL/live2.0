---
date: 2025-11-12
label: analysis
---

# Problem Nierealistycznych Cząsteczek - Analiza i Rozwiązanie

## 🔍 Objawy

- Cząsteczki są "dziwne" chemicznie
- Długie, rozciągnięte wiązania
- Niestabilne struktury (szybko się rozpadają)
- Połączenia między cząstkami które nie powinny tworzyć wiązań
- Density często < 0.3 (luźne struktury)

## 🎯 Główne Przyczyny

### 1. **Za Długie Maksymalne Odległości Wiązań**

**Problem w kodzie** (`binding.py` linia 341, 395):

```python
# Formation distance (stary kod):
if r <= PARTICLE_RADIUS_COMPILE * 6.8:  # = 3.4 jednostki
    # To jest BARDZO długie!

# Break distance:
if r > PARTICLE_RADIUS_COMPILE * 5.0:  # = 2.5 jednostki
    break_bond()
```

**Konsekwencje**:
- Cząstki mogą tworzyć wiązania na odległość **3.4 jednostki**
- To ~5-7x więcej niż realistyczne wiązania kowalencyjne (~0.5-1.5 Å)
- Powstają "pajęczyny" zamiast zwartych molekuł

### 2. **Za Niski Próg Prawdopodobieństwa Wiązania**

**Problem** (`binding.py` linia 344):

```python
# Stary kod:
if binding_probability > 0.005:  # 0.5% wystarczy!
    form_bond()

# Nowy kod (już w pliku):
if binding_probability > 0.15:  # 15% wymagane
```

**Konsekwencje**:
- Prawie każda para cząstek tworzy wiązania
- Brak selektywności chemicznej
- Powstają nierealistyczne połączenia

### 3. **Za Wysoki Próg Zrywania (Strain)**

**Problem** (`binding.py` linia 382):

```python
# Stary kod:
if strain > 3.0:  # 300% rozciągnięcie!
    break_bond()

# Nowy kod (już w pliku):
max_strain = 0.5  # 50% dla większości
max_strain = 0.3  # 30% dla kovalent
```

**Konsekwencje**:
- Wiązania mogą się rozciągać 3x bez zerwania
- Powstają "gumowe" cząsteczki zamiast sztywnych
- Nierealistyczne dynamika

### 4. **Parametry Sił Wiązań**

**Problem** (`binding.py` linia 524-527):

```python
self.bond_type_params = {
    0: {'k_spring': 2.0, 'rest_len': 1.0, ...},    # vdW
    1: {'k_spring': 500.0, 'rest_len': 0.8, ...},  # covalent - OK
    2: {'k_spring': 50.0, 'rest_len': 1.2, ...},   # H-bond
    3: {'k_spring': 100.0, 'rest_len': 0.9, ...}   # metallic
}
```

**Problem**:
- `rest_len` jest OK, ale...
- Formation distance (3.4) >> rest_len (0.8-1.2)
- Cząstki łączą się na odległość, potem "zjeżdżają" do rest_len
- To powoduje gwałtowne ruchy i niestabilność

## ✅ Kompletne Rozwiązanie

### Rozwiązanie 1: Poprawione Parametry Wiązań

Stwórz nowy plik konfiguracyjny `configs/realistic_chemistry.yaml`:

```yaml
mode: "open_chemistry"

# === CORE PARAMETERS ===
grid_width: 128
grid_height: 128
dt: 0.05  # Mniejszy timestep = bardziej stabilne
max_particles: 500
particle_radius: 0.5

# === REALISTIC BOND PARAMETERS ===
# These are calibrated for realistic chemistry

# Formation is more selective
binding_threshold: 0.15  # Wymagane 15% probability (nie 0.5%!)
unbinding_threshold: 0.12  # Stabilne wiązania

# Bond checking
bond_check_interval: 100  # Częstsze sprawdzanie
max_bond_length: 2.0  # MAX 2.0 jednostki (nie 3.4!)

# === REALISTIC DISTANCES ===
# In code, this translates to:
# - Formation distance: 1.2-1.8 units (depending on type)
# - Break distance: 2.0 units (hard limit)
# - Rest length: 0.8-1.2 units
# Result: bonds form near equilibrium, not far away!

# === PARTICLE DYNAMICS ===
initial_particle_count: 100
mass_range: [1.0, 16.0]  # H to O
charge_range: [-0.5, 0.5]

# === ENERGY ===
energy_transfer_rate: 0.01
energy_dissipation: 0.01
thermostat_enabled: true
thermostat_target_temp: 300.0
thermostat_alpha: 0.05

# === NOVELTY DETECTION ===
novelty_check_interval: 300  # Częstsza detekcja
min_cluster_size: 3
detect_novel_substances: true

# === VALIDATION ===
validate_every_n_steps: 200
enable_thermodynamic_validation: true
energy_tolerance: 0.002
momentum_tolerance: 0.0002

# === MUTATIONS ===
p_mut_base: 0.002
p_mut_gain: 20.0
attr_sigma: 0.1
```

### Rozwiązanie 2: Napraw Kod Wiązań

**WAŻNE**: Twój kod już ma pewne poprawki (linie 313-397), ale musisz się upewnić że są **aktywowane**.

Sprawdź w `backend/sim/core/binding.py`:

#### A. Formation Distance (linia ~341)

```python
# MUSI BYĆ:
max_formation_dist_covalent = PARTICLE_RADIUS_COMPILE * 2.0 * 1.2  # = 1.2 units
max_formation_dist_vdW = PARTICLE_RADIUS_COMPILE * 3.0 * 1.2      # = 1.8 units
max_formation_dist_hbond = PARTICLE_RADIUS_COMPILE * 2.5 * 1.2    # = 1.5 units

# NIE MOŻE BYĆ (stary kod):
# if r <= PARTICLE_RADIUS_COMPILE * 6.8:  # = 3.4 units - ZA DŁUGIE!
```

#### B. Binding Probability Threshold (linia ~344)

```python
# MUSI BYĆ:
if binding_probability > 0.15:  # 15% wymagane

# NIE MOŻE BYĆ (stary kod):
# if binding_probability > 0.005:  # 0.5% - ZA NISKI!
```

#### C. Strain Threshold (linia ~382-389)

```python
# MUSI BYĆ:
max_strain = 0.5  # default 50%
if bond_type == 1:  # covalent
    max_strain = 0.3  # 30% - REALISTYCZNE
elif bond_type == 0:  # vdW
    max_strain = 0.8  # 80%

# NIE MOŻE BYĆ (stary kod):
# if strain > 3.0:  # 300% - ZA WYSOKIE!
```

#### D. Max Distance (linia ~395-397)

```python
# MUSI BYĆ:
max_distance = PARTICLE_RADIUS_COMPILE * 4.0  # 2.0 units default
if bond_type == 1:  # covalent
    max_distance = PARTICLE_RADIUS_COMPILE * 3.0  # 1.5 units

# NIE MOŻE BYĆ (stary kod):
# if r > PARTICLE_RADIUS_COMPILE * 5.0:  # 2.5 units - ZA WYSOKIE!
```

### Rozwiązanie 3: Dodaj Walidację Chemiczną

Stwórz nowy plik `backend/sim/validation/chemistry_validator.py`:

```python
"""
Chemical realism validator
Checks if generated molecules make chemical sense
"""
import numpy as np
from typing import Dict, List, Tuple

class ChemistryValidator:
    """Validate chemical realism of clusters"""
    
    def __init__(self):
        # Realistic bond length ranges (in simulation units)
        # Based on literature: C-C = 1.54 Å, C-O = 1.43 Å, etc.
        self.realistic_bond_lengths = {
            'covalent': (0.7, 1.5),     # 0.7-1.5 units
            'vdW': (1.2, 2.0),           # van der Waals
            'h_bond': (1.0, 1.8),        # H-bonds
            'metallic': (0.8, 1.4)       # metallic
        }
        
        # Realistic density ranges
        self.realistic_density = {
            'min': 0.25,  # Very loose molecules
            'typical': 0.4,  # Typical organic molecules
            'max': 0.8   # Dense structures
        }
    
    def validate_cluster(self, cluster: Dict) -> Dict:
        """
        Validate if cluster is chemically realistic
        
        Returns:
            {
                'is_realistic': bool,
                'score': float (0-1),
                'issues': List[str],
                'warnings': List[str]
            }
        """
        issues = []
        warnings = []
        score = 1.0
        
        # 1. Check bond lengths
        if 'bonds' in cluster and 'positions' in cluster:
            bond_lengths = self._calculate_bond_lengths(cluster)
            
            for bond_len in bond_lengths:
                if bond_len > 2.0:
                    issues.append(f"Bond too long: {bond_len:.2f} units (max 2.0)")
                    score *= 0.5
                elif bond_len > 1.8:
                    warnings.append(f"Bond stretched: {bond_len:.2f} units")
                    score *= 0.9
        
        # 2. Check density
        density = cluster.get('properties', {}).get('graph_density', 0)
        if density < self.realistic_density['min']:
            issues.append(f"Density too low: {density:.3f} (min {self.realistic_density['min']})")
            score *= 0.7
        elif density < self.realistic_density['typical']:
            warnings.append(f"Density low: {density:.3f} (typical {self.realistic_density['typical']})")
            score *= 0.95
        
        # 3. Check cluster size vs bonds
        size = cluster.get('size', 0)
        bonds = cluster.get('properties', {}).get('bonds', 0)
        
        if size > 2:
            # For realistic molecules: bonds ≈ size - 1 (tree-like) to size * 1.5 (cyclic)
            min_bonds = size - 1
            max_bonds = int(size * 1.5)
            
            if bonds < min_bonds:
                issues.append(f"Too few bonds: {bonds} for {size} atoms (min {min_bonds})")
                score *= 0.6
            elif bonds > max_bonds:
                warnings.append(f"Many bonds: {bonds} for {size} atoms (max {max_bonds})")
                score *= 0.95
        
        # 4. Check bond length variance
        if 'bonds' in cluster and len(bond_lengths) > 1:
            variance = np.var(bond_lengths)
            if variance > 0.5:
                warnings.append(f"High bond length variance: {variance:.3f}")
                score *= 0.9
        
        is_realistic = len(issues) == 0 and score > 0.7
        
        return {
            'is_realistic': is_realistic,
            'score': score,
            'issues': issues,
            'warnings': warnings,
            'bond_lengths': bond_lengths if 'bonds' in cluster else [],
            'density': density
        }
    
    def _calculate_bond_lengths(self, cluster: Dict) -> List[float]:
        """Calculate all bond lengths in cluster"""
        positions = np.array(cluster['positions'])
        bonds = cluster.get('bonds', [])
        
        lengths = []
        for bond in bonds:
            i, j = bond[0], bond[1]
            if i < len(positions) and j < len(positions):
                pos_i = positions[i]
                pos_j = positions[j]
                length = np.linalg.norm(pos_i - pos_j)
                lengths.append(length)
        
        return lengths
```

### Rozwiązanie 4: Szybki Test

Stwórz skrypt testowy `test_realistic_chemistry.py`:

```python
"""
Test realistic chemistry parameters
"""
from backend.sim.config import SimulationConfig
from backend.sim.core.stepper import SimulationStepper
import taichi as ti

# Initialize Taichi
ti.init(arch=ti.cpu)

# Create config with realistic parameters
config = SimulationConfig(
    mode="open_chemistry",
    max_particles=100,
    initial_particle_count=50,
    binding_threshold=0.15,      # REALISTYCZNE (nie 0.5!)
    unbinding_threshold=0.12,
    particle_radius=0.5,
    dt=0.05,
    grid_width=128,
    grid_height=128,
    novelty_check_interval=300,
    min_cluster_size=3
)

# Create simulation
sim = SimulationStepper(config)

# Run 1000 steps
print("Running 1000 steps with realistic chemistry...")
for i in range(1000):
    sim.step()
    
    if i % 100 == 0:
        # Check clusters
        clusters = sim.binding.get_clusters(min_size=2)
        bonds = sim.binding.get_bonds()
        
        print(f"\nStep {i}:")
        print(f"  Clusters: {len(clusters)}")
        print(f"  Bonds: {len(bonds)}")
        
        if bonds:
            # Check bond lengths
            positions = sim.particles.position.to_numpy()
            bond_lengths = []
            for b in bonds[:10]:  # Check first 10
                i_idx, j_idx, strength = b
                pos_i = positions[i_idx]
                pos_j = positions[j_idx]
                length = ((pos_i[0] - pos_j[0])**2 + (pos_i[1] - pos_j[1])**2)**0.5
                bond_lengths.append(length)
            
            avg_length = sum(bond_lengths) / len(bond_lengths)
            max_length = max(bond_lengths)
            
            print(f"  Bond lengths: avg={avg_length:.2f}, max={max_length:.2f}")
            
            if max_length > 2.0:
                print(f"  ⚠️ WARNING: Bond too long! ({max_length:.2f} > 2.0)")
            elif avg_length < 1.5:
                print(f"  ✅ OK: Bonds are realistic")

print("\n✅ Test completed!")
```

## 🚀 Jak Zastosować

### Krok 1: Sprawdź Aktualny Kod

```powershell
# Sprawdź czy masz poprawki w binding.py
Get-Content backend\sim\core\binding.py | Select-String -Pattern "binding_probability > 0.15"

# Powinno pokazać linię ~344:
# if binding_probability > 0.15:  # INCREASED from 0.005
```

Jeśli **NIE** pokazuje tej linii, kod jest **STARY** i trzeba go naprawić.

### Krok 2: Użyj Realistycznej Konfiguracji

```powershell
# Skopiuj konfigurację
# (stwórz plik configs/realistic_chemistry.yaml z contentem powyżej)

# Restart backendu z nową konfiguracją
.\kill_backend.ps1

cd backend
python -m api.server --config ../configs/realistic_chemistry.yaml

# W drugim terminalu uruchom frontend
cd frontend
npm run dev
```

### Krok 3: Uruchom Test

```powershell
# Uruchom test (jeśli stworzysz test_realistic_chemistry.py)
python test_realistic_chemistry.py

# Powinno pokazać:
# ✅ OK: Bonds are realistic
# (nie: ⚠️ WARNING: Bond too long!)
```

### Krok 4: Monitoruj W Symulacji

W trakcie symulacji, sprawdzaj:

**Largest Connected Cluster** (panel po lewej):
- **Density** powinno być > 0.3 (lepiej > 0.4)
- **Avg Mass** powinno być sensowne (1-10)
- **Total Energy** nie powinno być 0.00 (jeśli jest 0, to problem z energią)

**Wizualizacja**:
- Wiązania powinny być **krótkie** (< 2 jednostki na ekranie)
- Cząsteczki powinny być **zwarte**, nie rozciągnięte
- Brak "pajęczyn" (długie cienkie wiązania)

## 📊 Oczekiwane Rezultaty

### Przed naprawą:
```
Cluster:
  Size: 9 particles
  Bonds: 8
  Density: 0.222  ← Niskie!
  Bond lengths: 1.2, 1.5, 2.8, 3.1  ← Za długie!
  Structure: Rozciągnięta "pajęczyna"
```

### Po naprawie:
```
Cluster:
  Size: 9 particles
  Bonds: 10-12
  Density: 0.45  ← OK!
  Bond lengths: 0.9, 1.1, 1.3, 1.4  ← Realistyczne!
  Structure: Zwarta molekuła
```

## 💡 Podsumowanie Zmian

| Parametr | Stara Wartość | Nowa Wartość | Efekt |
|----------|---------------|--------------|-------|
| Max formation distance | 3.4 units | 1.2-1.8 units | Krótsze wiązania |
| Binding probability threshold | 0.005 (0.5%) | 0.15 (15%) | Więcej selektywności |
| Max strain before break | 3.0 (300%) | 0.3-0.8 (30-80%) | Sztywniejsze wiązania |
| Max distance before break | 2.5 units | 1.5-2.0 units | Stabilniejsze struktury |

## ⚠️ Uwaga

Po zastosowaniu tych zmian:
- ✅ Cząsteczki będą bardziej realistyczne
- ✅ Density > 0.3-0.4
- ✅ Krótsze, stabilne wiązania
- ⚠️ Może być mniej klastrów (bo bardziej selektywne)
- ⚠️ Mniejsze struktury (bo wiązania tylko na krótki dystans)
- ✅ Ale za to chemicznie sensowne!

## 🎯 Priorytet Działań

1. **Najważniejsze**: Sprawdź czy `binding_probability > 0.15` jest w kodzie
2. **Ważne**: Użyj konfiguracji `realistic_chemistry.yaml`
3. **Pomocne**: Dodaj `ChemistryValidator` dla walidacji
4. **Opcjonalne**: Uruchom `test_realistic_chemistry.py` dla testu

---

**Pytanie**: Czy mam sprawdzić twój aktualny kod `binding.py` i naprawić go jeśli jest stary?

