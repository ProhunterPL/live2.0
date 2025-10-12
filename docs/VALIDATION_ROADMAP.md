# Live 2.0 - Validation Roadmap
**Mapa Drogowa Walidacji Naukowej**

*Ostatnia aktualizacja: 12 października 2025 (21:30 CET)*  
*Ostatnia zmiana: Completed Week 2 Day 1 - Physics Database Infrastructure*

---

## 📋 Executive Summary

Ten dokument łączy:
1. **SCIENTIFIC_INTEGRITY_VERIFICATION.md** - Status obecny (co działa)
2. **Live 2-plan walidacji naukowej.md** - Plan przyszłości (co zrobić)

**Teza**: Podstawy naukowe symulacji są solidne. Teraz budujemy na nich wiarygodność publikacyjną.

---

## ✅ PHASE 0: FOUNDATIONS (COMPLETED - Październik 2025)

### 🎯 Cel: Upewnić się że fundamenty fizyki działają poprawnie

#### ✅ Walidacja Termodynamiczna - AKTYWNA
**Status**: Zaimplementowana i działająca

- ✅ **ThermodynamicValidator** istnieje (`backend/sim/core/thermodynamics.py`)
  - Zachowanie energii: `E_after = E_before + E_injected - E_dissipated ± ε`
  - Tolerancja: 0.1% (1e-3)
  - Częstotliwość: Co 10,000 kroków
  
- ✅ **Zachowanie pędu**
  - Sprawdzane: `Σ(m·v) = const`
  - Tolerancja: 0.01% (1e-4)
  - Częstotliwość: Co 10,000 kroków

- ✅ **Rozkład Maxwella-Boltzmanna**
  - Weryfikacja rozkładu prędkości
  - Temperatura z energii kinetycznej: `T = m⟨v²⟩ / (2k_B)`
  - Sampling: 200 cząstek (wydajność vs dokładność)

- ✅ **II Zasada Termodynamiki**
  - Sprawdzane: `ΔS ≥ 0`
  - Entropia: Shannon (konfiguracyjna) + kinetyczna
  - Sampling: 200 cząstek

**Kod**:
```python
# backend/sim/core/stepper.py:300-323
if self.enable_validation and self.validator is not None:
    validation_results = self.validator.validate_essential_only(
        state_before, state_after, 
        energy_injected, energy_dissipated, 
        self.step_count
    )
```

#### ✅ Optymalizacje Wydajnościowe - COMPLETED
**Status**: Symulacja stabilna do 3000+ kroków

**Problem**: Freeze po ~1500 kroków z powodu:
- Wycieku pamięci (`energy_history` jako lista)
- O(n²) operacji każdy krok (`_attract_particles_for_bonding`)
- Zbyt częste kopiowanie GPU→CPU

**Rozwiązanie** (Październik 2025):
- ✅ `energy_history` → `deque(maxlen=1000)`
- ✅ Clustering assistance: co 50 kroków (było: każdy krok)
- ✅ Particle attraction: WYŁĄCZONE (O(n²) nieakceptowalne)
- ✅ Diagnostyka: co 500 kroków (było: 10)
- ✅ Garbage collection: co 500 kroków
- ✅ Sampling w metrics: 500 cząstek (było: 1000)

**Rezultat**:
- ✅ Symulacja działa płynnie ponad 3000+ kroków
- ✅ Zużycie RAM stabilne ~500-800 MB (było: 2-4 GB+)
- ✅ Wydajność nie degraduje się w czasie

**Dokumentacja**: `docs/MEMORY_PERFORMANCE_FIX.md`

#### ✅ Integralność Naukowa - VERIFIED
**Status**: Wszystkie kluczowe elementy fizyki zachowane

**Co zostało ZACHOWANE w 100%**:
1. ✅ Wszystkie prawa termodynamiki (I, II zasada)
2. ✅ Zachowanie energii i pędu
3. ✅ Rozkład Maxwella-Boltzmanna
4. ✅ System wiązań chemicznych
5. ✅ Potencjały międzycząsteczkowe (Lennard-Jones, Morse)
6. ✅ Dynamika ruchu (Euler/Verlet)
7. ✅ Novelty tracking
8. ✅ Adaptacyjny timestep z kontrolą błędu
9. ✅ Warunki brzegowe periodyczne
10. ✅ Thermal fluctuations

**Co zostało ZOPTYMALIZOWANE** (tylko częstotliwość):
- ⚠️ Wiązania: 100→150 kroków (50% wolniej, ale nadal szybko)
- ⚠️ Klastry: 200→300 kroków
- ⚠️ Mutacje: 200→300 kroków

**Co zostało WYŁĄCZONE** (tylko pomocnicze):
- ❌ `_attract_particles_for_bonding()` - O(n²) helper function
  - Główny system binding działa normalnie
  - Wiązania tworzą się przez potencjały

**Dokumentacja**: `docs/SCIENTIFIC_INTEGRITY_VERIFICATION.md`

#### ✅ Infrastruktura Techniczna
- ✅ Taichi GPU acceleration
- ✅ WebSocket real-time streaming
- ✅ Snapshot system
- ✅ Memory management
- ✅ Performance monitoring
- ✅ Metrics collection (deque-based, bounded)

---

## 📋 PHASE 1: VALIDATION SPRINT (4 Tygodnie - DO ZROBIENIA)

### 🎯 Cel: Osiągnąć publikowalną wiarygodność naukową

**Inspiracja**: "Physics-First, Emergence Second"
- Najpierw udowodnij że fizyka działa
- Potem pokaż emergencję chemii

---

### **Tydzień 1: Fundamenty Termodynamiczne** 🔥

#### Zadanie 1.1: Rozszerzyć ThermodynamicValidator
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 2 dni

**Co mamy**:
- ✅ Podstawowy validator działający
- ✅ 4 główne testy (energia, pęd, M-B, entropia)

**Co dodać**:
```python
# backend/sim/core/thermodynamics.py (rozszerzenie)

class ThermodynamicValidator:
    # ... istniejący kod ...
    
    def validate_virial_theorem(self, state):
        """
        Virial theorem: 2⟨T⟩ = -⟨V⟩ dla potencjałów ~ r^n
        Dobry test dla LJ potential
        """
        pass
    
    def validate_fluctuation_dissipation(self, trajectory):
        """
        Fluctuation-dissipation theorem
        Związek między fluktuacjami a odpowiedzią na perturbacje
        """
        pass
    
    def compute_heat_capacity(self, trajectory):
        """
        C_v = (⟨E²⟩ - ⟨E⟩²) / (k_B T²)
        """
        pass
```

**Deliverables**:
- [ ] `thermodynamics.py` z 3 nowymi metodami
- [ ] Unit testy dla każdej metody
- [ ] Notebook z przykładami: `notebooks/thermodynamics_extended.ipynb`

#### Zadanie 1.2: Continuous Validation Loop - UPGRADE
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

**Co mamy**:
- ✅ Walidacja co 10,000 kroków
- ✅ Logowanie wyników

**Co dodać**:
```python
# backend/sim/core/stepper.py (modyfikacja)

class SimulationStepper:
    def __init__(self, config):
        # ... istniejący kod ...
        
        # NOWE: Konfigurowalny validation
        self.validation_config = {
            'energy': {'enabled': True, 'interval': 10000},
            'momentum': {'enabled': True, 'interval': 10000},
            'maxwell_boltzmann': {'enabled': True, 'interval': 50000},
            'entropy': {'enabled': True, 'interval': 50000},
            'virial': {'enabled': False, 'interval': 100000},  # Nowe
        }
        
        # NOWE: Real-time alerts
        self.validation_alerts = []
        self.alert_threshold = {
            'energy': 0.01,  # 1% error triggers alert
            'momentum': 0.001,
        }
```

**Deliverables**:
- [ ] Konfigurowalny validation pipeline
- [ ] Real-time alerts system
- [ ] CSV export: `diagnostics/validation_log.csv`
- [ ] JSON export: `diagnostics/validation_summary.json`

#### Zadanie 1.3: Plots & Analysis
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

**Cel**: Wygenerować Figure 1 i Figure 2 do papera

```python
# scripts/analyze_thermodynamics.py

def plot_energy_conservation(validation_log, output_path):
    """
    Figure 1: Energy conservation over 10^6 steps
    - Panel A: E_total(t) ± expected range
    - Panel B: Relative error over time
    - Panel C: Cumulative drift
    """
    pass

def plot_maxwell_boltzmann_fit(validation_log, output_path):
    """
    Figure 2: Velocity distribution comparison
    - Histogram: empirical vs. theoretical M-B
    - Q-Q plot
    - Chi-square test results table
    """
    pass

def plot_entropy_evolution(validation_log, output_path):
    """
    Figure S1 (Supplementary): Entropy over time
    - ΔS distribution
    - Violations analysis
    """
    pass
```

**Deliverables**:
- [ ] `scripts/analyze_thermodynamics.py`
- [ ] Figure 1: `figures/fig1_energy_conservation.png` (300 DPI)
- [ ] Figure 2: `figures/fig2_maxwell_boltzmann.png` (300 DPI)
- [ ] Figure S1: `figures/figS1_entropy.png` (300 DPI)
- [ ] Data files: `data/results/thermodynamics/validation_10M_steps.csv`

#### Zadanie 1.4: Documentation Update
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

**Deliverables**:
- [ ] `docs/THERMODYNAMIC_VALIDATION.md` (complete guide)
- [ ] Mathematical derivations (LaTeX)
- [ ] Results summary with figures
- [ ] Supplementary materials draft

**Checkpoint Week 1**:
- [ ] Energy conservation < 0.1% drift? ✓ GO / ✗ NO-GO
- [ ] M-B distribution p > 0.05? ✓ GO / ✗ NO-GO
- [ ] Entropy increasing? ✓ GO / ✗ NO-GO

---

### **Tydzień 2: Parametry z Literatury** 📚

**Status**: ⚠️ KRYTYCZNA LUKA  
**Problem**: Obecnie parametry mogą być arbitralne ("z ręki")

#### Zadanie 2.1: Physics Database Schema
**Status**: ✅ COMPLETE (Oct 12, 2025)  
**Czas**: 1 dzień

```python
# backend/sim/core/physics_db.py

from dataclasses import dataclass
from typing import Optional, List

@dataclass
class Citation:
    doi: Optional[str]
    authors: List[str]
    title: str
    journal: str
    year: int
    url: Optional[str]
    notes: Optional[str] = None

@dataclass
class BondParameters:
    atom_pair: tuple  # ('C', 'C')
    bond_order: int   # 1/2/3
    
    # Morse: V = D_e*(1-exp(-a*(r-r_e)))^2
    D_e: float  # kJ/mol
    r_e: float  # Å
    a: float    # Å^-1
    
    # Spring (harmonic approx)
    k_spring: Optional[float]  # kJ/mol/Å²
    r_0: Optional[float]
    
    source: Citation
    confidence: str  # 'high'/'medium'/'low'
    method: str      # 'experimental'/'DFT'/'fitted'

@dataclass
class VanDerWaalsParameters:
    atom_type: str
    epsilon: float  # kJ/mol (well depth)
    sigma: float    # Å (zero-crossing)
    source: Citation
    method: str     # 'UFF'/'OPLS'/etc

class PhysicsDatabase:
    def __init__(self, db_path='data/physics_parameters.json'):
        self.bonds = {}
        self.vdw = {}
        self.load()
    
    def get_bond_parameters(self, atom_a, atom_b, order=1):
        """Get parameters with citation"""
        pass
    
    def get_vdw_parameters(self, atom_a, atom_b):
        """Lorentz-Berthelot combination rules"""
        pass
    
    def export_table_for_paper(self, output_path):
        """Generate LaTeX table (Table S1)"""
        pass
```

**Deliverables**:
- [x] `backend/sim/core/physics_db.py` ✅
- [x] JSON schema: `data/physics_parameters_schema.json` ✅
- [x] Validation script: `scripts/validate_parameters.py` ✅
- [x] Documentation: `docs/PHYSICS_DATABASE.md` ✅
- [x] Example database: `data/physics_parameters_example.json` ✅

#### Zadanie 2.2: Data Collection
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 3 dni (intensywne!)

**Źródła**:

**A) Bond Parameters**:
- NIST Chemistry WebBook (https://webbook.nist.gov/)
- CCCBDB (https://cccbdb.nist.gov/)
- Luo (2007) "Comprehensive Handbook of Chemical Bond Energies"

**B) Van der Waals**:
- UFF (Rappé et al. 1992) - doi:10.1021/ja00051a040
- OPLS force fields
- Experimental gas-phase data

```python
# scripts/collect_bond_parameters.py

def scrape_nist_bond_data():
    """Scrape bond data from NIST"""
    molecules = {
        'methane': '74-82-8',
        'formaldehyde': '50-00-0',
        'HCN': '74-90-8',
        'water': '7732-18-5',
        # ... 50+ prebiotic molecules
    }
    # Implementation...

def collect_literature_data():
    """Manual entry from papers"""
    luo_2007_data = [
        {
            'bond': ('C', 'C'), 'order': 1,
            'D_e': 348.0, 'r_e': 1.54,
            'source': 'doi:10.1201/9781420007282',
            'confidence': 'high',
            'method': 'experimental'
        },
        # ... 100+ entries
    ]
    return luo_2007_data
```

**Target**: 50+ bond types, 10 atom types (H, C, N, O, S, P, F, Cl, Br, I)

**Deliverables**:
- [ ] `data/physics_parameters.json` (complete database)
- [ ] `scripts/collect_bond_parameters.py`
- [ ] `scripts/collect_vdw_parameters.py`
- [ ] `data/raw/nist_bonds.csv`
- [ ] `data/raw/literature_params.csv`
- [ ] LaTeX table: `paper/tables/tableS1_parameters.tex`

#### Zadanie 2.3: Integration
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```python
# backend/sim/core/potentials.py (refactor)

class PotentialSystem:
    def __init__(self, config, physics_db):
        self.db = physics_db
        self.use_db = config.use_physics_db
    
    @ti.func
    def lennard_jones(self, type_i, type_j, r):
        """LJ with DB parameters"""
        if ti.static(self.use_db):
            eps, sig = self.db.get_vdw_parameters(type_i, type_j)
        else:
            # Fallback to config
            eps, sig = self.config.epsilon, self.config.sigma
        
        # ... standard LJ formula
```

**Deliverables**:
- [ ] Refactored `potentials.py`
- [ ] Config flag: `use_physics_db: true`
- [ ] Tests comparing old vs new (should be similar)
- [ ] Migration guide: `docs/MIGRATION_TO_DB_PARAMS.md`

**Checkpoint Week 2**:
- [ ] 30+ parameter sources documented? ✓ GO / ✗ NO-GO
- [ ] All parameters have DOIs? ✓ GO / ✗ NO-GO
- [ ] Database loading works? ✓ GO / ✗ NO-GO

---

### **Tydzień 3: Benchmark Reactions** 🧪

**Status**: ⚠️ KRYTYCZNA LUKA  
**Problem**: Brak walidacji przeciw znanym reakcjom chemicznym

#### Zadanie 3.1: Test Framework
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```python
# tests/benchmarks/test_known_reactions.py

import pytest
from backend.sim import Simulation
from backend.sim.analysis import ReactionDetector

class TestPreobioticReactions:
    """Benchmark against known prebiotic chemistry"""
    
    @pytest.mark.slow
    @pytest.mark.parametrize("seed", range(10))
    def test_formose_reaction(self, seed):
        """
        Formose: CH2O → sugars (autocatalytic)
        
        Expected (Breslow 1959):
        - Glycolaldehyde: 15-30% yield
        - Autocatalytic growth
        
        Validation:
        - Glycolaldehyde detected: YES
        - Yield in range: 15-30%
        - Rate increases over time: YES
        """
        sim = Simulation(config='configs/formose_test.yaml', seed=seed)
        sim.add_molecules('formaldehyde', count=1000)
        sim.add_catalyst('Ca(OH)2', concentration=0.01)
        
        trajectory = sim.run(steps=100000, record_every=100)
        
        detector = ReactionDetector(trajectory)
        products = detector.identify_products()
        
        assert 'glycolaldehyde' in products
        
        yield_ga = products['glycolaldehyde'].count / 1000
        assert 0.15 <= yield_ga <= 0.30, \
            f"Yield {yield_ga:.2%} outside [15%, 30%]"
        
        # Check autocatalysis
        rates = detector.compute_reaction_rates()
        is_autocatalytic = detector.test_autocatalysis(rates)
        assert is_autocatalytic
    
    @pytest.mark.slow
    def test_strecker_synthesis(self):
        """
        Strecker: RCHO + HCN + NH3 → amino acid
        Example: acetaldehyde → alanine
        
        Expected (Miller 1953): 5-15% yield
        """
        # Implementation...
    
    @pytest.mark.slow
    def test_HCN_polymerization(self):
        """
        HCN polymerization → oligomers → adenine
        
        Expected (Oró 1960):
        - Tetramer formation
        - Adenine (trace)
        """
        # Implementation...
```

**Deliverables**:
- [ ] `tests/benchmarks/test_formose.py`
- [ ] `tests/benchmarks/test_strecker.py`
- [ ] `tests/benchmarks/test_hcn_polymer.py`
- [ ] `tests/benchmarks/test_phosphorylation.py`
- [ ] `tests/benchmarks/test_detailed_balance.py`

#### Zadanie 3.2: Reference Data Collection
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```json
// data/benchmark_reactions.json

{
  "formose_reaction": {
    "reaction": "n·CH2O → sugars (autocatalytic)",
    "conditions": {
      "temperature": 298,
      "pH": 11.0,
      "catalyst": "Ca(OH)2"
    },
    "products": [
      {
        "name": "glycolaldehyde",
        "formula": "C2H4O2",
        "yield_range": [0.15, 0.30]
      }
    ],
    "sources": [
      {
        "doi": "10.1016/0040-4020(59)80055-X",
        "authors": ["Breslow"],
        "year": 1959
      }
    ]
  }
}
```

**Deliverables**:
- [ ] `data/benchmark_reactions.json`
- [ ] Literature review document
- [ ] Expected yields table

#### Zadanie 3.3: Analysis Tools
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 2 dni

```python
# backend/sim/analysis/reaction_detector.py

class ReactionDetector:
    """Analyze trajectories for reactions"""
    
    def identify_products(self) -> Dict[str, ProductInfo]:
        """Identify all molecules in trajectory"""
        pass
    
    def compute_reaction_rates(self) -> Dict[str, List[float]]:
        """d[A]/dt for each species"""
        pass
    
    def test_autocatalysis(self, rates, species=None) -> bool:
        """Test for autocatalytic behavior"""
        pass
    
    def compute_equilibrium_constant(self, reaction: str) -> float:
        """K_eq from concentrations"""
        pass
    
    def measure_rate_constants(self, reaction: str) -> Tuple[float, float]:
        """(k_forward, k_reverse)"""
        pass
```

**Deliverables**:
- [ ] `backend/sim/analysis/reaction_detector.py`
- [ ] `backend/sim/analysis/molecule_database.py`
- [ ] Unit tests
- [ ] Tutorial notebook

#### Zadanie 3.4: Comparison Plots
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```python
# scripts/plot_benchmark_comparisons.py

def plot_formose_validation(sim_results, lit_data, output):
    """
    Figure 3: Formose validation
    - A) Yields: sim vs literature
    - B) Kinetics: [glycolaldehyde](t)
    - C) Autocatalysis signature
    """
    pass

def generate_benchmark_table(all_results, output):
    """LaTeX table for paper"""
    pass
```

**Deliverables**:
- [ ] Figure 3: `figures/fig3_formose_validation.png`
- [ ] Figure 4: `figures/fig4_strecker_validation.png`
- [ ] Table 1: `paper/tables/table1_benchmarks.tex`
- [ ] Data: `data/results/benchmarks/formose_runs.csv`

**Checkpoint Week 3**:
- [ ] 3+ benchmarks reproduced? ✓ GO / ✗ NO-GO
- [ ] Yields within ±30% of literature? ✓ GO / ✗ NO-GO
- [ ] Structures confirmed? ✓ GO / ✗ NO-GO

---

### **Tydzień 4: PubChem Matcher v2** 🧬

**Status**: 🟡 PARTIAL (basic matcher exists)  
**Problem**: Obecny matcher używa heurystyk, potrzebny ML + walidacja

#### Zadanie 4.1: Atom Type Classifier (ML)
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 2 dni

```python
# matcher/ml/atom_classifier.py

from sklearn.ensemble import RandomForestClassifier

class AtomTypeClassifier:
    """ML classifier: node features → atom type"""
    
    def extract_features(self, node_info, graph):
        """
        Features:
        - degree, neighbor degrees
        - in ring, ring size
        - mass, energy
        - q_vector (6D charge)
        - valence estimate
        - electronegativity proxy
        """
        pass
    
    def train(self, training_data_path):
        """Train on 100k PubChem molecules"""
        pass
    
    def predict(self, node_info, graph, return_proba=False):
        """Predict atom type"""
        pass
```

**Training data**:
- 100k PubChem molecules (prebiotic subset)
- Features from graph topology + simulation attributes
- Target: H, C, N, O, S, P, F, Cl, Br, I

**Deliverables**:
- [ ] `matcher/ml/atom_classifier.py`
- [ ] `scripts/generate_training_data.py`
- [ ] `data/training/atom_features.pkl` (100k samples)
- [ ] Trained model: `models/atom_classifier.pkl`
- [ ] Accuracy > 85%
- [ ] Confusion matrix analysis

#### Zadanie 4.2: Multi-Metric Similarity
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```python
# matcher/similarity/multi_metric.py

class MultiMetricSimilarity:
    """Comprehensive similarity scoring"""
    
    def __init__(self, weights=None):
        self.weights = weights or {
            'topology': 0.25,      # Graph isomorphism
            'fingerprint': 0.35,   # Morgan/ECFP Tanimoto
            'energy': 0.15,        # Energy landscape
            'spectral': 0.15,      # Laplacian eigenvalues
            'geometric': 0.10      # 3D RMSD
        }
    
    def compute_similarity(self, cluster, pubchem_mol):
        """All metrics → weighted score"""
        pass
```

**Deliverables**:
- [ ] `matcher/similarity/multi_metric.py`
- [ ] 5 different similarity metrics
- [ ] Weighted combination
- [ ] Confidence scoring

#### Zadanie 4.3: Validation
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

```python
# matcher/validation/confidence.py

class MatchConfidenceEvaluator:
    """Evaluate match confidence"""
    
    def evaluate_match(self, cluster, match_result):
        """
        Returns:
        - confidence_score (0-1)
        - reliability ('high'/'medium'/'low')
        - warnings
        - validation_status
        """
        pass
    
    def is_chemically_plausible(self, cluster):
        """Check for valence violations, etc"""
        pass
```

**Validation set**: 50 clusters from known reactions

**Deliverables**:
- [ ] `matcher/validation/confidence.py`
- [ ] Validation report generator
- [ ] Accuracy on validation set > 80%

#### Zadanie 4.4: Integration
**Status**: ⏳ DO ZROBIENIA  
**Czas**: 1 dzień

**Deliverables**:
- [ ] `matcher/matcher_v2.py` (refactored)
- [ ] Updated documentation
- [ ] API for frontend
- [ ] Batch matching script

**Checkpoint Week 4**:
- [ ] Classifier accuracy > 85%? ✓ GO / ✗ NO-GO
- [ ] Match accuracy > 80%? ✓ GO / ✗ NO-GO
- [ ] Confidence scoring works? ✓ GO / ✗ NO-GO

---

## 🎓 PHASE 2: OPEN-ENDED EXPERIMENTS (Weeks 5-6 - FUTURE)

**Status**: 🔮 PRZYSZŁOŚĆ (po Phase 1)

### Cel: Generowanie nowych wyników dla publikacji

#### Scenariusze:
1. **Miller-Urey conditions**
   - CH₄, NH₃, H₂O, H₂
   - Electrical discharge (energy pulses)
   - 10⁷ steps × 10 runs

2. **Hydrothermal vent**
   - Alkaline (pH 9-11)
   - 50-150°C
   - H₂, H₂S, Fe²⁺

3. **Formamide-rich**
   - HCONH₂ as solvent
   - UV radiation
   - Mineral catalysts

#### Deliverables:
- [ ] 30 independent simulations
- [ ] Novel molecules catalog (100+)
- [ ] Top 20 PubChem matches
- [ ] DFT validation (top 5)

---

## 📊 PHASE 3: ANALYSIS & PAPER (Weeks 7-12 - FUTURE)

**Status**: 🔮 PRZYSZŁOŚĆ

### Week 7-8: Analysis & Visualization
- [ ] Reaction network construction
- [ ] Autocatalytic cycles (12+)
- [ ] Statistical analysis
- [ ] All paper figures (7+)
- [ ] Supplementary videos

### Week 9-10: Writing
- [ ] Introduction
- [ ] Methods
- [ ] Results
- [ ] Discussion
- [ ] Abstract + Conclusions
- [ ] Supplementary materials

### Week 11-12: Submission
- [ ] Internal review
- [ ] Revisions
- [ ] LaTeX formatting
- [ ] Submit to JCTC or Origins of Life
- [ ] ArXiv preprint
- [ ] GitHub release

---

## 📊 Progress Tracking

### Thermodynamic Validation
- [x] Energy conservation active
- [x] M-B distribution validated
- [x] Entropy monitoring active
- [ ] Extended tests (virial, etc)
- [ ] Figures generated
- [ ] Documentation complete

**Score**: 3/6 (50%)

### Parameter Database
- [x] Schema defined ✅
- [x] Infrastructure complete ✅
- [ ] Bond parameters (1/50+) - IN PROGRESS
- [ ] VDW parameters (1/10) - IN PROGRESS
- [ ] Citations collected (2/30+) - IN PROGRESS
- [ ] Integration complete
- [ ] LaTeX table generated

**Score**: 2/7 (29%) 📈 IMPROVING

### Benchmark Reactions
- [ ] Test framework
- [ ] Formose reaction
- [ ] Strecker synthesis
- [ ] HCN polymerization
- [ ] Phosphorylation
- [ ] Comparison plots

**Score**: 0/6 (0%) ⚠️ CRITICAL

### PubChem Matcher
- [x] Basic matcher exists
- [ ] ML classifier (0%)
- [ ] Multi-metric similarity
- [ ] Validation (0%)
- [ ] Integration
- [ ] Documentation

**Score**: 1/6 (17%)

### Overall Progress
**Phase 0 (Foundations)**: ✅ 100% DONE  
**Phase 1 (Validation Sprint)**: ⏳ 12% (Week 2 Day 1 complete)  
**Phase 2 (Experiments)**: ⏳ 0%  
**Phase 3 (Paper)**: ⏳ 0%  

**Total**: ✅ Phase 0 Complete, 📋 Phases 1-3 Planned

---

## 🎯 Critical Success Factors

### Must Have (dla publikacji):
- ✅ Thermodynamic validation (DONE)
- ⚠️ Literature parameters (TODO)
- ⚠️ Benchmark reactions (TODO)
- ✅ Statistical rigor framework (DONE)
- ✅ Open source (DONE)

### Nice to Have:
- 🎁 Experimental collaboration
- 🎁 Novel testable predictions
- 🎁 Video supplementary
- 🎁 Interactive demo

### Deal Breakers:
- ❌ Arbitrary parameters ← **MUST FIX (Week 2)**
- ❌ No validation vs real chemistry ← **MUST FIX (Week 3)**
- ❌ Thermodynamic violations ← **FIXED ✅**

---

## 📅 Timeline Summary

| Phase | Duration | Status | Priority |
|-------|----------|--------|----------|
| **Phase 0: Foundations** | Completed | ✅ DONE | - |
| **Phase 1: Validation Sprint** | 4 weeks | 📋 TODO | 🔴 HIGH |
| **Phase 2: Experiments** | 2 weeks | 🔮 FUTURE | 🟡 MEDIUM |
| **Phase 3: Paper** | 6 weeks | 🔮 FUTURE | 🟢 LOW |
| **Total to Publication** | 12 weeks | - | - |

---

## 🚀 Next Actions (Immediate)

### This Week:
1. [x] ✅ PhysicsDatabase infrastructure (Day 1 complete)
2. [ ] 🔄 Literature data collection (Days 2-3)
3. [ ] Integration with potentials.py (Day 4)

### Previously (Week 1):
1. [~] Extended thermodynamic tests (DEFERRED)
2. [~] Generate Figure 1 & 2 (DEFERRED)
3. [~] Update `THERMODYNAMIC_VALIDATION.md` (DEFERRED)

### This Month:
1. [ ] Complete Week 1-4 sprint
2. [ ] All validations passing
3. [ ] Ready for open-ended experiments

---

## 📚 References

### Completed Work:
- `docs/SCIENTIFIC_INTEGRITY_VERIFICATION.md` - Status po optymalizacjach
- `docs/MEMORY_PERFORMANCE_FIX.md` - Naprawa freezu
- `docs/RUNTIME_TIMER_FIX.md` - Timer czasu rzeczywistego
- `docs/MATCHER_BUTTON_FIX.md` - UI fixes

### Future Work:
- `docs/Live 2-plan walidacji naukowej.md` - Pełny plan (2805 linii)
- Paper target: JCTC lub Origins of Life and Evolution of Biospheres

---

## ✅ Checklist przed Publikacją

### Code
- [x] Tests passing
- [x] Code formatted
- [x] Documentation (partial)
- [ ] All validators implemented
- [ ] Parameter DB integrated
- [ ] Benchmark tests passing

### Data
- [ ] Raw data archived
- [ ] Parameters with citations
- [ ] Benchmark results
- [ ] Supplementary data

### Paper
- [ ] All figures (7+)
- [ ] All tables (5+)
- [ ] References complete
- [ ] Supplementary materials
- [ ] Word count < 6000

---

*Ten dokument jest żywym roadmapem. Aktualizuj po każdym milestone.*

**Last updated**: 12 października 2025 (21:30 CET)  
**Latest**: ✅ Week 2 Day 1 Complete - Physics Database Infrastructure  
**Next**: 📋 Days 2-3 - Literature Data Collection (50+ parameters)

