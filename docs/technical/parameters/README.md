# Parametry Naukowe

Dokumentacja dotycząca parametrów fizycznych i chemicznych w symulacji.

---

## 📄 Dokumenty

### [SCIENTIFIC_PARAMETERS_ANALYSIS.md](SCIENTIFIC_PARAMETERS_ANALYSIS.md)
Szczegółowa analiza parametrów naukowych używanych w symulacji:
- Parametry Lennard-Jones
- Parametry wiązań (Morse potential)
- Van der Waals radii
- Literaturowe źródła

### [QUICK_PARAMETER_UPDATE.md](QUICK_PARAMETER_UPDATE.md)
Szybki przewodnik aktualizacji parametrów:
- Jak zmienić parametry binding
- Jak dostosować siłę wiązań
- Jak zmodyfikować energia pulsów

---

## 🔬 Źródła Parametrów

### Physics Database (`data/physics_parameters.json`)
Centralna baza parametrów z cytacjami:
- **UFF Force Field** (Rappé et al. 1992)
- **Bond Energies** (Luo 2007)
- **Van der Waals** (Bondi 1964)

### Używane w Kodzie:
- `backend/sim/config.py` - domyślne parametry
- `backend/sim/core/binding.py` - parametry wiązań
- `backend/sim/core/potentials.py` - potencjały
- `backend/sim/core/spatial_hash.py` - LJ parameters

---

## 📊 Najważniejsze Parametry

### Binding:
- `binding_threshold`: 0.25 (frontend), 0.7 (backend default)
- `unbinding_threshold`: 0.15
- Zasięg wiązania: 3.25 Å (6.5× radius)

### Bond Types:
- vdW: k=2.0, strength=5.0
- Covalent: k=500.0, strength=100.0
- H-bond: k=50.0, strength=30.0
- Metallic: k=100.0, strength=50.0

### Lennard-Jones:
- σ (sigma): 3.4 Å (Carbon, UFF)
- ε (epsilon): 0.5 kJ/mol

---

## 🔗 Zobacz Też

- [Session 2024-10-22](../../sessions/2024-10-22-validation-parameters/) - Najnowsze zmiany parametrów
- [Physics Database](../../PHYSICS_DATABASE.md) - Dokumentacja bazy parametrów
- [Scientific Overview](../../SCIENTIFIC_OVERVIEW.md) - Przegląd naukowy

