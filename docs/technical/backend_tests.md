# Live 2.0 - Dokumentacja Testów

## 📁 Struktura Testów

Wszystkie testy są skoncentrowane w katalogu `backend/tests/` zgodnie z najlepszymi praktykami Python.

### 🧪 Testy Jednostkowe (Unit Tests)

#### `test_units.py`
- **TestMolecularGraph**: Testy funkcjonalności grafów molekularnych
  - `test_graph_creation` - Tworzenie grafów
  - `test_graph_hash_stability` - Stabilność hash grafu
  - `test_graph_complexity` - Obliczanie złożoności
- **TestSubstanceCatalog**: Testy katalogu substancji
  - `test_catalog_basic_operations` - Podstawowe operacje
  - `test_catalog_statistics` - Statystyki katalogu
- **TestConfiguration**: Testy konfiguracji
  - `test_simulation_config_defaults` - Domyślne ustawienia
  - `test_preset_config_defaults` - Konfiguracja preset
- **TestGraphCatalog**: Testy katalogu grafów
  - `test_graph_catalog_basic` - Podstawowe operacje

#### `test_core.py`
- **TestSimulationConfig**: Testy konfiguracji symulacji
- **TestRNG**: Testy generatora liczb losowych
- **TestGrid**: Testy siatki 2D
- **TestParticleSystem**: Testy systemu cząstek
- **TestPotentialSystem**: Testy systemu potencjałów
- **TestBindingSystem**: Testy systemu wiązań
- **TestMolecularGraph**: Testy grafów molekularnych
- **TestSubstanceCatalog**: Testy katalogu substancji
- **TestMetricsCollector**: Testy kolektora metryk
- **TestNoveltyTracker**: Testy trackera nowości
- **TestComplexityAnalyzer**: Testy analizatora złożoności
- **TestEnergySystem**: Testy systemu energii
- **TestPresetPrebioticSimulator**: Testy symulatora preset

### 🔄 Testy Property-Based

#### `test_property_based.py`
- **TestInvariants**: Testy inwariantów
  - `test_energy_conservation_no_impulses` - Zachowanie energii
  - `test_particle_count_conservation` - Zachowanie liczby cząstek
  - `test_mass_conservation` - Zachowanie masy
- **TestLocality**: Testy lokalności
  - `test_force_cutoff_distance` - Odcięcie sił
  - `test_neighbor_list_locality` - Lokalność list sąsiadów
- **TestNumericalStability**: Testy stabilności numerycznej
  - `test_adaptive_timestep_bounds` - Granice adaptacyjnego kroku
  - `test_velocity_clamping` - Ograniczenie prędkości
  - `test_particle_position_bounds` - Granice pozycji cząstek
- **TestGraphInvariants**: Testy inwariantów grafów
  - `test_graph_hash_stability` - Stabilność hash
  - `test_catalog_determinism` - Determinizm katalogu

### 💾 Testy Snapshotów

#### `test_snapshots.py`
- **TestSnapshotManager**: Testy zarządzania snapshotami
  - `test_snapshot_creation` - Tworzenie snapshotów
  - `test_snapshot_loading` - Ładowanie snapshotów
  - `test_snapshot_validation` - Walidacja snapshotów
  - `test_snapshot_deletion` - Usuwanie snapshotów
- **TestSnapshotAPI**: Testy API snapshotów
  - `test_snapshot_api_operations` - Operacje API

### 🌐 Testy API

#### `test_api.py`
- **TestSimulationAPI**: Testy API symulacji
  - `test_create_simulation_open_chemistry` - Tworzenie symulacji OC
  - `test_create_simulation_preset_prebiotic` - Tworzenie symulacji PP
  - `test_get_simulation_status` - Status symulacji
  - `test_start_simulation` - Uruchamianie symulacji
  - `test_pause_simulation` - Pauzowanie symulacji
  - `test_resume_simulation` - Wznawianie symulacji
  - `test_stop_simulation` - Zatrzymywanie symulacji
  - `test_reset_simulation` - Resetowanie symulacji
  - `test_get_novel_substances` - Pobieranie nowych substancji
  - `test_get_metrics` - Pobieranie metryk
  - `test_save_snapshot` - Zapisywanie snapshotów
  - `test_load_snapshot` - Ładowanie snapshotów
- **TestWebSocketAPI**: Testy WebSocket
  - `test_websocket_connection` - Połączenie WebSocket
  - `test_websocket_invalid_simulation` - Nieprawidłowa symulacja
  - `test_websocket_data_format` - Format danych
- **TestErrorHandling**: Testy obsługi błędów
- **TestConcurrentSimulations**: Testy współbieżnych symulacji

### ⚡ Testy Wydajności

#### `test_performance.py`
- **TestPerformanceRequirements**: Testy wymagań wydajności
  - `test_particle_creation_performance` - Wydajność tworzenia cząstek
  - `test_binding_system_performance` - Wydajność systemu wiązań
  - `test_graph_computation_performance` - Wydajność obliczeń grafów
  - `test_catalog_performance` - Wydajność katalogu

#### `test_performance_integration.py`
- Testy wydajności integracyjne (10,000 cząstek @ 60 FPS)

#### `phase_0_performance_test.py`
- Testy wydajności Phase 0

### 🔄 Testy Długiego Biegu

#### `test_stability_24h.py`
- Testy stabilności 24-godzinnej
- Testy zachowania energii i stabilności numerycznej

### 🔗 Testy Łączności

#### `simple_connectivity_test.py`
- Testy podstawowej łączności systemu

### 📊 Testy Snapshotów z Obrazami

#### `test_snapshot_with_images.py`
- Testy generowania snapshotów z obrazami wizualizacji

### 🛠️ Narzędzia Testowe

#### `profile_simulation.py`
- Profilowanie wydajności symulacji

#### `conftest.py`
- Konfiguracja pytest z inicjalizacją Taichi
- Fixtures dla testów

## 🚀 Uruchamianie Testów

### Wszystkie testy
```bash
cd backend
python -m pytest tests/ -v
```

### Testy jednostkowe
```bash
python -m pytest tests/test_units.py -v
```

### Testy snapshotów
```bash
python -m pytest tests/test_snapshots.py -v
```

### Testy property-based
```bash
python -m pytest tests/test_property_based.py -v
```

### Testy API (wymagają uruchomionego serwera)
```bash
python -m pytest tests/test_api.py -v
```

### Testy wydajności
```bash
python -m pytest tests/test_performance.py -v
```

## 📋 Status Testów

### ✅ Przechodzące (15/15)
- Testy jednostkowe: 8/8
- Testy snapshotów: 5/5
- Testy property-based (grafy): 2/2

### ⚠️ Wymagające Implementacji
- Testy property-based (inwarianty fizyczne)
- Testy wydajności (N=200k cząstek)
- Testy długiego biegu (8-24h)
- Testy API (wymagają serwera)

## 🔧 Konfiguracja

### Taichi Initialization
Testy używają `conftest.py` do inicjalizacji Taichi z architekturą CPU dla stabilności.

### Mocki
Testy używają mocków dla:
- Obiektów symulacji w testach snapshotów
- API calls w testach jednostkowych
- WebSocket connections

## 📈 Metryki Testów

- **Pokrycie**: Podstawowe funkcjonalności przetestowane
- **Stabilność**: Testy przechodzą deterministycznie
- **Wydajność**: Testy jednostkowe < 1s, integracyjne < 10s
- **Maintainability**: Struktura modularna, łatwa rozbudowa

## 🎯 Kryteria Akceptacji v1

### ✅ Spełnione
- System działa w Trybie B (testy jednostkowe przechodzą)
- Snapshot/restore odtwarza stan (testy przechodzą)
- Testy bazowe przechodzą (15/15)

### ⚠️ Wymagające Weryfikacji
- Generuje nowe substancje (wymaga testów długiego biegu)
- Novelty>0 w dłuższym biegu (wymaga testów długiego biegu)
- Wydajność akceptowalna na GPU (wymaga testów wydajności)

## 📚 Dokumentacja Powiązana

- `TEST_SUMMARY.md` - Szczegółowe podsumowanie zgodnie z planem v1
- `pytest.ini` - Konfiguracja pytest
- `conftest.py` - Fixtures i konfiguracja testów

---
*Dokumentacja wygenerowana automatycznie - aktualizacja: 2024*
