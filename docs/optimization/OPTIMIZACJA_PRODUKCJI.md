# Optymalizacja Wydajności Produkcji

## 🔍 Analiza Problemów

### Obecna Sytuacja:
- **Step wykonuje się szybko**: 6.0ms (dobrze!)
- **Ale między krokami jest opóźnienie**: ~6 sekund między step 1000 a 1100
- **Używa**: `SimulationStepper` (nie Hybrid)
- **Backend**: CPU (dobrze dla chemii)

### Główne Problemy:

1. **`get_visualization_data()` jest wolne** (402ms z testów)
   - Bonds/Clusters: 339ms (84% czasu!)
   - Wywoływane co 0.1s (10 FPS) w broadcast loop
   - Blokuje symulację

2. **Nie używa Hybrid Mode**
   - Chemia wykonywana synchronicznie
   - Blokuje symulację podczas analizy

3. **Broadcast loop może być zbyt częsty**
   - 10 FPS może być za dużo dla wolnej wizualizacji

---

## ✅ ROZWIĄZANIA

### 1. Użyj HybridSimulationStepper (NAJWAŻNIEJSZE!)

**Zmiana w `backend/api/server.py`:**

```python
# Zamiast:
from sim.core.stepper import SimulationStepper

# Użyj:
from sim.core.hybrid_stepper import HybridSimulationStepper

# W create_simulation:
simulation = HybridSimulationStepper(config)  # Zamiast SimulationStepper(config)
```

**Zalety:**
- ✅ Chemia wykonywana w tle (nie blokuje symulacji)
- ✅ Symulacja może działać szybciej
- ✅ Pełna analiza chemiczna bez spowalniania

### 2. Zwiększ Interwał Broadcast

**W `broadcast_loop` (linia ~906):**

```python
# Zamiast:
await asyncio.sleep(0.1)  # 10 FPS

# Użyj:
await asyncio.sleep(0.2)  # 5 FPS - mniej obciążenia
# LUB
await asyncio.sleep(0.5)  # 2 FPS - dla bardzo wolnej wizualizacji
```

**Zalety:**
- ✅ Mniej wywołań `get_visualization_data()`
- ✅ Mniej obciążenia CPU
- ✅ Szybsza symulacja

### 3. Zoptymalizuj Konfigurację

**W `create_simulation` dodaj optymalizacje:**

```python
config.chemistry_snapshot_interval = 200  # Zamiast 100 - rzadsze analizy
config.metrics_update_interval = 1000  # Rzadsze aktualizacje metryk
config.enable_diagnostics = False  # Wyłącz jeśli nie potrzebujesz
```

### 4. Cache Wizualizacji

**Dodaj cache dla wolnych operacji:**

```python
# W broadcast_loop, przed get_visualization_data:
last_viz_step = getattr(simulation, '_last_viz_step', -1)
if simulation.step_count - last_viz_step < 10:  # Cache przez 10 kroków
    data = getattr(simulation, '_cached_viz_data', None)
    if data:
        # Użyj cache
        continue

# Po get_visualization_data:
simulation._last_viz_step = simulation.step_count
simulation._cached_viz_data = data
```

---

## 🚀 IMPLEMENTACJA

### Krok 1: Zmień na HybridSimulationStepper

```python
# backend/api/server.py, linia ~95
from sim.core.hybrid_stepper import HybridSimulationStepper

# Linia ~325
simulation = HybridSimulationStepper(config)  # Zmiana tutaj
```

### Krok 2: Zwiększ Interwał Broadcast

```python
# backend/api/server.py, linia ~906
await asyncio.sleep(0.2)  # 5 FPS zamiast 10 FPS
```

### Krok 3: Dodaj Optymalizacje Konfiguracji

```python
# backend/api/server.py, w create_simulation (po linii ~320)
if request.mode == "open_chemistry":
    config = SimulationConfig(**request.config)
    # Optymalizacje
    config.chemistry_snapshot_interval = getattr(config, 'chemistry_snapshot_interval', 200)
    config.metrics_update_interval = getattr(config, 'metrics_update_interval', 1000)
```

---

## 📊 Oczekiwane Ulepszenia

### Przed:
- Step: 6ms ✅
- Opóźnienie między krokami: ~6s ❌
- Wizualizacja: 402ms (339ms Bonds/Clusters) ❌

### Po:
- Step: 6ms ✅ (bez zmian)
- Opóźnienie: <100ms ✅ (chemia w tle)
- Wizualizacja: rzadsze wywołania (5 FPS zamiast 10) ✅

**Szacowany speedup: 10-60x dla symulacji!**

---

## 🔧 Szybka Naprawa (Minimalne Zmiany)

Jeśli chcesz szybką poprawę bez dużych zmian:

1. **Zwiększ interwał broadcast:**
   ```python
   await asyncio.sleep(0.5)  # 2 FPS zamiast 10 FPS
   ```

2. **Dodaj warunek dla wolnej wizualizacji:**
   ```python
   # Tylko wywołuj get_visualization_data jeśli minęło >0.5s
   last_viz_time = getattr(simulation, '_last_viz_time', 0)
   if time.time() - last_viz_time < 0.5:
       continue
   simulation._last_viz_time = time.time()
   ```

---

## ✅ REKOMENDACJA FINALNA

**Najlepsze rozwiązanie:**
1. ✅ Użyj `HybridSimulationStepper` 
2. ✅ Zwiększ interwał broadcast do 0.2-0.5s
3. ✅ Zoptymalizuj konfigurację (większe interwały)

**To da największy speedup przy minimalnych zmianach!**

