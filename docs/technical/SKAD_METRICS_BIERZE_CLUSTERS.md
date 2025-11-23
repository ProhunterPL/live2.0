# Skąd Metrics Bierze Liczbę Clusters?

## 🔍 Analiza

### Endpoint Metrics:
```
GET /simulation/{simulation_id}/metrics
```

**Kod:** `backend/api/server.py` linia 491-500
```python
metrics = simulation.aggregator.get_aggregated_stats()
return {"metrics": metrics}
```

### Ścieżka Danych:

1. **`aggregator.get_aggregated_stats()`** → zwraca `self.aggregated_stats`
2. **`update_aggregated_stats()`** → aktualizuje stats z `metrics_collector.get_current_metrics()`
3. **`get_current_metrics()`** → zwraca `cluster_count` z `self.cluster_count[None]`
4. **`update_metrics()`** → ustawia `self.metrics.cluster_count[None]`

---

## ⚠️ PROBLEM: Uproszczone Przybliżenie!

### W `SimulationStepper` (`backend/sim/core/stepper.py` linia 1066-1094):

```python
# Update cluster metrics manually - OPTIMIZED with simplified algorithm
cluster_count = 0
try:
    max_check = min(particle_count, 1000)
    if max_check > 0:
        bond_active = self.binding.bond_active.to_numpy()[:max_check, :max_check]
        particles_active = self.particles.active.to_numpy()[:max_check]
        
        # Count particles that have at least one bond
        particles_with_bonds = np.sum(particles_active * np.any(bond_active, axis=1))
        
        # Simple approximation: particles_with_bonds / 2 (assuming average 2 particles per cluster)
        cluster_count = max(1, int(particles_with_bonds / 2))  # ⚠️ PRZYBLIŻENIE!
    else:
        cluster_count = 1
    
    self.metrics.cluster_count[None] = cluster_count
```

**To jest tylko przybliżenie!** Zakłada średnio 2 cząstki na klaster.

### Jeśli masz 499.0 clusters:
- `particles_with_bonds ≈ 998` (499 * 2)
- **To nie jest rzeczywista liczba klastrów!**

---

## ✅ Rozwiązanie: HybridSimulationStepper

### W `HybridSimulationStepper` (`backend/sim/core/hybrid_stepper.py`):

**CPU Worker używa NetworkX** (linia 233-254):
```python
def _detect_clusters_cpu(self, bonds: List[tuple], n_particles: int) -> List[List[int]]:
    """CPU-based cluster detection using NetworkX"""
    if not bonds:
        return []
    
    # Build graph
    G = nx.Graph()
    for i, j, strength in bonds:
        G.add_edge(i, j, weight=strength)
    
    # Find connected components (clusters) - PRAWDZIWA ANALIZA!
    clusters = []
    for component in nx.connected_components(G):
        cluster = list(component)
        if len(cluster) >= self.config.min_cluster_size:
            clusters.append(cluster)
    
    return clusters
```

**Aktualizacja metrics** (linia 427-428):
```python
# Update cluster count
if 'n_clusters' in cpu_metrics:
    self.metrics.cluster_count[None] = int(cpu_metrics['n_clusters'])  # ✅ PRAWDZIWA LICZBA!
```

**To jest prawdziwa liczba klastrów** z analizy grafowej (connected components)!

---

## 📊 Porównanie

| Metoda | Dokładność | Wydajność | Używane przez |
|--------|------------|-----------|---------------|
| **SimulationStepper** | ❌ Przybliżenie (`particles_with_bonds / 2`) | ✅ Szybkie | `SimulationStepper` |
| **HybridSimulationStepper** | ✅ Prawdziwa (NetworkX connected components) | ✅ Szybkie (CPU) | `HybridSimulationStepper` |
| **BindingSystem.get_cluster_stats()** | ✅ Prawdziwa | ⚠️ Wolniejsze | Używane w `get_simulation_state()` |

---

## 🔧 Jak Naprawić (jeśli używasz SimulationStepper)

### Opcja 1: Użyj HybridSimulationStepper (REKOMENDOWANE)
✅ Już zaimplementowane w `server.py` - używa HybridSimulationStepper!

### Opcja 2: Napraw SimulationStepper
Zmień w `backend/sim/core/stepper.py` linia 1083:

```python
# Zamiast:
cluster_count = max(1, int(particles_with_bonds / 2))

# Użyj:
cluster_stats = self.binding.get_cluster_stats()
cluster_count = cluster_stats['num_clusters']
```

**Uwaga:** To będzie wolniejsze, bo `get_cluster_stats()` wykonuje pełną analizę grafową.

---

## ✅ Sprawdzenie

### Jeśli używasz HybridSimulationStepper:
- ✅ Powinieneś mieć **prawdziwą** liczbę klastrów
- ✅ Sprawdź logi: `"Using HybridSimulationStepper for sim_..."`
- ✅ Liczba powinna być dokładna (z NetworkX)

### Jeśli używasz SimulationStepper:
- ⚠️ Masz **przybliżenie** (`particles_with_bonds / 2`)
- ⚠️ Liczba może być niedokładna
- ⚠️ Rozważ przejście na HybridSimulationStepper

---

## 📝 Podsumowanie

**Skąd metrics bierze liczbę clusters:**

1. **SimulationStepper**: Uproszczone przybliżenie `particles_with_bonds / 2` ❌
2. **HybridSimulationStepper**: Prawdziwa liczba z NetworkX connected components ✅
3. **BindingSystem.get_cluster_stats()**: Prawdziwa liczba, ale wolniejsza ⚠️

**Twoja wartość 499.0:**
- Jeśli używasz HybridSimulationStepper → ✅ Prawdziwa liczba
- Jeśli używasz SimulationStepper → ❌ Przybliżenie (rzeczywista liczba może być inna)

---

*Sprawdź logi serwera aby zobaczyć który stepper jest używany!*

