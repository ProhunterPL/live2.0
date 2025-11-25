# SYSTEM CRASH REPORT
# ====================
**Data:** 15 października 2025, 21:27  
**Severity:** CRITICAL ⚠️

## 🔴 INCIDENT

Test `phase2_quick_fix_test.yaml` spowodował **restart systemu**.

### Konfiguracja która spowodowała crash:
- **Particles:** 500 (1775 atoms estimated)
- **GPU:** CUDA (RTX 5070)
- **Test:** 10,000 kroków
- **Physics:** Reactions disabled, ale forces ENABLED

### Prawdopodobna przyczyna:

**GPU Driver Hang/Crash:**
1. CUDA kernel compilation dla force computation
2. O(n²) algorithm z 1775 atoms = ~1,574,000 par
3. GPU driver timeout (TDR - Timeout Detection and Recovery)
4. Windows forced system reset

## 🔍 ANALIZA

### Problem #1: GPU Watchdog Timer
Windows ma "TDR" (Timeout Detection and Recovery):
- Jeśli GPU kernel nie odpowiada przez 2-5 sekund
- System uznaje że GPU się zawiesił
- **Resetuje driver (lub cały system)**

### Problem #2: Kernel Compilation Timeout
CUDA kompiluje kernele przy pierwszym użyciu:
- `compute_forces_kernel` dla 1775 atoms jest OGROMNY
- Kompilacja może zająć 10-30+ sekund
- Windows TDR timeout = 2-5 sekund
- **Crash podczas kompilacji!**

### Problem #3: O(n²) na GPU
- 1775 atoms = 1,574,000 par
- Każda para: distance, LJ force, Coulomb force, etc.
- GPU próbuje wykonać to w jednym kernel call
- **Zbyt duże obciążenie → hang → crash**

## ✅ ROZWIĄZANIE

### Natychmiastowe działania:

1. **WYŁĄCZ GPU** - używaj TYLKO CPU:
```python
ti.init(arch=ti.cpu, cpu_max_num_threads=28)
```

2. **DRASTYCZNIE zmniejsz particles** do 50-100:
```yaml
n_particles: 50  # Not 500!
```

3. **Krótkie testy** (100 steps max):
```yaml
max_steps: 100
```

### Długoterminowe rozwiązanie:

**GPU jest NIEUŻYWALNE bez spatial hashing!**

Obecny O(n²) algorithm:
- ✅ Działa na CPU (wolno ale bezpiecznie)
- ❌ **Crashuje GPU** (kernels zbyt duże)

Po implementacji spatial hashing (O(n)):
- ✅ Działa na CPU (szybko)
- ✅ **Będzie działać na GPU** (małe kernels)

## 📊 BEZPIECZNE LIMITY

### CPU (bezpieczne):
- Max particles: **2000-5000** (wolno ale działa)
- Algorytm: O(n²) akceptowalny
- Crash risk: **NISKI**

### GPU bez spatial hashing (NIEBEZPIECZNE):
- Max particles: **<100** (powyżej = crash risk)
- Algorytm: O(n²) **NIE działa na GPU!**
- Crash risk: **WYSOKI** ⚠️

### GPU ze spatial hashing (docelowe):
- Max particles: **10,000-100,000** (bardzo szybko)
- Algorytm: O(n) = GPU friendly
- Crash risk: **NISKI**

## 🎯 REKOMENDACJE

### NA TERAZ:

1. ✅ **Używaj TYLKO CPU**
2. ✅ **Max 100 particles** dla testów
3. ✅ **Max 1000 steps** dla testów
4. ❌ **NIE używaj GPU** dopóki nie będzie spatial hashing

### PLAN:

**Etap 1: Bezpieczne testy (CPU only)**
```yaml
simulation:
  n_particles: 100
  max_steps: 1000
physics:
  enable_reactions: false
```
**Expected:** 10-20 steps/s (BEZPIECZNIE)

**Etap 2: Implementuj Spatial Hashing (4h pracy)**
- Zmień O(n²) → O(n)
- Test na CPU z 2000 particles
- Sprawdź czy działa

**Etap 3: Re-enable GPU (po spatial hashing)**
- GPU będzie bezpieczne z O(n)
- Expected: 500-2000 steps/s
- **Brak crash risk**

## ⚠️ OSTRZEŻENIA

### NIE RÓB:
- ❌ GPU + >100 particles
- ❌ GPU + O(n²) forces
- ❌ CUDA bez spatial hashing
- ❌ Tests >1000 steps bez sprawdzenia

### RÓB:
- ✅ CPU dla wszystkich testów
- ✅ <100 particles dla GPU tests
- ✅ Implementuj spatial hashing ASAP
- ✅ Krótkie testy (100-1000 steps)

---

**Wniosek:** GPU + O(n²) algorithm = **SYSTEM CRASH**. 

Musisz:
1. Używać CPU (bezpiecznie)
2. Albo zaimplementować spatial hashing (4h)
3. **Potem** możesz używać GPU bezpiecznie


