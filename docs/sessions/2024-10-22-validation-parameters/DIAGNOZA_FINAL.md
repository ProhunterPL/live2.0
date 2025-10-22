# 🔴 DIAGNOZA: CRASH NA KROKU 63000

## GŁÓWNA PRZYCZYNA: dt=0.035 → 25% ENERGY DRIFT!

### Frontend wysyła w `App.tsx:182`:
```typescript
dt: 0.035,  // Zwiększone z 0.01 ❌ TO JEST PROBLEM!
```

### Co się dzieje:
1. **dt=0.035 to 7× więcej niż default (0.005)**
2. **Każdy krok generuje ogromne błędy numeryczne**
3. **Po 63000 krokach błąd kumulatywny = 25% drift**
4. **Python runtime/Windows zabija proces (out of bounds, NaN, inf)**

---

## 🎯 NATYCHMIASTOWA NAPRAWA

### ZMIANA 1: Napraw timestep w frontend
```typescript
// frontend/src/App.tsx linia 182
dt: 0.005,  // NAPRAWIONE: było 0.035 (7× za dużo!)
```

### ZMIANA 2: Wyłącz kosztowną walidację
```typescript
// frontend/src/App.tsx - dodaj po linii 210
enable_thermodynamic_validation: false,  // Wyłącz dla produkcji
validate_every_n_steps: 10000,  // Waliduj rzadziej
```

---

## 🔬 PROBLEM 2: MAŁE KLASTRY

### Obecna config frontend (App.tsx:198-199):
```typescript
binding_threshold: 0.25,      // Za niski! (default: 0.7)
unbinding_threshold: 0.15,    // OK dla stability
pulse_amplitude: 8.0,         // Wysoka energia - OK!
```

### Problem jest w backend/sim/core/binding.py:
- **Zasięg wiązania: 2× radius** (linia 310) - za mały!
- **Próg probability: 0.6** (linia 315) - za wysoki!
- **Siła wiązań kowalencyjnych: 20** (linia 519) - za słaba!

---

## ✅ PLAN NAPRAWY

### KROK 1: Frontend (natychmiastowe)
Edytuj `frontend/src/App.tsx`:
```typescript
// Linia 182
dt: 0.005,  // było: 0.035

// Po linii 210 dodaj:
enable_thermodynamic_validation: false,
validate_every_n_steps: 10000,
```

### KROK 2: Backend (większe klastry)
Edytuj `backend/sim/core/binding.py`:

**A. Linia 310 - zwiększ zasięg:**
```python
if r <= PARTICLE_RADIUS_COMPILE * 3.5:  # było: 2.0
```

**B. Linia 315 - zmniejsz próg:**
```python
if binding_probability > 0.35:  # było: 0.6
```

**C. Linia 519 - zwiększ siłę wiązań:**
```python
1: {'k_spring': 50.0, 'rest_len': 0.8, 'strength': 100.0},  # było: k=10, str=20
```

**D. Linia 329 - zmniejsz mass_ratio:**
```python
if mass_ratio > 0.5:  # było: 0.7 - pozwala O-H bonds!
```

---

## 📊 OCZEKIWANE REZULTATY

### Po naprawie dt:
- ✅ Symulacja stabilna przez >200k steps
- ✅ Energy drift <5%
- ✅ Brak crash

### Po naprawie binding:
- ✅ Klastry 8-15 atomów
- ✅ Stabilne cząsteczki organiczne
- ✅ NH3, H2O, CH4, formaldehyd, HCN

---

## 🚀 CZY MAM WPROWADZIĆ TE ZMIANY?

Powiedz:
1. **"TAK" - wprowadź wszystkie zmiany**
2. **"TYLKO FRONTEND" - naprawa crash (stabilność)**
3. **"TYLKO BACKEND" - większe klastry**
4. **"CUSTOM" - wybiorę co zmienić**

