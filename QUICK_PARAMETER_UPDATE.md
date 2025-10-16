# 🚀 Szybki Przewodnik - Nowe Parametry Naukowe

## ✅ Co Zostało Zmienione?

Parametry symulacji zostały **zaktualizowane na podstawie literatury naukowej** (Miller-Urey 1953, energie wiązań chemicznych):

```diff
# backend/sim/config.py

- binding_threshold: 0.6      →  0.45  (łatwiejsze wiązanie)
- theta_break: 1.0             →  1.5   (stabilniejsze klastry)  
- pulse_amplitude: 2.5         →  1.8   (łagodniejsze wyładowania)
- pulse_every: 50              →  100   (rzadsze impulsy)
```

---

## 🎯 Oczekiwane Efekty

Po restarcie backendu i uruchomieniu nowej symulacji:

### Przed (stare parametry):
- ❌ Max cluster size: **3 atomy**
- ❌ Novelty rate: **0**
- ❌ Buttony matchera: **niewidoczne**
- ❌ Struktury H-H-H: **chemicznie niemożliwe**

### Po (nowe parametry):
- ✅ Cluster size: **4-10 atomów** (realistyczne!)
- ✅ Novelty rate: **> 0.1** (aktywna ewolucja)
- ✅ Buttony matchera: **widoczne od ~500 kroków**
- ✅ Struktury H-O-H: **chemicznie poprawne** (woda!)

---

## 📝 Instrukcja Krok po Kroku

### 1️⃣ Zrestartuj Backend

```powershell
# Zatrzymaj (Ctrl+C w terminalu backendu)

# Uruchom ponownie
.\start_backend.ps1
```

### 2️⃣ Utwórz Nową Symulację

W interfejsie:
1. Kliknij **Stop** (czerwony kwadrat) ⏹️
2. Kliknij **Reset** (ikona odświeżania) 🔄
3. Kliknij **Start** (zielony play) ▶️

### 3️⃣ Obserwuj Rezultaty

**Po ~500-1000 kroków** powinieneś zobaczyć:

```
🔬 PubChem Matcher
Emergence Rate: (0.15)  ← > 0 !
💚 Active

Recent Discoveries
[Match All (5)] ← Button pojawił się!

┌─────────────────────┐
│ SUB_abc12...        │
│ Size: 5             │  ← Większe klastry!
│ Complexity: 0.42    │
│ [Download 🔽]       │  ← Można zapisać!
└─────────────────────┘
```

### 4️⃣ Testuj Matcher

1. Kliknij **Download (🔽)** przy substancji
2. Sprawdź plik w `matches/cluster_*.mol`
3. Zweryfikuj: środkowy atom **NIE jest H**

```mol
  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  ...  ✅ Stopień 1
    1.2990    0.7500    0.0000 O   0  0  ...  ✅ Stopień 2 (nie H!)
    2.5981   -0.0000    0.0000 C   0  0  ...  ✅ Stopień 3
    3.8971    0.7500    0.0000 O   0  0  ...  ✅ Stopień 2
    5.1962   -0.0000    0.0000 H   0  0  ...  ✅ Stopień 1
```

---

## 📊 Monitoring Postępów

### Metryki do Obserwowania:

| Metryka | Cel | Co Oznacza |
|---------|-----|------------|
| **Novelty Rate** | > 0.05 | Nowe struktury się pojawiają |
| **Total Discovered** | Rośnie | Katalog się wypełnia |
| **Novel Substances** | > 0 | Buttony matchera widoczne |
| **Max Cluster Size** | 4-10 | Realistyczne molekuły |

### Troubleshooting:

**Problem**: Novelty rate dalej 0 po 1000 kroków
- ✓ Sprawdź czy backend się zrestartował
- ✓ Sprawdź czy utworzono NOWĄ symulację (nie kontynuowano starej)
- ✓ Poczekaj dłużej (2000+ kroków)

**Problem**: Cluster size max 3
- ✓ Sprawdź `backend/sim/config.py` czy zmiany są obecne
- ✓ Sprawdź logi backendu: `backend/logs.txt`

**Problem**: Buttony matchera nie pojawiają się
- ✓ Sprawdź konsolę przeglądarki (F12) - błędy API?
- ✓ Sprawdź endpoint: `http://localhost:8000/simulation/{id}/novel-substances`

---

## 🔬 Podstawy Naukowe

### Źródła Wartości Parametrów:

1. **binding_threshold = 0.45**
   - Oparte na: Energie wiązań vdW (2-10 kJ/mol), H-bond (10-40 kJ/mol)
   - Źródło: Steiner (2002), Stone (2013), Luo (2007)

2. **theta_break = 1.5**
   - Oparte na: Energie aktywacji dysocjacji peptydów (80-100 kJ/mol)
   - Źródło: Radzicka & Wolfenden (1996)

3. **pulse_amplitude = 1.8**
   - Oparte na: Energia wyładowań w eksperymencie Miller-Urey (1-10 eV)
   - Źródło: Miller (1953), Stribling & Miller (1987)

**Pełna analiza z 20+ referencjami**: `SCIENTIFIC_PARAMETERS_ANALYSIS.md`

---

## 📖 Dodatkowe Zasoby

- **`MATCHER_FIX_SUMMARY.md`**: Pełne podsumowanie wszystkich poprawek
- **`SCIENTIFIC_PARAMETERS_ANALYSIS.md`**: Szczegółowa analiza naukowa
- **`backend/sim/config.py`**: Plik z parametrami (już zaktualizowany!)

---

## ❓ FAQ

**Q: Czy muszę zmieniać coś w kodzie?**  
A: **NIE!** Wszystkie zmiany są już w `backend/sim/config.py`. Wystarczy restart backendu.

**Q: Czy stara symulacja będzie działać z nowymi parametrami?**  
A: **NIE** - musisz utworzyć NOWĄ symulację. Stara używa starych wartości.

**Q: Czy mogę wrócić do starych parametrów?**  
A: Tak, zmień wartości w `backend/sim/config.py` i zrestartuj backend.

**Q: Co jeśli klastry są za duże/za małe?**  
A: Zobacz `SCIENTIFIC_PARAMETERS_ANALYSIS.md` sekcja "Opcja 2: Agresywna" dla ekstremalnych wartości.

---

**Powodzenia! 🧪🔬**

Jeśli coś nie działa, sprawdź logi:
- Backend: `backend/logs.txt`
- Frontend: Konsola przeglądarki (F12)

