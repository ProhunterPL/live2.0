# Porównanie skryptów monitorujących symulacje Phase 2B

## 🎯 Cel: Sprawdzenie aktywnych symulacji i wykrywanie stuck

---

## 📊 Ranking skryptów (od najlepszego)

### 🥇 **1. `check_real_progress.py`** ⭐⭐⭐⭐⭐
**Najlepszy do: wykrywania stuck i szacowania rzeczywistego postępu**

**Zalety:**
- ✅ **Porównuje postęp między uruchomieniami** (cache w `.progress_cache.json`)
- ✅ **Szacuje rzeczywisty postęp** na podstawie CPU usage (uwzględnia log buffering)
- ✅ **Wykrywa stuck** - porównuje CPU time vs czas od ostatniego logu
- ✅ **Szacuje czas zakończenia** na podstawie tempa postępu
- ✅ **Wykrywa completed** - sprawdza `results.json`
- ✅ **Szczegółowa analiza** - pokazuje czy proces faktycznie pracuje mimo starych logów

**Wady:**
- ⚠️ Wymaga co najmniej 2 uruchomień (pierwsze tworzy cache)
- ⚠️ Szacunki mogą być niedokładne przy bardzo starych logach

**Użycie:**
```bash
python3 aws_test/scripts/check_real_progress.py [results_dir]
```

**Przykładowy output:**
```
🔍 miller_urey_extended/run_5 (PID: 12345)
  📊 Last logged step: 438,000/500,000 (87.6%)
  ⏰ Log age: 240.0 minutes (4.00 hours)
  💻 CPU usage: 450%
  📈 Estimated current step: ~475,000 (95.0%)
  ⏱️  Estimated time remaining: ~2.1 hours
  ✅ Process is actively computing (using multiple cores)
```

**Kiedy używać:**
- Gdy chcesz wiedzieć **rzeczywisty postęp** mimo log buffering
- Gdy chcesz **wykryć stuck** (brak postępu między uruchomieniami)
- Gdy chcesz **oszacować czas zakończenia**

---

### 🥈 **2. `check_actual_progress.py`** ⭐⭐⭐⭐
**Najlepszy do: kompleksowej weryfikacji stanu wszystkich symulacji**

**Zalety:**
- ✅ **Sprawdza stan procesów** (R/D/S/Z) - wykrywa zombie/stuck
- ✅ **Sprawdza aktywność plików** (mtime, rozmiary) - wykrywa czy pliki się zmieniają
- ✅ **Sprawdza logi** - pokazuje ostatnie wpisy
- ✅ **Sprawdza rozmiary katalogów** - wykrywa wzrost danych
- ✅ **Kompleksowy przegląd** - wszystko w jednym miejscu

**Wady:**
- ⚠️ Nie porównuje postępu między uruchomieniami
- ⚠️ Nie szacuje rzeczywistego postępu przy log buffering

**Użycie:**
```bash
python3 aws_test/scripts/check_actual_progress.py --results-dir ~/live2.0/results/phase2b_additional
```

**Przykładowy output:**
```
🔍 PROCESS STATE CHECK
📊 PID 12345 (run_5):
   State: R | CPU: 450% | Memory: 12% | Running: 12:34:56
   ✅ Process is actively running
   🧵 Threads: 8

📁 FILE ACTIVITY CHECK
🔍 miller_urey_extended/run_5:
   📄 simulation.log: 1,234,567 bytes, modified 240.0 min ago
   📸 Snapshots: 10 files, latest 240.0 min ago
   💾 Checkpoints: 5 files, latest 120.0 min ago
```

**Kiedy używać:**
- Gdy chcesz **kompleksowy przegląd** wszystkich symulacji
- Gdy chcesz **wykryć zombie procesy** lub problemy z I/O
- Gdy chcesz **sprawdzić aktywność plików** (czy snapshots się tworzą)

---

### 🥉 **3. `check_if_simulation_stuck.sh`** ⭐⭐⭐
**Najlepszy do: szybkiej weryfikacji pojedynczego runu**

**Zalety:**
- ✅ **Szybki** - sprawdza jeden run
- ✅ **Sprawdza completed** - wykrywa `results.json`
- ✅ **Sprawdza proces** - wykrywa czy proces działa
- ✅ **Sprawdza CPU usage** - wykrywa stuck (niskie CPU)
- ✅ **Szacuje czas zakończenia** - na podstawie ostatniego stepu

**Wady:**
- ⚠️ Tylko jeden run na raz
- ⚠️ Nie uwzględnia log buffering (może fałszywie wykryć stuck)

**Użycie:**
```bash
bash aws_test/scripts/check_if_simulation_stuck.sh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_5
```

**Kiedy używać:**
- Gdy chcesz **szybko sprawdzić jeden konkretny run**
- Gdy chcesz **wykryć completed** (sprawdza `results.json`)
- Gdy chcesz **sprawdzić czy proces działa**

---

### **4. `check_simulation_status.sh`** ⭐⭐⭐
**Najlepszy do: szybkiego przeglądu wszystkich runów**

**Zalety:**
- ✅ **Przegląd wszystkich scenariuszy** - miller_urey, hydrothermal, formamide
- ✅ **Pokazuje completed vs running** - szybki status
- ✅ **Pokazuje ostatni step** - podstawowy postęp

**Wady:**
- ⚠️ Nie wykrywa stuck (nie porównuje postępu)
- ⚠️ Nie uwzględnia log buffering

**Użycie:**
```bash
bash aws_test/scripts/check_simulation_status.sh [results_dir]
```

**Kiedy używać:**
- Gdy chcesz **szybki przegląd** wszystkich runów
- Gdy chcesz **sprawdzić completed** (ile runów skończonych)
- Gdy chcesz **podstawowy status** bez szczegółów

---

### **5. `monitor_by_filesize.py`** ⭐⭐
**Najlepszy do: monitorowania zmian plików między uruchomieniami**

**Zalety:**
- ✅ **Wykrywa zmiany rozmiarów** - logi, checkpoints, snapshots
- ✅ **Cache'uje stan** - porównuje między uruchomieniami
- ✅ **Wykrywa nowe pliki** - checkpoints/snapshots

**Wady:**
- ⚠️ Wymaga co najmniej 2 uruchomień
- ⚠️ Nie sprawdza procesów (tylko pliki)
- ⚠️ Może nie wykryć stuck jeśli pliki się nie zmieniają (log buffering)

**Użycie:**
```bash
python3 aws_test/scripts/monitor_by_filesize.py --results-dir ~/live2.0/results/phase2b_additional
```

**Kiedy używać:**
- Gdy chcesz **monitorować zmiany plików** między uruchomieniami
- Gdy chcesz **wykryć nowe checkpoints/snapshots**
- Gdy chcesz **sprawdzić czy pliki się zmieniają** (przy log buffering)

---

### **6. `check_process_details.sh`** ⭐⭐
**Najlepszy do: szczegółowej analizy procesów**

**Zalety:**
- ✅ **Pełna linia komend** - pokazuje wszystkie parametry
- ✅ **Szczegóły procesu** - PID, state, CPU, memory
- ✅ **Wykrywa stuck** - sprawdza wiek logu (>24h)

**Wady:**
- ⚠️ Tylko procesy (nie sprawdza plików)
- ⚠️ Prosty wykrywacz stuck (tylko wiek logu)

**Użycie:**
```bash
bash aws_test/scripts/check_process_details.sh
```

**Kiedy używać:**
- Gdy chcesz **zobaczyć pełną linię komend** procesów
- Gdy chcesz **sprawdzić szczegóły procesów** (PID, state, CPU)
- Gdy chcesz **szybko wykryć bardzo stare logi** (>24h)

---

## 🎯 Rekomendacja: Który skrypt użyć?

### **Dla codziennego monitorowania (wykrywanie stuck):**
```bash
# Najlepszy: check_real_progress.py
python3 aws_test/scripts/check_real_progress.py ~/live2.0/results/phase2b_additional
```
**Dlaczego:** Porównuje postęp między uruchomieniami, szacuje rzeczywisty postęp mimo log buffering, wykrywa stuck.

### **Dla kompleksowej weryfikacji (wszystkie symulacje):**
```bash
# Najlepszy: check_actual_progress.py
python3 aws_test/scripts/check_actual_progress.py --results-dir ~/live2.0/results/phase2b_additional
```
**Dlaczego:** Sprawdza procesy, pliki, logi, rozmiary - wszystko w jednym miejscu.

### **Dla szybkiego sprawdzenia (jeden run):**
```bash
# Najlepszy: check_if_simulation_stuck.sh
bash aws_test/scripts/check_if_simulation_stuck.sh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_5
```
**Dlaczego:** Szybki, sprawdza completed, proces, CPU usage.

### **Dla szybkiego przeglądu (wszystkie runy):**
```bash
# Najlepszy: check_simulation_status.sh
bash aws_test/scripts/check_simulation_status.sh ~/live2.0/results/phase2b_additional
```
**Dlaczego:** Szybki przegląd wszystkich scenariuszy, pokazuje completed vs running.

---

## 🔍 Jak wykryć stuck?

### **Metoda 1: Porównanie postępu (najlepsza)**
```bash
# Uruchom 2 razy z odstępem 1-2h
python3 aws_test/scripts/check_real_progress.py
# Jeśli postęp = 0 między uruchomieniami → STUCK
```

### **Metoda 2: Stan procesu + CPU**
```bash
python3 aws_test/scripts/check_actual_progress.py
# Jeśli:
# - State = D (I/O wait) + CPU < 10% → STUCK
# - State = Z (zombie) → CRASHED
# - State = R (running) + CPU > 100% → OK (może być log buffering)
```

### **Metoda 3: Wiek logu + CPU**
```bash
bash aws_test/scripts/check_if_simulation_stuck.sh <run_dir>
# Jeśli:
# - Log > 24h + CPU < 10% → STUCK
# - Log > 24h + CPU > 100% → Log buffering (OK, ale sprawdź postęp)
```

---

## 📋 Checklist: Czy symulacja jest stuck?

- [ ] **Proces istnieje?** (`ps aux | grep run_phase2_full`)
- [ ] **Proces state = R?** (running, nie D/S/Z)
- [ ] **CPU usage > 100%?** (aktywna praca)
- [ ] **Log się zmienia?** (mtime < 1h lub rozmiar rośnie)
- [ ] **Postęp między uruchomieniami?** (check_real_progress.py)
- [ ] **Snapshots się tworzą?** (nowe pliki w snapshots/)
- [ ] **Checkpoints się tworzą?** (nowe pliki w checkpoints/)

**Jeśli wszystkie = NIE → STUCK**

---

## 🚀 Szybkie komendy

```bash
# 1. Sprawdź wszystkie aktywne symulacje (najlepsze)
python3 aws_test/scripts/check_real_progress.py ~/live2.0/results/phase2b_additional

# 2. Kompleksowa weryfikacja (wszystko)
python3 aws_test/scripts/check_actual_progress.py --results-dir ~/live2.0/results/phase2b_additional

# 3. Szybki przegląd (wszystkie runy)
bash aws_test/scripts/check_simulation_status.sh ~/live2.0/results/phase2b_additional

# 4. Sprawdź jeden konkretny run
bash aws_test/scripts/check_if_simulation_stuck.sh ~/live2.0/results/phase2b_additional/miller_urey_extended/run_5

# 5. Sprawdź szczegóły procesów
bash aws_test/scripts/check_process_details.sh
```

---

## 📝 Uwagi

1. **Log buffering:** Stare logi (np. 4h) + wysokie CPU (>100%) = **normalne** (log buffering). Użyj `check_real_progress.py` do szacowania rzeczywistego postępu.

2. **Cache files:** `check_real_progress.py` i `monitor_by_filesize.py` tworzą cache files (`.progress_cache.json`, `~/.phase2b_filesize_cache.json`). Można je bezpiecznie usunąć (zostaną odtworzone).

3. **Completed detection:** Wszystkie skrypty sprawdzają `results.json` - jeśli istnieje, symulacja jest completed.

4. **Stuck detection:** Najlepsze wyniki daje **porównanie postępu** między uruchomieniami (`check_real_progress.py`).

---

**Ostatnia aktualizacja:** 2025-11-17

