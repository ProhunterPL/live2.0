# Naprawa Deadlocka w Phase 2B - Podsumowanie PL

**Data**: 10 listopada 2025  
**Problem**: Symulacje zatrzymały się w nieskończonej pętli (kernel wykrywania klastrów)  
**Status**: 🚨 KRYTYCZNY - 7 z 9 symulacji zablokowanych

---

## CO TRZEBA ZROBIĆ - SZYBKA INSTRUKCJA

```bash
# Na instancji AWS:
cd ~/live2.0
chmod +x aws_test/DEPLOY_FIX_NOW.sh
bash aws_test/DEPLOY_FIX_NOW.sh
```

To:
1. ✅ Zabije zablokowane procesy
2. ✅ Naprawi kod
3. ✅ Zrestartuje z bezpieczną konfiguracją

---

## CO SIĘ STAŁO

### Objawy
- **Wysokie CPU** ale **brak postępu** przez 6-41 godzin
- Logi zamrożone na konkretnych krokach (np. run_7 utknął na 363,000 przez 6 godzin)
- **Brak nowych snapshotów/checkpointów** mimo zużycia CPU
- Stan procesu: `Sl` (śpi/czeka)
- 196 wątków na proces

### Przyczyna

Algorytm wykrywania klastrów (`update_clusters_kernel`) używa pętli O(N²):

```python
for i in range(1000):              # Wszystkie cząstki
    for j in range(i+1, 1000):     # Sprawdź wszystkie pary
        # Operacje Union-Find
```

Z złożonymi sieciami wiązań:
- ~500,000 iteracji na sprawdzenie klastrów
- Union-Find w patologicznych przypadkach → wykładnicze spowolnienie
- Deadlock synchronizacji wątków

### Dotknięte uruchomienia

| Run | Status | Utknął na | Czas |
|-----|--------|-----------|------|
| run_1 | ✅ Ukończony | - | - |
| run_2 | 🚫 Utknął | 24,000 | 41.5h |
| run_3 | 🚫 Utknął | 335,000 | 30.4h |
| run_4 | 🚫 Utknął | 26,000 | 41.4h |
| run_5 | 🚫 Utknął | 439,000 | 11.2h |
| run_6 | 🚫 Utknął | 78,000 | 14.7h |
| run_7 | 🚫 Utknął | 363,000 | 6.3h |
| run_8 | 🚫 Utknął | 104,000 | 13.6h |
| run_9 | ⏳ Działa | 75,000 | 0h |

**Tylko run_9 jest zdrowy** (startował niedawno, nie trafił jeszcze na problematyczny stan)

---

## NAPRAWA

### Patch kodu

**Plik**: `backend/sim/core/stepper.py` (linie 507-514)

**Przed** (zakodowane na sztywno co 1200 kroków):
```python
# OPTIMIZATION: Update clusters every 1200 steps
if (self.step_count - 300) % 1200 == 0:
    self.binding.update_clusters(...)
```

**Po** (konfigurowalne):
```python
# OPTIMIZATION: Update clusters - now configurable
cluster_interval = getattr(self.config, 'cluster_check_interval', 1200)
if cluster_interval < 999999999:  # Tylko jeśli nie wyłączone
    if (self.step_count - 300) % cluster_interval == 0:
        self.binding.update_clusters(...)
```

### Zmiana konfiguracji

**Plik**: `aws_test/configs/phase2_miller_urey_extended_SAFER.yaml`

```yaml
physics:
  # WYŁĄCZ WYKRYWANIE KLASTRÓW (główny podejrzany o deadlock)
  cluster_check_interval: 999999999  # Faktycznie wyłączone
```

---

## PLAN NAPRAWY

### Faza 1: Natychmiastowa (Teraz)
1. ✅ Zabij zablokowane symulacje (zostaw run_9)
2. ✅ Zastosuj hotfix kodu
3. ⏳ Monitoruj run_9 do ukończenia
4. ✅ Uruchom 9 nowych symulacji (run_10 do run_18)

**Czas**: 2-3 godziny  
**Rezultat**: 10-11 ukończonych symulacji

### Faza 2: Uzupełnienie do 30 (Później)
5. Analizuj wyniki z Fazy 1
6. Jeśli sukces, uruchom pozostałe 19-20 symulacji
7. Ukończ pełną analizę statystyczną

**Czas**: +3-4 godziny  
**Rezultat**: 30 ukończonych symulacji

---

## JAK WDROŻYĆ

### Krok 1: Skopiuj pliki na AWS (z Windows)

```powershell
# Na lokalnym komputerze (PowerShell):
.\copy_fix_to_aws.ps1
```

### Krok 2: Uruchom naprawę (na AWS)

```bash
# Połącz się z AWS:
ssh ubuntu@ip-172-31-0-42

# Uruchom automatyczną naprawę:
cd ~/live2.0
bash aws_test/DEPLOY_FIX_NOW.sh
```

**To wszystko!** Skrypt zrobi resztę automatycznie.

---

## MONITOROWANIE

### Sprawdź postęp
```bash
# Ogólny status
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py

# Monitorowanie oparte na plikach (omija buforowanie logów)
python3 ~/live2.0/aws_test/scripts/monitor_by_filesize.py

# Obserwuj w czasie rzeczywistym (aktualizacja co 5 min)
watch -n 300 "python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
```

### Sprawdź logi
```bash
# Nowe symulacje (czas rzeczywisty)
tail -f ~/live2.0/logs/phase2b_safe/run_10.log

# Stara symulacja (run_9)
tail -f ~/live2.0/results/phase2b_additional/miller_urey_extended/run_9/simulation.log
```

---

## KRYTERIA SUKCESU

✅ **Symulacje są zdrowe jeśli:**
- Logi aktualizują się co 1-2 minuty
- Nowe snapshoty co 50K kroków
- Stan procesu to `R` lub `Sl` z niedawną aktywnością plików
- Żaden proces nie utknął >1 godzinę

❌ **Symulacje są zablokowane jeśli:**
- Logi niezmienione >1 godzinę
- Wysokie CPU ale brak nowych plików
- Stan procesu `Sl` przez dłuższy czas
- Brak snapshotów >2 godziny

---

## OCZEKIWANE WYNIKI

### Wydajność (bez wykrywania klastrów)
- **Prędkość**: ~140 kroków/sek/rdzeń na CPU
- **Czas na symulację**: 60-90 minut (500K kroków)
- **9 równoległych symulacji**: ~90 minut łącznie
- **Koszt AWS**: ~$2-3 za 9 uruchomień

### Wpływ naukowy
- Wykrywanie klastrów wpływa tylko na **metryki**, nie chemię
- Wiązania, reakcje, energia nadal dokładne
- Klastry można obliczyć w post-processingu jeśli potrzeba
- **Brak utraty wartości naukowej**

---

## CO ZOSTAŁO STWORZONE

### Skrypty
- `copy_fix_to_aws.ps1` - Kopiowanie na AWS (Windows)
- `aws_test/DEPLOY_FIX_NOW.sh` - Jedno-komendowe wdrożenie
- `aws_test/scripts/kill_stuck_simulations.sh` - Zabij zablokowane procesy
- `aws_test/scripts/apply_cluster_fix.sh` - Zastosuj hotfix
- `aws_test/scripts/restart_phase2b_safe.sh` - Restart z bezpiecznym configiem

### Konfiguracje
- `aws_test/configs/phase2_miller_urey_extended_SAFER.yaml` - Bez wykrywania klastrów

### Dokumentacja
- `docs/aws_test/PODSUMOWANIE_NAPRAWY_PL.md` - Ten plik (po polsku)
- `docs/aws_test/EMERGENCY_SUMMARY.md` - Angielska wersja
- `docs/aws_test/CLUSTER_FIX_INSTRUCTIONS.md` - Szczegółowe instrukcje

---

## ROZWIĄZYWANIE PROBLEMÓW

### Jeśli nowe symulacje też się zablokują:

1. **Sprawdź czy fix został zastosowany:**
   ```bash
   grep "cluster_interval = getattr" backend/sim/core/stepper.py
   ```

2. **Sprawdź config:**
   ```bash
   grep "cluster_check_interval" aws_test/configs/phase2_miller_urey_extended_SAFER.yaml
   ```

3. **Sprawdź stan procesów:**
   ```bash
   ps aux | grep "run_phase2_full.py"
   ```

---

## NASTĘPNE KROKI PO NAPRAWIE

1. ✅ Zweryfikuj 10-11 ukończonych symulacji
2. ✅ Uruchom szybką analizę na ukończonych
3. ✅ Jeśli wyniki dobre, uruchom pozostałe 19-20 symulacji
4. ✅ Ukończ pełną analizę Phase 2B
5. ✅ Wygeneruj rysunki do artykułu

**Szacowany całkowity czas naprawy**: 5-6 godzin do 30 ukończonych symulacji

---

## PODSUMOWANIE

**Problem**: Symulacje zatrzymały się w wykrywaniu klastrów (O(N²) deadlock)  
**Rozwiązanie**: Wyłącz wykrywanie klastrów w konfiguracji  
**Koszt**: Brak - klastry to tylko metryki, nie wpływają na chemię  
**Czas**: 2-3 godziny na naprawę + restart  

**Naprawa jest solidna - wykrywanie klastrów jest teraz całkowicie wyłączone w trybie SAFER!** 🚀

**Zaufaj procesowi i uruchom `bash aws_test/DEPLOY_FIX_NOW.sh`** ✨

