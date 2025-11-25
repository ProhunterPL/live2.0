# FAQ - Problemy z Klastrami i PubChem Matcher

## ❓ Dlaczego PubChem Matcher pokazuje "Waiting for clusters..."?

**Odpowiedź**: Detekcja novel substances działa co **700 kroków** (domyślnie). Jeśli jesteś na kroku 32000, ostatnia detekcja była przy kroku 31500, a następna będzie przy 32200.

**Rozwiązanie**:
- Poczekaj do następnego interwalu (sprawdź: `python force_cluster_detection.py`)
- Lub użyj konfiguracji z `novelty_check_interval: 200`

---

## ❓ Dlaczego metryki pokazują 498 klastrów, ale PubChem Matcher pokazuje 0?

**Odpowiedź**: Metryki używają **uproszczonego przybliżenia**:
```python
cluster_count = particles_with_bonds / 2  # Zakłada średnio 2 cząstki = 1 klaster
```

To nie jest prawdziwa liczba! PubChem Matcher używa **rzeczywistej detekcji** klastrów z katalogu.

**Rozwiązanie**:
- Metryki to tylko estymacja (dla wydajności)
- Prawdziwa liczba jest w PubChem Matcher (po detekcji)
- Możesz naprawić metryki (patrz: `ROZWIAZANIE_KROK_PO_KROKU.md`)

---

## ❓ Dlaczego klastry na wizualizacji wyglądają "dziwnie" (długie wiązania)?

**Odpowiedź**: Maksymalna długość wiązania to **6.8 * radius = 3.4 jednostki**. To pozwala na bardzo rozciągnięte wiązania.

**Rozwiązanie**:
- Zmniejsz max bond length do `4.0 * radius = 2.0 jednostki`
- Edytuj `backend/sim/core/binding.py` linia 313
- Lub użyj konfiguracji `configs/fast_cluster_detection.yaml`

---

## ❓ Jak często są odświeżane klastry na fronten

dzie?

**Odpowiedź**: 
- Klastry są pobierane co **500 kroków** (cache dla wydajności)
- Novel substances są wykrywane co **700 kroków** (domyślnie)
- Frontend otrzymuje dane co **~0.067s** (15 FPS)

**Rozwiązanie**:
- Zmniejsz cache interval do 200 kroków (patrz: `ROZWIAZANIE_KROK_PO_KROKU.md`)
- Zmniejsz novelty_check_interval do 200 kroków

---

## ❓ Dlaczego nie widzę małych klastrów (2 cząstki)?

**Odpowiedź**: `min_cluster_size = 3` (domyślnie). Klastry <3 cząstek **nie są rejestrowane** jako novel substances.

**Rozwiązanie**:
```yaml
# W config YAML:
min_cluster_size: 2  # Zamiast 3
```

**Uwaga**: Małe klastry są mniej stabilne i ciekawe chemicznie.

---

## ❓ Klastry szybko się rozpadają - jak je ustabilizować?

**Odpowiedź**: Wiązania zrywają się gdy:
- Dystans > `5.0 * radius = 2.5 jednostki`
- Strain > 300% (rozciągnięcie)
- Age > max_age (starzenie)

**Rozwiązanie**:
```yaml
# W config YAML:
unbinding_threshold: 0.12  # Zamiast 0.18 (trudniej się zrywają)

# Lub w kodzie (binding.py), zwiększ max_age:
max_age = 20000.0  # Zamiast 10000.0 (wiązania żyją dłużej)
```

---

## ❓ Jak wymusi natychmiastową detekcję klastrów?

**Odpowiedź**: Nie ma bezpośredniego API do wymuszenia. Ale możesz:

**Opcja 1**: Poczekaj do następnego interwalu
```powershell
python force_cluster_detection.py  # Pokaże kiedy będzie następna detekcja
```

**Opcja 2**: Zrestartuj z mniejszym interwałem
```yaml
novelty_check_interval: 100  # Detekcja co 100 kroków
```

**Opcja 3**: Użyj API do ręcznej detekcji (TODO - do implementacji)

---

## ❓ Jak sprawdzić czy detekcja działa?

**Odpowiedź**: Sprawdź logi:
```powershell
cat logs\logs.txt | Select-String "detect_novel_substances"
```

Powinno być:
```
INFO - Detecting novel substances at step 700
INFO - Detecting novel substances at step 1400
INFO - Detecting novel substances at step 2100
...
```

Jeśli nie ma logów:
- `detect_novel_substances: true` w config
- Backend działa
- Symulacja jest uruchomiona (started)

---

## ❓ PubChem Matcher pokazuje klastry, ale są "mało prawdopodobne"

**Odpowiedź**: To może być spowodowane:
1. **Długimi wiązaniami** (>2.0 jednostki) - nierealistyczne
2. **Niską gęstością** klastrów (density < 0.3)
3. **Niestabilną strukturą** (klaster zaraz się rozpadnie)

**Rozwiązanie**:
- Zmniejsz max bond length do 4.0 * radius
- Zwiększ stabilność wiązań (unbinding_threshold)
- Filtruj klastry po density (>0.3)

---

## ❓ Jak działa PubChem Matcher?

**Odpowiedź**: 
1. Symulacja wykrywa klastry (connected components)
2. Co `novelty_check_interval` kroków, klastry są dodawane do katalogu
3. Frontend pobiera novel substances z API (`/novel-substances`)
4. PubChem Matcher pokazuje je w liście
5. Kliknięcie "Match All" wywołuje dopasowanie do PubChem

**Krok szczegółowy**:
```
Step N (N % 700 == 0):
  → detect_novel_substances()
  → binding.get_clusters() [prawdziwa detekcja]
  → catalog.add_substance() [jeśli ≥3 cząstki]
  → Frontend: API /novel-substances
  → PubChem Matcher: pokazuje listę
```

---

## ❓ Czy mogę wyłączyć cache dla testowania?

**Odpowiedź**: Tak, ale będzie bardzo wolne!

```python
# W stepper.py:
if self.step_count % 1 == 0:  # Co krok (zamiast 500)
    clusters = self.binding.get_clusters()
```

**Uwaga**: FPS spadnie z ~5 do ~0.5 (10x wolniej)

**Lepsze rozwiązanie**: Zmniejsz do 100-200 kroków (kompromis)

---

## ❓ Jak zoptymalizować dla maksymalnej detekcji klastrów?

**Odpowiedź**: Użyj tej konfiguracji:
```yaml
# Szybka detekcja
novelty_check_interval: 100
min_cluster_size: 2

# Stabilne wiązania
binding_threshold: 0.5
unbinding_threshold: 0.12

# Częste odświeżanie (wolniejsze!)
# W kodzie: cache co 100 kroków zamiast 500
```

**Koszt**: ~20-30% wolniejsza symulacja

---

## ❓ Czy długość wiązań 6.8 jest błędna?

**Odpowiedź**: Nie, to **zamierzone**! 6.8 * 0.5 = 3.4 Å to zasięg van der Waals (z literatury).

**Ale**: Większość klastrów powinna mieć wiązania kowalencyjne (~1.5 Å), nie vdW.

**Rozwiązanie**:
- Zmniejsz do 4.0 * radius dla bardziej stabilnych struktur
- Lub zostaw 6.8 i zaakceptuj luźne klastry

---

## ❓ Co znaczy "Density: 0.222" w klastrze?

**Odpowiedź**: 
```
Density = liczba_wiązań / max_możliwych_wiązań
        = bonds / (n_particles * (n_particles - 1) / 2)
```

Dla klastra 9 cząstek, 8 wiązań:
```
Density = 8 / (9 * 8 / 2) = 8 / 36 = 0.222
```

**Interpretacja**:
- Density = 1.0: pełny graf (każdy z każdym)
- Density = 0.5: połowa możliwych wiązań
- Density = 0.222: rzadki graf (22% wiązań)

**Dla chemii**:
- Density < 0.3: luźna struktura, może niestabilna
- Density 0.3-0.6: typowe molekuły
- Density > 0.6: gęste klastry (rzadkie w chemii)

---

## ❓ Skrypty diagnostyczne nie działają (błąd połączenia)

**Odpowiedź**: Backend nie jest uruchomiony.

**Rozwiązanie**:
```powershell
# Uruchom backend
.\\start_backend.ps1

# Sprawdź czy działa
curl http://localhost:8000/health

# Następnie uruchom diagnostykę
python force_cluster_detection.py
```

---

## ❓ Gdzie są logi?

**Odpowiedź**: `logs/logs.txt` (w katalogu projektu)

```powershell
# Zobacz wszystkie logi
cat logs\logs.txt

# Filtruj po słowie kluczowym
cat logs\logs.txt | Select-String "novel"
cat logs\logs.txt | Select-String "cluster"
cat logs\logs.txt | Select-String "ERROR"
```

---

## ❓ Szybkie rozwiązanie - jeden komenda?

**Odpowiedź**:
```powershell
# Diagnostyka (powie co robić)
python force_cluster_detection.py

# Jeśli mówi "poczekaj X kroków":
# → Poczekaj i odśwież frontend (Ctrl+R)

# Jeśli mówi "brak novel substances, zmień config":
# → Użyj configs/fast_cluster_detection.yaml
cd backend
python -m api.server --config ../configs/fast_cluster_detection.yaml
```

---

## 🎓 Więcej informacji

- **Szczegółowa analiza**: `docs/CLUSTER_DETECTION_ISSUE.md`
- **Rozwiązanie krok po kroku**: `ROZWIAZANIE_KROK_PO_KROKU.md`
- **Podsumowanie (PL)**: `PODSUMOWANIE_PROBLEM_KLASTROW.md`
- **Skrypty**:
  - `check_real_clusters.py` - stan klastrów
  - `force_cluster_detection.py` - timing detekcji
- **Konfiguracja**: `configs/fast_cluster_detection.yaml`

