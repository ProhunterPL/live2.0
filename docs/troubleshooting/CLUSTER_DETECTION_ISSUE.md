# Problem z Detekcją Klastrów - Analiza i Rozwiązania

## 🔍 Diagnoza Problemu

Twoja symulacja na kroku 32k pokazuje rozbieżności w liczeniu klastrów:

### Symptomy:
1. **Metryki pokazują**: 498 klastrów
2. **PubChem Matcher pokazuje**: "Waiting for clusters..." (0 substancji)
3. **Wizualizacja**: Dziwne, rozciągnięte klastry z długimi wiązaniami

### Przyczyny:

#### 1. **Metryki używają uproszczonego przybliżenia**

W `backend/sim/core/stepper.py:1083`:
```python
# Simple approximation: particles_with_bonds / 2 (assuming average 2 particles per cluster)
cluster_count = max(1, int(particles_with_bonds / 2))
```

**To nie jest prawdziwa liczba klastrów!** To tylko estymacja zakładająca że średnio 2 cząstki = 1 klaster.

#### 2. **Detekcja novel substances działa co 700 kroków**

W `backend/sim/config.py:49`:
```python
novelty_check_interval: int = Field(default=700, gt=0)
```

Przy kroku 32000:
- Ostatnia detekcja: krok `32000 - (32000 % 700) = 31500`
- Następna detekcja: krok `32200`

#### 3. **Klastry są cache'owane co 500 kroków**

W `backend/sim/core/stepper.py:1339-1346`:
```python
if self.step_count % 500 == 0:
    bonds = self.binding.get_bonds()
    clusters = self.binding.get_clusters()
    # Cache the data for intermediate steps
```

Wizualizacja może pokazywać stare dane między tymi interwałami.

#### 4. **Novel substances wymagają ≥3 cząstek**

W `backend/sim/config.py:48`:
```python
min_cluster_size: int = Field(default=3, ge=1)
```

Klastry 2-cząsteczkowe (które widzisz na zdjęciach) **nie są** rejestrowane jako novel substances.

#### 5. **Długie wiązania są dozwolone**

W `backend/sim/core/binding.py`:
```python
PARTICLE_RADIUS_COMPILE = 0.5

# Wiązania formują się przy:
if r <= PARTICLE_RADIUS_COMPILE * 6.8:  # = 3.4
    # ...

# Wiązania zrywają się przy:
if r > PARTICLE_RADIUS_COMPILE * 5.0:  # = 2.5
    result = 1  # break bond
```

To pozwala na **bardzo rozciągnięte** wiązania (do 3.4 jednostki), co jest nierealistyczne chemicznie.

## 🔧 Rozwiązania

### Rozwiązanie 1: Sprawdź czy detekcja działa

Uruchom skrypt diagnostyczny:
```powershell
python check_real_clusters.py
```

To pokaże:
- Rzeczywisty stan katalogu (ile novel substances)
- Kiedy była ostatnia detekcja
- Kiedy będzie następna
- Statystyki wiązań

### Rozwiązanie 2: Zmniejsz interwały detekcji (tymczasowe)

Edytuj aktualną konfigurację symulacji lub stwórz nową z:

```yaml
novelty_check_interval: 200  # Zamiast 700
min_cluster_size: 3          # Zachowaj (sensowne klastry)
```

**UWAGA**: Częstsza detekcja = niższy FPS

### Rozwiązanie 3: Napraw metryki (zalecane długoterminowo)

Zamień uproszczone liczenie w `backend/sim/core/stepper.py:1083` na:

```python
# BEFORE (wrong):
cluster_count = max(1, int(particles_with_bonds / 2))

# AFTER (correct):
# Use real cluster detection from binding system
real_clusters = self.binding.get_clusters(min_size=2)
cluster_count = len(real_clusters)
```

**Kompromis**: To będzie wolniejsze, ale dokładne.

### Rozwiązanie 4: Zmniejsz dozwoloną długość wiązań

W `backend/sim/core/binding.py:313` zmień:

```python
# BEFORE:
if r <= PARTICLE_RADIUS_COMPILE * 6.8:  # = 3.4 (bardzo długie!)

# AFTER:
if r <= PARTICLE_RADIUS_COMPILE * 4.0:  # = 2.0 (bardziej realistyczne)
```

To sprawi że wiązania będą krótsze i bardziej stabilne chemicznie.

### Rozwiązanie 5: Zmniejsz cache'owanie (jeśli masz wydajność)

W `backend/sim/core/stepper.py:1339` zmień:

```python
# BEFORE:
if self.step_count % 500 == 0:

# AFTER:
if self.step_count % 100 == 0:  # Częstsze odświeżanie
```

## 📊 Co Sprawdzić Teraz

1. **Poczekaj na następny interval novelty detection**:
   - Jeśli jesteś na kroku 32000, następna detekcja będzie przy kroku 32200
   - Sprawdź czy wtedy pojawią się novel substances

2. **Sprawdź logi backendu**:
   ```powershell
   cat logs\logs.txt | Select-String "detect_novel_substances"
   ```

3. **Sprawdź czy klastry są rzeczywiście ≥3 cząstki**:
   - Na twoich zdjęciach widzę klastry 2-cząsteczkowe
   - Te **NIE SĄ** rejestrowane jako novel substances
   - Większe klastry (9 cząstek na drugim zdjęciu) **POWINNY** być wykryte

4. **Sprawdź density klastrów**:
   - Na pierwszym zdjęciu: 3 cząstki, 2 wiązania
   - Density = 2 / (3 * 2 / 2) = 2/3 = 0.667 ✅ (OK)
   - Na drugim: 9 cząstek, 8 wiązań
   - Density = 8 / (9 * 8 / 2) = 8/36 = 0.222 ⚠️ (niska)

## 🎯 Najprawdopodobniejsza Przyczyna

**Klastry są wykrywane, ale:**
1. **Cache** - frontend pokazuje stare dane (sprzed max 500 kroków)
2. **Timing** - ostatnia novelty detection była 500 kroków temu
3. **Rozmiar** - większość klastrów ma <3 cząstki (nie są rejestrowane)
4. **Niestabilność** - długie wiązania powodują że klastry szybko się rozpadają

## 💡 Zalecana Akcja

**Krótkoterminowo:**
```powershell
# 1. Sprawdź diagnostykę (gdy backend działa)
python check_real_clusters.py

# 2. Poczekaj na krok 32200 (następna novelty detection)

# 3. Sprawdź czy wtedy pojawią się novel substances w PubChem Matcher
```

**Długoterminowo:**
1. Zmniejsz `novelty_check_interval` do 200-300
2. Zmniejsz `max_bond_length` z 6.8 do 4.0
3. Zamień metryki na prawdziwe liczenie (zamiast przybliżenia)
4. Zmniejsz cache interval dla klastrów z 500 do 100-200

## 📁 Pliki do Edycji

1. **backend/sim/config.py** - zmień `novelty_check_interval`
2. **backend/sim/core/binding.py** - zmień max bond length
3. **backend/sim/core/stepper.py** - zmień metryki i cache intervals

## ⚠️ Ważne Uwagi

- Częstsza detekcja = niższy FPS (kompromis)
- Prawdziwe liczenie klastrów = wolniejsze ale dokładne
- Krótsze wiązania = bardziej realistyczne chemicznie
- Cache'owanie jest potrzebne dla wydajności

