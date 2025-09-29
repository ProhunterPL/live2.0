# Live 2.0 --- v1 (Open-Ended Prebiotic Simulator, 2D, GPU)

**Dokument wykonawczy (MD) do Cursor AI / agenta budującego system**

> Cel: zbudować **otwartą** symulację 2D zdolną do generowania
> **nowych** cząstek („pierwiastków"), **nowych substancji**
> (klastrów/grafów) i **nowych reakcji** (przebudów grafów) bez
> zamykania świata w sztywnej liście gatunków. Dwa tryby: **Preset
> Prebiotic (A)** dla szybkich dem i walidacji oraz **Open Chemistry
> (B)** dla otwartej emergencji (domyślny do badań).

------------------------------------------------------------------------

## 0) Założenia i zasady projektowe

-   **Otwartość**: brak statycznej listy cząstek i reakcji w trybie B.
    Cząstki mają **ciągłe atrybuty**, wiązania powstają z
    **potencjałów**. Substancje identyfikujemy **po grafie** wiązań, nie
    po nazwie.\
-   **Lokalność**: tylko interakcje krótkiego zasięgu, listy sąsiadów,
    periodyczne brzegi (torus).\
-   **Inwarianty**: kontrolujemy (konfigurowalne) zachowanie energii (w
    przybliżeniu), masy i wybranych składowych wektorów „ładunków" `q⃗`.\
-   **Stabilność numeryczna**: adaptacyjny krok, ograniczenia gęstości,
    wygaszanie energii.\
-   **Skalowalność**: GPU-first (Taichi/CUDA) + binarny streaming do
    React.\
-   **Obserwowalność**: metryki „novelty" (tempo pojawiania się
    **nowych** substancji), złożoność grafów, rozkłady atrybutów.\
-   **Powtarzalność**: seedy RNG, snapshoty (checkpointy) stanu.

### 0.1 Obserwacje po weryfikacji v1 i szybkie usprawnienia

-   **Jedno źródło prawdy dla API i planu** – unikać duplikacji endpointów
    oraz niespójnych kontraktów; każda funkcja powinna mieć „owner'a” w
    planie i w kodzie.
-   **Trzymajmy się pętli iteracyjnej**: `Plan → Implementacja → Pomiar →
    Retrospekcja`. Każda iteracja kończy się krótkim raportem i audytem
    metryk (w repo `reports/` + wpis w changelogu).
-   **Priorytetyzacja przez "safety rails"**: zanim rozszerzymy fizykę,
    domykamy stabilność, monitoring oraz regresje API. Bez działających
    testów i obserwowalności nie przechodzimy do kolejnych faz.

### 0.2 Rytm pracy (operational cadence)

1. **Poniedziałek** – aktualizacja roadmapy (ten dokument) o
   zrealizowane punkty i nowe odkrycia z ostatniego tygodnia.
2. **Wtorek–Czwartek** – implementacja sprintowa + codzienny log
   (skrócony changelog) w `reports/devlog.md`.
3. **Piątek** – `stability run` (12h), analiza KPI, aktualizacja backlogu
   ryzyk.
4. **Stały kanał feedbacku** – wszystkie nowe pomysły / dema trafiają do
   sekcji „Discovery Backlog” poniżej.

### 0.3 Discovery Backlog (utrzymywać max 10 pozycji)

1. Light-weight serializer dla klastrów (odchudzi stream o 30–40%).
2. Narzędzia do inspekcji konfiguracji (diff seeda / parametrów między
   snapshotami).
3. Minimalna wizualizacja 3D (ortogonalna projekcja) jako proof-of-concept
   dla fazy 2.4.
4. Integracja z `wandb` / `mlflow` dla logów eksperymentalnych.
5. Tryb headless (CLI) do batchowych runów na klastrze.
6. Mechanizm walidacji kontraktów API (schemat JSON) – spięty z CI.
7. Automatyczne generowanie changelogu na podstawie commitów/PR (konwencje
   Conventional Commits).
8. Samocertyfikacja wersji (semver + plik `VERSION`).
9. Lekka biblioteka JS do dekodowania strumienia (wydzielić z frontendu).
10. Analiza `novelty` w czasie rzeczywistym z adaptacyjnym progiem alertów.

------------------------------------------------------------------------

## 1) Zakres v1 (MVP)

**Definition of Done v1 (MVP)** – wszystkie punkty poniżej muszą być
odhaczone, zanim przejdziemy do fazy 1:

1. **Siatka 2D** H×W (domyślnie 256×256), periodyczne brzegi oraz
   aktualizowany spatial hash (✅ w repo; utrzymywać testy regresyjne).
2. **Tryby uruchomienia**: - **A) Preset Prebiotic**: ciągłe stężenia,
   podstawowe reakcje, walidacja strumienia. - **B) Open Chemistry**
   *(domyślny)*: cząsteczki off-lattice, potencjały, wiązania,
   katalogowanie nowych substancji.
3. **Źródła energii**: impulsy czasoprzestrzenne + zanikanie;
   parametryzowane z poziomu configu i snapshotów.
4. **Transmutacje**: rzadkie mutacje atrybutów cząstek w obszarach
   wysokiej energii → **narodziny nowych klas „pierwiastków"** (obecnie
   TODO – zaplanowane w sprincie 0.3).
5. **Detektory nowości**: logowanie pojawień się nowych substancji (ID
   grafu), rozmiaru, czasu życia; pipeline `metrics → aggregator → API →
   frontend`.
6. **Frontend** (React + TS): heatmapy (gęstość/energia), licznik
   nowości, podgląd klastrów jako grafów (SVG), sterowanie start/pause,
   wybór trybu. Zadbajmy o fallback „headless” (CLI) – wymóg QA.
7. **API** (FastAPI + WebSocket): spójne endpointy REST + binarny stream
   (msgpack) z docelowym 15 FPS. Kontrakt opisany w `backend/api/protocol.md`.
8. **Snapshoty**: zapis/odczyt stanu (particles, grafy, parametry, seed)
   + test regresyjny snapshotów.
9. **Testy**: jednostkowe, property-based (inwarianty energii/masy),
   wydajnościowe, long-run (novelty\>0). Coverage min. 70% przy
   zakończeniu v1.

**Poza zakresem v1** (opcjonalne w v1.1+): błony/kompartmenty, gradienty
pH, formalne bilanse termodynamiczne, ewolucja potencjałów środowiska w
czasie globalnie (na razie tylko mutacje lokalne).

------------------------------------------------------------------------

## 2) Architektura systemu

### 2.1 Monorepo

    live2/
      backend/
        sim/
          core/
            grid.py            
            rng.py             
            fields_ca.py        
            particles.py        
            potentials.py       
            binding.py          
            graphs.py           
            catalog.py          
            metrics.py          
            energy.py           
            stepper.py          
          io/
            snapshot.py        
            schema.md          
          config.py            
        api/
          server.py            
          protocol.md          
        requirements.txt
      frontend/
        src/
          App.tsx
          components/
            HeatmapCanvas.tsx
            GraphPreview.tsx
            Controls.tsx
          lib/
            ws.ts
          index.css
        package.json
        vite.config.ts
      docker/
        backend.Dockerfile
        frontend.Dockerfile
      README.md

### 2.2 Przepływ danych

-   **Backend** (GPU) wykonuje kroki symulacji i co `n_vis` kroków:
    -   publikuje przez **WS** ramkę `Float32Array` (np. gęstość,
        energia) + paczkę meta (np. licznik nowości),
    -   co `n_log` kroków wysyła listę **NOWYCH** substancji (ID,
        rozmiar, cechy) od poprzedniego ticku.\
-   **Frontend** odbiera strumień, renderuje heatmapy na `<canvas>`, a
    nowo wykryte substancje pokazuje jako **mini-grafy** w panelu
    „Nowości".

------------------------------------------------------------------------

## 3) Model symulacji

### 3.1 Tryb A (Preset Prebiotic --- walidacja)

-   **Pola stężeń** `conc[S, H, W]`, dyfuzja: dyskretny Laplasjan
    (4-sąsiedztwo).\
-   **Reakcje**: kilka równań masowo-akcyjnych (np. HCN→NH₂CHO +
    by-products), modyfikowane energią/katalizą.\
-   **Energia**: plamy energii podnoszą efektywną stałą reakcji.\
-   **Wizualizacja**: heatmapa NH₂CHO, HCN, energia.

> *Uwaga*: Tryb A służy tylko do weryfikacji „rurek" (GPU→WS→React) i
> sanity-check; nie ogranicza Trybu B.

### 3.2 Tryb B (Open Chemistry --- docelowy v1)

... (tu cała reszta planu jak w wiadomości powyżej) ...
