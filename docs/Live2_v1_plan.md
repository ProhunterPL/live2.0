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

------------------------------------------------------------------------

## 1) Zakres v1 (MVP)

**Musi zawierać**: 1. **Siatka 2D** H×W (domyślnie 256×256), periodyczne
brzegi.\
2. **Dwa tryby uruchomienia**: - **A) Preset Prebiotic**: ciągłe
stężenia kilku „chemikaliów" + parę reakcji (np. HCN→NH₂CHO) do
weryfikacji pipeline'u i wizualizacji. - **B) Open Chemistry**
*(domyślny)*: cząsteczki off-lattice, potencjały, tworzenie/zrywanie
wiązań, **rejestr nowych substancji** (hash grafu).\
3. **Źródła energii**: impulsy czasoprzestrzenne (UV/błyskawice),
rozchodzące się plamy energii z zanikiem.\
4. **Transmutacje**: rzadkie mutacje atrybutów cząstek w obszarach
wysokiej energii → **narodziny nowych klas „pierwiastków"**.\
5. **Detektory nowości**: liczenie i logowanie pojawień się **nowych
substancji** (ID grafu), ich rozmiaru, czasu życia; metryka novelty
vs. historia.\
6. **Frontend** (React + TS): heatmapy (gęstość/energia), licznik
nowości, podgląd kilku wylosowanych klastrów jako grafów (SVG),
pauza/wznowienie, wybór trybu.\
7. **API** (FastAPI + WebSocket): stream ramek (binarny), endpointy
sterujące, zrzuty snapshotów.\
8. **Snapshoty**: zapis/odczyt stanu (particles, grafy, parametry,
seed).\
9. **Testy**: jednostkowe (hash grafu, katalog), property-based
(inwarianty), testy wydajności (profiling), testy długiego biegu
(utrzymanie novelty\>0 w czasie).

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
