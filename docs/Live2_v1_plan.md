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


3.2.1 Cząstki (SOA)

Każda cząstka ma atrybuty (ciągłe, konfigurowalne długości wektorów):

x⃗: pozycja (ciągła) 2D
v⃗: prędkość 2D
m: masa (skalar nieujemny)
q⃗: wektor „ładunków” (np. 2–6 składowych)
w⃗: wektor „walencji” / kierunkowości (np. 2–4 „orbitali” z orientacjami)
r0, ε: parametry potencjału podstawowego (zasięg, głębokość)
id_elem: identyfikator klasy „pierwiastka” (powstaje endogennie; patrz transmutacje)
cluster_id: identyfikator aktualnego klastra (-1 jeśli wolna)

3.2.2 Potencjały i wiązania
	•	Potencjał bazowy (np. Lennard-Jones-like) + terminy kierunkowe dopasowujące orientacje w⃗ (preferencje kątowe).
	•	Terminy wielociałowe (słabe) – stabilizują określone geometrie (np. trójkąty, łańcuchy).
	•	Binding rule: jeśli ΔE_bind < −θ_bind i lokalne ograniczenia walencji nie są przekroczone → utwórz krawędź.
	•	Breaking rule: jeśli ΔE_break > θ_break (np. silny impuls energii) → przerwij krawędź.

3.2.3 Listy sąsiadów
	•	Siatka przestrzenna / binning → O(N) budowa list na GPU, promień wyszukiwania r_cut.
	•	Aktualizacja co k kroków lub adaptacyjnie (na podstawie maksymalnego przesunięcia cząstek).

3.2.4 Klastrowanie i grafy
	•	Każdy klaster (połączony spójny komponent wiązań) ma graf: do limitu N_max węzłów (np. 12–20 w v1).
	•	Graf normalizujemy (kanoniczna etykietyzacja) i liczymy hash WL → ID substancji.
	•	Rejestr substancji (Catalog):
	•	jeśli hash nieznany → nowa substancja (zapis atrybutów: rozmiar, rozkład q⃗, średnie parametry).
	•	przechowujemy historię: czas pierwszego pojawienia, liczebność w bieżącym ticku, średni czas życia.

3.2.5 Transmutacje (narodziny „pierwiastków”)
	•	Jeżeli energy(x⃗,t) > E*, to z prawdopodobieństwem p_mut modyfikuj (mało!) atrybuty cząstki (q⃗, w⃗, r0, ε, m).
	•	Identyfikacja pierwiastka: cząstki z „bliskimi” atrybutami grupujemy (np. przez siatkowanie w przestrzeni atrybutów) → dostają ten sam id_elem.
	•	Konserwacja inwariantów (opcjonalnie): np. sum( wybranych składowych q⃗) ≈ const (w granicach tolerancji).

3.2.6 Energia
	•	Impulsy: co T_pulse kroków losowe plamy (promień R, amplituda A), tłumione współczynnikiem λ_decay.
	•	Efekt: (a) wzrost temperatury lokalnej → większe v⃗, (b) modyfikacja skutecznych θ_bind/θ_break, (c) wzrost p_mut.

3.2.7 Integracja ruchu
	•	Prosty Verlet / symplectic Euler z ograniczeniem prędkości (clamp) i adaptacyjnym dt (na podstawie maks. sił).

⸻

4) Metryki, logowanie, „novelty”

Per tick:
	•	Novelty(t): liczba NOWYCH ID substancji wykrytych w tym kroku (i/lub EMA z okna).
	•	Dynamika substancji: top-K najczęstszych, ich średni rozmiar, średni czas życia, modularność grafów.
	•	Spektrum złożoności: histogram rozmiarów grafów, cykli, rozkład stopni.
	•	Energia: suma energii potencjalnej/kinetycznej (przybliżona), rozkład po siatce.
	•	Transmutacje: liczba mutacji cząstek oraz „narodziny” nowych id_elem.

Zdarzenia specjalne (event log):
	•	(SUBSTANCE_NEW) {id, time, signature (rozmiar; cechy), seed}
	•	(ELEMENT_NEW) {id_elem, time, centroid atrybutów}
	•	(PEAK_ACTIVITY) {bbox, time, stats}

⸻

5) Interfejsy i protokoły

5.1 REST (FastAPI)
	•	GET /status → { mode: "A"|"B", step, H, W, seeds, metrics_summary }
	•	POST /control body: {action: "pause"|"resume"|"reset"|"snapshot"|"restore", payload?}
	•	POST /mode body: { mode: "A"|"B" } (przełącza tryb)
	•	GET /snapshot/{id} → plik (binarny)
	•	POST /params body: { ... } (zmiana wybranych parametrów runtime — whitelist)

5.2 WebSocket /ws
	•	FrameType 1: HEATMAP — nagłówek <u32 H, u32 W, u32 channel> + Float32Array(H*W)
	•	channels: 0=energy, 1=density (liczba cząstek w komórce), >=2 zarezerwowane
	•	FrameType 2: METRICS — JSON (krótki): {step, novelty, new_substances:[{id,size,lifespan0}], transmutations:int}
	•	FrameType 3: SUBSTANCE_SAMPLES — list(a) grafów do podglądu: {id, nodes:[{elem, attr...}], edges:[[i,j],...]} (ograniczyć do małych N, np. ≤12)

Wysyłamy naprzemiennie: HEATMAP co n_vis kroków i METRICS przy każdym.

⸻

6) Frontend (React + TS)
	•	HeatmapCanvas: renderuje Float32Array(H*W) → <canvas> (normalizacja lokalna).
	•	GraphPreview: SVG dla kilku nowych substancji (ID + graf + krótkie cechy).
	•	Controls: przyciski Pause/Resume/Reset, select Mode A/B; suwak intensywności impulsów (A) i częstości (T_pulse).
	•	Panel metryk: licznik novelty (ostatnie N), wykres liniowy (Chart.js), top-K substancji.

⸻

7) Konfiguracja (YAML)

grid:
  H: 256
  W: 256
  periodic: true

time:
  dt_init: 0.25
  dt_min: 0.05
  dt_max: 1.0
  adapt_strength: 0.6

open_chemistry:
  neighbor_radius: 3.0
  rebuild_neighbors_every: 10
  theta_bind: 0.8
  theta_break: 1.4
  vmax: 4.5
  clamp_density_per_cell: 64

energy:
  pulse_every: 120
  pulse_radius: 18
  pulse_amplitude: 2.5
  decay: 0.92
  p_mut_base: 1e-5
  E_star: 1.2
  p_mut_gain: 6.0

catalog:
  max_graph_size: 12
  wl_iterations: 3
  dedup_tol: 0.05

visualization:
  every_n_steps: 2
  metrics_every: 1
  sample_graphs_top_k: 6

rng:
  master_seed: 1337


⸻

😎 Algorytmy (pseudokod)

8.1 Krok symulacji (Tryb B)

for step in 1..:
  if step % rebuild_neighbors_every == 0:
    neighbors = build_neighbor_lists(particles, r_cut)

  // energia
  energy.apply_pulses(step)
  energy.decay()

  // siły międzycząsteczkowe
  F = compute_forces(particles, neighbors, potentials, energy_local)

  // integracja ruchu (adapt dt)
  dt = adapt_dt(F, dt_prev)
  integrate(particles, F, dt)

  // wiązania (binding/breaking)
  for each neighbor pair (i,j):
    dE = binding_energy_delta(i,j, potentials, energy_local)
    if dE < -theta_bind and valence_ok(i,j): make_bond(i,j)
    if dE >  theta_break: break_bond(i,j)

  // transmutacje (wysoka energia)
  for each particle in high_energy_regions:
    maybe_mutate_attributes(particle, p_mut(E), invariants)

  // klastrowanie i katalogowanie
  clusters = connected_components(bond_graph)
  new_substances = []
  for c in clusters with size <= max_graph_size:
    G = canonical_graph(c)
    id = wl_hash(G)
    if !catalog.contains(id): catalog.add(id, G); new_substances.append(id)

  // metryki i streaming
  novelty = |new_substances|
  push_ws_heatmap(...)
  push_ws_metrics(novelty, new_substances, ...)


⸻

9) Testy i kryteria akceptacji

9.1 Jednostkowe
	•	graphs.hash: identyczne grafy → identyczny hash; permutacja węzłów nie zmienia hash.
	•	catalog: dodanie nowej substancji rejestruje ją dokładnie raz.
	•	binding/breaking: progi θ działają deterministycznie w warunkach testowych.

9.2 Property-based
	•	Inwarianty: przy braku impulsów energii i mutacji, energia/cząstki nie „eksplodują” (limity).
	•	Lokalność: siła = 0 dla par poza r_cut.
	•	Stabilność numeryczna: dt adaptuje się do zakresu [dt_min, dt_max].

9.3 Wydajność
	•	N=200k cząstek na GPU: utrzymanie ≥ 10 steps/s w konfiguracji domyślnej (orientacyjnie; dopuszczalne 5–10 w zależności od GPU).
	•	Budowa list sąsiadów O(N), koszt kroków ~ liniowy.

9.4 Długi bieg (8–24h)
	•	Novelty(t) nie spada do zera po krótkim czasie (utrzymuje się > 0 w oknach), co wskazuje na otwartość.
	•	Snapshot+restore odtwarza stan deterministycznie (dla tych samych seedów).

⸻

10) Plan prac (kamienie milowe)

M1 — Skeleton & Tryb A (2–3 dni)
	•	FastAPI + WS, React canvas, GPU init (Taichi), siatka 2D, dyfuzja + 1 reakcja, impulsy energii, Pause/Resume.

M2 — Open Chemistry Core (5–7 dni)
	•	particles.py (SOA), neighbors, integracja ruchu, potencjały, binding/breaking, adapt dt, ograniczenia gęstości.

M3 — Katalog substancji (3–4 dni)
	•	Grafy klastrów (≤N_max), WL-hash, rejestr, metryki novelty, event log, WS METRICS.

M4 — UI „Nowości” (2 dni)
	•	Panel metryk + GraphPreview (SVG), konfiguracja trybów i parametrów (min. energia, częstość impulsów).

M5 — Snapshot/Restore + Testy (3–4 dni)
	•	Snapshot (particles, bonds, seeds, params), property tests, testy wydajności.

M6 — Raport v1 i preset experiment packs (2 dni)
	•	Zestaw presetów (YAML), skrypty uruchomieniowe, dokumentacja uruchomienia i eksperymentów.

⸻

11) Uruchomienie lokalne

Backend

conda create -n live2 python=3.11 -y
conda activate live2
pip install taichi fastapi uvicorn[standard] numpy pydantic msgpack
uvicorn api.server:app --port 8008 --reload

Frontend

cd frontend
npm i
npm run dev
# http://localhost:5173


⸻

12) Deploy (AWS — szkic)
	•	EC2 GPU (g5/g6) z AMI CUDA; docker compose dla backend/frontend.
	•	NGINX: terminacja TLS, reverse proxy do WS.
	•	S3: snapshoty + logi; CloudWatch: metryki.
	•	(Opcja) ECR + ECS/EKS dla skalowania eksperymentów; worker’y GPU.

⸻

13) Ryzyka i mitigacje
	•	Eksplozja złożoności: limit max_graph_size, clamp gęstości, adapt dt.
	•	Brak emergencji: regulacja energii (A, T_pulse), p_mut, kształtu potencjałów.
	•	Wydajność: ograniczyć r_cut, rebuild co k kroków, kompresja ramek, zmniejszyć H×W.
	•	Zabetonowanie: unikać twardych list; potencjały i atrybuty parametryzowane, transmutacje aktywne.

⸻

14) Zadania dla agenta (Cursor AI) — checklist
	1.	Utwórz monorepo wg struktury z pkt. 2.1.
	2.	Zaimplementuj Tryb A (fields_ca, energy, stepper, WS/React); potwierdź wizualizację.
	3.	Dodaj Tryb B: particles, neighbors, potentials, binding, graphs, catalog, metrics.
	4.	Implementuj transmutacje atrybutów warunkowane energią.
	5.	Dodaj METRICS WS + panel React (novelty + top-K substancji).
	6.	Wykonaj snapshot/restore i testy (unit + property).
	7.	Dostarcz preset YAML + skrypty uruchomieniowe.
	8.	Profiling i krótkie strojenie parametrów, raport v1.

⸻

15) Kryteria „DONE” (akceptacja v1)
	•	System działa w Trybie B, generuje nowe substancje (ID pojawiają się w czasie), novelty>0 w dłuższym biegu.
	•	Frontend pokazuje heatmapy oraz mini-grafy nowo odkrytych substancji.
	•	Snapshot/restore odtwarza stan; metryki wysyłane regularnie; sterowanie pauzą/trybem działa.
	•	Testy bazowe przechodzą; wydajność akceptowalna na lokalnym GPU (≥ ~5–10 steps/s przy N≈100–200k cząstek, zależnie od GPU i ustawień).

⸻

16) Dalsze kierunki (v1.1+)
	•	Kompartmenty (pęcherzyki/błony) jako nowy typ struktury stabilizującej reakcje.
	•	Różne „klimaty” przestrzenne (gradienty parametrów potencjałów).
	•	Replikatory symboliczne (łańcuchy z regułami kopiowania) jako prototyp „genów”.
	•	Search of novelty sterujący seedami i parametrami (automatyczna eksploracja).

⸻

Uwagi końcowe:
	•	Ten projekt jest eksperymentem naukowym — oczekujemy nieprzewidywalności.
	•	Priorytetem jest otwartość: minimalna liczba twardych reguł, maksimum parametryzacji i lokalnych mechanizmów.
	•	v1 ma dowieść, że system potrafi wytwarzać nowość (substancje i „pierwiastki”) bez dopisywania list reakcji.
