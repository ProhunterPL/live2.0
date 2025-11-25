---
date: 2025-11-23
label: manual
---
# ARCHIVE POLICY – Live 2.0  
**Bezpieczne zasady sprzątania, porządkowania i wersjonowania repozytorium**

Projekt Live 2.0 jest długoterminowym projektem naukowym (Phase 2B → Phase 3 → publikacja).  
Repo zawiera wiele eksperymentów, wyników, starych wersji kodu i dokumentacji.  

Celem tego dokumentu jest ustalenie **jasnych zasad archiwizacji**, aby:

- niczego **nie utracić**,  
- zachować **odtwarzalność** eksperymentów,  
- utrzymać **porządek** w repo,  
- umożliwić efektywną pracę agentom AI i developerom  
- nie uszkodzić **Phase 2B** ani materiałów do publikacji.

---

# 📁 1. Główna idea archiwizacji

## 🧩 **NIC nie usuwamy.**  
W Live 2.0 *nic* nie jest kasowane bez dodatkowej decyzji naukowej.

Zamiast usuwania stosujemy:

- **przenoszenie (`git mv`)** do dedykowanej przestrzeni `archive/`
- **wskazanie pliku „kanonicznego”** (ten właściwy, aktualny)
- **oznaczenie plików jako legacy** (te stare, zachowane dla historii)
- **adnotację**, dlaczego dany plik został zarchiwizowany

---

# 📂 2. Struktura katalogu `archive/`

Suggestowana struktura — jeśli katalog jeszcze nie istnieje, utwórz go:

archive/
│
├── one_off_scripts/ # jednorazowe, debug, testy chwilowe, kopie
│
├── old_docs/ # dokumenty zastąpione nowymi w docs/
│
├── experiments/ # prototypy, próbne wersje algorytmów, alternatywne pipeline’y
│
├── tmp_results/ # wyniki nieużywane do publikacji / analizy
│
└── deprecated/ # kod/konfiguracje zastąpione finalną wersją


> Uwaga: nie mieszamy tych kategorii.  
> Każdy z katalogów ma jasne przeznaczenie.

---

# 🟥 3. Czego **nie wolno** archiwizować (READ-ONLY)

Te elementy są krytyczne naukowo i muszą zostać na miejscu:

### 🔬 CORE SIMULATION ENGINE  


backend/sim/core/**
backend/sim/chemistry/**
backend/sim/config.py
backend/sim/molecule_extractor.py


### 🔁 PHASE 2B RUNTIME  


scripts/run_phase2_full.py
scripts/phase2.py
aws_test/scripts/**
aws_test/configs/phase2_*.yaml
aws_test/configs/SUPER_FAST.yaml


### 📚 SCIENTIFIC DOCUMENTATION  


docs/phase2b/**
docs/technical/**
docs/VALIDATION_ROADMAP.md
docs/SCIENTIFIC_VALIDITY_ANALYSIS.md
docs/NAVIGATION_GUIDE_LIVE2.md


### 📦 RESULTS  
(tylko do archiwów offline, nigdy usuwać/przenosić w repo)


results/**
phase2b_results/**
run/snapshots/**
run_*/checkpoints/**


**Nigdy nie przenosimy powyższych elementów do `archive/`.**  
**Nigdy nie usuwamy ich.**  
**Każda zmiana = potwierdzenie od Michała.**

---

# 🟦 4. Co można archiwizować

## 4.1 **Jednorazowe skrypty** (idealne do `archive/one_off_scripts/`)
- skrypty debugujące, robione ad hoc:
  - np. `test_molecule_extraction_x.py`
  - `debug_cluster_issue.py`
  - `old_monitoring_2025_10_12.py`
- stare workflowy do analizy, nieużywane
- eksperymentalne skrypty porzucone po wprowadzeniu wersji finalnych

> Zasada: jeśli skrypt był użyty raz i faza eksperymentu zakończona → archiwizujemy.

---

## 4.2 **Dokumenty zastąpione nowszymi wersjami** (`archive/old_docs/`)

Przykłady:

- starsze wersje raportów, np. `report_v1.md` po stworzeniu `report_v3.md`
- zduplikowane dokumenty, np. różne analizy tych samych wyników
- stare TODO / sesje notatek, których treść została przeniesiona do właściwego `docs/`

---

## 4.3 **Eksperymenty, prototypy, alternatywne algorytmy** (`archive/experiments/`)

Przykłady:

- alternatywny cluster detector (stara wersja)
- alternatywne steppers
- “test ml approach”, “experimental chemical rules”
- prototypy do Phase 3, które nie są częścią Phase 2B

> Zasada: jeżeli coś *mogłoby* być wartościowe historycznie lub porównawczo → archiwizujemy, nie kasujemy.

---

## 4.4 **Wyniki i artefakty nieużywane w publikacji** (`archive/tmp_results/`)

- częściowe wyniki runów
- przerwane symulacje
- dane referencyjne, które są za duże na repo
- fragmentaryczne snapshoty

**Uwaga:**  
Jeśli wyniki są duże → przenosimy poza repo (np. lokalnie / cloud storage)  
i wstawiamy w repo jedynie link / README o lokalizacji.

---

# 🟩 5. Proces archiwizacji krok po kroku

## Krok 1 – AI lub człowiek identyfikuje pliki do archiwizacji  
Przykład listy:

- `scripts/debug_old_chemical_rules.py`
- `docs/session_notes_2025_10_04.md`
- `scripts/new_cluster_algorithm_v2_backup.py`

## Krok 2 – Agent generuje **plan archiwizacji**
Lista:

1. Dokładne ścieżki plików  
2. Proponowane nowe lokalizacje  
3. Uzasadnienie kategorii (one-off, old-docs, experiments, results)

**Agent nie rusza plików bez akceptacji!**

## Krok 3 – Po akceptacji wykonywane są operacje:


git mv scripts/debug_old_chemical_rules.py archive/one_off_scripts/
git mv docs/session_notes_2025_10_04.md archive/old_docs/


## Krok 4 – Agent tworzy plik `ARCHIVE_LOG.md`  
Format wpisu:


[2025-11-23] Archived: scripts/debug_old_chemical_rules.py

Reason: One-off debugging script, replaced by new chemistry validation in Phase 2B.
Moved to: archive/one_off_scripts/


---

# 🟧 6. Zasady dla agentów AI (Cursor, GPT, itp.)

1. **Zawsze pokazujesz PLAN przed przenoszeniem.**
2. **Nigdy nie archiwizujesz CORE / Phase2B / Scientific Docs / Results.**
3. **Jeśli istnieje nowszy dokument → starszy trafia do `archive/old_docs/`.**
4. **Jeśli istnieje zduplikowany kod → wybierasz kanoniczny i archiwizujesz resztę.**
5. **Jeśli plik ma charakter jednorazowy lub testowy → do `archive/one_off_scripts/`.**
6. **Jeśli plik to prototyp algorytmu → do `archive/experiments/`.**
7. **Nie tworzymy nowego bałaganu: nowe prototypy → od razu do `archive/experiments/`.**

---

# 🧭 7. Quick Reference

| Typ pliku | Gdzie przenieść |
|-----------|------------------|
| Debug scripts | `archive/one_off_scripts/` |
| Stare dokumenty | `archive/old_docs/` |
| Prototypy | `archive/experiments/` |
| Tymczasowe wyniki | `archive/tmp_results/` |
| Duże dane | poza repo + README z lokalizacją |
| CORE / Phase2B | **ZAKAZ ARCHIWIZACJI** |

---

# 🏁 8. Final Notes

- Zasady te są obowiązujące do zakończenia **Phase 3 (paper)**.
- Po publikacji można wprowadzić wersję “cleanup v2”, ale tylko po backupie.
- To repo to kombinacja **kod + nauka + HPC** → zachowujemy całą historię.

---

**Maintainer:** Michał Klawikowski  
**Last Updated:** 2025-11-23  
**Purpose:** Bezpieczne utrzymanie repo i praca agentów AI bez ryzyka utraty pracy naukowej