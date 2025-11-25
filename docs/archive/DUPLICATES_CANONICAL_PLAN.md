---
date: 2025-11-23
label: plan
---

# PLAN ARCHIWIZACJI DUPLIKATÓW - Wersje Kanoniczne

**Data:** 2025-11-23  
**Podstawa:** `docs/STRUCTURE_DEEP_ANALYSIS.md`  
**Zasady:** Zgodnie z `docs/ARCHIVE_POLICY.md` sekcja 6.4 - "Jeśli istnieje zduplikowany kod → wybierasz kanoniczny i archiwizujesz resztę"

---

## ⚠️ READ-ONLY ZONES (NIE DOTYKAMY)

Zgodnie z `ARCHIVE_POLICY.md`, **NIE PROponujemy** archiwizacji:
- `backend/sim/core/**`
- `backend/sim/chemistry/**`
- `scripts/run_phase2_full.py`
- `aws_test/configs/**`
- `docs/phase2b/**`
- `docs/technical/**`
- `results/**`

---

## 📋 1. DUPLIKATY `check_status*.ps1` (16 plików)

### 1.1 Analiza duplikatów

**Lokalizacje:**
- `backend/tests/check_status.ps1` do `check_status8.ps1` (8 plików)
- `tests/check_status.ps1` do `check_status8.ps1` (8 plików)

**Daty modyfikacji:**
- `backend/tests/check_status.ps1`: 2025-09-27 22:21:35
- `tests/check_status.ps1`: 2025-10-06 19:45:41

**Obserwacja:** `tests/check_status.ps1` jest nowszy (6.10.2025 vs 27.09.2025)

### 1.2 Wersje kanoniczne

| Plik | Status | Uzasadnienie |
|------|--------|--------------|
| `tests/check_status.ps1` | ✅ **KANONICZNY** | Nowsza data modyfikacji (6.10.2025), lokalizacja w głównym katalogu `tests/` |
| `backend/tests/check_status.ps1` | ⚠️ **DO ARCHIWUM** | Starsza wersja (27.09.2025), duplikat |

### 1.3 Wersje iteracyjne (2-8)

**Obserwacja:** Wersje 2-8 to iteracyjne wersje rozwojowe. `check_status.ps1` (wersja 1) jest prawdopodobnie finalną wersją.

| Plik | Status | Uzasadnienie |
|------|--------|--------------|
| `backend/tests/check_status2.ps1` do `check_status8.ps1` | ⚠️ **DO ARCHIWUM** | Iteracyjne wersje rozwojowe, starsza lokalizacja |
| `tests/check_status2.ps1` do `check_status8.ps1` | ⚠️ **DO ARCHIWUM** | Iteracyjne wersje rozwojowe (wersja 1 jest kanoniczna) |

### 1.4 Plan archiwizacji

| Ścieżka źródłowa | Nowa lokalizacja | Uzasadnienie |
|------------------|------------------|--------------|
| `backend/tests/check_status.ps1` | `archive/deprecated/backend_tests_check_status.ps1` | Starsza wersja, zastąpiona przez `tests/check_status.ps1` |
| `backend/tests/check_status2.ps1` | `archive/one_off_scripts/backend_tests_check_status2.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status3.ps1` | `archive/one_off_scripts/backend_tests_check_status3.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status4.ps1` | `archive/one_off_scripts/backend_tests_check_status4.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status5.ps1` | `archive/one_off_scripts/backend_tests_check_status5.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status6.ps1` | `archive/one_off_scripts/backend_tests_check_status6.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status7.ps1` | `archive/one_off_scripts/backend_tests_check_status7.ps1` | Iteracyjna wersja rozwojowa |
| `backend/tests/check_status8.ps1` | `archive/one_off_scripts/backend_tests_check_status8.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status2.ps1` | `archive/one_off_scripts/tests_check_status2.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status3.ps1` | `archive/one_off_scripts/tests_check_status3.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status4.ps1` | `archive/one_off_scripts/tests_check_status4.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status5.ps1` | `archive/one_off_scripts/tests_check_status5.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status6.ps1` | `archive/one_off_scripts/tests_check_status6.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status7.ps1` | `archive/one_off_scripts/tests_check_status7.ps1` | Iteracyjna wersja rozwojowa |
| `tests/check_status8.ps1` | `archive/one_off_scripts/tests_check_status8.ps1` | Iteracyjna wersja rozwojowa |

**Podsumowanie:**
- ✅ **KANONICZNY:** `tests/check_status.ps1` (zostaje)
- ⚠️ **DO ARCHIWUM:** 15 plików (1 deprecated + 14 one_off_scripts)

---

## 📋 2. DUPLIKATY `matcher.py` vs `matcher_v2.py`

### 2.1 Analiza duplikatów

**Lokalizacja:** `matcher/`  
**Pliki:**
- `matcher/matcher.py`
- `matcher/matcher_v2.py`

**Daty modyfikacji:**
- `matcher.py`: 2025-10-22 17:03:46
- `matcher_v2.py`: 2025-10-13 17:02:51

**Obserwacja:** `matcher.py` jest nowszy (22.10.2025 vs 13.10.2025), ale nazwa `matcher_v2.py` sugeruje że to nowsza wersja. Wymaga weryfikacji zawartości.

### 2.2 Wersje kanoniczne (DO WERYFIKACJI)

| Plik | Status | Uzasadnienie |
|------|--------|--------------|
| `matcher/matcher.py` | ⚠️ **DO WERYFIKACJI** | Nowsza data modyfikacji (22.10.2025), ale nazwa sugeruje że to stara wersja |
| `matcher/matcher_v2.py` | ⚠️ **DO WERYFIKACJI** | Starsza data modyfikacji (13.10.2025), ale nazwa sugeruje że to nowsza wersja |

**Uwaga:** Wymaga weryfikacji zawartości plików lub sprawdzenia w kodzie, która wersja jest używana.

### 2.3 Plan archiwizacji (DO POTWIERDZENIA)

**Opcja A:** Jeśli `matcher_v2.py` jest kanoniczny (nowsza wersja):
| Ścieżka źródłowa | Nowa lokalizacja | Uzasadnienie |
|------------------|------------------|--------------|
| `matcher/matcher.py` | `archive/deprecated/matcher/matcher.py` | Stara wersja, zastąpiona przez `matcher_v2.py` |

**Opcja B:** Jeśli `matcher.py` jest kanoniczny (nowsza wersja):
| Ścieżka źródłowa | Nowa lokalizacja | Uzasadnienie |
|------------------|------------------|--------------|
| `matcher/matcher_v2.py` | `archive/deprecated/matcher/matcher_v2.py` | Eksperymentalna wersja v2, zastąpiona przez `matcher.py` |

**Status:** ⚠️ **DO POTWIERDZENIA** - wymaga weryfikacji która wersja jest aktualnie używana

---

## 📋 3. DUPLIKATY `start_backend*.ps1`

### 3.1 Analiza duplikatów

**Lokalizacja:** Root  
**Pliki:**
- `start_backend.ps1`
- `start_backend_simple.ps1`

**Daty modyfikacji:**
- `start_backend.ps1`: 2025-09-30 22:34:14
- `start_backend_simple.ps1`: 2025-10-06 19:45:41

**Obserwacja:** `start_backend_simple.ps1` jest nowszy (6.10.2025), ale nazwa "simple" sugeruje że to uproszczona wersja.

### 3.2 Wersje kanoniczne

| Plik | Status | Uzasadnienie |
|------|--------|--------------|
| `start_backend.ps1` | ✅ **KANONICZNY** | Główna wersja, nazwa bez sufiksu "simple" sugeruje że to wersja produkcyjna |
| `start_backend_simple.ps1` | ⚠️ **DO ARCHIWUM** | Wersja "simple" - prawdopodobnie uproszczona wersja testowa lub jednorazowa |

### 3.3 Plan archiwizacji

| Ścieżka źródłowa | Nowa lokalizacja | Uzasadnienie |
|------------------|------------------|--------------|
| `start_backend_simple.ps1` | `archive/deprecated/start_backend_simple.ps1` | Wersja "simple" - uproszczona wersja, zastąpiona przez `start_backend.ps1` |

**Status:** ✅ **DO ARCHIWUM** (deprecated)

---

## 📋 4. DUPLIKAT `live2.0/backend/` vs `backend/`

### 4.1 Analiza duplikatów

**Lokalizacje:**
- `live2.0/backend/` (root) - zagnieżdżony katalog
- `backend/` (root) - główny katalog

**Obserwacja:** `live2.0/backend/` wygląda na błąd agenta - zagnieżdżony katalog projektu w samym projekcie.

### 4.2 Wersje kanoniczne

| Katalog | Status | Uzasadnienie |
|---------|--------|--------------|
| `backend/` (root) | ✅ **KANONICZNY** | Główny katalog projektu, standardowa lokalizacja |
| `live2.0/backend/` (root) | ⚠️ **DO ARCHIWUM** | Prawdopodobny błąd agenta, zagnieżdżony duplikat |

### 4.3 Plan archiwizacji (DO POTWIERDZENIA)

| Ścieżka źródłowa | Nowa lokalizacja | Uzasadnienie | Status |
|------------------|------------------|--------------|--------|
| `live2.0/` (root) | `archive/deprecated/live2.0/` | Zagnieżdżony katalog projektu - prawdopodobny błąd agenta. Zawiera duplikat `backend/`. **WYMAGA WERYFIKACJI** czy `live2.0/backend/snapshots/` zawiera unikalne debug snapshots | ⚠️ DO POTWIERDZENIA |

**Uwaga:** Przed archiwizacją wymaga weryfikacji:
1. Czy `live2.0/backend/snapshots/` zawiera unikalne dane
2. Czy `live2.0/backend/` to rzeczywiście duplikat czy ma unikalne pliki

**Status:** ⚠️ **DO POTWIERDZENIA** (wysoki priorytet)

---

## 📋 5. DUPLIKATY `fix_taichi_version.*`

### 5.1 Status

**Obserwacja:** Oba pliki zostały już zarchiwizowane w poprzedniej operacji:
- `fix_taichi_version.ps1` → `archive/one_off_scripts/`
- `fix_taichi_version.sh` → `archive/one_off_scripts/`

**Status:** ✅ **JUŻ ZARCHIWIZOWANE** - nie wymaga dalszych działań

---

## 📊 PODSUMOWANIE PLANU

### Wersje kanoniczne (zostają):
1. ✅ `tests/check_status.ps1` - kanoniczna wersja check_status
2. ✅ `start_backend.ps1` - kanoniczna wersja start_backend
3. ⚠️ `matcher/matcher.py` lub `matcher/matcher_v2.py` - **DO WERYFIKACJI**
4. ✅ `backend/` (root) - kanoniczny katalog backendu

### Do archiwizacji:

#### archive/deprecated/ (zastąpione wersje):
- `backend/tests/check_status.ps1` - starsza wersja, zastąpiona przez `tests/check_status.ps1`
- `start_backend_simple.ps1` - wersja "simple", zastąpiona przez `start_backend.ps1`
- `live2.0/` (root) - **DO POTWIERDZENIA** - zagnieżdżony duplikat backendu
- `matcher/matcher.py` LUB `matcher/matcher_v2.py` - **DO WERYFIKACJI** - jedna wersja jest zastąpiona

#### archive/one_off_scripts/ (iteracyjne wersje rozwojowe):
- `backend/tests/check_status2.ps1` do `check_status8.ps1` (7 plików)
- `tests/check_status2.ps1` do `check_status8.ps1` (7 plików)

### Statystyki:
- ✅ **KANONICZNE (zostają):** 3-4 elementy (w zależności od weryfikacji matcher)
- ⚠️ **DO ARCHIWUM (deprecated):** 2-3 elementy (w zależności od weryfikacji)
- ⚠️ **DO ARCHIWUM (one_off_scripts):** 14 plików (iteracyjne wersje)
- ⚠️ **DO POTWIERDZENIA:** 2 elementy (live2.0/, matcher)

---

## ⚠️ WAŻNE UWAGI

1. **Weryfikacja przed archiwizacją:**
   - `live2.0/backend/snapshots/` - sprawdzić czy zawiera unikalne debug snapshots
   - `matcher/matcher.py` vs `matcher_v2.py` - sprawdzić która wersja jest używana w kodzie

2. **Zasada kanoniczna (ARCHIVE_POLICY 6.4):**
   - Wybieramy kanoniczną wersję na podstawie: daty modyfikacji, nazwy, użycia w kodzie
   - Archiwizujemy resztę jako deprecated lub one_off_scripts

3. **Git mv:**
   - Wszystkie operacje będą wykonane przez `git mv` aby zachować historię

4. **ARCHIVE_LOG.md:**
   - Po wykonaniu planu zostanie zaktualizowany log archiwizacji

---

## ❓ CZY ZATWIERDZASZ PLAN?

**Plan zawiera:**
- ✅ **3-4 wersje kanoniczne** (zostają na miejscu)
- ⚠️ **16-17 elementów do archiwizacji** (2-3 deprecated + 14 one_off_scripts)
- ⚠️ **2 elementy DO POTWIERDZENIA** przed archiwizacją

**Proponowane działanie:**
1. **Zweryfikuj** `live2.0/backend/snapshots/` i `matcher/*.py`
2. **Zatwierdź** wersje kanoniczne
3. **Zatwierdź** plan archiwizacji
4. **Wykonanie** `git mv` dla zatwierdzonych elementów

**Czy zatwierdzasz plan?**

