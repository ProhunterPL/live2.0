---
date: 2025-11-23
label: guide
---

# Archive Directory

Ten katalog zawiera **zarchiwizowane pliki** z projektu Live 2.0, które nie są już aktywnie używane, ale są zachowane dla:
- historii projektu,
- możliwości odtworzenia wyników,
- referencji w przyszłości.

## 📁 Struktura

- **`one_off_scripts/`** - jednorazowe skrypty, debug, testy chwilowe, kopie
- **`old_docs/`** - dokumenty zastąpione nowymi w `docs/`
- **`experiments/`** - prototypy, próbne wersje algorytmów, alternatywne pipeline'y
- **`tmp_results/`** - wyniki nieużywane do publikacji / analizy
- **`deprecated/`** - kod/konfiguracje zastąpione finalną wersją

## ⚠️ Zasady

1. **Nigdy nie usuwamy** plików z archiwum bez wyraźnej zgody
2. Używamy **`git mv`** do przenoszenia plików tutaj (zachowuje historię)
3. Przed przeniesieniem do archiwum, agent musi pokazać plan użytkownikowi
4. **Nie przenosimy** plików z read-only zones (CORE, Phase 2B, RESULTS)

## 🔍 Jak używać

Gdy chcesz przenieść plik do archiwum:
```bash
git mv <ścieżka_pliku> archive/<odpowiedni_podkatalog>/
```

---

*Zobacz `docs/NAVIGATION_GUIDE.md` sekcja 3 (ARCHIVE POLICY) dla szczegółów.*

