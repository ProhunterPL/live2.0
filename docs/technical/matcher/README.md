# Matcher - Identyfikacja Molekuł

Dokumentacja systemu rozpoznawania molekuł (matcher).

---

## 📄 Dokumenty

### [PUBCHEM_MATCHER_FIX.md](PUBCHEM_MATCHER_FIX.md)
Naprawa integracji z PubChem:
- Problem: Matcher nie znajdował molekuł w PubChem
- Rozwiązanie: Ulepszone zapytania i fallback
- Nowy algorytm wyszukiwania

### [MATCHER_FIX_SUMMARY.md](MATCHER_FIX_SUMMARY.md)
Podsumowanie napraw:
- Historia zmian w matcher
- Aktualne możliwości
- Przykłady użycia

---

## 🔍 Jak Działa Matcher

### 1. Detekcja Klastrów
Symulacja wykrywa klastry cząsteczek:
- Minimum 2 atomy
- Stabilne przez określony czas
- Graph representation (bonds between atoms)

### 2. Konwersja do Formatu Chemicznego
Klaster → Molekuła:
- Generacja SMILES
- Generacja InChI
- Tworzenie plików .mol, .xyz

### 3. Wyszukiwanie w Bazach
- **PubChem** - główna baza
- **Local database** - cache'owane wyniki
- **Similarity search** - podobne molekuły

### 4. Wizualizacja
Wyniki:
- `matches/cluster_*.png` - wizualizacja
- `matches/cluster_*.mol` - struktura
- `matches/cluster_*.xyz` - geometria

---

## 🛠️ Pliki Źródłowe

### Główne Moduły:
- `matcher/matcher.py` - główny matcher
- `matcher/matcher_v2.py` - ulepszona wersja
- `matcher/chem.py` - konwersje chemiczne
- `matcher/similarity.py` - podobieństwo molekuł

### Frontend:
- `frontend/src/components/NoveltyPanel.tsx` - UI dla matcher

---

## 📊 Przykłady Wykrytych Molekuł

Z sesji 2024-10-22:
- **NH₃** (Amoniak) - 4 atomy
- **N₃H₃** (Cykliczny trimer azotu) - 6 atomów

Oczekiwane po optymalizacji:
- Glikol (C₂H₆O₂)
- Formamid (CH₃NO)
- Mocznik (CH₄N₂O)
- HCN
- Formaldehyd

---

## 🔗 Zobacz Też

- [Matcher V2](../../MATCHER_V2.md) - Dokumentacja v2
- [Novelty Panel](../../README_MATCHER.md) - UI documentation

