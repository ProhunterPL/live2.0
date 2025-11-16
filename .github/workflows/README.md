# GitHub Actions CI/CD

Ten folder zawiera konfiguracje GitHub Actions dla automatycznego testowania kodu.

## 🚀 Workflows

### `ci-tests.yml` - Główny pipeline testów

**Triggery:**
- Push do brancha `main`
- Pull Request do brancha `main`

**Jobs:**

1. **test** (Matrix: Python 3.10, 3.11)
   - Uruchamia szybkie testy jednostkowe
   - Pomija testy oznaczone jako `slow`
   - Timeout: 30 minut
   - Używa Taichi w trybie CPU

2. **integration-tests** (tylko po merge do main)
   - Uruchamia testy integracyjne
   - Timeout: 60 minut
   - Wykonuje się tylko na pushach do main

3. **code-quality**
   - Black (formatowanie kodu)
   - isort (sortowanie importów)
   - mypy (sprawdzanie typów - nie blokuje)

## 📋 Wymagania

- Python 3.10 lub 3.11
- Ubuntu latest (GitHub-hosted runner)
- Wszystkie zależności z `requirements.txt`

## 🔧 Lokalne testowanie

### Uruchomienie testów jak w CI:

```bash
# Backend unit tests (szybkie)
cd backend
pytest tests/ -v -m "not slow" --tb=short --color=yes

# Root tests (bez długich testów stabilności)
pytest tests/ -v -k "not stability and not 24h" --tb=short --color=yes
```

### Sprawdzenie formatowania:

```bash
# Formatowanie
black --check backend/ scripts/ matcher/

# Sortowanie importów
isort --check-only backend/ scripts/ matcher/

# Typowanie
mypy backend/sim/ --ignore-missing-imports
```

### Naprawienie formatowania:

```bash
# Auto-fix
black backend/ scripts/ matcher/
isort backend/ scripts/ matcher/
```

## 🎯 Statusy testów

Po każdym pushu do main lub utworzeniu PR, GitHub automatycznie:

1. ✅ Uruchomi wszystkie testy jednostkowe
2. ✅ Sprawdzi jakość kodu
3. ✅ (Po merge) Uruchomi testy integracyjne

## 🐛 Troubleshooting

### Testy failują lokalnie ale nie w CI
- Sprawdź czy używasz `TI_ARCH=cpu` dla Taichi
- Upewnij się że masz te same wersje zależności: `pip install -r requirements.txt`

### Testy failują w CI ale nie lokalnie
- CI używa Ubuntu, lokalnie może być Windows/Mac
- CI nie ma GPU, sprawdź czy test wymaga CPU mode
- Sprawdź logi w zakładce "Actions" w GitHub

### RDKit installation issues
- RDKit może mieć problemy na niektórych systemach
- W CI używamy standardowego `rdkit` package
- Jeśli wystąpią problemy, można stworzyć `requirements-ci.txt` bez RDKit

## 📊 Badge statusu

Możesz dodać badge do README.md:

```markdown
![CI Tests](https://github.com/TWOJ-USERNAME/live2.0/workflows/CI%20Tests/badge.svg)
```

## 🔐 Secrets i konfiguracja

Obecnie workflow nie wymaga żadnych secrets. Jeśli w przyszłości będziesz potrzebować:

1. Przejdź do Settings → Secrets and variables → Actions
2. Dodaj nowy secret
3. Użyj w workflow: `${{ secrets.SECRET_NAME }}`

## ⚡ Optymalizacja

Obecna konfiguracja:
- Używa cache dla pip dependencies
- Uruchamia testy na matrix (Python 3.10 i 3.11)
- Zapisuje artefakty z logów testów (7 dni retention)
- Timeout zabezpiecza przed zablokowanymi testami

## 📝 Dodatkowe zasoby

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [pytest Documentation](https://docs.pytest.org/)
- [Taichi CI/CD](https://docs.taichi-lang.org/docs/hello_world#cpu-and-gpu)

