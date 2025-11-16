# 🚀 CI/CD Quick Reference

## Szybkie Komendy

### 🧪 Testy Lokalne

```bash
# Uruchom wszystko (jak w CI)
bash .github/scripts/run_local_ci.sh                # Linux/Mac
.\.github\scripts\run_local_ci.ps1                  # Windows

# Backend testy (szybkie)
cd backend && pytest tests/ -v -m "not slow"

# Root testy (bez długich)
pytest tests/ -v -k "not stability and not 24h"

# Konkretny test
pytest tests/test_config.py::test_default_config -v
```

### 🎨 Formatowanie

```bash
# Check (nie zmienia plików)
black --check backend/ scripts/ matcher/
isort --check-only backend/ scripts/ matcher/

# Fix (naprawia automatycznie)
black backend/ scripts/ matcher/
isort backend/ scripts/ matcher/
```

### 📋 Test Markers

```bash
pytest -m unit              # Tylko unit tests
pytest -m integration       # Tylko integration tests
pytest -m "not slow"        # Bez slow tests
pytest -m "unit and not slow"  # Unit bez slow
```

### 🔍 Debugging

```bash
pytest -v                   # Verbose
pytest -s                   # Pokaż printy
pytest -x                   # Stop na pierwszym błędzie
pytest -l                   # Pokaż local variables
pytest --tb=long            # Pełny traceback
pytest --pdb                # Debugger przy failure
```

## Workflow Status

| Status | Znaczenie | Akcja |
|--------|-----------|-------|
| ✅ passing | Wszystko OK | Możesz merge |
| ❌ failing | Coś nie działa | Sprawdź logi |
| 🟡 pending | W trakcie | Poczekaj |
| ⚪ skipped | Pominięty | Normalnie |

## Markery w Testach

```python
@pytest.mark.unit           # Test jednostkowy
@pytest.mark.integration    # Test integracyjny  
@pytest.mark.slow           # Test długi (pominięty w CI)
@pytest.mark.skip           # Zawsze pomiń
@pytest.mark.skipif(...)    # Pomiń warunkowo
```

## CI Jobs

| Job | Kiedy | Czas | Co robi |
|-----|-------|------|---------|
| **test** | PR + push main | 10-15 min | Unit tests (Py 3.10, 3.11) |
| **integration-tests** | tylko push main | 30-45 min | Pełne testy |
| **code-quality** | PR + push main | 3-5 min | Black, isort, mypy |

## Troubleshooting Quick Fixes

### "Test failed in CI but passes locally"
```bash
export TI_ARCH=cpu
pip install --force-reinstall -r requirements.txt
pytest tests/ -v -m "not slow"
```

### "Import error: No module named 'sim'"
```bash
export PYTHONPATH=$(pwd)/backend  # dla backend tests
export PYTHONPATH=$(pwd)          # dla root tests
```

### "Black/isort failures"
```bash
black backend/ scripts/ matcher/
isort backend/ scripts/ matcher/
git add . && git commit -m "style: fix formatting"
```

### "Taichi error in CI"
```bash
# CI zawsze używa CPU mode
export TI_ARCH=cpu
pytest tests/
```

## Git Workflow z CI

```bash
# 1. Nowy branch
git checkout -b feature/my-feature

# 2. Zmiany
# ... edytuj pliki ...

# 3. Test lokalnie
bash .github/scripts/run_local_ci.sh

# 4. Commit
git add .
git commit -m "feat: add new feature"

# 5. Push
git push origin feature/my-feature

# 6. Otwórz PR na GitHub
# CI uruchomi się automatycznie

# 7. Po przejściu testów - merge do main
```

## Conventional Commits

| Prefix | Użycie | Przykład |
|--------|--------|----------|
| `feat:` | Nowa funkcjonalność | `feat: add particle collision` |
| `fix:` | Naprawa błędu | `fix: correct energy calculation` |
| `docs:` | Dokumentacja | `docs: update API reference` |
| `style:` | Formatowanie | `style: fix indentation` |
| `refactor:` | Refactoring | `refactor: simplify binding logic` |
| `test:` | Testy | `test: add integration tests` |
| `perf:` | Wydajność | `perf: optimize grid search` |
| `chore:` | Maintenance | `chore: update dependencies` |

## Przydatne Linki

- **Actions Tab:** `github.com/USER/live2.0/actions`
- **Workflow File:** `.github/workflows/ci-tests.yml`
- **Full Guide:** `docs/CI_CD_GUIDE.md`
- **Workflow README:** `.github/workflows/README.md`

## Environment Variables

```bash
# Taichi CPU mode (required for CI)
export TI_ARCH=cpu

# Python path
export PYTHONPATH=$(pwd)
export PYTHONPATH=$(pwd)/backend

# Pytest verbosity
export PYTEST_CURRENT_TEST=1
```

## Badge dla README

```markdown
![CI Tests](https://github.com/USER/live2.0/workflows/CI%20Tests/badge.svg)
```

## Szybka Diagnoza

```bash
# Sprawdź co failuje
pytest tests/ -v --tb=short

# Uruchom tylko failures z ostatniego razu
pytest --lf

# Uruchom failures i kolejny test
pytest --lf --ff

# Coverage report
pytest --cov=backend/sim --cov-report=term-missing

# Timing (które testy są najwolniejsze)
pytest --durations=10
```

## Pre-commit Hook (Opcjonalny)

Stwórz `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Auto-format przed commitem
black backend/ scripts/ matcher/
isort backend/ scripts/ matcher/
git add -u
```

```bash
chmod +x .git/hooks/pre-commit
```

---

**Szybka pomoc:** `pytest --help` | `pytest --markers` | `pytest --fixtures`

