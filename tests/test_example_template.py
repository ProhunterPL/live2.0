"""
Przykładowy template dla testów w Live 2.0

Ten plik pokazuje jak pisać różne typy testów.
Możesz go skopiować i zmodyfikować dla swoich potrzeb.
"""

import pytest
import numpy as np

# ============================================================================
# SIMPLE UNIT TESTS - Podstawowe testy jednostkowe
# ============================================================================

def test_basic_assertion():
    """Najprostszy możliwy test"""
    assert 1 + 1 == 2


def test_with_variables():
    """Test z użyciem zmiennych"""
    expected = 42
    actual = 40 + 2
    assert actual == expected


def test_numpy_arrays():
    """Test z użyciem numpy"""
    arr = np.array([1, 2, 3])
    assert len(arr) == 3
    assert arr.sum() == 6
    assert arr.dtype == np.int_


# ============================================================================
# FIXTURES - Współdzielone dane dla testów
# ============================================================================

@pytest.fixture
def sample_config():
    """Fixture zwracający przykładową konfigurację"""
    return {
        "n_particles": 100,
        "box_size": 50.0,
        "dt": 0.01
    }


@pytest.fixture
def sample_particles():
    """Fixture zwracający przykładowe dane cząstek"""
    positions = np.random.rand(10, 3) * 50.0
    velocities = np.random.randn(10, 3)
    return {
        "positions": positions,
        "velocities": velocities
    }


def test_with_fixture(sample_config):
    """Test używający fixture"""
    assert sample_config["n_particles"] == 100
    assert sample_config["dt"] > 0


def test_with_multiple_fixtures(sample_config, sample_particles):
    """Test używający wielu fixtures"""
    n_particles = sample_config["n_particles"]
    positions = sample_particles["positions"]
    assert len(positions) <= n_particles


# ============================================================================
# PARAMETRIZED TESTS - Uruchamianie tego samego testu z różnymi danymi
# ============================================================================

@pytest.mark.parametrize("input,expected", [
    (1, 2),
    (2, 4),
    (3, 6),
    (10, 20),
])
def test_double_function(input, expected):
    """Test z wieloma przypadkami testowymi"""
    assert input * 2 == expected


@pytest.mark.parametrize("element,expected_valence", [
    (1, 1),   # H - wodór
    (6, 4),   # C - węgiel
    (7, 3),   # N - azot
    (8, 2),   # O - tlen
])
def test_element_valence(element, expected_valence):
    """Przykład testu parametryzowanego dla chemii"""
    # To jest tylko przykład - faktyczna implementacja byłaby inna
    valences = {1: 1, 6: 4, 7: 3, 8: 2}
    assert valences.get(element) == expected_valence


# ============================================================================
# MARKERS - Kategoryzacja testów
# ============================================================================

@pytest.mark.unit
def test_unit_example():
    """Test jednostkowy (unit)"""
    assert True


@pytest.mark.integration
def test_integration_example():
    """Test integracyjny - testuje wiele komponentów razem"""
    # Ten test uruchomi się tylko w jobie integration-tests
    assert True


@pytest.mark.slow
def test_slow_simulation():
    """Test oznaczony jako slow - pominięty w szybkim CI"""
    # Symulacja długotrwałego procesu
    for i in range(1000):
        _ = i ** 2
    assert True


# ============================================================================
# EXCEPTION TESTING - Testowanie czy funkcja rzuca wyjątki
# ============================================================================

def divide(a, b):
    """Helper function for testing exceptions"""
    if b == 0:
        raise ValueError("Cannot divide by zero")
    return a / b


def test_exception_raised():
    """Test sprawdzający czy wyjątek jest rzucany"""
    with pytest.raises(ValueError):
        divide(10, 0)


def test_exception_message():
    """Test sprawdzający treść wyjątku"""
    with pytest.raises(ValueError, match="Cannot divide by zero"):
        divide(10, 0)


def test_no_exception():
    """Test sprawdzający że NIE ma wyjątku"""
    result = divide(10, 2)
    assert result == 5.0


# ============================================================================
# APPROXIMATE COMPARISONS - Porównania z tolerancją
# ============================================================================

def test_approximate_float():
    """Test z przybliżonym porównaniem float"""
    result = 0.1 + 0.2
    expected = 0.3
    assert result == pytest.approx(expected)


def test_approximate_array():
    """Test z przybliżonym porównaniem numpy array"""
    result = np.array([0.1 + 0.2, 0.3 + 0.4])
    expected = np.array([0.3, 0.7])
    np.testing.assert_allclose(result, expected)


# ============================================================================
# CONDITIONAL TESTS - Testy warunkowe
# ============================================================================

@pytest.mark.skipif(
    np.__version__ < "1.20.0",
    reason="Requires numpy >= 1.20.0"
)
def test_requires_new_numpy():
    """Test który wymaga nowszej wersji numpy"""
    assert True


try:
    import taichi as ti
    HAS_TAICHI = True
except ImportError:
    HAS_TAICHI = False


@pytest.mark.skipif(not HAS_TAICHI, reason="Requires Taichi")
def test_requires_taichi():
    """Test który wymaga Taichi"""
    assert HAS_TAICHI


# ============================================================================
# CLEANUP - Sprzątanie po testach
# ============================================================================

@pytest.fixture
def temp_file(tmp_path):
    """Fixture tworzący tymczasowy plik (automatycznie usuwany)"""
    file_path = tmp_path / "test_data.txt"
    file_path.write_text("test content")
    yield file_path
    # Cleanup happens automatically with tmp_path


def test_with_temp_file(temp_file):
    """Test używający tymczasowego pliku"""
    content = temp_file.read_text()
    assert content == "test content"


# ============================================================================
# MOCKING - Symulowanie obiektów (zaawansowane)
# ============================================================================

class MockSimulation:
    """Przykładowa klasa do mockowania"""
    def __init__(self):
        self.steps = 0
    
    def step(self):
        self.steps += 1
        return self.steps


@pytest.fixture
def mock_sim():
    """Fixture zwracający mock symulacji"""
    return MockSimulation()


def test_with_mock(mock_sim):
    """Test używający mocka"""
    assert mock_sim.steps == 0
    mock_sim.step()
    assert mock_sim.steps == 1


# ============================================================================
# INTEGRATION TEST EXAMPLE
# ============================================================================

@pytest.mark.integration
@pytest.mark.slow
def test_full_simulation_pipeline():
    """
    Przykład pełnego testu integracyjnego.
    Ten test jest oznaczony jako slow i integration - uruchomi się tylko
    w dedykowanym jobie CI po merge do main.
    """
    # 1. Setup
    config = {
        "n_particles": 10,
        "steps": 100,
        "dt": 0.01
    }
    
    # 2. Run simulation (tu byłaby faktyczna symulacja)
    results = []
    for step in range(config["steps"]):
        results.append(step)
    
    # 3. Verify
    assert len(results) == config["steps"]
    assert results[0] == 0
    assert results[-1] == config["steps"] - 1


# ============================================================================
# PROPERTY-BASED TESTING (z hypothesis)
# ============================================================================

from hypothesis import given, strategies as st


@given(st.integers(min_value=0, max_value=1000))
def test_double_is_even(n):
    """Test property-based: podwojona liczba jest parzysta"""
    assert (n * 2) % 2 == 0


@given(st.floats(min_value=-100, max_value=100, allow_nan=False))
def test_absolute_value_positive(x):
    """Test property-based: wartość bezwzględna jest nieujemna"""
    assert abs(x) >= 0


# ============================================================================
# HELPFUL PYTEST COMMANDS
# ============================================================================

"""
Przydatne komendy pytest:

# Uruchom wszystkie testy
pytest

# Uruchom ten plik
pytest tests/test_example_template.py

# Uruchom konkretny test
pytest tests/test_example_template.py::test_basic_assertion

# Uruchom z verbose output
pytest -v

# Uruchom tylko testy unit
pytest -m unit

# Uruchom wszystko oprócz slow
pytest -m "not slow"

# Uruchom testy zawierające "fixture" w nazwie
pytest -k fixture

# Pokaż print statements
pytest -s

# Zatrzymaj przy pierwszym failurze
pytest -x

# Pokaż local variables przy failure
pytest -l

# Pokaż pełny traceback
pytest --tb=long

# Uruchom równolegle (wymaga pytest-xdist)
pytest -n auto

# Generuj coverage report
pytest --cov=backend/sim --cov-report=html
"""

