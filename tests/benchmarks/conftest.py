"""
Pytest configuration for benchmark tests
=========================================

Provides fixtures and configuration for benchmark reaction tests.
"""

import pytest
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def pytest_configure(config):
    """Register custom markers"""
    config.addinivalue_line(
        "markers", "benchmark: mark test as benchmark reaction test (slow)"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow (> 1 minute)"
    )


@pytest.fixture(scope="session")
def literature_data():
    """
    Fixture providing literature data for all benchmark reactions
    
    Returns dict with formose, strecker, hcn_polymerization data
    """
    from test_formose import FORMOSE_LITERATURE
    from test_strecker import STRECKER_LITERATURE
    from test_hcn_polymerization import HCN_POLYMERIZATION_LITERATURE
    
    return {
        'formose': FORMOSE_LITERATURE,
        'strecker': STRECKER_LITERATURE,
        'hcn_polymerization': HCN_POLYMERIZATION_LITERATURE
    }


@pytest.fixture
def tolerance_config():
    """
    Fixture providing tolerance configuration for yield comparisons
    
    Benchmark reactions have inherent variability, so we allow ±30% tolerance
    """
    return {
        'yield_tolerance': 0.30,     # ±30% for yields
        'rate_tolerance': 0.50,      # ±50% for rate constants
        'detection_tolerance': 2.0   # 2x for detection times
    }


@pytest.fixture
def simulation_config_benchmark():
    """
    Fixture providing simulation configuration for benchmark tests
    
    Optimized for benchmark reactions (longer runs, more particles)
    """
    config = {
        'max_particles': 5000,
        'max_steps': 100000,
        'dt': 0.005,
        'temperature': 298,
        'seed': 42,
        'enable_diagnostics': True,
        'enable_validation': True
    }
    return config


@pytest.fixture
def reaction_analyzer():
    """
    Fixture providing ReactionAnalyzer instance
    
    Will be implemented in Week 3.3
    """
    # Placeholder
    class MockReactionAnalyzer:
        def detect_products(self, trajectory):
            return []
        
        def compute_yields(self, products, initial_reactants):
            return {}
        
        def test_autocatalysis(self, concentration_trajectory):
            return False
    
    return MockReactionAnalyzer()


def pytest_collection_modifyitems(config, items):
    """
    Automatically mark slow tests
    
    Tests with 'slow' in name or marked as 'benchmark' are automatically slow
    """
    for item in items:
        if "slow" in item.nodeid.lower() or "benchmark" in item.keywords:
            item.add_marker(pytest.mark.slow)


@pytest.fixture
def benchmark_report_dir(tmp_path):
    """
    Fixture providing directory for benchmark reports
    
    Each test run gets a unique directory for results
    """
    report_dir = tmp_path / "benchmark_results"
    report_dir.mkdir(exist_ok=True)
    return report_dir

