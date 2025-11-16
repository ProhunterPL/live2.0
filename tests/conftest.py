"""
Pytest configuration and fixtures for Live 2.0 root tests
Handles Taichi initialization and test setup
"""

import pytest
import taichi as ti
import os
import sys

# Add project root to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

def pytest_configure(config):
    """Configure pytest with Taichi initialization"""
    # Initialize Taichi for tests - use CPU for consistency in CI
    try:
        ti.init(arch=ti.cpu, debug=False, cpu_max_num_threads=1)
        print("Taichi initialized for tests (CPU)")
    except Exception as e:
        print(f"Failed to initialize Taichi: {e}")
        # Try fallback initialization
        ti.init(arch=ti.cpu, debug=False)
        print("Taichi initialized with fallback settings")

def pytest_sessionstart(session):
    """Called after the Session object has been created"""
    # Ensure Taichi is initialized before any tests run
    pytest_configure(session.config)

@pytest.fixture(scope="session", autouse=True)
def taichi_session():
    """Session-scoped fixture to ensure Taichi is initialized"""
    try:
        ti.init(arch=ti.cpu, debug=False, cpu_max_num_threads=1)
    except:
        ti.init(arch=ti.cpu, debug=False)
    yield
    # Cleanup after all tests
    try:
        ti.reset()
    except:
        pass

@pytest.fixture(autouse=True)
def taichi_test():
    """Test-scoped fixture for Taichi cleanup between tests"""
    yield
    # Note: NOT resetting Taichi between tests to avoid threading issues
    # Taichi doesn't handle resets well in multi-threaded environments (FastAPI)
    # Only reset at session end via taichi_session fixture
    pass

