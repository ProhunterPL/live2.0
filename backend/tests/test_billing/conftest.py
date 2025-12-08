"""
Pytest configuration for billing tests.
"""

import pytest
import os
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

# Use SQLite file-based for tests (to avoid threading issues with FastAPI TestClient)
import tempfile
import os

# Create temporary database file
_test_db_file = tempfile.NamedTemporaryFile(delete=False, suffix='.db')
_test_db_file.close()
TEST_DATABASE_URL = f"sqlite:///{_test_db_file.name}"

@pytest.fixture(scope="session")
def test_engine():
    """Create test database engine."""
    # Use file-based SQLite to avoid threading issues
    engine = create_engine(
        TEST_DATABASE_URL, 
        echo=False,
        connect_args={"check_same_thread": False}
    )
    
    # Create tables
    from backend.billing.models import Base
    Base.metadata.create_all(bind=engine)
    
    yield engine
    
    # Cleanup
    Base.metadata.drop_all(bind=engine)
    # Remove temporary database file
    try:
        os.unlink(_test_db_file.name)
    except:
        pass


@pytest.fixture
def db(test_engine):
    """Create database session for testing."""
    # Ensure tables exist
    from backend.billing.models import Base
    Base.metadata.create_all(bind=test_engine)
    
    SessionLocal = sessionmaker(bind=test_engine)
    db = SessionLocal()
    try:
        yield db
    finally:
        db.rollback()
        db.close()

