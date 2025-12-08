"""
Integration tests for billing routes.
"""

import pytest
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from backend.api.v1.main import app
from backend.billing.models import User
from backend.billing.auth import create_user
# Use conftest db fixture


@pytest.fixture
def client(db):
    """Create test client."""
    from fastapi import FastAPI
    from backend.billing.routes import auth as billing_auth
    from unittest.mock import patch
    
    # Create test app with billing routes
    test_app = FastAPI()
    
    # Override get_db dependency to use our test db
    def override_get_db():
        yield db
    
    # Override dependency before including router
    from backend.billing.database import get_db
    test_app.dependency_overrides[get_db] = override_get_db
    
    # Mock STRIPE_SECRET_KEY to avoid PaymentProcessor errors
    with patch('backend.billing.routes.auth.STRIPE_SECRET_KEY', ''):
        test_app.include_router(billing_auth.router)
        yield TestClient(test_app)


# Use conftest db fixture


def test_register_user(client, db):
    """Test user registration."""
    import uuid
    unique_email = f"test_register_{uuid.uuid4().hex[:8]}@example.com"
    response = client.post(
        "/auth/register",  # Router has prefix="/auth"
        json={
            "email": unique_email,
            "password": "test_password",
            "tier": "research"
        }
    )
    
    if response.status_code != 200:
        print(f"Response status: {response.status_code}")
        print(f"Response body: {response.text}")
    assert response.status_code == 200, f"Expected 200, got {response.status_code}: {response.text}"
    data = response.json()
    assert "access_token" in data
    assert "api_key" in data
    assert data["user"]["email"] == unique_email
    assert data["user"]["tier"] == "research"


def test_register_duplicate_email(client, db):
    """Test registration with duplicate email."""
    import uuid
    unique_email = f"test_duplicate_{uuid.uuid4().hex[:8]}@example.com"
    
    # Register first user
    client.post(
        "/auth/register",  # Router has prefix="/auth"
        json={
            "email": unique_email,
            "password": "password1",
            "tier": "hobby"
        }
    )
    
    # Try to register again with same email
    response = client.post(
        "/auth/register",  # Router has prefix="/auth"
        json={
            "email": unique_email,
            "password": "password2",
            "tier": "hobby"
        }
    )
    
    assert response.status_code == 400
    assert "already registered" in response.json()["detail"].lower()


def test_login(client, db):
    """Test user login."""
    import uuid
    unique_email = f"test_login_{uuid.uuid4().hex[:8]}@example.com"
    
    # Register user first
    register_response = client.post(
        "/auth/register",  # Router has prefix="/auth"
        json={
            "email": unique_email,
            "password": "test_password",
            "tier": "hobby"
        }
    )
    
    # Login
    response = client.post(
        "/auth/login",  # Router has prefix="/auth"
        json={
            "email": unique_email,
            "password": "test_password"
        }
    )
    
    assert response.status_code == 200
    data = response.json()
    assert "access_token" in data
    assert "api_key" in data
    assert data["user"]["email"] == unique_email


def test_login_invalid_credentials(client, db):
    """Test login with invalid credentials."""
    response = client.post(
        "/auth/login",  # Router has prefix="/auth"
        json={
            "email": "nonexistent@example.com",
            "password": "wrong_password"
        }
    )
    
    assert response.status_code == 401
    assert "invalid" in response.json()["detail"].lower()

