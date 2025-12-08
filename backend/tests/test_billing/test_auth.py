"""
Tests for billing authentication.
"""

import pytest
from sqlalchemy.orm import Session
from uuid import UUID

from backend.billing.auth import (
    hash_password,
    verify_password,
    generate_jwt_token,
    verify_jwt_token,
    create_user,
    authenticate_user
)
from backend.billing.models import User


def test_hash_password():
    """Test password hashing."""
    password = "test_pass"  # Shorter password to avoid bcrypt 72-byte limit
    hashed = hash_password(password)
    
    assert hashed != password
    assert len(hashed) > 50  # bcrypt hash is long
    # bcrypt format can be $2b$ or $2a$
    assert hashed.startswith("$2")  # bcrypt format


def test_verify_password():
    """Test password verification."""
    password = "test_pass"  # Shorter password
    hashed = hash_password(password)
    
    assert verify_password(password, hashed) is True
    assert verify_password("wrong_pass", hashed) is False


def test_generate_jwt_token():
    """Test JWT token generation."""
    user_id = UUID("12345678-1234-5678-1234-567812345678")
    token = generate_jwt_token(user_id)
    
    assert isinstance(token, str)
    assert len(token) > 50  # JWT tokens are long


def test_verify_jwt_token():
    """Test JWT token verification."""
    user_id = UUID("12345678-1234-5678-1234-567812345678")
    token = generate_jwt_token(user_id)
    
    verified_id = verify_jwt_token(token)
    assert verified_id == user_id
    
    # Test invalid token
    invalid_id = verify_jwt_token("invalid_token")
    assert invalid_id is None


def test_create_user(db: Session):
    """Test user creation."""
    user = create_user(
        db=db,
        email="test_user@example.com",
        password="test_pass",  # Shorter password
        tier="research"
    )
    
    assert user.id is not None
    assert user.email == "test_user@example.com"
    assert user.tier == "research"
    assert user.subscription_status == "trial"
    assert user.api_key.startswith("sk_live_")
    assert len(user.api_key) > 20
    
    # Verify password is hashed
    assert user.password_hash != "test_pass"
    assert verify_password("test_pass", user.password_hash) is True


def test_create_user_duplicate_email(db: Session):
    """Test user creation with duplicate email."""
    create_user(
        db=db,
        email="test_duplicate@example.com",
        password="pass1",
        tier="hobby"
    )
    
    # Try to create again with same email
    with pytest.raises(ValueError, match="Email already registered"):
        create_user(
            db=db,
            email="test_duplicate@example.com",
            password="pass2",
            tier="hobby"
        )


def test_authenticate_user(db: Session):
    """Test user authentication."""
    user = create_user(
        db=db,
        email="test_auth@example.com",
        password="correct_pass",
        tier="hobby"
    )
    
    # Correct credentials
    authenticated = authenticate_user(db, "test_auth@example.com", "correct_pass")
    assert authenticated is not None
    assert authenticated.id == user.id
    
    # Wrong password
    authenticated = authenticate_user(db, "test_auth@example.com", "wrong_pass")
    assert authenticated is None
    
    # Non-existent user
    authenticated = authenticate_user(db, "nonexistent@example.com", "password")
    assert authenticated is None

