"""
Authentication utilities: password hashing and JWT tokens.
"""

import logging
from datetime import datetime, timedelta
from typing import Optional
from jose import JWTError, jwt
import bcrypt
from sqlalchemy.orm import Session
import uuid

from backend.billing.config import JWT_SECRET_KEY, JWT_ALGORITHM, JWT_EXPIRATION_HOURS, BCRYPT_ROUNDS
from backend.billing.models import User
from backend.billing.api_key import generate_api_key

logger = logging.getLogger(__name__)


def hash_password(password: str) -> str:
    """
    Hash password using bcrypt.
    
    Args:
        password: Plain text password
    
    Returns:
        Hashed password
    """
    # Convert password to bytes
    password_bytes = password.encode('utf-8')
    # Generate salt and hash
    salt = bcrypt.gensalt(rounds=BCRYPT_ROUNDS)
    hashed = bcrypt.hashpw(password_bytes, salt)
    # Return as string
    return hashed.decode('utf-8')


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """
    Verify password against hash.
    
    Args:
        plain_password: Plain text password
        hashed_password: Hashed password
    
    Returns:
        True if password matches
    """
    try:
        password_bytes = plain_password.encode('utf-8')
        hashed_bytes = hashed_password.encode('utf-8')
        return bcrypt.checkpw(password_bytes, hashed_bytes)
    except Exception:
        return False


def generate_jwt_token(user_id: uuid.UUID) -> str:
    """
    Generate JWT token for user.
    
    Args:
        user_id: User ID
    
    Returns:
        JWT token string
    """
    expiration = datetime.utcnow() + timedelta(hours=JWT_EXPIRATION_HOURS)
    
    payload = {
        "sub": str(user_id),  # Subject (user ID)
        "exp": expiration,
        "iat": datetime.utcnow()
    }
    
    token = jwt.encode(payload, JWT_SECRET_KEY, algorithm=JWT_ALGORITHM)
    return token


def verify_jwt_token(token: str) -> Optional[uuid.UUID]:
    """
    Verify JWT token and return user ID.
    
    Args:
        token: JWT token string
    
    Returns:
        User ID if valid, None otherwise
    """
    try:
        payload = jwt.decode(token, JWT_SECRET_KEY, algorithms=[JWT_ALGORITHM])
        user_id = payload.get("sub")
        if user_id:
            return uuid.UUID(user_id)
    except JWTError as e:
        logger.warning(f"JWT verification failed: {e}")
        return None
    
    return None


def create_user(
    db: Session,
    email: str,
    password: str,
    tier: str = "hobby"
) -> User:
    """
    Create new user with hashed password and API key.
    
    Args:
        db: Database session
        email: User email
        password: Plain text password
        tier: Tier name (default: hobby)
    
    Returns:
        Created User object
    
    Raises:
        ValueError: If email already exists
    """
    # Check if user exists
    existing_user = db.query(User).filter(User.email == email).first()
    if existing_user:
        raise ValueError(f"Email already registered: {email}")
    
    # Hash password
    password_hash = hash_password(password)
    
    # Generate API key
    api_key = generate_api_key()
    
    # Ensure API key is unique
    while db.query(User).filter(User.api_key == api_key).first():
        api_key = generate_api_key()
    
    # Create user
    user = User(
        email=email,
        password_hash=password_hash,
        api_key=api_key,
        tier=tier,
        subscription_status="trial"  # Start with trial
    )
    
    db.add(user)
    db.commit()
    db.refresh(user)
    
    logger.info(f"Created user: {user.email} (tier: {tier})")
    return user


def authenticate_user(
    db: Session,
    email: str,
    password: str
) -> Optional[User]:
    """
    Authenticate user by email and password.
    
    Args:
        db: Database session
        email: User email
        password: Plain text password
    
    Returns:
        User object if authenticated, None otherwise
    """
    user = db.query(User).filter(User.email == email).first()
    if not user:
        return None
    
    if not verify_password(password, user.password_hash):
        return None
    
    return user

