"""
SQLAlchemy models for billing module.
"""

from sqlalchemy import Column, String, Integer, Boolean, DateTime, Date, ForeignKey, UniqueConstraint
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from datetime import datetime, date
import uuid

Base = declarative_base()


class User(Base):
    """User model."""
    __tablename__ = "users"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    email = Column(String(255), unique=True, nullable=False, index=True)
    password_hash = Column(String(255), nullable=False)
    api_key = Column(String(64), unique=True, nullable=False, index=True)
    
    tier = Column(String(20), nullable=False, default="hobby")  # hobby, research, pro, enterprise
    subscription_status = Column(String(20), nullable=False, default="trial")  # active, cancelled, expired, trial
    
    stripe_customer_id = Column(String(255), nullable=True, index=True)
    stripe_subscription_id = Column(String(255), nullable=True, index=True)
    
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated_at = Column(DateTime, nullable=False, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    subscriptions = relationship("Subscription", back_populates="user", cascade="all, delete-orphan")
    usage_records = relationship("Usage", back_populates="user", cascade="all, delete-orphan")
    
    def __repr__(self):
        return f"<User(id={self.id}, email={self.email}, tier={self.tier})>"


class Subscription(Base):
    """Subscription model."""
    __tablename__ = "subscriptions"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False, index=True)
    
    tier = Column(String(20), nullable=False)  # hobby, research, pro, enterprise
    status = Column(String(20), nullable=False)  # active, cancelled, expired, trialing
    
    current_period_start = Column(Date, nullable=False)
    current_period_end = Column(Date, nullable=False)
    cancel_at_period_end = Column(Boolean, nullable=False, default=False)
    
    stripe_subscription_id = Column(String(255), unique=True, nullable=True, index=True)
    
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated_at = Column(DateTime, nullable=False, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    user = relationship("User", back_populates="subscriptions")
    
    def __repr__(self):
        return f"<Subscription(id={self.id}, user_id={self.user_id}, tier={self.tier}, status={self.status})>"


class Usage(Base):
    """Monthly usage aggregates."""
    __tablename__ = "usage"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False, index=True)
    
    date = Column(Date, nullable=False, index=True)  # YYYY-MM-01 (first day of month)
    reactions_count = Column(Integer, nullable=False, default=0)
    api_calls_count = Column(Integer, nullable=False, default=0)
    tier = Column(String(20), nullable=False)  # Tier at time of usage
    
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated_at = Column(DateTime, nullable=False, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Unique constraint: one record per user per month
    __table_args__ = (
        UniqueConstraint('user_id', 'date', name='uq_user_date'),
    )
    
    # Relationships
    user = relationship("User", back_populates="usage_records")
    
    def __repr__(self):
        return f"<Usage(id={self.id}, user_id={self.user_id}, date={self.date}, reactions={self.reactions_count})>"


# Helper function for API v1 integration
def get_user_by_api_key(api_key: str, db_session=None):
    """
    Get user by API key (for API v1 integration).
    
    Args:
        api_key: API key string
        db_session: Optional database session (if None, creates new)
    
    Returns:
        User object or None
    """
    if db_session is None:
        from backend.billing.database import SessionLocal
        db = SessionLocal()
        try:
            user = db.query(User).filter(User.api_key == api_key).first()
            return user
        finally:
            db.close()
    else:
        return db_session.query(User).filter(User.api_key == api_key).first()

