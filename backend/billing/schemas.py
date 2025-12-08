"""
Pydantic schemas for billing module.
"""

from pydantic import BaseModel, EmailStr
from typing import Optional, Dict, List
from datetime import date
from uuid import UUID


# Request schemas
class RegisterRequest(BaseModel):
    """User registration request."""
    email: EmailStr
    password: str
    tier: str = "hobby"  # Default tier


class LoginRequest(BaseModel):
    """User login request."""
    email: EmailStr
    password: str


class UpgradeRequest(BaseModel):
    """Subscription upgrade request."""
    tier: str


class CancelRequest(BaseModel):
    """Subscription cancellation request."""
    cancel_at_period_end: bool = True


# Response schemas
class UserResponse(BaseModel):
    """User response."""
    id: UUID
    email: str
    tier: str
    subscription_status: str
    api_key: str
    
    class Config:
        from_attributes = True


class AuthResponse(BaseModel):
    """Authentication response."""
    access_token: str
    api_key: str
    user: UserResponse


class SubscriptionResponse(BaseModel):
    """Subscription response."""
    tier: str
    status: str
    current_period_start: date
    current_period_end: date
    usage: Dict
    
    class Config:
        from_attributes = True


class UsageMonthResponse(BaseModel):
    """Usage for a single month."""
    month: str  # YYYY-MM
    reactions: int
    reactions_quota: Optional[int]
    api_calls: int
    api_calls_quota: Optional[int]
    percentage_used: float


class UsageResponse(BaseModel):
    """Usage response."""
    current_month: UsageMonthResponse
    history: List[UsageMonthResponse]

