"""
Authentication routes for billing module.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
import logging

from backend.billing.schemas import RegisterRequest, LoginRequest, AuthResponse, UserResponse
from backend.billing.auth import create_user, authenticate_user, generate_jwt_token
from backend.billing.database import get_db
from backend.billing.subscriptions import SubscriptionManager
from backend.billing.payments import PaymentProcessor
from backend.billing.config import STRIPE_SECRET_KEY

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/auth", tags=["auth"])


@router.post("/register", response_model=AuthResponse)
async def register(
    request: RegisterRequest,
    db: Session = Depends(get_db)
):
    """
    Register new user.
    
    Creates user account, generates API key, and creates trial subscription.
    """
    try:
        # Create user
        user = create_user(
            db=db,
            email=request.email,
            password=request.password,
            tier=request.tier
        )
        
        # Create trial subscription (without Stripe payment)
        # Only create subscription if Stripe is configured
        if STRIPE_SECRET_KEY:
            try:
                payment_processor = PaymentProcessor(STRIPE_SECRET_KEY)
                subscription_manager = SubscriptionManager(db, payment_processor)
                subscription_manager.create_subscription(
                    user_id=user.id,
                    tier=request.tier,
                    payment_method_id=None  # Trial without payment
                )
            except Exception as e:
                # If subscription creation fails, continue anyway (trial mode)
                logger.warning(f"Failed to create subscription: {e}")
        # If no Stripe key, user starts with trial status (already set in create_user)
        
        # Generate JWT token
        access_token = generate_jwt_token(user.id)
        
        return AuthResponse(
            access_token=access_token,
            api_key=user.api_key,
            user=UserResponse(
                id=user.id,
                email=user.email,
                tier=user.tier,
                subscription_status=user.subscription_status,
                api_key=user.api_key
            )
        )
    
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to register user: {str(e)}"
        )


@router.post("/login", response_model=AuthResponse)
async def login(
    request: LoginRequest,
    db: Session = Depends(get_db)
):
    """
    Login user.
    
    Returns JWT token and API key.
    """
    user = authenticate_user(db, request.email, request.password)
    
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid email or password"
        )
    
    # Generate JWT token
    access_token = generate_jwt_token(user.id)
    
    return AuthResponse(
        access_token=access_token,
        api_key=user.api_key,
        user=UserResponse(
            id=user.id,
            email=user.email,
            tier=user.tier,
            subscription_status=user.subscription_status,
            api_key=user.api_key
        )
    )

