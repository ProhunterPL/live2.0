"""
Configuration for billing module.
"""

import os
from pathlib import Path

# Database configuration
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    "postgresql://live2:password@localhost:5432/live2_billing"
)

# Stripe configuration
STRIPE_SECRET_KEY = os.getenv("STRIPE_SECRET_KEY", "")
STRIPE_PUBLISHABLE_KEY = os.getenv("STRIPE_PUBLISHABLE_KEY", "")
STRIPE_WEBHOOK_SECRET = os.getenv("STRIPE_WEBHOOK_SECRET", "")

# Stripe Price IDs (must be created in Stripe dashboard)
STRIPE_PRICE_IDS = {
    "hobby": os.getenv("STRIPE_PRICE_ID_HOBBY", "price_hobby_monthly"),
    "research": os.getenv("STRIPE_PRICE_ID_RESEARCH", "price_research_monthly"),
    "pro": os.getenv("STRIPE_PRICE_ID_PRO", "price_pro_monthly"),
    "enterprise": None  # Custom pricing
}

# JWT configuration
JWT_SECRET_KEY = os.getenv("JWT_SECRET_KEY", "change-this-secret-key-in-production")
JWT_ALGORITHM = "HS256"
JWT_EXPIRATION_HOURS = 24

# Password hashing
BCRYPT_ROUNDS = 12

# API Key configuration
API_KEY_PREFIX = "sk_live_"
API_KEY_LENGTH = 32  # hex characters (16 bytes)

# Redis configuration (reuse from API v1)
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_USERNAME = os.getenv("REDIS_USERNAME", None)  # None for local Redis, "default" for Redis Labs
REDIS_PASSWORD = os.getenv("REDIS_PASSWORD", None)  # None for local Redis, password for Redis Labs
REDIS_DB_USAGE = 0  # Same as API v1
REDIS_DECODE_RESPONSES = os.getenv("REDIS_DECODE_RESPONSES", "True").lower() == "true"

