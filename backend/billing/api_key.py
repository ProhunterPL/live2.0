"""
API key generation for users.
"""

import secrets
from backend.billing.config import API_KEY_PREFIX, API_KEY_LENGTH


def generate_api_key() -> str:
    """
    Generate cryptographically secure API key.
    
    Format: sk_live_{32_hex_characters}
    
    Returns:
        API key string
    """
    random_hex = secrets.token_hex(API_KEY_LENGTH // 2)  # token_hex returns 2x length
    return f"{API_KEY_PREFIX}{random_hex}"

