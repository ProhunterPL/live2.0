#!/usr/bin/env python3
"""
Full monetization system test - checks all components.
"""

import sys
import os
from dotenv import load_dotenv

load_dotenv()

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def test_redis():
    """Test Redis connection."""
    print("[1/4] Testing Redis...")
    try:
        import redis
        redis_host = os.getenv("REDIS_HOST", "localhost")
        redis_port = int(os.getenv("REDIS_PORT", 6379))
        redis_username = os.getenv("REDIS_USERNAME", None)
        redis_password = os.getenv("REDIS_PASSWORD", None)
        
        redis_kwargs = {
            "host": redis_host,
            "port": redis_port,
            "db": 0,
            "decode_responses": True
        }
        if redis_username:
            redis_kwargs["username"] = redis_username
        if redis_password:
            redis_kwargs["password"] = redis_password
        
        client = redis.Redis(**redis_kwargs)
        client.ping()
        print("  OK - Redis connection working")
        return True
    except Exception as e:
        print(f"  FAILED - {e}")
        return False

def test_database():
    """Test database connection."""
    print("[2/4] Testing Database...")
    try:
        from sqlalchemy import create_engine, text
        from sqlalchemy.pool import NullPool
        
        database_url = os.getenv("DATABASE_URL")
        if not database_url:
            print("  FAILED - DATABASE_URL not set")
            return False
        
        engine = create_engine(database_url, poolclass=NullPool, echo=False)
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        print("  OK - Database connection working")
        return True
    except Exception as e:
        print(f"  FAILED - {e}")
        return False

def test_stripe_config():
    """Test Stripe configuration."""
    print("[3/4] Testing Stripe Configuration...")
    stripe_secret = os.getenv("STRIPE_SECRET_KEY", "").strip()
    stripe_webhook = os.getenv("STRIPE_WEBHOOK_SECRET", "").strip()
    price_hobby = os.getenv("STRIPE_PRICE_ID_HOBBY", "").strip()
    price_research = os.getenv("STRIPE_PRICE_ID_RESEARCH", "").strip()
    price_pro = os.getenv("STRIPE_PRICE_ID_PRO", "").strip()
    
    issues = []
    
    if not stripe_secret:
        issues.append("STRIPE_SECRET_KEY not set")
    if not stripe_webhook:
        issues.append("STRIPE_WEBHOOK_SECRET not set")
    elif not stripe_webhook.startswith("whsec_"):
        issues.append(f"STRIPE_WEBHOOK_SECRET format invalid (should start with 'whsec_', got: {stripe_webhook[:20]}...)")
    if not price_hobby:
        issues.append("STRIPE_PRICE_ID_HOBBY not set")
    elif not price_hobby.startswith("price_"):
        # Accept both price_ and prod_ (user might have Product IDs)
        if not price_hobby.startswith("prod_"):
            issues.append(f"STRIPE_PRICE_ID_HOBBY format invalid (should be 'price_...' or 'prod_...', got: {price_hobby[:10]}...)")
    if not price_research:
        issues.append("STRIPE_PRICE_ID_RESEARCH not set")
    elif not price_research.startswith("price_"):
        if not price_research.startswith("prod_"):
            issues.append(f"STRIPE_PRICE_ID_RESEARCH format invalid (should be 'price_...' or 'prod_...', got: {price_research[:10]}...)")
    if not price_pro:
        issues.append("STRIPE_PRICE_ID_PRO not set")
    elif not price_pro.startswith("price_"):
        if not price_pro.startswith("prod_"):
            issues.append(f"STRIPE_PRICE_ID_PRO format invalid (should be 'price_...' or 'prod_...', got: {price_pro[:10]}...)")
    
    if issues:
        for issue in issues:
            print(f"  WARNING - {issue}")
        return False
    
    print("  OK - Stripe configuration complete")
    return True

def test_jwt_config():
    """Test JWT configuration."""
    print("[4/4] Testing JWT Configuration...")
    jwt_secret = os.getenv("JWT_SECRET_KEY")
    if not jwt_secret or jwt_secret == "change-this-secret-key-in-production":
        print("  WARNING - JWT_SECRET_KEY not changed from default")
        return False
    
    print("  OK - JWT secret configured")
    return True

if __name__ == "__main__":
    print("=" * 60)
    print("Monetization System Test")
    print("=" * 60)
    print()
    
    results = []
    results.append(("Redis", test_redis()))
    results.append(("Database", test_database()))
    results.append(("Stripe Config", test_stripe_config()))
    results.append(("JWT Config", test_jwt_config()))
    
    print()
    print("=" * 60)
    print("Results:")
    print("=" * 60)
    
    all_ok = True
    for name, result in results:
        status = "OK" if result else "FAILED/WARNING"
        symbol = "[OK]" if result else "[FAIL]"
        print(f"  {symbol} {name}: {status}")
        if not result:
            all_ok = False
    
    print()
    if all_ok:
        print("SUCCESS: All components ready!")
        sys.exit(0)
    else:
        print("WARNING: Some components need attention (see above)")
        sys.exit(1)

