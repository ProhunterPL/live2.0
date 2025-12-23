#!/usr/bin/env python3
"""
Quick test script to verify Redis connection for billing module.

Usage:
    python scripts/test_redis_billing.py
"""

import sys
import os

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from dotenv import load_dotenv
load_dotenv()

def test_redis_connection():
    """Test Redis connection using billing module dependencies."""
    print("=" * 60)
    print("Testing Redis Connection for Billing Module")
    print("=" * 60)
    
    try:
        # Import dependencies directly to avoid circular imports
        import importlib.util
        deps_path = os.path.join(project_root, 'backend', 'api', 'v1', 'dependencies.py')
        deps_spec = importlib.util.spec_from_file_location("api_deps", deps_path)
        deps_module = importlib.util.module_from_spec(deps_spec)
        deps_spec.loader.exec_module(deps_module)
        get_redis_usage = deps_module.get_redis_usage
        print("\n[1/3] Importing Redis dependency... OK")
        
        redis_client = get_redis_usage()
        print("[2/3] Getting Redis client... OK")
        
        # Test ping
        try:
            redis_client.ping()
            print("[3/3] Redis ping()... OK")
            print("\n✅ Redis connection successful!")
            
            # Test basic operations
            print("\n" + "=" * 60)
            print("Testing Redis Operations")
            print("=" * 60)
            
            test_key = "test:billing:connection"
            test_value = "test_value_123"
            
            # Set
            redis_client.setex(test_key, 60, test_value)
            print(f"[SET] {test_key} = {test_value}... OK")
            
            # Get
            retrieved = redis_client.get(test_key)
            if retrieved == test_value:
                print(f"[GET] {test_key} = {retrieved}... OK")
            else:
                print(f"[GET] Expected {test_value}, got {retrieved}... FAILED")
                return False
            
            # Delete
            redis_client.delete(test_key)
            print(f"[DEL] {test_key}... OK")
            
            # Test usage tracking key format
            usage_key = "usage:test-user-id:2025-12:reactions"
            redis_client.incrby(usage_key, 1)
            redis_client.expire(usage_key, 3600)
            count = redis_client.get(usage_key)
            print(f"[USAGE] Increment usage key... OK (count: {count})")
            redis_client.delete(usage_key)
            
            print("\n✅ All Redis operations successful!")
            return True
            
        except Exception as e:
            if "Redis not available" in str(e) or "not available" in str(e).lower():
                print(f"[3/3] Redis ping()... FAILED (Redis not available)")
                print(f"\n❌ Redis connection failed: {e}")
                print("\n💡 Troubleshooting:")
                print("   1. Check if Redis is running: redis-cli ping")
                print("   2. Check REDIS_HOST and REDIS_PORT in .env")
                print("   3. For Redis Labs, check REDIS_USERNAME and REDIS_PASSWORD")
            else:
                print(f"[3/3] Redis ping()... FAILED")
                print(f"\n❌ Redis connection failed: {e}")
            return False
            
    except ImportError as e:
        print(f"\n❌ Import failed: {e}")
        print("💡 Make sure you're in the project root and dependencies are installed")
        return False
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_billing_config():
    """Test billing configuration."""
    print("\n" + "=" * 60)
    print("Billing Configuration Check")
    print("=" * 60)
    
    # Import config directly to avoid circular imports
    import importlib.util
    config_path = os.path.join(project_root, 'backend', 'billing', 'config.py')
    config_spec = importlib.util.spec_from_file_location("billing_config", config_path)
    config_module = importlib.util.module_from_spec(config_spec)
    config_spec.loader.exec_module(config_module)
    
    REDIS_HOST = config_module.REDIS_HOST
    REDIS_PORT = config_module.REDIS_PORT
    REDIS_USERNAME = config_module.REDIS_USERNAME
    REDIS_PASSWORD = config_module.REDIS_PASSWORD
    DATABASE_URL = config_module.DATABASE_URL
    STRIPE_SECRET_KEY = config_module.STRIPE_SECRET_KEY
    JWT_SECRET_KEY = config_module.JWT_SECRET_KEY
    
    print(f"REDIS_HOST: {REDIS_HOST}")
    print(f"REDIS_PORT: {REDIS_PORT}")
    print(f"REDIS_USERNAME: {'***' if REDIS_USERNAME else 'None (local Redis)'}")
    print(f"REDIS_PASSWORD: {'***' if REDIS_PASSWORD else 'None (local Redis)'}")
    print(f"DATABASE_URL: {DATABASE_URL[:30]}..." if len(DATABASE_URL) > 30 else f"DATABASE_URL: {DATABASE_URL}")
    print(f"STRIPE_SECRET_KEY: {'Set' if STRIPE_SECRET_KEY else 'Not set (optional for testing)'}")
    print(f"JWT_SECRET_KEY: {'Set' if JWT_SECRET_KEY != 'change-this-secret-key-in-production' else 'Using default (change in production!)'}")


if __name__ == "__main__":
    test_billing_config()
    success = test_redis_connection()
    
    print("\n" + "=" * 60)
    if success:
        print("✅ All tests passed!")
        sys.exit(0)
    else:
        print("❌ Some tests failed. Check the output above.")
        sys.exit(1)

