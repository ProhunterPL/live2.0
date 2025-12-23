#!/usr/bin/env python3
"""
Simple Redis connection test for billing module.
"""

import sys
import os
from dotenv import load_dotenv

load_dotenv()

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def test_redis_direct():
    """Test Redis connection directly."""
    print("=" * 60)
    print("Redis Connection Test (Direct)")
    print("=" * 60)
    
    try:
        import redis
        
        # Get config from env
        redis_host = os.getenv("REDIS_HOST", "localhost")
        redis_port = int(os.getenv("REDIS_PORT", 6379))
        redis_username = os.getenv("REDIS_USERNAME", None)
        redis_password = os.getenv("REDIS_PASSWORD", None)
        redis_db = 0  # Usage DB
        decode_responses = os.getenv("REDIS_DECODE_RESPONSES", "True").lower() == "true"
        
        print(f"\nConfiguration:")
        print(f"  Host: {redis_host}")
        print(f"  Port: {redis_port}")
        print(f"  Username: {redis_username or 'None (local)'}")
        print(f"  Password: {'***' if redis_password else 'None (local)'}")
        print(f"  DB: {redis_db}")
        
        # Create client
        redis_kwargs = {
            "host": redis_host,
            "port": redis_port,
            "db": redis_db,
            "decode_responses": decode_responses
        }
        if redis_username:
            redis_kwargs["username"] = redis_username
        if redis_password:
            redis_kwargs["password"] = redis_password
        
        print(f"\n[1/4] Creating Redis client...")
        client = redis.Redis(**redis_kwargs)
        
        print(f"[2/4] Testing connection (ping)...")
        result = client.ping()
        if result:
            print("  OK - Redis is responding")
        else:
            print("  FAILED - No response")
            return False
        
        print(f"[3/4] Testing write operation...")
        test_key = "test:billing:connection_check"
        test_value = "test_123"
        client.setex(test_key, 10, test_value)
        print("  OK - Write successful")
        
        print(f"[4/4] Testing read operation...")
        retrieved = client.get(test_key)
        if retrieved == test_value:
            print("  OK - Read successful")
        else:
            print(f"  FAILED - Expected '{test_value}', got '{retrieved}'")
            return False
        
        # Cleanup
        client.delete(test_key)
        
        # Test usage key format
        print(f"\n[5/5] Testing usage tracking format...")
        usage_key = "usage:test-user-123:2025-12:reactions"
        client.incrby(usage_key, 5)
        client.expire(usage_key, 60)
        count = int(client.get(usage_key) or 0)
        if count == 5:
            print("  OK - Usage tracking format works")
        else:
            print(f"  FAILED - Expected 5, got {count}")
            return False
        client.delete(usage_key)
        
        print("\n" + "=" * 60)
        print("SUCCESS: Redis connection and operations working!")
        print("=" * 60)
        return True
        
    except redis.ConnectionError as e:
        print(f"\nFAILED: Redis connection error: {e}")
        print("\nTroubleshooting:")
        print("  1. Check if Redis is running")
        print("  2. Verify REDIS_HOST and REDIS_PORT in .env")
        print("  3. For Redis Labs, check REDIS_USERNAME and REDIS_PASSWORD")
        return False
    except redis.AuthenticationError as e:
        print(f"\nFAILED: Redis authentication error: {e}")
        print("\nTroubleshooting:")
        print("  1. Check REDIS_USERNAME and REDIS_PASSWORD in .env")
        return False
    except ImportError:
        print("\nFAILED: redis module not installed")
        print("Install with: pip install redis")
        return False
    except Exception as e:
        print(f"\nFAILED: Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_redis_direct()
    sys.exit(0 if success else 1)

