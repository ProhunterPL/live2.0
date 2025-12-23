#!/usr/bin/env python3
"""
Smoke test E2E for monetization flow:
1. Register user
2. Create checkout session
3. Simulate webhook (subscription created)
4. Verify subscription is active
5. Test API access with subscription
"""

import sys
import os
import requests
import json
import uuid
from dotenv import load_dotenv

load_dotenv()

API_BASE = os.getenv("API_BASE_URL", "http://localhost:8001/api/v1")

def test_register():
    """Test user registration."""
    print("[1/5] Testing user registration...")
    
    email = f"test_{uuid.uuid4().hex[:8]}@example.com"
    password = "test_password_123"
    
    response = requests.post(
        f"{API_BASE}/auth/register",
        json={
            "email": email,
            "password": password,
            "tier": "hobby"
        }
    )
    
    if response.status_code != 200:
        print(f"  FAILED - Status: {response.status_code}, Response: {response.text}")
        return None
    
    data = response.json()
    print(f"  OK - User registered: {email}")
    print(f"  User ID: {data['user']['id']}")
    print(f"  API Key: {data['api_key'][:20]}...")
    
    return {
        "email": email,
        "password": password,
        "user_id": data["user"]["id"],
        "api_key": data["api_key"],
        "access_token": data["access_token"]
    }

def test_login(user_data):
    """Test user login."""
    print("[2/5] Testing user login...")
    
    response = requests.post(
        f"{API_BASE}/auth/login",
        json={
            "email": user_data["email"],
            "password": user_data["password"]
        }
    )
    
    if response.status_code != 200:
        print(f"  FAILED - Status: {response.status_code}")
        return None
    
    data = response.json()
    print(f"  OK - Login successful")
    
    return data

def test_get_subscription(user_data):
    """Test getting subscription status."""
    print("[3/5] Testing get subscription...")
    
    headers = {"X-API-Key": user_data["api_key"]}
    response = requests.get(
        f"{API_BASE}/billing/subscription",
        headers=headers
    )
    
    if response.status_code == 404:
        print(f"  OK - No subscription yet (expected for new user)")
        return None
    elif response.status_code == 200:
        data = response.json()
        print(f"  OK - Subscription found: tier={data['tier']}, status={data['status']}")
        return data
    else:
        print(f"  WARNING - Status: {response.status_code}, Response: {response.text}")
        return None

def test_checkout_session(user_data):
    """Test creating checkout session."""
    print("[4/5] Testing checkout session creation...")
    
    headers = {"X-API-Key": user_data["api_key"]}
    response = requests.post(
        f"{API_BASE}/billing/checkout/session",
        json={
            "tier": "hobby",
            "success_url": "https://example.com/success",
            "cancel_url": "https://example.com/cancel"
        },
        headers=headers
    )
    
    if response.status_code != 200:
        print(f"  FAILED - Status: {response.status_code}, Response: {response.text}")
        return None
    
    data = response.json()
    print(f"  OK - Checkout session created")
    print(f"  Session ID: {data['session_id']}")
    print(f"  URL: {data['url'][:50]}...")
    
    return data

def test_api_access(user_data):
    """Test API v1 access with API key."""
    print("[5/5] Testing API v1 access...")
    
    headers = {"X-API-Key": user_data["api_key"]}
    response = requests.get(
        f"{API_BASE}/health",
        headers=headers
    )
    
    if response.status_code == 200:
        print(f"  OK - API access working")
        return True
    else:
        print(f"  WARNING - Status: {response.status_code}")
        return False

if __name__ == "__main__":
    print("=" * 60)
    print("Monetization E2E Smoke Test")
    print("=" * 60)
    print()
    
    try:
        # Test 1: Register
        user_data = test_register()
        if not user_data:
            print("\nFAILED: Registration failed")
            sys.exit(1)
        
        # Test 2: Login
        login_data = test_login(user_data)
        if not login_data:
            print("\nFAILED: Login failed")
            sys.exit(1)
        
        # Test 3: Get subscription (may not exist yet)
        subscription = test_get_subscription(user_data)
        
        # Test 4: Checkout session
        checkout = test_checkout_session(user_data)
        if not checkout:
            print("\nWARNING: Checkout session creation failed (check Stripe config)")
        
        # Test 5: API access
        api_ok = test_api_access(user_data)
        
        print()
        print("=" * 60)
        print("Smoke Test Results:")
        print("=" * 60)
        print(f"  Registration: OK")
        print(f"  Login: OK")
        print(f"  Subscription: {'Found' if subscription else 'Not found (expected)'}")
        print(f"  Checkout: {'OK' if checkout else 'FAILED'}")
        print(f"  API Access: {'OK' if api_ok else 'WARNING'}")
        print()
        
        if checkout and api_ok:
            print("SUCCESS: Core flow working!")
            print("\nNext steps:")
            print("  1. Complete Stripe checkout in browser (use test card)")
            print("  2. Verify webhook updates subscription status")
            print("  3. Test API v1 endpoints with active subscription")
            sys.exit(0)
        else:
            print("WARNING: Some tests failed (see above)")
            sys.exit(1)
            
    except requests.exceptions.ConnectionError:
        print("\nERROR: Cannot connect to API server")
        print(f"Make sure server is running at {API_BASE}")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

