"""
Load test scenario for Legally subscription endpoints.

Based on MONETIZATION_STATUS.md requirements:
- Concurrent users: 100, 500, 1000, 5000
- Ramp-up time: 60 seconds
- Duration: 5 minutes per level
- Target RPS: 1000 requests/second

Run with:
    locust -f legally_subscription.py --host=http://localhost:8000 --users=500 --spawn-rate=10 --run-time=5m
"""

import os
import random
from locust import HttpUser, task, between

class LegallySubscriptionUser(HttpUser):
    """
    Load test scenario for Legally subscription endpoints.
    """
    wait_time = between(1, 3)  # Wait 1-3 seconds between requests
    host = os.getenv("LOAD_TEST_HOST", "http://localhost:8000")
    
    def on_start(self):
        """Setup: Register and login"""
        self.email = f"test_{random.randint(1000, 9999)}@example.com"
        self.password = "test123"
        
        # Register
        try:
            self.client.post("/api/auth/register", json={
                "email": self.email,
                "password": self.password
            })
        except:
            pass  # May already exist
        
        # Login and get API key
        try:
            response = self.client.post("/api/auth/login", json={
                "email": self.email,
                "password": self.password
            })
            if response.status_code == 200:
                data = response.json()
                self.api_key = data.get("api_key", "test_key")
            else:
                self.api_key = "test_key"
        except:
            self.api_key = "test_key"
        
        self.headers = {"X-API-Key": self.api_key}
    
    @task(3)
    def check_subscription_status(self):
        """Check subscription status (most common)"""
        self.client.get(
            "/api/subscription/status",
            headers=self.headers,
            name="/api/subscription/status"
        )
    
    @task(2)
    def get_usage_stats(self):
        """Get usage statistics"""
        self.client.get(
            "/api/usage/stats",
            headers=self.headers,
            name="/api/usage/stats"
        )
    
    @task(1)
    def create_subscription(self):
        """Create subscription (less frequent)"""
        self.client.post(
            "/api/subscription/create",
            json={"tier": "starter"},
            headers=self.headers,
            name="/api/subscription/create"
        )
    
    @task(1)
    def upgrade_subscription(self):
        """Upgrade subscription"""
        self.client.post(
            "/api/subscription/upgrade",
            json={"tier": "professional"},
            headers=self.headers,
            name="/api/subscription/upgrade"
        )
