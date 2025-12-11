"""
Main Locust configuration file.

Run with:
    locust -f locustfile.py --host=http://localhost:8000

Or use specific scenario:
    locust -f legally_subscription.py --host=http://localhost:8000
    locust -f live2_api.py --host=http://localhost:8000
"""

import os
from locust import HttpUser, task, between

class APIUser(HttpUser):
    """
    Base load test user.
    """
    wait_time = between(1, 3)
    host = os.getenv("LOAD_TEST_HOST", "http://localhost:8000")
    
    def on_start(self):
        """Setup: Get API key"""
        self.api_key = os.getenv("TEST_API_KEY", "test_key")
        self.headers = {"X-API-Key": self.api_key}
    
    @task
    def health_check(self):
        """Health check endpoint"""
        self.client.get("/health", headers=self.headers)
