"""
Load test scenario for Live 2.0 API v1 endpoints.

Based on MONETIZATION_STATUS.md requirements:
- Concurrent users: 50, 200, 500
- Ramp-up time: 120 seconds (simulations are time-consuming)
- Duration: 10 minutes per level
- Target RPS: 200 requests/second

Run with:
    locust -f live2_api.py --host=http://localhost:8000 --users=200 --spawn-rate=5 --run-time=10m
"""

import os
from locust import HttpUser, task, between

class Live2APIUser(HttpUser):
    """
    Load test scenario for Live 2.0 API v1.
    """
    wait_time = between(2, 5)  # Longer wait for simulation endpoints
    host = os.getenv("LOAD_TEST_HOST", "http://localhost:8000")
    
    def on_start(self):
        """Setup: Get API key"""
        # Assume API key from environment or test user
        self.api_key = os.getenv("TEST_API_KEY", "test_key")
        self.headers = {"X-API-Key": self.api_key}
    
    @task(5)
    def get_molecules(self):
        """Get molecules (most common)"""
        self.client.get(
            "/api/v1/molecules?limit=100",
            headers=self.headers,
            name="/api/v1/molecules"
        )
    
    @task(3)
    def get_reactions(self):
        """Get reactions"""
        self.client.get(
            "/api/v1/reactions?limit=50",
            headers=self.headers,
            name="/api/v1/reactions"
        )
    
    @task(2)
    def generate_dataset(self):
        """Generate dataset (async job)"""
        self.client.post(
            "/api/v1/generate_dataset",
            json={
                "dataset_type": "reaction_trajectories",
                "params": {"runs": ["run1"]}
            },
            headers=self.headers,
            name="/api/v1/generate_dataset"
        )
    
    @task(1)
    def run_simulation(self):
        """Run simulation (longest operation)"""
        self.client.post(
            "/api/v1/run_simulation",
            json={
                "initial_conditions": {},
                "steps": 1000
            },
            headers=self.headers,
            name="/api/v1/run_simulation"
        )
    
    @task(1)
    def get_job_status(self):
        """Get job status"""
        self.client.get(
            "/api/v1/jobs/test_job_id",
            headers=self.headers,
            name="/api/v1/jobs/{job_id}"
        )
