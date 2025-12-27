"""
Configuration for API v1.
"""

import os
from pathlib import Path
from typing import Optional

# Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_USERNAME = os.getenv("REDIS_USERNAME", None)  # None for local Redis, "default" for Redis Labs
REDIS_PASSWORD = os.getenv("REDIS_PASSWORD", None)  # None for local Redis, password for Redis Labs
REDIS_DECODE_RESPONSES = os.getenv("REDIS_DECODE_RESPONSES", "True").lower() == "true"
REDIS_DB_JOBS = 1  # Separate DB for jobs
REDIS_DB_USAGE = 0  # DB for usage tracking

# Storage configuration
STORAGE_TYPE = os.getenv("STORAGE_TYPE", "local")  # "local" | "s3" | "minio"
STORAGE_BASE_PATH = Path(os.getenv("STORAGE_BASE_PATH", "datasets/api"))
S3_BUCKET = os.getenv("AWS_S3_BUCKET")  # Changed from S3_BUCKET to AWS_S3_BUCKET
S3_REGION = os.getenv("S3_REGION", "us-east-1")
S3_ACCESS_KEY_ID = os.getenv("S3_ACCESS_KEY_ID")
S3_SECRET_ACCESS_KEY = os.getenv("S3_SECRET_ACCESS_KEY")
S3_ENDPOINT_URL = os.getenv("S3_ENDPOINT_URL")  # For S3-compatible services (Supabase, MinIO)

# Results directory
BASE_RESULTS_DIR = os.getenv(
    "BASE_RESULTS_DIR",
    str(Path(__file__).parent.parent.parent.parent / "results" / "phase2b_additional")
)

# Job configuration
JOB_TTL_DAYS = 7  # Job data retention (days)
JOB_QUEUE_NAME = "job_queue"

# Rate limiting configuration
RATE_LIMIT_WINDOW_SECONDS = 60  # Per-minute rate limiting window

