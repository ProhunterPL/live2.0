"""
Test script for AWS Batch integration.

Tests job submission, status checking, and cancellation.
"""

import os
import sys
import asyncio
import json
from pathlib import Path

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    project_root = Path(__file__).parent.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
        print(f"Loaded environment variables from {env_file}")
    else:
        print(f"Warning: .env file not found at {env_file}")
except ImportError:
    print("Warning: python-dotenv not installed. Using system environment variables only.")

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Check if boto3 is available before importing
try:
    import boto3
    BOTO3_AVAILABLE = True
except ImportError:
    BOTO3_AVAILABLE = False

# Import directly to avoid circular imports (only if boto3 available)
if BOTO3_AVAILABLE:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "aws_batch",
        project_root / "backend" / "api" / "v1" / "aws_batch.py"
    )
    aws_batch_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(aws_batch_module)
    AWSBatchClient = aws_batch_module.AWSBatchClient


async def test_batch_integration():
    """Test AWS Batch integration."""
    print("=" * 60)
    print("Testing AWS Batch Integration")
    print("=" * 60)
    print()
    
    # Initialize client
    print("1. Initializing AWS Batch client...")
    try:
        batch_client = AWSBatchClient()
        print(f"   ✓ Client initialized: queue={batch_client.job_queue}, definition={batch_client.job_definition}")
    except Exception as e:
        print(f"   ✗ Failed to initialize client: {e}")
        return False
    
    print()
    
    # Test job submission
    print("2. Testing job submission...")
    test_job_id = f"test-{os.urandom(8).hex()}"
    test_params = {
        "simulation_config": {
            "mode": "open_chemistry",
            "config": {
                "num_particles": 100,
                "temperature": 300,
                "steps": 1000
            }
        }
    }
    
    try:
        batch_job_id = batch_client.submit_job(
            job_id=test_job_id,
            job_params=test_params,
            user_id="test-user",
            environment_vars={
                "TEST_MODE": "true"
            }
        )
        print(f"   ✓ Job submitted: {batch_job_id}")
    except Exception as e:
        print(f"   ✗ Failed to submit job: {e}")
        return False
    
    print()
    
    # Test status check
    print("3. Testing status check...")
    try:
        status = batch_client.get_job_status(batch_job_id)
        print(f"   ✓ Status retrieved: {status['status']} ({status['batch_status']})")
        print(f"     Started at: {status.get('started_at', 'N/A')}")
    except Exception as e:
        print(f"   ✗ Failed to get status: {e}")
        return False
    
    print()
    
    # Test cancellation
    print("4. Testing job cancellation...")
    try:
        cancelled = batch_client.cancel_job(batch_job_id, "Test cancellation")
        if cancelled:
            print(f"   ✓ Job cancelled successfully")
        else:
            print(f"   ✗ Cancellation returned False")
            return False
    except Exception as e:
        print(f"   ✗ Failed to cancel job: {e}")
        return False
    
    print()
    
    # Test status after cancellation
    print("5. Checking status after cancellation...")
    try:
        status = batch_client.get_job_status(batch_job_id)
        print(f"   ✓ Final status: {status['status']} ({status['batch_status']})")
    except Exception as e:
        print(f"   ✗ Failed to get final status: {e}")
        return False
    
    print()
    print("=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    return True


if __name__ == "__main__":
    # Check environment variables
    required_vars = ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_REGION"]
    missing = [var for var in required_vars if not os.getenv(var)]
    
    if missing:
        print(f"Error: Missing required environment variables: {', '.join(missing)}")
        print("Please set them in your .env file or environment")
        sys.exit(1)
    
    # Check if boto3 is available
    if not BOTO3_AVAILABLE:
        print("=" * 60)
        print("Warning: boto3 not installed locally")
        print("=" * 60)
        print("Environment variables loaded successfully from .env:")
        print(f"  AWS_ACCESS_KEY_ID: {'OK' if os.getenv('AWS_ACCESS_KEY_ID') else 'MISSING'}")
        print(f"  AWS_SECRET_ACCESS_KEY: {'OK' if os.getenv('AWS_SECRET_ACCESS_KEY') else 'MISSING'}")
        print(f"  AWS_REGION: {os.getenv('AWS_REGION', 'NOT SET')}")
        print()
        print("This test will work on DO Droplet where boto3 is installed.")
        print("To test locally, install: pip install boto3")
        sys.exit(0)
    
    success = asyncio.run(test_batch_integration())
    sys.exit(0 if success else 1)

