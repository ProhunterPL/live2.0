"""
AWS Batch integration for job execution.

Handles submission, status checking, and cancellation of AWS Batch jobs.
"""

import os
import logging
from typing import Dict, Optional, List
import boto3
from botocore.exceptions import ClientError

logger = logging.getLogger(__name__)


class AWSBatchClient:
    """Client for AWS Batch operations."""
    
    def __init__(
        self,
        region: Optional[str] = None,
        job_queue: Optional[str] = None,
        job_definition: Optional[str] = None,
        access_key_id: Optional[str] = None,
        secret_access_key: Optional[str] = None
    ):
        """
        Initialize AWS Batch client.
        
        Args:
            region: AWS region (default: from env or eu-central-1)
            job_queue: Batch job queue name (default: from env or live2-job-queue)
            job_definition: Batch job definition name (default: from env or live2-simulation)
            access_key_id: AWS access key (default: from env)
            secret_access_key: AWS secret key (default: from env)
        """
        self.region = region or os.getenv("AWS_REGION", "eu-central-1")
        self.job_queue = job_queue or os.getenv("AWS_BATCH_JOB_QUEUE", "live2-job-queue")
        self.job_definition = job_definition or os.getenv("AWS_BATCH_JOB_DEFINITION", "live2-simulation")
        
        # Initialize boto3 client
        session_kwargs = {}
        if access_key_id or os.getenv("AWS_ACCESS_KEY_ID"):
            session_kwargs["aws_access_key_id"] = access_key_id or os.getenv("AWS_ACCESS_KEY_ID")
        if secret_access_key or os.getenv("AWS_SECRET_ACCESS_KEY"):
            session_kwargs["aws_secret_access_key"] = secret_access_key or os.getenv("AWS_SECRET_ACCESS_KEY")
        
        if session_kwargs:
            session = boto3.Session(**session_kwargs)
            self.batch_client = session.client("batch", region_name=self.region)
        else:
            # Use default credentials (IAM role, env vars, etc.)
            self.batch_client = boto3.client("batch", region_name=self.region)
        
        logger.info(f"AWS Batch client initialized: queue={self.job_queue}, definition={self.job_definition}")
    
    def submit_job(
        self,
        job_id: str,
        job_params: Dict,
        user_id: str,
        environment_vars: Optional[Dict[str, str]] = None
    ) -> str:
        """
        Submit job to AWS Batch.
        
        Args:
            job_id: Unique job ID (used as Batch job name)
            job_params: Job parameters (passed as environment variables)
            user_id: User ID (for S3 path)
            environment_vars: Additional environment variables
        
        Returns:
            AWS Batch job ID
        """
        # Prepare environment variables
        env = [
            {"name": "JOB_ID", "value": job_id},
            {"name": "USER_ID", "value": user_id},
            {"name": "JOB_PARAMS", "value": str(job_params)}  # JSON string
        ]
        
        # Add custom environment variables
        if environment_vars:
            for key, value in environment_vars.items():
                env.append({"name": key, "value": str(value)})
        
        # Prepare job parameters
        job_params_dict = {
            "jobName": f"live2-{job_id[:32]}",  # Batch job names must be <= 128 chars
            "jobQueue": self.job_queue,
            "jobDefinition": self.job_definition,
            "containerOverrides": {
                "environment": env
            }
        }
        
        try:
            response = self.batch_client.submit_job(**job_params_dict)
            batch_job_id = response["jobId"]
            logger.info(f"Submitted job {job_id} to AWS Batch: {batch_job_id}")
            return batch_job_id
        except ClientError as e:
            logger.error(f"Failed to submit job {job_id} to AWS Batch: {e}")
            raise Exception(f"AWS Batch submission failed: {str(e)}")
    
    def get_job_status(self, batch_job_id: str) -> Dict:
        """
        Get AWS Batch job status.
        
        Args:
            batch_job_id: AWS Batch job ID
        
        Returns:
            Dict with status, statusReason, startedAt, stoppedAt, etc.
        """
        try:
            response = self.batch_client.describe_jobs(jobs=[batch_job_id])
            if not response.get("jobs"):
                raise ValueError(f"Job {batch_job_id} not found")
            
            job = response["jobs"][0]
            
            # Map AWS Batch status to our status
            batch_status = job.get("status", "UNKNOWN")
            status_map = {
                "SUBMITTED": "queued",
                "PENDING": "queued",
                "RUNNABLE": "queued",
                "RUNNING": "running",
                "SUCCEEDED": "completed",
                "FAILED": "failed",
                "CANCELLED": "cancelled"
            }
            
            status = status_map.get(batch_status, "unknown")
            
            return {
                "status": status,
                "batch_status": batch_status,
                "status_reason": job.get("statusReason"),
                "started_at": job.get("startedAt"),
                "stopped_at": job.get("stoppedAt"),
                "created_at": job.get("createdAt"),
                "exit_code": job.get("container", {}).get("exitCode"),
                "log_stream_name": job.get("container", {}).get("logStreamName")
            }
        except ClientError as e:
            logger.error(f"Failed to get status for job {batch_job_id}: {e}")
            raise Exception(f"AWS Batch status check failed: {str(e)}")
    
    def cancel_job(self, batch_job_id: str, reason: str = "User requested cancellation") -> bool:
        """
        Cancel AWS Batch job.
        
        Args:
            batch_job_id: AWS Batch job ID
            reason: Cancellation reason
        
        Returns:
            True if successful
        """
        try:
            self.batch_client.cancel_job(
                jobId=batch_job_id,
                reason=reason
            )
            logger.info(f"Cancelled AWS Batch job {batch_job_id}: {reason}")
            return True
        except ClientError as e:
            logger.error(f"Failed to cancel job {batch_job_id}: {e}")
            raise Exception(f"AWS Batch cancellation failed: {str(e)}")
    
    def list_job_artifacts(self, job_id: str, user_id: str, bucket: str = "live2-artifacts") -> List[Dict]:
        """
        List S3 artifacts for a job.
        
        Args:
            job_id: Job ID
            user_id: User ID
            bucket: S3 bucket name
        
        Returns:
            List of artifact dicts with keys: key, size, last_modified
        """
        try:
            s3_client = self.batch_client._client_config.region_name  # Get region
            s3 = boto3.client("s3", region_name=self.region)
            
            prefix = f"prod/{user_id}/{job_id}/"
            response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
            
            artifacts = []
            if "Contents" in response:
                for obj in response["Contents"]:
                    artifacts.append({
                        "key": obj["Key"],
                        "size": obj["Size"],
                        "last_modified": obj["LastModified"].isoformat()
                    })
            
            return artifacts
        except ClientError as e:
            logger.error(f"Failed to list artifacts for job {job_id}: {e}")
            return []
    
    def generate_presigned_url(self, bucket: str, key: str, expiration: int = 3600) -> str:
        """
        Generate presigned URL for S3 object.
        
        Args:
            bucket: S3 bucket name
            key: S3 object key
            expiration: URL expiration in seconds (default: 1 hour)
        
        Returns:
            Presigned URL
        """
        try:
            # Use same credentials as batch client
            session_kwargs = {}
            if os.getenv("AWS_ACCESS_KEY_ID"):
                session_kwargs["aws_access_key_id"] = os.getenv("AWS_ACCESS_KEY_ID")
            if os.getenv("AWS_SECRET_ACCESS_KEY"):
                session_kwargs["aws_secret_access_key"] = os.getenv("AWS_SECRET_ACCESS_KEY")
            
            if session_kwargs:
                session = boto3.Session(**session_kwargs)
                s3 = session.client("s3", region_name=self.region)
            else:
                s3 = boto3.client("s3", region_name=self.region)
            
            url = s3.generate_presigned_url(
                "get_object",
                Params={"Bucket": bucket, "Key": key},
                ExpiresIn=expiration
            )
            return url
        except ClientError as e:
            logger.error(f"Failed to generate presigned URL for {key}: {e}")
            raise Exception(f"Failed to generate presigned URL: {str(e)}")

