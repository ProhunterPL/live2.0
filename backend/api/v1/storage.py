"""
Storage manager for API v1 results.

Handles file storage (local, S3, MinIO).
"""

from typing import Optional
from pathlib import Path
import os
import logging
import shutil

logger = logging.getLogger(__name__)

# Try to import boto3 for S3 support
try:
    import boto3
    from botocore.exceptions import ClientError
    BOTO3_AVAILABLE = True
except ImportError:
    BOTO3_AVAILABLE = False
    logger.debug("boto3 not available (optional). Install with: pip install boto3")


class StorageManager:
    """
    Storage manager for API results (S3/MinIO).
    
    Placeholder implementation - can be extended with boto3 for S3.
    """
    
    def __init__(
        self,
        storage_type: str = "local",  # "local" | "s3" | "minio"
        base_path: str = "datasets/api",
        s3_bucket: Optional[str] = None,
        s3_region: Optional[str] = None,
        s3_access_key_id: Optional[str] = None,
        s3_secret_access_key: Optional[str] = None,
        s3_endpoint_url: Optional[str] = None
    ):
        """
        Initialize storage manager.
        
        Args:
            storage_type: Storage backend type
            base_path: Base path for local storage
            s3_bucket: S3 bucket name (if using S3)
            s3_region: S3 region (if using S3)
            s3_access_key_id: S3 access key ID
            s3_secret_access_key: S3 secret access key
            s3_endpoint_url: S3 endpoint URL (for S3-compatible services)
        """
        self.storage_type = storage_type
        self.base_path = Path(base_path)
        if storage_type == "local":
            self.base_path.mkdir(parents=True, exist_ok=True)
        self.s3_bucket = s3_bucket
        self.s3_region = s3_region
        self.s3_access_key_id = s3_access_key_id
        self.s3_secret_access_key = s3_secret_access_key
        self.s3_endpoint_url = s3_endpoint_url
        
        # Initialize S3 client if using S3
        self.s3_client = None
        if storage_type == "s3":
            if not BOTO3_AVAILABLE:
                raise ImportError("boto3 required for S3 storage. Install with: pip install boto3")
            if not s3_access_key_id or not s3_secret_access_key:
                raise ValueError("S3 credentials required for S3 storage")
            if not s3_bucket:
                raise ValueError("S3 bucket name required for S3 storage")
            
            # Create S3 client
            s3_config = {
                "aws_access_key_id": s3_access_key_id,
                "aws_secret_access_key": s3_secret_access_key,
                "region_name": s3_region or "us-east-1"
            }
            
            # Use custom endpoint for S3-compatible services (Supabase, MinIO)
            if s3_endpoint_url:
                s3_config["endpoint_url"] = s3_endpoint_url
            
            self.s3_client = boto3.client("s3", **s3_config)
    
    async def upload_file(self, local_path: str, job_id: str) -> str:
        """
        Upload file to storage and return URL.
        
        Args:
            local_path: Local file path
            job_id: Job ID (for naming)
        
        Returns:
            URL to access file
        """
        if self.storage_type == "local":
            # Copy to storage directory
            file_path = Path(local_path)
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {local_path}")
            
            dest_path = self.base_path / f"{job_id}{file_path.suffix}"
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            
            shutil.copy2(file_path, dest_path)
            logger.info(f"Uploaded file to local storage: {dest_path}")
            
            return f"/api/v1/download/{job_id}"
        
        elif self.storage_type == "s3":
            if not self.s3_client:
                raise RuntimeError("S3 client not initialized")
            
            file_path = Path(local_path)
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {local_path}")
            
            # S3 object key (path in bucket)
            s3_key = f"api/{job_id}{file_path.suffix}"
            
            try:
                # Upload file to S3
                self.s3_client.upload_file(
                    str(file_path),
                    self.s3_bucket,
                    s3_key,
                    ExtraArgs={"ContentType": self._get_content_type(file_path.suffix)}
                )
                logger.info(f"Uploaded file to S3: s3://{self.s3_bucket}/{s3_key}")
                
                # Generate presigned URL directly using the same key
                try:
                    url = self.s3_client.generate_presigned_url(
                        "get_object",
                        Params={"Bucket": self.s3_bucket, "Key": s3_key},
                        ExpiresIn=3600
                    )
                    return url
                except ClientError as e:
                    logger.error(f"Failed to generate presigned URL: {e}")
                    # Return S3 URL even if presigned URL generation fails
                    return f"s3://{self.s3_bucket}/{s3_key}"
            except ClientError as e:
                logger.error(f"Failed to upload file to S3: {e}")
                raise RuntimeError(f"S3 upload failed: {e}")
        
        else:
            raise ValueError(f"Unknown storage type: {self.storage_type}")
    
    def get_presigned_url(self, job_id: str, expires_in: int = 3600) -> str:
        """
        Get presigned URL for file download.
        
        Args:
            job_id: Job ID
            expires_in: URL expiration in seconds
        
        Returns:
            Presigned URL
        """
        if self.storage_type == "local":
            return f"/api/v1/download/{job_id}"
        
        elif self.storage_type == "s3":
            if not self.s3_client:
                raise RuntimeError("S3 client not initialized")
            
            # Try to find file extension by checking common formats
            s3_key = None
            for ext in [".parquet", ".json", ".graphml", ".csv", ".txt"]:
                test_key = f"api/{job_id}{ext}"
                try:
                    # Check if object exists
                    self.s3_client.head_object(Bucket=self.s3_bucket, Key=test_key)
                    s3_key = test_key
                    break
                except ClientError as e:
                    # Check if it's a 404 (not found) or other error
                    error_code = e.response.get("Error", {}).get("Code", "")
                    if error_code != "404":
                        logger.warning(f"Error checking S3 object {test_key}: {e}")
                    continue
            
            if not s3_key:
                raise FileNotFoundError(f"File not found in S3 for job_id: {job_id}")
            
            try:
                # Generate presigned URL
                url = self.s3_client.generate_presigned_url(
                    "get_object",
                    Params={"Bucket": self.s3_bucket, "Key": s3_key},
                    ExpiresIn=expires_in
                )
                return url
            except ClientError as e:
                logger.error(f"Failed to generate presigned URL: {e}")
                raise RuntimeError(f"Failed to generate presigned URL: {e}")
        
        else:
            raise ValueError(f"Unknown storage type: {self.storage_type}")
    
    def get_file_path(self, job_id: str) -> Optional[Path]:
        """
        Get local file path for job (local storage only).
        
        Args:
            job_id: Job ID
        
        Returns:
            Path to file or None if not found
        """
        if self.storage_type != "local":
            return None
        
        # Try to find file by job_id (check common extensions)
        for ext in [".parquet", ".json", ".graphml", ".csv"]:
            file_path = self.base_path / f"{job_id}{ext}"
            if file_path.exists():
                return file_path
        
        return None
    
    def _get_content_type(self, suffix: str) -> str:
        """Get content type for file extension."""
        content_types = {
            ".parquet": "application/parquet",
            ".json": "application/json",
            ".graphml": "application/xml",
            ".csv": "text/csv",
            ".txt": "text/plain"
        }
        return content_types.get(suffix.lower(), "application/octet-stream")

