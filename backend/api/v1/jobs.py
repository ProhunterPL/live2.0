"""
Async job processor for API v1.

Handles long-running operations like dataset generation.
"""

from typing import Dict, Optional, List, TYPE_CHECKING
from enum import Enum
import uuid
import asyncio
import logging
from datetime import datetime
import json
from fastapi import HTTPException, status

if TYPE_CHECKING:
    import redis

try:
    import redis
except ImportError:
    redis = None  # type: ignore

from backend.dataset_export import DatasetExporter

logger = logging.getLogger(__name__)


class JobStatus(str, Enum):
    """Job status enum."""
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobProcessor:
    """
    Async job processor for long-running operations.
    
    Uses Redis for job queue and status tracking.
    """
    
    def __init__(
        self, 
        redis_client: Optional["redis.Redis"] = None,
        base_results_dir: str = "results/phase2b_additional",
        storage_manager: Optional["StorageManager"] = None
    ):
        """
        Initialize job processor.
        
        Args:
            redis_client: Redis client
            base_results_dir: Base directory with simulation results
            storage_manager: Storage manager for results
        """
        self.redis = redis_client
        self.base_results_dir = base_results_dir
        self.storage = storage_manager
        self.exporter = DatasetExporter(base_results_dir)
        self._processor_task: Optional[asyncio.Task] = None
        self._start_processor()
    
    def _start_processor(self):
        """Start background job processor task."""
        if self._processor_task is None or self._processor_task.done():
            self._processor_task = asyncio.create_task(self._process_jobs())
            logger.info("Started job processor background task")
    
    async def create_job(
        self,
        job_type: str,
        params: Dict,
        user_id: str,
        webhook_url: Optional[str] = None
    ) -> str:
        """
        Create async job and return job_id.
        
        Args:
            job_type: "generate_dataset" or "run_simulation"
            params: Job parameters
            user_id: User ID who created the job
            webhook_url: Optional webhook URL for notifications
        
        Returns:
            Job ID string
        """
        job_id = f"job_{uuid.uuid4().hex[:12]}"
        
        job_data = {
            "job_id": job_id,
            "job_type": job_type,
            "params": params,
            "user_id": user_id,
            "status": JobStatus.QUEUED.value,
            "progress": 0,
            "created_at": datetime.utcnow().isoformat(),
            "webhook_url": webhook_url,
            "error": None,
            "result_url": None
        }
        
        # Store job in Redis
        job_key = f"job:{job_id}"
        self.redis.setex(
            job_key,
            86400 * 7,  # 7 days TTL
            json.dumps(job_data)
        )
        
        # Queue job for processing
        self.redis.lpush("job_queue", job_id)
        
        logger.info(f"Created job {job_id} of type {job_type}")
        return job_id
    
    async def get_job_status(self, job_id: str) -> Dict:
        """
        Get job status.
        
        Args:
            job_id: Job ID
        
        Returns:
            Job status dict
        
        Raises:
            HTTPException 404: Job not found
        """
        job_key = f"job:{job_id}"
        job_data = self.redis.get(job_key)
        
        if not job_data:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Job {job_id} not found"
            )
        
        return json.loads(job_data)
    
    async def list_user_jobs(self, user_id: str, limit: int = 20) -> List[Dict]:
        """
        List all jobs for a user.
        
        Args:
            user_id: User ID
            limit: Maximum number of jobs to return
        
        Returns:
            List of job status dicts
        """
        jobs = []
        
        # Scan all job keys (job:*)
        cursor = 0
        while True:
            cursor, keys = self.redis.scan(cursor, match="job:*", count=100)
            
            for key in keys:
                try:
                    job_data = self.redis.get(key)
                    if job_data:
                        job_dict = json.loads(job_data)
                        # Filter by user_id
                        if job_dict.get("user_id") == user_id:
                            jobs.append(job_dict)
                except (json.JSONDecodeError, KeyError) as e:
                    logger.warning(f"Failed to parse job data from {key}: {e}")
                    continue
            
            if cursor == 0:
                break
        
        # Sort by created_at (newest first)
        jobs.sort(key=lambda x: x.get("created_at", ""), reverse=True)
        
        # Limit results
        return jobs[:limit]
    
    async def _process_jobs(self):
        """Background task to process queued jobs."""
        while True:
            try:
                # Pop job from queue (blocking)
                result = self.redis.brpop("job_queue", timeout=1)
                if not result:
                    continue
                
                job_id = result[1]  # brpop returns (key, value)
                
                # Get job data
                job_key = f"job:{job_id}"
                job_data_str = self.redis.get(job_key)
                if not job_data_str:
                    logger.warning(f"Job {job_id} not found in Redis")
                    continue
                
                job_data = json.loads(job_data_str)
                
                # Update status to running
                job_data["status"] = JobStatus.RUNNING.value
                self.redis.setex(job_key, 86400 * 7, json.dumps(job_data))
                
                # Process job
                try:
                    if job_data["job_type"] == "generate_dataset":
                        result = await self._process_dataset_job(job_data)
                    elif job_data["job_type"] == "run_simulation":
                        result = await self._process_simulation_job(job_data)
                    else:
                        raise ValueError(f"Unknown job type: {job_data['job_type']}")
                    
                    # Update job as completed
                    job_data["status"] = JobStatus.COMPLETED.value
                    job_data["progress"] = 100
                    job_data["result_url"] = result["url"]
                    job_data["completed_at"] = datetime.utcnow().isoformat()
                    
                except Exception as e:
                    logger.error(f"Job {job_id} failed: {e}", exc_info=True)
                    job_data["status"] = JobStatus.FAILED.value
                    job_data["error"] = str(e)
                
                # Save updated job data
                self.redis.setex(job_key, 86400 * 7, json.dumps(job_data))
                
                # Send webhook if configured
                if job_data.get("webhook_url"):
                    await self._send_webhook(job_data)
            
            except Exception as e:
                logger.error(f"Error processing jobs: {e}", exc_info=True)
                await asyncio.sleep(1)
    
    async def _process_dataset_job(self, job_data: Dict) -> Dict:
        """
        Process dataset generation job.
        
        Args:
            job_data: Job data dict
        
        Returns:
            Dict with result URL
        """
        params = job_data["params"]
        dataset_type = params["dataset_type"]
        
        # Create progress callback
        def progress_callback(current: int, total: int, message: str):
            job_key = f"job:{job_data['job_id']}"
            try:
                job_data_local = json.loads(self.redis.get(job_key) or "{}")
                job_data_local["progress"] = int(100 * current / total) if total > 0 else 0
                self.redis.setex(job_key, 86400 * 7, json.dumps(job_data_local))
            except Exception as e:
                logger.warning(f"Failed to update progress: {e}")
        
        # Run export in thread pool (blocking operation)
        loop = asyncio.get_event_loop()
        
        try:
            if dataset_type == "reaction_trajectories":
                output_path = await loop.run_in_executor(
                    None,
                    lambda: self.exporter.export_reaction_trajectories(
                        runs=params["runs"],
                        output_format=params.get("output_format", "parquet"),
                        filters=params.get("filters"),
                        progress_callback=progress_callback
                    )
                )
            elif dataset_type == "autocatalysis_network":
                output_path = await loop.run_in_executor(
                    None,
                    lambda: self.exporter.export_autocatalysis_network(
                        runs=params["runs"],
                        output_format=params.get("output_format", "json"),
                        include_metrics=params.get("include_metrics", True)
                    )
                )
            elif dataset_type == "novel_molecules":
                output_path = await loop.run_in_executor(
                    None,
                    lambda: self.exporter.export_novel_molecules(
                        runs=params["runs"],
                        novelty_threshold=params.get("novelty_threshold", 0.7),
                        limit=params.get("limit"),
                        output_format=params.get("output_format", "json"),
                        include_graphs=params.get("include_graphs", True)
                    )
                )
            else:
                raise ValueError(f"Unknown dataset type: {dataset_type}")
            
            # Upload to storage
            if self.storage:
                result_url = await self.storage.upload_file(output_path, job_data["job_id"])
            else:
                # Fallback: return local path
                result_url = f"/api/v1/download/{job_data['job_id']}"
            
            return {"url": result_url}
        
        except Exception as e:
            logger.error(f"Dataset job {job_data['job_id']} failed: {e}", exc_info=True)
            raise
    
    async def _process_simulation_job(self, job_data: Dict) -> Dict:
        """
        Process simulation job (placeholder - wymaga integracji z symulacją).
        
        Args:
            job_data: Job data dict
        
        Returns:
            Dict with result URL
        """
        # TODO: Implement simulation job processing
        # This requires integration with backend/sim/core/stepper.py
        raise NotImplementedError("Simulation jobs not yet implemented")
    
    async def _send_webhook(self, job_data: Dict):
        """Send webhook notification (placeholder)."""
        # TODO: Implement webhook sending with httpx
        logger.info(f"Webhook notification for job {job_data['job_id']} (not implemented)")

