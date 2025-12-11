"""
Status page API endpoints.
"""

from fastapi import APIRouter, Depends, Query
from typing import Dict, Optional
from datetime import datetime, date

from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.uptime import UptimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker

# Optional import (may fail if sqlalchemy not installed)
try:
    from backend.monitoring.metrics.sla import SLACalculator
except ImportError:
    SLACalculator = None  # type: ignore

router = APIRouter(prefix="/status", tags=["status"])


# Dependencies (singleton pattern)
_response_time_tracker: Optional[ResponseTimeTracker] = None
_uptime_tracker: Optional[UptimeTracker] = None
_error_rate_tracker: Optional[ErrorRateTracker] = None
_sla_calculator: Optional[SLACalculator] = None


def get_response_time_tracker() -> ResponseTimeTracker:
    """Get response time tracker instance."""
    global _response_time_tracker
    if _response_time_tracker is None:
        _response_time_tracker = ResponseTimeTracker()
    return _response_time_tracker


def get_uptime_tracker() -> UptimeTracker:
    """Get uptime tracker instance."""
    global _uptime_tracker
    if _uptime_tracker is None:
        _uptime_tracker = UptimeTracker()
    return _uptime_tracker


def get_error_rate_tracker() -> ErrorRateTracker:
    """Get error rate tracker instance."""
    global _error_rate_tracker
    if _error_rate_tracker is None:
        _error_rate_tracker = ErrorRateTracker()
    return _error_rate_tracker


def get_sla_calculator() -> Optional[SLACalculator]:
    """Get SLA calculator instance."""
    if SLACalculator is None:
        return None
    
    global _sla_calculator
    if _sla_calculator is None:
        # Try to get DB, but allow None if billing module not available
        db = None
        try:
            from backend.billing.database import get_db
            db = next(get_db())
        except:
            pass
        
        _sla_calculator = SLACalculator(
            db=db,
            response_time_tracker=get_response_time_tracker(),
            uptime_tracker=get_uptime_tracker(),
            error_rate_tracker=get_error_rate_tracker()
        )
    return _sla_calculator


@router.get("/health")
async def health_check():
    """Public health check endpoint."""
    return {
        "status": "operational",
        "timestamp": datetime.utcnow().isoformat()
    }


@router.get("/metrics")
async def get_metrics(
    project: str = Query("live2", description="Project name (live2 or legally)"),
    tier: str = Query("hobby", description="Tier name")
):
    """
    Get current metrics for status page.
    
    Returns:
        Dict with current metrics (uptime, response time, error rate)
    """
    response_time_tracker = get_response_time_tracker()
    uptime_tracker = get_uptime_tracker()
    error_rate_tracker = get_error_rate_tracker()
    
    # Get current metrics
    if project == "live2":
        main_endpoint = "/api/v1/molecules"
    else:
        main_endpoint = "/api/subscription/status"
    
    percentiles = response_time_tracker.get_percentiles(main_endpoint, tier)
    uptime = uptime_tracker.get_current_uptime(tier)
    error_rate = error_rate_tracker.get_current_error_rate(tier)
    
    return {
        "project": project,
        "tier": tier,
        "uptime": uptime,
        "response_time": {
            "p50": percentiles.get("p50", 0),
            "p95": percentiles.get("p95", 0),
            "p99": percentiles.get("p99", 0)
        },
        "error_rate": error_rate,
        "timestamp": datetime.utcnow().isoformat()
    }


@router.get("/sla")
async def get_sla_compliance(
    project: str = Query("live2", description="Project name"),
    tier: str = Query("hobby", description="Tier name"),
    month: Optional[str] = Query(None, description="Month in YYYY-MM format")
):
    """
    Get SLA compliance status.
    
    Args:
        project: Project name
        tier: Tier name
        month: Month in YYYY-MM format (default: current month)
    
    Returns:
        SLA compliance status
    """
    sla_calculator = get_sla_calculator()
    
    if sla_calculator is None:
        from fastapi import HTTPException
        raise HTTPException(
            status_code=503,
            detail="SLA calculator not available. SQLAlchemy may not be installed."
        )
    
    if month:
        try:
            month_date = datetime.strptime(month, "%Y-%m").date().replace(day=1)
        except ValueError:
            return {"error": "Invalid month format. Use YYYY-MM"}
    else:
        month_date = None
    
    compliance = sla_calculator.check_sla_compliance(project, tier, month_date)
    return compliance
