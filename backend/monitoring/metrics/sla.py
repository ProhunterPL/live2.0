"""
SLA compliance calculation per tier.
"""

import logging
from typing import Dict, List, Optional
from datetime import datetime, date, timedelta
from enum import Enum
from sqlalchemy.orm import Session

from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.uptime import UptimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker

logger = logging.getLogger(__name__)


class SLAStatus(str, Enum):
    """SLA compliance status."""
    COMPLIANT = "compliant"
    VIOLATED = "violated"
    WARNING = "warning"  # Close to violation


class SLACalculator:
    """
    Calculates SLA compliance per tier.
    
    Compares actual metrics against SLA targets from MONETIZATION_STATUS.md.
    """
    
    # SLA targets from MONETIZATION_STATUS.md
    SLA_TARGETS = {
        "legally": {
            "free": {
                "uptime": 95.0,
                "response_time_p95": 2000.0,  # 2s
                "error_rate": 1.0  # 1%
            },
            "starter": {
                "uptime": 99.0,
                "response_time_p95": 1000.0,  # 1s
                "error_rate": 0.5
            },
            "professional": {
                "uptime": 99.5,
                "response_time_p95": 500.0,  # 500ms
                "error_rate": 0.1
            },
            "law_firm": {
                "uptime": 99.9,
                "response_time_p95": 300.0,  # 300ms
                "error_rate": 0.1
            },
            "enterprise": {
                "uptime": 99.95,
                "response_time_p95": 200.0,  # 200ms
                "error_rate": 0.05
            }
        },
        "live2": {
            "hobby": {
                "uptime": 95.0,
                "response_time_p95": 10000.0,  # 10s for simulations
                "error_rate": 1.0
            },
            "research": {
                "uptime": 99.0,
                "response_time_p95": 5000.0,  # 5s
                "error_rate": 0.5
            },
            "pro": {
                "uptime": 99.5,
                "response_time_p95": 3000.0,  # 3s
                "error_rate": 0.1
            },
            "enterprise": {
                "uptime": 99.9,
                "response_time_p95": 2000.0,  # 2s
                "error_rate": 0.1
            }
        }
    }
    
    def __init__(
        self,
        db: Optional[Session],
        response_time_tracker: ResponseTimeTracker,
        uptime_tracker: UptimeTracker,
        error_rate_tracker: ErrorRateTracker
    ):
        """
        Initialize SLA calculator.
        
        Args:
            db: Database session (optional, for storing SLA records)
            response_time_tracker: Response time tracker
            uptime_tracker: Uptime tracker
            error_rate_tracker: Error rate tracker
        """
        self.db = db
        self.response_time_tracker = response_time_tracker
        self.uptime_tracker = uptime_tracker
        self.error_rate_tracker = error_rate_tracker
    
    def check_sla_compliance(
        self,
        project: str,  # "legally" or "live2"
        tier: str,
        month: Optional[date] = None
    ) -> Dict:
        """
        Check SLA compliance for tier in month.
        
        Args:
            project: Project name ("legally" or "live2")
            tier: Tier name
            month: Month to check (default: current month)
        
        Returns:
            Dict with compliance status and details
        """
        if month is None:
            month = date.today().replace(day=1)
        
        targets = self.SLA_TARGETS.get(project, {}).get(tier)
        if not targets:
            return {
                "status": SLAStatus.COMPLIANT.value,
                "message": "No SLA targets defined for tier",
                "project": project,
                "tier": tier,
                "month": month.isoformat()
            }
        
        # Get actual metrics
        end_date = (month + timedelta(days=32)).replace(day=1) - timedelta(days=1)
        uptime = self.uptime_tracker.calculate_uptime(tier, month, end_date)
        
        # Get response time (p95) for main endpoints
        if project == "live2":
            main_endpoint = "/api/v1/molecules"
        else:
            main_endpoint = "/api/subscription/status"
        
        response_times = self.response_time_tracker.get_percentiles(
            main_endpoint,
            tier,
            [95]
        )
        response_time_p95 = response_times.get("p95", 0.0)
        
        # Get error rate
        error_rate = self.error_rate_tracker.get_error_rate(tier, month)
        
        # Check compliance
        violations = []
        warnings = []
        
        # Uptime check
        if uptime < targets["uptime"]:
            violations.append(f"Uptime {uptime:.2f}% < {targets['uptime']}%")
        elif uptime < targets["uptime"] + 0.5:  # Within 0.5% of target
            warnings.append(f"Uptime {uptime:.2f}% close to target {targets['uptime']}%")
        
        # Response time check
        if response_time_p95 > targets["response_time_p95"]:
            violations.append(
                f"Response time p95 {response_time_p95:.0f}ms > {targets['response_time_p95']:.0f}ms"
            )
        elif response_time_p95 > targets["response_time_p95"] * 0.95:  # Within 5% of target (was 10%)
            warnings.append(
                f"Response time p95 {response_time_p95:.0f}ms close to target {targets['response_time_p95']:.0f}ms"
            )
        
        # Error rate check
        if error_rate > targets["error_rate"]:
            violations.append(
                f"Error rate {error_rate:.2f}% > {targets['error_rate']}%"
            )
        elif error_rate > targets["error_rate"] * 0.95:  # Within 5% of target (was 10%)
            warnings.append(
                f"Error rate {error_rate:.2f}% close to target {targets['error_rate']}%"
            )
        
        # Determine status
        if violations:
            status = SLAStatus.VIOLATED
        elif warnings:
            status = SLAStatus.WARNING
        else:
            status = SLAStatus.COMPLIANT
        
        return {
            "status": status.value,
            "month": month.isoformat(),
            "project": project,
            "tier": tier,
            "metrics": {
                "uptime": uptime,
                "uptime_target": targets["uptime"],
                "response_time_p95": response_time_p95,
                "response_time_p95_target": targets["response_time_p95"],
                "error_rate": error_rate,
                "error_rate_target": targets["error_rate"]
            },
            "violations": violations,
            "warnings": warnings,
            "compliant": len(violations) == 0
        }
