"""
Monitoring module for Live 2.0 and Legally.
"""

from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.uptime import UptimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker

# Optional imports (may fail if sqlalchemy not installed)
try:
    from backend.monitoring.metrics.sla import SLACalculator, SLAStatus
    _sla_available = True
except ImportError:
    SLACalculator = None  # type: ignore
    SLAStatus = None  # type: ignore
    _sla_available = False

__all__ = [
    "ResponseTimeTracker",
    "UptimeTracker",
    "ErrorRateTracker",
]

if _sla_available:
    __all__.extend(["SLACalculator", "SLAStatus"])
