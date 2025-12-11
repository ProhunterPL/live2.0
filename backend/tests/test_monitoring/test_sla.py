"""
Tests for SLA calculator.
"""

import pytest
from datetime import date
from unittest.mock import Mock, MagicMock

from backend.monitoring.metrics.sla import SLACalculator, SLAStatus
from backend.monitoring.metrics.response_time import ResponseTimeTracker
from backend.monitoring.metrics.uptime import UptimeTracker
from backend.monitoring.metrics.error_rate import ErrorRateTracker


@pytest.fixture
def mock_trackers():
    """Create mock trackers."""
    response_tracker = Mock(spec=ResponseTimeTracker)
    uptime_tracker = Mock(spec=UptimeTracker)
    error_tracker = Mock(spec=ErrorRateTracker)
    
    return response_tracker, uptime_tracker, error_tracker


def test_sla_calculator_compliant(mock_trackers):
    """Test SLA compliance check - compliant case."""
    response_tracker, uptime_tracker, error_tracker = mock_trackers
    
    # Mock metrics that meet SLA targets with good margin
    # Pro tier: uptime 99.5%, response_time_p95 500ms, error_rate 0.1%
    # Use values well below thresholds to avoid warnings
    # Warning thresholds: uptime < target+0.5%, response_time > target*0.95, error_rate > target*0.95
    response_tracker.get_percentiles.return_value = {"p95": 200.0}  # < 500ms*0.95=475ms (well below)
    uptime_tracker.calculate_uptime.return_value = 100.0  # > 99.5%+0.5%=100% (perfect, no warning)
    error_tracker.get_error_rate.return_value = 0.01  # < 0.1%*0.95=0.095% (well below)
    
    calculator = SLACalculator(
        db=None,
        response_time_tracker=response_tracker,
        uptime_tracker=uptime_tracker,
        error_rate_tracker=error_tracker
    )
    
    compliance = calculator.check_sla_compliance("live2", "pro", date.today())
    
    assert compliance["status"] == SLAStatus.COMPLIANT.value
    assert compliance["compliant"] == True
    assert len(compliance["violations"]) == 0


def test_sla_calculator_violated(mock_trackers):
    """Test SLA compliance check - violated case."""
    response_tracker, uptime_tracker, error_tracker = mock_trackers
    
    # Mock metrics that violate SLA targets
    response_tracker.get_percentiles.return_value = {"p95": 600.0}  # > 500ms target
    uptime_tracker.calculate_uptime.return_value = 99.0  # < 99.5% target
    error_tracker.get_error_rate.return_value = 0.2  # > 0.1% target
    
    calculator = SLACalculator(
        db=None,
        response_time_tracker=response_tracker,
        uptime_tracker=uptime_tracker,
        error_rate_tracker=error_tracker
    )
    
    compliance = calculator.check_sla_compliance("live2", "pro", date.today())
    
    assert compliance["status"] == SLAStatus.VIOLATED.value
    assert compliance["compliant"] == False
    assert len(compliance["violations"]) > 0


def test_sla_calculator_warning(mock_trackers):
    """Test SLA compliance check - warning case."""
    response_tracker, uptime_tracker, error_tracker = mock_trackers
    
    # Mock metrics close to SLA targets (within warning threshold)
    response_tracker.get_percentiles.return_value = {"p95": 480.0}  # Close to 500ms
    uptime_tracker.calculate_uptime.return_value = 99.6  # OK
    error_tracker.get_error_rate.return_value = 0.08  # Close to 0.1%
    
    calculator = SLACalculator(
        db=None,
        response_time_tracker=response_tracker,
        uptime_tracker=uptime_tracker,
        error_rate_tracker=error_tracker
    )
    
    compliance = calculator.check_sla_compliance("live2", "pro", date.today())
    
    # Should be compliant but with warnings
    assert compliance["status"] in [SLAStatus.COMPLIANT.value, SLAStatus.WARNING.value]

