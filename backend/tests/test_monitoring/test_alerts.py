"""
Tests for alerts and remediation.
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from datetime import date

from backend.monitoring.alerts.notifier import AlertNotifier, Alert
from backend.monitoring.alerts.remediation import RemediationEngine
from backend.monitoring.metrics.sla import SLACalculator


@pytest.fixture
def mock_smtp():
    """Mock SMTP server."""
    with patch("smtplib.SMTP") as mock_smtp:
        mock_server = MagicMock()
        mock_smtp.return_value.__enter__.return_value = mock_server
        yield mock_server


def test_alert_notifier_send_email(mock_smtp):
    """Test email alert sending."""
    notifier = AlertNotifier(
        smtp_host="smtp.test.com",
        smtp_port=587,
        smtp_user="test@test.com",
        smtp_password="password",
        admin_emails=["admin@test.com"]
    )
    
    alert = Alert(
        severity="critical",
        title="Test Alert",
        message="Test message",
        project="live2",
        tier="hobby",
        violations=["Test violation"],
        timestamp=date.today()
    )
    
    result = notifier.send_alert(alert)
    
    # Should attempt to send email
    assert mock_smtp.starttls.called
    assert mock_smtp.login.called
    assert mock_smtp.send_message.called


def test_remediation_engine_check_compliant():
    """Test remediation engine with compliant SLA."""
    mock_calc = Mock(spec=SLACalculator)
    mock_calc.check_sla_compliance.return_value = {
        "status": "compliant",
        "compliant": True,
        "violations": []
    }
    
    engine = RemediationEngine(
        db=None,
        sla_calculator=mock_calc,
        alert_notifier=None
    )
    
    result = engine.check_and_remediate("live2", "hobby")
    
    assert result["remediated"] == False
    assert result["reason"] == "SLA compliant"


def test_remediation_engine_check_violated():
    """Test remediation engine with violated SLA."""
    mock_calc = Mock(spec=SLACalculator)
    mock_calc.check_sla_compliance.return_value = {
        "status": "violated",
        "compliant": False,
        "violations": ["Uptime below target"],
        "metrics": {
            "uptime": 94.0,
            "uptime_target": 95.0
        }
    }
    
    engine = RemediationEngine(
        db=None,
        sla_calculator=mock_calc,
        alert_notifier=None
    )
    
    result = engine.check_and_remediate("live2", "hobby")
    
    assert result["remediated"] == True  # Credit calculated
    assert "credit_amount" in result

