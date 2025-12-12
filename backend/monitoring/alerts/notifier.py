"""
Alert notification system for SLA violations and system issues.
"""

import os
import logging
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, List, Optional
from datetime import datetime
from dataclasses import dataclass

logger = logging.getLogger(__name__)


def getenv_strip_quotes(key: str, default: Optional[str] = None) -> Optional[str]:
    """
    Get environment variable and strip quotes if present.
    
    This handles cases where values in .env files are quoted to escape special characters.
    For example: SMTP_PASSWORD="your-password#with#hash" will return "your-password#with#hash" (without quotes).
    
    Args:
        key: Environment variable key
        default: Default value if not found
    
    Returns:
        Environment variable value with quotes stripped, or default
    """
    value = os.getenv(key, default)
    if value and isinstance(value, str):
        # Strip surrounding quotes (single or double)
        if (value.startswith('"') and value.endswith('"')) or \
           (value.startswith("'") and value.endswith("'")):
            value = value[1:-1]
    return value


@dataclass
class Alert:
    """Alert data structure."""
    severity: str  # "critical", "warning", "info"
    title: str
    message: str
    project: str
    tier: str
    violations: List[str]
    timestamp: datetime


class AlertNotifier:
    """
    Sends alert notifications via email and other channels.
    """
    
    def __init__(
        self,
        smtp_host: Optional[str] = None,
        smtp_port: int = 587,
        smtp_user: Optional[str] = None,
        smtp_password: Optional[str] = None,
        from_email: Optional[str] = None,
        admin_emails: Optional[List[str]] = None
    ):
        """
        Initialize alert notifier.
        
        Args:
            smtp_host: SMTP server host
            smtp_port: SMTP server port
            smtp_user: SMTP username
            smtp_password: SMTP password
            from_email: From email address
            admin_emails: List of admin emails to notify
        """
        self.smtp_host = smtp_host or os.getenv("SMTP_HOST", "smtp.gmail.com")
        self.smtp_port = int(smtp_port or os.getenv("SMTP_PORT", "587"))
        self.smtp_user = smtp_user or getenv_strip_quotes("SMTP_USER")
        self.smtp_password = smtp_password or getenv_strip_quotes("SMTP_PASSWORD")
        self.from_email = from_email or getenv_strip_quotes("ALERT_FROM_EMAIL", "alerts@live2.ai")
        
        # Handle admin_emails - can be list or comma-separated string
        if admin_emails:
            if isinstance(admin_emails, list):
                self.admin_emails = admin_emails
            else:
                self.admin_emails = [e.strip() for e in admin_emails.split(",") if e.strip()]
        else:
            admin_emails_str = getenv_strip_quotes("ALERT_ADMIN_EMAILS", "")
            self.admin_emails = [e.strip() for e in admin_emails_str.split(",") if e.strip()] if admin_emails_str else []
        self.admin_emails = [e.strip() for e in self.admin_emails if e.strip()]
    
    def send_sla_violation_alert(
        self,
        project: str,
        tier: str,
        violations: List[str],
        metrics: Dict,
        month: str
    ) -> bool:
        """
        Send SLA violation alert.
        
        Args:
            project: Project name (live2 or legally)
            tier: Tier name
            violations: List of violation descriptions
            metrics: Current metrics dict
            month: Month string (YYYY-MM)
        
        Returns:
            True if sent successfully
        """
        severity = "critical" if len(violations) > 0 else "warning"
        
        alert = Alert(
            severity=severity,
            title=f"SLA Violation: {project} - {tier} tier",
            message=self._format_sla_violation_message(project, tier, violations, metrics, month),
            project=project,
            tier=tier,
            violations=violations,
            timestamp=datetime.utcnow()
        )
        
        return self.send_alert(alert)
    
    def send_alert(self, alert: Alert) -> bool:
        """
        Send alert notification.
        
        Args:
            alert: Alert object
        
        Returns:
            True if sent successfully
        """
        if not self.admin_emails:
            logger.warning("No admin emails configured, skipping alert")
            return False
        
        try:
            # Send email
            return self._send_email(alert)
        except Exception as e:
            logger.error(f"Failed to send alert: {e}")
            return False
    
    def _send_email(self, alert: Alert) -> bool:
        """
        Send email alert.
        
        Args:
            alert: Alert object
        
        Returns:
            True if sent successfully
        """
        if not self.smtp_user or not self.smtp_password:
            logger.warning("SMTP credentials not configured, skipping email alert")
            return False
        
        try:
            msg = MIMEMultipart()
            msg['From'] = self.from_email
            msg['To'] = ", ".join(self.admin_emails)
            msg['Subject'] = f"[{alert.severity.upper()}] {alert.title}"
            
            # HTML body
            html_body = self._format_email_html(alert)
            msg.attach(MIMEText(html_body, 'html'))
            
            # Plain text body
            text_body = self._format_email_text(alert)
            msg.attach(MIMEText(text_body, 'plain'))
            
            # Send email
            with smtplib.SMTP(self.smtp_host, self.smtp_port) as server:
                server.starttls()
                server.login(self.smtp_user, self.smtp_password)
                server.send_message(msg)
            
            logger.info(f"Alert sent successfully: {alert.title}")
            return True
        except Exception as e:
            logger.error(f"Failed to send email alert: {e}")
            return False
    
    def _format_sla_violation_message(
        self,
        project: str,
        tier: str,
        violations: List[str],
        metrics: Dict,
        month: str
    ) -> str:
        """Format SLA violation message."""
        lines = [
            f"SLA Violation detected for {project} - {tier} tier",
            f"Month: {month}",
            "",
            "Violations:",
        ]
        for v in violations:
            lines.append(f"  - {v}")
        
        lines.extend([
            "",
            "Current Metrics:",
            f"  - Uptime: {metrics.get('uptime', 0):.2f}% (target: {metrics.get('uptime_target', 0)}%)",
            f"  - Response Time (p95): {metrics.get('response_time_p95', 0):.0f}ms (target: {metrics.get('response_time_p95_target', 0)}ms)",
            f"  - Error Rate: {metrics.get('error_rate', 0):.2f}% (target: {metrics.get('error_rate_target', 0)}%)",
        ])
        
        return "\n".join(lines)
    
    def _format_email_html(self, alert: Alert) -> str:
        """Format email as HTML."""
        severity_color = {
            "critical": "#ef4444",
            "warning": "#f59e0b",
            "info": "#3b82f6"
        }.get(alert.severity, "#666")
        
        return f"""
        <html>
        <body>
            <h2 style="color: {severity_color};">{alert.title}</h2>
            <p><strong>Severity:</strong> {alert.severity.upper()}</p>
            <p><strong>Project:</strong> {alert.project}</p>
            <p><strong>Tier:</strong> {alert.tier}</p>
            <p><strong>Time:</strong> {alert.timestamp.isoformat()}</p>
            <hr>
            <pre style="background: #f5f5f5; padding: 10px; border-radius: 4px;">{alert.message}</pre>
        </body>
        </html>
        """
    
    def _format_email_text(self, alert: Alert) -> str:
        """Format email as plain text."""
        return f"""
{alert.title}

Severity: {alert.severity.upper()}
Project: {alert.project}
Tier: {alert.tier}
Time: {alert.timestamp.isoformat()}

{alert.message}
        """.strip()

