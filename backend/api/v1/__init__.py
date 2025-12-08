"""
API v1 for Live 2.0 Synthetic Data as a Service.

Provides REST API endpoints for:
- Dataset generation (async jobs)
- Molecule and reaction queries
- Reaction prediction
- Job status tracking
"""

from backend.api.v1.main import app

__all__ = ["app"]

