"""
FastAPI dependencies for billing module.

Provides a single source of truth for "current user" resolution:
- Authorization: Bearer <JWT> (billing JWT)
- X-API-Key: <api_key> (API v1 compatible)
"""

from __future__ import annotations

import logging
from typing import Optional

from fastapi import Depends, HTTPException, Request, status
from sqlalchemy.orm import Session

from backend.billing.auth import verify_jwt_token
from backend.billing.database import get_db
from backend.billing.models import User, get_user_by_api_key

logger = logging.getLogger(__name__)


def _get_bearer_token(request: Request) -> Optional[str]:
    auth = request.headers.get("Authorization", "")
    if not auth:
        return None
    parts = auth.split(" ", 1)
    if len(parts) != 2:
        return None
    scheme, token = parts[0].strip().lower(), parts[1].strip()
    if scheme != "bearer" or not token:
        return None
    return token


async def get_current_user(
    request: Request,
    db: Session = Depends(get_db),
) -> User:
    """
    Resolve current user from request headers.

    Accepted:
    - Authorization: Bearer <JWT>
    - X-API-Key: <api_key>
    """
    # 1) JWT
    token = _get_bearer_token(request)
    if token:
        user_id = verify_jwt_token(token)
        if not user_id:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token")
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="User not found")
        return user

    # 2) API key
    api_key = request.headers.get("X-API-Key")
    if api_key:
        user = get_user_by_api_key(api_key, db_session=db)
        if not user:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid API key")
        return user

    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Missing authentication (Authorization: Bearer or X-API-Key)",
    )


