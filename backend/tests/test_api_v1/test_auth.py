"""
Tests for API v1 authentication.
"""

import pytest
from fastapi import Request, HTTPException
from unittest.mock import Mock, patch

from backend.api.v1.auth import verify_api_key, get_user_by_api_key, User


def test_get_user_by_api_key_not_found():
    """Test get_user_by_api_key returns None for invalid key."""
    user = get_user_by_api_key("invalid_key")
    assert user is None


@pytest.mark.asyncio
async def test_verify_api_key_missing():
    """Test verify_api_key raises 401 for missing API key."""
    request = Mock(spec=Request)
    request.headers.get.return_value = None
    
    with pytest.raises(HTTPException) as exc_info:
        await verify_api_key(request)
    
    assert exc_info.value.status_code == 401
    assert "Missing API key" in exc_info.value.detail


@pytest.mark.asyncio
async def test_verify_api_key_invalid():
    """Test verify_api_key raises 401 for invalid API key."""
    request = Mock(spec=Request)
    request.headers.get.return_value = "invalid_key"
    
    with patch('backend.api.v1.auth.get_user_by_api_key', return_value=None):
        with pytest.raises(HTTPException) as exc_info:
            await verify_api_key(request)
        
        assert exc_info.value.status_code == 401
        assert "Invalid or expired" in exc_info.value.detail


@pytest.mark.asyncio
async def test_verify_api_key_valid():
    """Test verify_api_key returns User for valid API key."""
    request = Mock(spec=Request)
    request.headers.get.return_value = "valid_key"
    
    mock_user = User(
        id="user_123",
        api_key="valid_key",
        tier="research",
        subscription_status="active"
    )
    
    with patch('backend.api.v1.auth.get_user_by_api_key', return_value=mock_user):
        user = await verify_api_key(request)
        assert user.id == "user_123"
        assert user.tier == "research"
        assert user.subscription_status == "active"


@pytest.mark.asyncio
async def test_verify_api_key_inactive():
    """Test verify_api_key raises 401 for inactive subscription."""
    request = Mock(spec=Request)
    request.headers.get.return_value = "expired_key"
    
    mock_user = User(
        id="user_123",
        api_key="expired_key",
        tier="research",
        subscription_status="expired"
    )
    
    with patch('backend.api.v1.auth.get_user_by_api_key', return_value=mock_user):
        with pytest.raises(HTTPException) as exc_info:
            await verify_api_key(request)
        
        assert exc_info.value.status_code == 401

