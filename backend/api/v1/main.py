"""
FastAPI app for API v1.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.api.v1.routes import datasets, simulations, molecules, reactions, predictions, jobs

# Import billing routes
try:
    from backend.billing.routes import auth as billing_auth, subscription, usage, webhooks, checkout
    BILLING_AVAILABLE = True
except ImportError:
    BILLING_AVAILABLE = False

# Create FastAPI app
app = FastAPI(
    title="Live 2.0 API v1",
    version="1.0.0",
    description="Synthetic Data as a Service API for Live 2.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify actual origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(datasets.router)
app.include_router(simulations.router)
app.include_router(molecules.router)
app.include_router(reactions.router)
app.include_router(predictions.router)
app.include_router(jobs.router)

# Include billing routes (if available)
if BILLING_AVAILABLE:
    app.include_router(billing_auth.router)
    app.include_router(subscription.router)
    app.include_router(usage.router)
    app.include_router(webhooks.router)
    app.include_router(checkout.router)


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "ok", "version": "1.0.0"}

