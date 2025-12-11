# Live 2.0 & Legally Status Page

React frontend for displaying system status and SLA compliance.

## Setup

```bash
cd frontend/status
npm install
```

## Development

```bash
npm run dev
```

Opens at http://localhost:3001

## Build

```bash
npm run build
```

Outputs to `dist/` directory.

## Deployment

Build output can be served statically or mounted to FastAPI:

```python
# In backend/api/server.py
from fastapi.staticfiles import StaticFiles

app.mount("/status", StaticFiles(directory="frontend/status/dist", html=True), name="status")
```

## Environment Variables

Create `.env` file:

```
VITE_API_BASE=http://localhost:8000
```
