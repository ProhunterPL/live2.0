# Live 2.0 API Protocol Documentation

## Overview
This document describes the REST API and WebSocket protocol for the Live 2.0 simulation system.

## Base URL
```
http://localhost:8000
```

## REST API Endpoints

### Simulation Management

#### Create Simulation
```http
POST /simulation/create
Content-Type: application/json

{
  "config": {
    "grid_height": 256,
    "grid_width": 256,
    "mode": "open_chemistry",
    "max_particles": 10000,
    "dt": 0.01
  },
  "mode": "open_chemistry"
}
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation created successfully",
  "simulation_id": "sim_1234567890"
}
```

#### Get Simulation Status
```http
GET /simulation/{simulation_id}/status
```

**Response:**
```json
{
  "simulation_id": "sim_1234567890",
  "is_running": true,
  "is_paused": false,
  "current_time": 45.2,
  "step_count": 4520,
  "particle_count": 150,
  "novelty_rate": 0.15,
  "health_score": 0.85
}
```

#### Start Simulation
```http
POST /simulation/{simulation_id}/start
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation started"
}
```

#### Pause Simulation
```http
POST /simulation/{simulation_id}/pause
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation paused"
}
```

#### Resume Simulation
```http
POST /simulation/{simulation_id}/resume
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation resumed"
}
```

#### Stop Simulation
```http
POST /simulation/{simulation_id}/stop
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation stopped"
}
```

#### Reset Simulation
```http
POST /simulation/{simulation_id}/reset
```

**Response:**
```json
{
  "success": true,
  "message": "Simulation reset"
}
```

### Data Access

#### Get Novel Substances
```http
GET /simulation/{simulation_id}/novel-substances?count=10
```

**Response:**
```json
{
  "substances": [
    {
      "id": "SUB_abc12345_1234567890",
      "timestamp": 45.2,
      "size": 5,
      "complexity": 12.3,
      "properties": {
        "density": 0.4,
        "diameter": 3.2
      }
    }
  ]
}
```

#### Get Metrics
```http
GET /simulation/{simulation_id}/metrics
```

**Response:**
```json
{
  "metrics": {
    "particle_count": 150,
    "total_energy": 1250.5,
    "bond_count": 75,
    "novelty_rate": 0.15,
    "health_score": 0.85
  }
}
```

### Snapshot Management

#### Save Snapshot
```http
POST /simulation/{simulation_id}/snapshot/save
Content-Type: application/json

{
  "filename": "my_snapshot.json"
}
```

**Response:**
```json
{
  "success": true,
  "filename": "my_snapshot.json"
}
```

#### Load Snapshot
```http
POST /simulation/{simulation_id}/snapshot/load
Content-Type: application/json

{
  "filename": "my_snapshot.json"
}
```

**Response:**
```json
{
  "success": true,
  "message": "Snapshot loaded"
}
```

## WebSocket Protocol

### Connection
```javascript
const ws = new WebSocket('ws://localhost:8000/simulation/{simulation_id}/stream');
```

### Data Format
The WebSocket sends binary data encoded with msgpack containing:

```json
{
  "particles": {
    "positions": [[x1, y1], [x2, y2], ...],
    "attributes": [[mass1, charge_x1, charge_y1, charge_z1], ...],
    "active_mask": [1, 1, 0, 1, ...]
  },
  "energy_field": [[e11, e12, ...], [e21, e22, ...], ...],
  "bonds": [[i1, j1, strength1], [i2, j2, strength2], ...],
  "clusters": [[p1, p2, p3], [p4, p5], ...],
  "metrics": {
    "particle_count": 150,
    "bond_count": 75,
    "novelty_rate": 0.15,
    "health_score": 0.85
  }
}
```

### Client Implementation Example
```javascript
import msgpack from 'msgpack-lite';

const ws = new WebSocket('ws://localhost:8000/simulation/sim_123/stream');

ws.onopen = function() {
  console.log('Connected to simulation stream');
};

ws.onmessage = function(event) {
  // Decode binary data
  const data = msgpack.decode(event.data);
  
  // Update visualization
  updateParticles(data.particles);
  updateEnergyField(data.energy_field);
  updateBonds(data.bonds);
  updateMetrics(data.metrics);
};

ws.onclose = function() {
  console.log('Disconnected from simulation stream');
};

ws.onerror = function(error) {
  console.error('WebSocket error:', error);
};
```

## Error Handling

### HTTP Error Responses
```json
{
  "detail": "Simulation not found",
  "status_code": 404
}
```

### WebSocket Error Codes
- `1008`: Policy violation (e.g., simulation not found)
- `1000`: Normal closure
- `1001`: Going away
- `1002`: Protocol error
- `1003`: Unsupported data

## Rate Limiting

### WebSocket Streaming
- Maximum 10 connections per simulation
- 10 FPS update rate
- Automatic disconnection after 5 minutes of inactivity

### REST API
- 100 requests per minute per IP
- 1000 requests per hour per IP

## Authentication

Currently, the API does not require authentication. In production, consider implementing:
- API key authentication
- JWT tokens
- OAuth 2.0

## CORS Configuration

The API allows cross-origin requests from any origin. In production, restrict to specific domains:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://yourdomain.com"],
    allow_credentials=True,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)
```

## Performance Considerations

### WebSocket Optimization
- Binary encoding reduces bandwidth by ~60%
- Automatic client disconnection handling
- Efficient data structures for real-time updates

### Memory Management
- Automatic cleanup of inactive simulations
- Snapshot compression
- Garbage collection of old data

### Scalability
- Multiple simulation instances
- Load balancing support
- Horizontal scaling capabilities

## Monitoring and Logging

### Metrics
- Connection count
- Message throughput
- Error rates
- Response times

### Logging
- Request/response logging
- Error logging
- Performance metrics
- Security events

## Security Best Practices

### Input Validation
- Validate all input parameters
- Sanitize file paths
- Check data types and ranges

### Rate Limiting
- Implement rate limiting
- Prevent abuse
- Monitor for suspicious activity

### Error Handling
- Don't expose internal errors
- Log errors securely
- Provide meaningful error messages

## Testing

### Unit Tests
```python
def test_create_simulation():
    response = client.post("/simulation/create", json={
        "config": {"grid_height": 256, "grid_width": 256},
        "mode": "open_chemistry"
    })
    assert response.status_code == 200
    assert response.json()["success"] == True
```

### Integration Tests
```python
def test_websocket_streaming():
    with client.websocket_connect("/simulation/sim_123/stream") as websocket:
        data = websocket.receive_bytes()
        assert len(data) > 0
```

### Load Testing
- Test with multiple concurrent connections
- Measure performance under load
- Identify bottlenecks

## Deployment

### Docker
```dockerfile
FROM python:3.9-slim
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
CMD ["uvicorn", "api.server:app", "--host", "0.0.0.0", "--port", "8000"]
```

### Environment Variables
```bash
LIVE2_HOST=0.0.0.0
LIVE2_PORT=8000
LIVE2_WORKERS=4
LIVE2_LOG_LEVEL=info
```

### Production Configuration
- Use production ASGI server (e.g., Gunicorn)
- Enable HTTPS
- Configure reverse proxy
- Set up monitoring
- Implement logging
- Configure backup strategies
