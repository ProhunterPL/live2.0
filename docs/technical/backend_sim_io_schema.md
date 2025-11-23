# Live 2.0 Data Schema Documentation

## Overview
This document describes the data formats and schemas used in the Live 2.0 simulation system.

## WebSocket Streaming Format

### Binary Data Format (msgpack)
The WebSocket stream sends binary data encoded with msgpack containing:

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
    "health_score": 0.85,
    ...
  }
}
```

## REST API Data Formats

### Simulation Configuration

#### Open Chemistry Mode
```json
{
  "mode": "open_chemistry",
  "config": {
    "grid_height": 256,
    "grid_width": 256,
    "dt": 0.01,
    "max_time": 1000.0,
    "energy_decay": 0.95,
    "energy_threshold": 0.1,
    "max_particles": 10000,
    "particle_radius": 0.5,
    "binding_threshold": 0.8,
    "unbinding_threshold": 0.2,
    "novelty_window": 100,
    "min_cluster_size": 2,
    "vis_frequency": 10,
    "log_frequency": 100,
    "seed": 42
  }
}
```

#### Preset Prebiotic Mode
```json
{
  "mode": "preset_prebiotic",
  "config": {
    "species": {
      "HCN": 0.1,
      "NH2CHO": 0.0,
      "H2O": 0.5
    },
    "reaction_rates": {
      "HCN_to_NH2CHO": 0.01
    },
    "diffusion_coeffs": {
      "HCN": 0.1,
      "NH2CHO": 0.05,
      "H2O": 0.2
    }
  }
}
```

### Simulation Status
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

### Novel Substance
```json
{
  "id": "SUB_abc12345_1234567890",
  "timestamp": 45.2,
  "size": 5,
  "complexity": 12.3,
  "properties": {
    "density": 0.4,
    "diameter": 3.2,
    "clustering_coeff": 0.6
  },
  "graph": {
    "particles": [1, 2, 3, 4, 5],
    "bonds": [[1, 2], [2, 3], [3, 4], [4, 5]],
    "particle_attributes": {
      "1": [1.0, 0.5, -0.3, 0.1],
      "2": [1.2, -0.2, 0.4, -0.1],
      ...
    }
  }
}
```

## Snapshot Format

### Snapshot File Structure
```json
{
  "config": { /* Simulation configuration */ },
  "current_time": 45.2,
  "step_count": 4520,
  "is_running": false,
  "is_paused": false,
  "mode": "open_chemistry",
  "particles": {
    "positions": [[x1, y1], [x2, y2], ...],
    "attributes": [[mass1, charge_x1, charge_y1, charge_z1], ...],
    "active_mask": [1, 1, 0, 1, ...]
  },
  "energy_field": [[e11, e12, ...], [e21, e22, ...], ...],
  "bonds": [[i1, j1, strength1], [i2, j2, strength2], ...],
  "clusters": [[p1, p2, p3], [p4, p5], ...],
  "novel_substances": [ /* Array of novel substances */ ],
  "metrics": { /* Aggregated metrics */ },
  "catalog_stats": { /* Substance catalog statistics */ },
  "timestamp": 1234567890.123
}
```

## Data Types

### Particle Attributes
- **mass**: float (positive)
- **charge_x**: float (can be negative)
- **charge_y**: float (can be negative)
- **charge_z**: float (can be negative)

### Position
- **x**: float (0 to grid_width)
- **y**: float (0 to grid_height)

### Bond
- **particle_i**: int (particle index)
- **particle_j**: int (particle index)
- **strength**: float (0.0 to 1.0)

### Energy Field
- 2D array of floats representing energy density at each grid point

### Metrics
```json
{
  "particle_count": 150,
  "total_energy": 1250.5,
  "total_mass": 180.3,
  "bond_count": 75,
  "cluster_count": 12,
  "energy_field_sum": 5000.0,
  "energy_field_max": 25.5,
  "energy_field_mean": 0.076,
  "novelty_rate": 0.15,
  "discovery_rate": 2.3,
  "total_substances": 45,
  "health_score": 0.85
}
```

## Error Formats

### API Error Response
```json
{
  "detail": "Error message",
  "status_code": 400,
  "timestamp": "2023-12-01T12:00:00Z"
}
```

### WebSocket Error
```json
{
  "type": "error",
  "message": "Error message",
  "code": 1008,
  "timestamp": "2023-12-01T12:00:00Z"
}
```

## Validation Rules

### Particle Data
- Positions must be within grid bounds
- Mass must be positive
- Active mask must be 0 or 1

### Bond Data
- Particle indices must be valid
- Strength must be between 0 and 1
- No self-bonds (i != j)

### Energy Field
- Must match grid dimensions
- Values must be non-negative
- Must be finite numbers

### Configuration
- Grid dimensions must be positive
- Time step must be positive
- Particle limits must be reasonable
- Thresholds must be between 0 and 1

## Performance Considerations

### WebSocket Streaming
- Binary format reduces bandwidth by ~60% compared to JSON
- 10 FPS update rate for smooth visualization
- Automatic client disconnection handling

### Snapshot Management
- Compressed JSON format
- Metadata tracking for quick listing
- Automatic cleanup of old snapshots

### Memory Usage
- Particle data: ~100 bytes per particle
- Energy field: ~1MB for 256x256 grid
- Bond matrix: ~400MB for 10k particles (sparse)

## Security Considerations

- Input validation on all API endpoints
- Rate limiting on WebSocket connections
- Snapshot file access restrictions
- Configuration parameter bounds checking
