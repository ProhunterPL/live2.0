"""
FastAPI server for Live 2.0 simulation
Provides WebSocket streaming and REST API endpoints
"""

import asyncio
import json
import time
import logging
from typing import Dict, List, Optional, Any
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import uvicorn
import numpy as np
import msgpack
import taichi as ti

# Setup logging to file with rotation (max 5MB)
from logging.handlers import RotatingFileHandler

# Create rotating file handler with max 5MB per file, keep 3 backup files
file_handler = RotatingFileHandler(
    'logs.txt', 
    maxBytes=5*1024*1024,  # 5MB
    backupCount=3
)

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        file_handler,
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Initialize Taichi: prefer CUDA, then Vulkan, then CPU
try:
    if hasattr(ti, 'cuda') and ti.cuda.is_available():
        ti.init(arch=ti.cuda)
    elif hasattr(ti, 'vulkan'):
        ti.init(arch=ti.vulkan)
    else:
        ti.init(arch=ti.cpu)
except Exception:
    ti.init(arch=ti.cpu)

from sim.config import SimulationConfig, PresetPrebioticConfig, OpenChemistryConfig
from sim.core.stepper import SimulationStepper
from sim.io.snapshot import SnapshotManager

# Pydantic models for API
class SimulationRequest(BaseModel):
    config: Dict[str, Any]
    mode: str = "open_chemistry"

class SimulationResponse(BaseModel):
    success: bool
    message: str
    simulation_id: Optional[str] = None

class SimulationStatus(BaseModel):
    simulation_id: str
    is_running: bool
    is_paused: bool
    current_time: float
    step_count: int
    particle_count: int
    novelty_rate: float
    health_score: float

class NovelSubstance(BaseModel):
    id: str
    timestamp: float
    size: int
    complexity: float
    properties: Dict[str, Any]

class Live2Server:
    """Main server class for Live 2.0 simulation"""
    
    def __init__(self):
        self.app = FastAPI(title="Live 2.0 Simulation API", version="1.0.0")
        self.setup_cors()
        self.setup_routes()
        
        # Simulation management
        self.simulations: Dict[str, SimulationStepper] = {}
        self.active_connections: Dict[str, List[WebSocket]] = {}
        self.snapshot_manager = SnapshotManager()
        
        # WebSocket broadcasting and simulation loops
        self.broadcast_tasks: Dict[str, asyncio.Task] = {}
        self.simulation_tasks: Dict[str, asyncio.Task] = {}
        
        # Default configuration
        self.default_config = SimulationConfig()
    
    def setup_cors(self):
        """Setup CORS middleware"""
        self.app.add_middleware(
            CORSMiddleware,
            allow_origins=["*"],  # In production, specify actual origins
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )
    
    def setup_routes(self):
        """Setup API routes"""
        
        @self.app.get("/")
        async def root():
            return {"message": "Live 2.0 Simulation API", "version": "1.0.0"}
        
        @self.app.post("/simulation/create", response_model=SimulationResponse)
        async def create_simulation(request: SimulationRequest):
            """Create a new simulation"""
            try:
                # Enforce single simulation: stop and cleanup existing ones
                self.stop_all_simulations()
                # Create simulation ID
                simulation_id = f"sim_{int(time.time() * 1000)}"
                
                # Parse configuration
                if request.mode == "preset_prebiotic":
                    config = PresetPrebioticConfig(**request.config)
                elif request.mode == "open_chemistry":
                    config = SimulationConfig(**request.config)
                else:
                    raise ValueError(f"Unknown mode: {request.mode}")
                
                # Create simulation
                simulation = SimulationStepper(config)
                self.simulations[simulation_id] = simulation
                self.active_connections[simulation_id] = []
                
                return SimulationResponse(
                    success=True,
                    message="Simulation created successfully",
                    simulation_id=simulation_id
                )
            except Exception as e:
                raise HTTPException(status_code=400, detail=str(e))
        
        @self.app.get("/simulation/{simulation_id}/status", response_model=SimulationStatus)
        async def get_simulation_status(simulation_id: str):
            """Get simulation status"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            state = simulation.get_simulation_state()
            
            return SimulationStatus(
                simulation_id=simulation_id,
                is_running=state['is_running'],
                is_paused=state['is_paused'],
                current_time=state['current_time'],
                step_count=state['step_count'],
                particle_count=state['particle_count'],
                novelty_rate=state['novelty_rate'],
                health_score=state['health_score']
            )
        
        @self.app.post("/simulation/{simulation_id}/start")
        async def start_simulation(simulation_id: str):
            """Start simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.start()
            # Ensure simulation loop is running in background (non-blocking)
            self.start_simulation_loop(simulation_id)
            
            return {"success": True, "message": "Simulation started"}
        
        @self.app.post("/simulation/{simulation_id}/pause")
        async def pause_simulation(simulation_id: str):
            """Pause simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.pause()
            
            return {"success": True, "message": "Simulation paused"}
        
        @self.app.post("/simulation/{simulation_id}/resume")
        async def resume_simulation(simulation_id: str):
            """Resume simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.resume()
            
            return {"success": True, "message": "Simulation resumed"}
        
        @self.app.post("/simulation/{simulation_id}/stop")
        async def stop_simulation(simulation_id: str):
            """Stop simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.stop()
            
            # Stop broadcasting
            if simulation_id in self.broadcast_tasks:
                self.broadcast_tasks[simulation_id].cancel()
                del self.broadcast_tasks[simulation_id]
            # Stop simulation loop and cleanup resources
            if simulation_id in self.simulation_tasks:
                self.simulation_tasks[simulation_id].cancel()
                del self.simulation_tasks[simulation_id]
            self.cleanup_simulation(simulation_id)
            
            return {"success": True, "message": "Simulation stopped"}
        
        @self.app.post("/simulation/{simulation_id}/reset")
        async def reset_simulation(simulation_id: str):
            """Reset simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.reset()
            
            return {"success": True, "message": "Simulation reset"}
        
        @self.app.get("/simulation/{simulation_id}/novel-substances")
        async def get_novel_substances(simulation_id: str, count: int = 10):
            """Get recent novel substances"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            substances = simulation.get_novel_substances(count)
            
            return {"substances": substances}
        
        @self.app.post("/simulation/{simulation_id}/snapshot/save")
        async def save_snapshot(simulation_id: str, filename: str = None):
            """Save simulation snapshot"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            
            if filename is None:
                filename = f"snapshot_{simulation_id}_{int(time.time())}.json"
            
            simulation.save_snapshot(filename)
            
            return {"success": True, "filename": filename}
        
        @self.app.post("/simulation/{simulation_id}/snapshot/load")
        async def load_snapshot(simulation_id: str, filename: str):
            """Load simulation snapshot"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            simulation.load_snapshot(filename)
            
            return {"success": True, "message": "Snapshot loaded"}
        
        @self.app.get("/simulation/{simulation_id}/metrics")
        async def get_metrics(simulation_id: str):
            """Get simulation metrics"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            metrics = simulation.aggregator.get_aggregated_stats()
            
            return {"metrics": metrics}
        
        @self.app.get("/simulation/{simulation_id}/novel-substances")
        async def get_novel_substances(simulation_id: str, limit: int = 10):
            """Get novel substances discovered in simulation"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            
            # Get novel substances from simulation state or catalog
            try:
                if hasattr(simulation, 'catalog') and hasattr(simulation.catalog, 'get_novel_substances'):
                    substances = simulation.catalog.get_novel_substances(limit)
                else:
                    # Fallback: return empty list or mock data
                    substances = []
                
                return {"substances": substances}
            except Exception as e:
                logger.warning(f"Failed to get novel substances for simulation {simulation_id}: {e}")
                return {"substances": []}
        
        @self.app.websocket("/simulation/{simulation_id}/stream")
        async def websocket_endpoint(websocket: WebSocket, simulation_id: str):
            """WebSocket endpoint for real-time data streaming"""
            logger.info(f"WebSocket connection attempt for simulation {simulation_id}")
            
            # Validate simulation exists before accepting connection
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                try:
                    await websocket.close(code=1008, reason="Simulation not found")
                except Exception as e:
                    logger.warning(f"Failed to close WebSocket for simulation {simulation_id}: {e}")
                return
            
            await websocket.accept()
            logger.info(f"WebSocket accepted for simulation {simulation_id}")
            
            # Add connection to active connections
            if simulation_id not in self.active_connections:
                self.active_connections[simulation_id] = []
            self.active_connections[simulation_id].append(websocket)
            logger.info(f"Added WebSocket connection for simulation {simulation_id}")
            
            # Don't auto-start simulation - wait for explicit start command from frontend
            simulation = self.simulations[simulation_id]
            logger.info(f"WebSocket connected for simulation {simulation_id} (not auto-started)")
            
            # Start broadcasting if not already started
            if simulation_id not in self.broadcast_tasks:
                logger.info(f"Starting broadcast task for simulation {simulation_id}")
                self.broadcast_tasks[simulation_id] = asyncio.create_task(
                    self.broadcast_simulation_data(simulation_id)
                )
            
            try:
                # Passive keep-alive; disconnections are handled in broadcast on send failure
                while True:
                    await asyncio.sleep(30)
            except WebSocketDisconnect:
                # Remove connection
                if websocket in self.active_connections[simulation_id]:
                    self.active_connections[simulation_id].remove(websocket)
                
                # Stop broadcasting if no connections
                if not self.active_connections[simulation_id]:
                    if simulation_id in self.broadcast_tasks:
                        self.broadcast_tasks[simulation_id].cancel()
                        del self.broadcast_tasks[simulation_id]
    
    async def broadcast_simulation_data(self, simulation_id: str):
        """Broadcast simulation data to connected clients"""
        logger.info(f"Starting broadcast for simulation {simulation_id}")
        
        # Check if simulation exists
        if simulation_id not in self.simulations:
            logger.error(f"Simulation {simulation_id} not found for broadcast")
            return
            
        simulation = self.simulations[simulation_id]
        
        while True:
            try:
                # Check if simulation still exists
                if simulation_id not in self.simulations:
                    logger.info(f"Simulation {simulation_id} removed, stopping broadcast")
                    break
                # Skip frames to reduce bandwidth based on vis_frequency
                vis_freq = getattr(simulation.config, 'vis_frequency', 10)
                if vis_freq is None or vis_freq <= 0:
                    vis_freq = 10
                if simulation.step_count % int(vis_freq) != 0:
                    await asyncio.sleep(0.01)
                    continue

                # Get visualization data
                logger.debug(f"Getting visualization data for simulation {simulation_id}")
                data = simulation.get_visualization_data()

                # Lightweight downsampling for large fields
                def _downsample_2x(arr2d):
                    try:
                        import numpy as _np
                        a = _np.asarray(arr2d)
                        h, w = a.shape
                        if h >= 2 and w >= 2:
                            a = a[:h - (h % 2), :w - (w % 2)]
                            a = a.reshape(a.shape[0]//2, 2, a.shape[1]//2, 2).mean(axis=(1, 3))
                        return a.tolist()
                    except Exception:
                        return arr2d

                if isinstance(data.get('energy_field'), list) and len(data['energy_field']) >= 256:
                    data['energy_field'] = _downsample_2x(data['energy_field'])
                if isinstance(data.get('concentration_view'), list) and len(data['concentration_view']) >= 256:
                    data['concentration_view'] = _downsample_2x(data['concentration_view'])
                
                # Convert numpy types to native Python types for serialization
                def convert_numpy_types(obj):
                    import numpy as np
                    if isinstance(obj, np.ndarray):
                        return obj.tolist()
                    elif isinstance(obj, dict):
                        return {k: convert_numpy_types(v) for k, v in obj.items()}
                    elif isinstance(obj, list):
                        return [convert_numpy_types(item) for item in obj]
                    elif isinstance(obj, tuple):
                        return tuple(convert_numpy_types(item) for item in obj)
                    elif hasattr(obj, 'dtype'):  # Any numpy scalar
                        if np.issubdtype(obj.dtype, np.integer):
                            return int(obj)
                        elif np.issubdtype(obj.dtype, np.floating):
                            return float(obj)
                        elif np.issubdtype(obj.dtype, np.bool_):
                            return bool(obj)
                        else:
                            return obj.item() if hasattr(obj, 'item') else str(obj)
                    elif hasattr(obj, 'item'):  # Handle numpy scalars
                        return obj.item()
                    else:
                        return obj
                
                # Convert numpy types before serialization
                serializable_data = convert_numpy_types(data)
                binary_data = msgpack.packb(serializable_data)
                
                # Send to all connected clients
                disconnected = []
                if simulation_id in self.active_connections:
                    for websocket in self.active_connections[simulation_id]:
                        try:
                            await websocket.send_bytes(binary_data)
                        except:
                            disconnected.append(websocket)
                
                # Remove disconnected clients
                if simulation_id in self.active_connections:
                    for websocket in disconnected:
                        if websocket in self.active_connections[simulation_id]:
                            self.active_connections[simulation_id].remove(websocket)
                
                # Wait before next broadcast (about 15 FPS for better performance)
                await asyncio.sleep(0.067)
                
            except Exception as e:
                logger.error(f"BROADCAST ERROR: {e}")
                import traceback
                logger.error(f"BROADCAST TRACEBACK: {traceback.format_exc()}")
                # Stop broadcasting if simulation is gone
                if simulation_id not in self.simulations:
                    logger.info(f"Simulation {simulation_id} removed, stopping broadcast due to error")
                    break
                # Wait a bit before retrying
                await asyncio.sleep(1)
    
    async def run_simulation_loop(self, simulation_id: str):
        """Run simulation loop for a specific simulation"""
        simulation = self.simulations[simulation_id]
        
        while simulation.is_running:
            if not simulation.is_paused:
                simulation.step()
            
            await asyncio.sleep(0.02)  # 50 FPS simulation for better performance
    
    def start_simulation_loop(self, simulation_id: str):
        """Start simulation loop"""
        if simulation_id in self.simulations:
            if simulation_id not in self.simulation_tasks or self.simulation_tasks[simulation_id].done():
                self.simulation_tasks[simulation_id] = asyncio.create_task(self.run_simulation_loop(simulation_id))
    
    def cleanup_simulation(self, simulation_id: str):
        """Clean up simulation resources"""
        if simulation_id in self.simulations:
            del self.simulations[simulation_id]
        
        if simulation_id in self.active_connections:
            del self.active_connections[simulation_id]
        
        if simulation_id in self.broadcast_tasks:
            self.broadcast_tasks[simulation_id].cancel()
            del self.broadcast_tasks[simulation_id]
        if simulation_id in self.simulation_tasks:
            self.simulation_tasks[simulation_id].cancel()
            del self.simulation_tasks[simulation_id]

    def stop_all_simulations(self):
        """Stop and cleanup all running simulations (single-sim policy)"""
        for sim_id in list(self.simulations.keys()):
            try:
                self.simulations[sim_id].stop()
            except Exception:
                pass
            if sim_id in self.broadcast_tasks:
                try:
                    self.broadcast_tasks[sim_id].cancel()
                except Exception:
                    pass
                self.broadcast_tasks.pop(sim_id, None)
            if sim_id in self.simulation_tasks:
                try:
                    self.simulation_tasks[sim_id].cancel()
                except Exception:
                    pass
                self.simulation_tasks.pop(sim_id, None)
            self.active_connections.pop(sim_id, None)
            self.simulations.pop(sim_id, None)
        
        # Clear any remaining broadcast tasks
        for sim_id in list(self.broadcast_tasks.keys()):
            try:
                self.broadcast_tasks[sim_id].cancel()
            except Exception:
                pass
            self.broadcast_tasks.pop(sim_id, None)
    
    def cleanup_inactive_simulations(self):
        """Clean up simulations that have no active connections"""
        inactive_sims = []
        for sim_id, connections in self.active_connections.items():
            if not connections:  # No active connections
                inactive_sims.append(sim_id)
        
        for sim_id in inactive_sims:
            logger.info(f"Cleaning up inactive simulation {sim_id}")
            # Stop broadcast task
            if sim_id in self.broadcast_tasks:
                task = self.broadcast_tasks[sim_id]
                if not task.done():
                    task.cancel()
                self.broadcast_tasks.pop(sim_id, None)
            
            # Remove from active connections
            self.active_connections.pop(sim_id, None)
            
            # Remove simulation
            self.simulations.pop(sim_id, None)

# Global server instance
server = Live2Server()
app = server.app

if __name__ == "__main__":
    import asyncio
    import sys
    
    # Fix for Windows asyncio issues
    if sys.platform == "win32":
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
    
    uvicorn.run(app, host="localhost", port=8000)
