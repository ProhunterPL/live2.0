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

# Setup logging to file
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs.txt'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Initialize Taichi first
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
        
        # WebSocket broadcasting
        self.broadcast_tasks: Dict[str, asyncio.Task] = {}
        
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
        
        @self.app.websocket("/simulation/{simulation_id}/stream")
        async def websocket_endpoint(websocket: WebSocket, simulation_id: str):
            """WebSocket endpoint for real-time data streaming"""
            logger.info(f"WebSocket connection attempt for simulation {simulation_id}")
            await websocket.accept()
            logger.info(f"WebSocket accepted for simulation {simulation_id}")
            
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                await websocket.close(code=1008, reason="Simulation not found")
                return
            
            # Add connection to active connections
            self.active_connections[simulation_id].append(websocket)
            logger.info(f"Added WebSocket connection for simulation {simulation_id}")
            
            # Ensure simulation is running when a client connects
            simulation = self.simulations[simulation_id]
            if not simulation.is_running:
                logger.info(f"Starting simulation {simulation_id}")
                simulation.start()
            
            # Start broadcasting if not already started
            if simulation_id not in self.broadcast_tasks:
                logger.info(f"Starting broadcast task for simulation {simulation_id}")
                self.broadcast_tasks[simulation_id] = asyncio.create_task(
                    self.broadcast_simulation_data(simulation_id)
                )
            
            try:
                while True:
                    # Keep connection alive
                    await websocket.receive_text()
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
        simulation = self.simulations[simulation_id]
        
        while True:
            try:
                # Advance simulation to ensure time progresses while streaming
                logger.info(f"BROADCAST: About to step simulation {simulation_id} - current_time={simulation.current_time:.6f}, step_count={simulation.step_count}")
                simulation.step()
                logger.info(f"BROADCAST: Step completed for simulation {simulation_id} - current_time={simulation.current_time:.6f}, step_count={simulation.step_count}")
                
                # Get visualization data
                logger.debug(f"Getting visualization data for simulation {simulation_id}")
                data = simulation.get_visualization_data()
                
                # Serialize data
                if len(data['particles']['positions']) > 0:
                    # Use msgpack for binary data
                    binary_data = msgpack.packb(data)
                    
                    # Send to all connected clients
                    disconnected = []
                    for websocket in self.active_connections[simulation_id]:
                        try:
                            await websocket.send_bytes(binary_data)
                        except:
                            disconnected.append(websocket)
                    
                    # Remove disconnected clients
                    for websocket in disconnected:
                        self.active_connections[simulation_id].remove(websocket)
                
                # Wait before next broadcast (about 30 FPS)
                await asyncio.sleep(0.033)
                
            except Exception as e:
                logger.error(f"BROADCAST ERROR: {e}")
                import traceback
                logger.error(f"BROADCAST TRACEBACK: {traceback.format_exc()}")
                # continue broadcasting despite errors
    
    async def run_simulation_loop(self, simulation_id: str):
        """Run simulation loop for a specific simulation"""
        simulation = self.simulations[simulation_id]
        
        while simulation.is_running:
            if not simulation.is_paused:
                simulation.step()
            
            await asyncio.sleep(0.01)  # 100 FPS simulation
    
    def start_simulation_loop(self, simulation_id: str):
        """Start simulation loop"""
        if simulation_id in self.simulations:
            asyncio.create_task(self.run_simulation_loop(simulation_id))
    
    def cleanup_simulation(self, simulation_id: str):
        """Clean up simulation resources"""
        if simulation_id in self.simulations:
            del self.simulations[simulation_id]
        
        if simulation_id in self.active_connections:
            del self.active_connections[simulation_id]
        
        if simulation_id in self.broadcast_tasks:
            self.broadcast_tasks[simulation_id].cancel()
            del self.broadcast_tasks[simulation_id]

# Global server instance
server = Live2Server()
app = server.app

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
