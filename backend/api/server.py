"""
FastAPI server for Live 2.0 simulation
Provides WebSocket streaming and REST API endpoints
"""

import asyncio
import json
import time
import logging
from typing import Dict, List, Optional, Any
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException, Request
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
    level=logging.INFO,  # Changed from DEBUG to INFO
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        file_handler,
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Set level for all loggers
logging.getLogger().setLevel(logging.INFO)  # Changed from DEBUG to INFO
logging.getLogger('sim.core.stepper').setLevel(logging.INFO)  # Changed from DEBUG to INFO

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
        # Snapshot management - use absolute path
        import os
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        snapshots_path = os.path.join(project_root, "snapshots")
        self.snapshot_manager = SnapshotManager(snapshot_dir=snapshots_path)
        
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
        
        @self.app.get("/health")
        async def health_check():
            """Health check endpoint with memory usage info"""
            import psutil
            memory = psutil.virtual_memory()
            
            # Get simulation count
            active_simulations = len([s for s in self.simulations.values() if s.is_running])
            
            return {
                "status": "healthy",
                "memory_percent": memory.percent,
                "active_simulations": active_simulations,
                "total_simulations": len(self.simulations)
            }
        
        @self.app.get("/simulations/active")
        async def get_active_simulations():
            """Get list of active simulation IDs"""
            return {
                "simulations": list(self.simulations.keys()),
                "count": len(self.simulations)
            }
        
        @self.app.get("/simulation/current")
        async def get_current_simulation():
            """Get current simulation ID (single simulation policy)"""
            if self.simulations:
                # Return the most recent simulation
                sim_id = max(self.simulations.keys())
                return {
                    "simulation_id": sim_id,
                    "exists": True
                }
            else:
                return {
                    "simulation_id": None,
                    "exists": False
                }
        
        @self.app.post("/simulations/cleanup")
        async def cleanup_simulations():
            """Clean up old simulations"""
            try:
                old_count = len(self.simulations)
                self.cleanup_old_simulations(max_simulations=3)
                new_count = len(self.simulations)
                return {
                    "success": True,
                    "message": f"Cleaned up {old_count - new_count} simulations",
                    "remaining": new_count
                }
            except Exception as e:
                logger.error(f"Failed to cleanup simulations: {e}")
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.get("/simulations/debug")
        async def debug_simulations():
            """Debug simulation tasks"""
            debug_info = {
                "simulations": list(self.simulations.keys()),
                "simulation_tasks": {},
                "broadcast_tasks": {},
                "active_connections": {}
            }
            
            for sim_id in self.simulations.keys():
                # Check simulation tasks
                if sim_id in self.simulation_tasks:
                    task = self.simulation_tasks[sim_id]
                    debug_info["simulation_tasks"][sim_id] = {
                        "exists": True,
                        "done": task.done(),
                        "cancelled": task.cancelled()
                    }
                else:
                    debug_info["simulation_tasks"][sim_id] = {"exists": False}
                
                # Check broadcast tasks
                if sim_id in self.broadcast_tasks:
                    task = self.broadcast_tasks[sim_id]
                    debug_info["broadcast_tasks"][sim_id] = {
                        "exists": True,
                        "done": task.done(),
                        "cancelled": task.cancelled()
                    }
                else:
                    debug_info["broadcast_tasks"][sim_id] = {"exists": False}
                
                # Check active connections
                debug_info["active_connections"][sim_id] = len(self.active_connections.get(sim_id, []))
            
            return debug_info
        
        @self.app.post("/simulation/create", response_model=SimulationResponse)
        async def create_simulation(request: SimulationRequest):
            """Create a new simulation"""
            try:
                # DEBUG: Print received request
                # Debug prints removed for performance
                
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
                
                logger.info(f"Created simulation {simulation_id}")
                
                return SimulationResponse(
                    success=True,
                    message="Simulation created successfully",
                    simulation_id=simulation_id
                )
            except Exception as e:
                logger.error(f"Failed to create simulation: {e}")
                import traceback
                traceback.print_exc()
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
            # Debug print removed for performance
            logger.info(f"START SIMULATION called for {simulation_id}")
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                raise HTTPException(status_code=404, detail="Simulation not found")

            simulation = self.simulations[simulation_id]
            logger.info(f"Calling simulation.start() for {simulation_id}")
            simulation.start()
            logger.info(f"Simulation started, calling start_simulation_loop")
            # Ensure simulation loop is running in background (non-blocking)
            self.start_simulation_loop(simulation_id)
            logger.info(f"start_simulation_loop called for {simulation_id}")
            
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
        async def save_snapshot(simulation_id: str, request: Request):
            """Save simulation snapshot with optional image generation"""
            # Debug print removed for performance
            
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            # Get parameters from query string
            query_params = request.query_params
            filename = query_params.get("filename")
            save_images = query_params.get("save_images", "true").lower() == "true"
            
            # Debug print removed for performance
            
            if filename is None:
                filename = f"snapshot_{simulation_id}_{int(time.time())}.json"
            
            simulation = self.simulations[simulation_id]
            # Debug print removed for performance
            
            try:
                saved_filename = simulation.save_snapshot(filename, save_images=save_images)
                logger.info(f"Snapshot saved: {saved_filename}")
                return {"success": True, "filename": saved_filename, "images_generated": save_images}
            except Exception as e:
                logger.error(f"Snapshot save failed: {e}")
                import traceback
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=f"Failed to save snapshot: {str(e)}")
        
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
        
        @self.app.get("/simulation/{simulation_id}/performance")
        async def get_performance_metrics(simulation_id: str):
            """Get simulation performance metrics"""
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            performance_metrics = simulation.performance_monitor.get_performance_metrics(simulation.step_count)
            warnings = simulation.performance_monitor.check_performance_warnings()
            
            return {
                "performance": performance_metrics,
                "warnings": warnings,
                "status": "ok"
            }
        
        @self.app.get("/simulation/{simulation_id}/novel-substances")
        async def get_novel_substances(simulation_id: str, count: int = 10):
            """Get novel substances discovered in simulation"""
            logger.info(f"Getting novel substances for simulation {simulation_id}, count={count}")
            
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            logger.info(f"Simulation found: {type(simulation)}")
            
            # Get novel substances from simulation state or catalog
            try:
                if hasattr(simulation, 'catalog'):
                    logger.info(f"Simulation has catalog: {type(simulation.catalog)}")
                    if hasattr(simulation.catalog, 'get_recent_substances'):
                        logger.info(f"Catalog has get_recent_substances method")
                        substances = simulation.catalog.get_recent_substances(count)
                        logger.info(f"Retrieved {len(substances)} substances")
                        
                        # Convert to expected format
                        formatted_substances = []
                        for substance in substances:
                            try:
                                # Convert numpy types to Python types for JSON serialization
                                complexity = float(substance.graph.get_complexity())
                                density = float(substance.graph.get_density())
                                
                                # Handle properties safely
                                properties = substance.properties if hasattr(substance, 'properties') else {}
                                avg_mass = float(properties.get('avg_mass', 1.0)) if 'avg_mass' in properties else 1.0
                                avg_charge = properties.get('avg_charge', [0.0, 0.0, 0.0])
                                
                                # Convert charge to Python floats if it's numpy array
                                if hasattr(avg_charge, 'tolist'):
                                    avg_charge = avg_charge.tolist()
                                elif isinstance(avg_charge, np.ndarray):
                                    avg_charge = [float(x) for x in avg_charge]
                                
                                formatted_substances.append({
                                    "id": str(substance.id),
                                    "timestamp": float(substance.last_seen),
                                    "size": int(substance.graph.get_node_count()),
                                    "complexity": complexity,
                                    "properties": {
                                        "mass": avg_mass,
                                        "charge": avg_charge,
                                        "bonds": int(substance.graph.get_edge_count()),
                                        "graph_density": density
                                    }
                                })
                            except Exception as e:
                                logger.warning(f"Error formatting substance {substance.id}: {e}")
                                continue
                        
                        substances = formatted_substances
                        logger.info(f"Formatted {len(substances)} substances")
                    else:
                        logger.warning(f"Catalog does not have get_recent_substances method")
                        substances = []
                else:
                    logger.warning(f"Simulation does not have catalog attribute")
                    substances = []
                
                logger.info(f"Returning {len(substances)} substances")
                return {"substances": substances}
            except Exception as e:
                logger.error(f"Failed to get novel substances for simulation {simulation_id}: {e}", exc_info=True)
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
            else:
                # Check if broadcast task is still running
                task = self.broadcast_tasks[simulation_id]
                if task.done():
                    logger.info(f"Broadcast task finished for simulation {simulation_id}, restarting")
                    self.broadcast_tasks[simulation_id] = asyncio.create_task(
                        self.broadcast_simulation_data(simulation_id)
                    )
                else:
                    logger.debug(f"Broadcast task already exists and running for simulation {simulation_id}")
            
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
                        logger.info(f"Cancelling broadcast task for simulation {simulation_id} - no connections")
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
                logger.debug(f"Broadcast loop iteration for {simulation_id}")
                # Check if simulation still exists
                if simulation_id not in self.simulations:
                    logger.info(f"Simulation {simulation_id} removed, stopping broadcast")
                    break
                # OPTIMIZED: Send data less frequently to avoid overwhelming slow get_visualization_data
                # Wait for simulation to progress (don't spam if simulation is slow)
                await asyncio.sleep(0.1)  # 10 FPS broadcast - FASTER for better responsiveness

                # Get visualization data with timing
                t0 = time.time()
                logger.debug(f"Calling get_visualization_data for {simulation_id} at step {simulation.step_count}")
                
                # Start broadcast timing
                simulation.performance_monitor.start_broadcast_timing()
                
                data = simulation.get_visualization_data()
                t1 = time.time()
                logger.debug(f"get_visualization_data returned in {t1-t0:.2f}s")
                if (t1 - t0) > 0.5:
                    logger.warning(f"get_visualization_data took {t1-t0:.2f}s - TOO SLOW!")

                # OPTIMIZED: Downsample 4x for energy field to reduce bandwidth by 16x
                def _downsample_4x(arr):
                    """Downsample NumPy array by 4x (16x fewer pixels) and convert to list"""
                    try:
                        if not isinstance(arr, np.ndarray):
                            arr = np.asarray(arr)
                        h, w = arr.shape
                        if h >= 4 and w >= 4:
                            # Trim to multiple of 4
                            h_trim = h - (h % 4)
                            w_trim = w - (w % 4)
                            arr = arr[:h_trim, :w_trim]
                            # Downsample by 4x using mean pooling
                            arr = arr.reshape(h_trim//4, 4, w_trim//4, 4).mean(axis=(1, 3))
                        return arr.tolist()
                    except Exception as e:
                        logger.warning(f"Downsample failed: {e}, returning original")
                        return arr.tolist() if isinstance(arr, np.ndarray) else arr

                # Downsample energy field (128x128 -> 32x32 = 16x reduction)
                t2 = time.time()
                if isinstance(data.get('energy_field'), np.ndarray):
                    data['energy_field'] = _downsample_4x(data['energy_field'])
                elif isinstance(data.get('energy_field'), list):
                    data['energy_field'] = _downsample_4x(np.asarray(data['energy_field']))
                t3 = time.time()
                if (t3 - t2) > 0.5:
                    logger.warning(f"Energy field downsample took {t3-t2:.2f}s")
                
                # Downsample concentration view if present
                if isinstance(data.get('concentration_view'), np.ndarray):
                    data['concentration_view'] = _downsample_4x(data['concentration_view'])
                elif isinstance(data.get('concentration_view'), list):
                    data['concentration_view'] = _downsample_4x(np.asarray(data['concentration_view']))
                
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
                t4 = time.time()
                serializable_data = convert_numpy_types(data)
                t5 = time.time()
                if (t5 - t4) > 0.5:
                    logger.warning(f"convert_numpy_types took {t5-t4:.2f}s")
                
                t6 = time.time()
                binary_data = msgpack.packb(serializable_data)
                t7 = time.time()
                if (t7 - t6) > 0.5:
                    logger.warning(f"msgpack.packb took {t7-t6:.2f}s")
                
                # Send to all connected clients
                disconnected = []
                if simulation_id in self.active_connections:
                    client_count = len(self.active_connections[simulation_id])
                    # logger.debug(f"BROADCAST: Sending data to {client_count} clients, step {simulation.step_count}")
                    for websocket in self.active_connections[simulation_id]:
                        try:
                            await websocket.send_bytes(binary_data)
                        except:
                            disconnected.append(websocket)
                    
                    # End broadcast timing
                    broadcast_time = simulation.performance_monitor.end_broadcast_timing()
                    if broadcast_time > 0.05:  # Log if broadcast takes more than 50ms
                        logger.warning(f"Slow broadcast: {broadcast_time*1000:.1f}ms")
                        
                else:
                    pass  # No active connections
                
                # Remove disconnected clients
                if simulation_id in self.active_connections:
                    for websocket in disconnected:
                        if websocket in self.active_connections[simulation_id]:
                            self.active_connections[simulation_id].remove(websocket)
                
                # Wait before next broadcast (about 10 FPS for better performance)
                await asyncio.sleep(0.1)  # 10 FPS broadcasting - CONSISTENT with above
                
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
        logger.info(f"Starting simulation loop for {simulation_id}")
        simulation = self.simulations[simulation_id]
        logger.info(f"Starting simulation loop for {simulation_id}")
        
        step_count = 0
        loop_count = 0
        while simulation.is_running:
            loop_count += 1
            if loop_count <= 10:  # Debug first 10 loops
                logger.info(f"SIMULATION LOOP {loop_count}: is_running={simulation.is_running}, is_paused={simulation.is_paused}")
            
            if not simulation.is_paused:
                try:
                    simulation.step()
                    step_count += 1
                    # Log every 100 steps to avoid spam
                    if simulation.step_count % 100 == 0:
                        logger.info(f"Simulation {simulation_id}: step {simulation.step_count}, time {simulation.current_time:.3f}")
                except Exception as e:
                    logger.error(f"Error in simulation step for {simulation_id}: {e}")
                    import traceback
                    traceback.print_exc()
                    break
            
            await asyncio.sleep(0.01)  # 100 FPS simulation - MUCH FASTER for better performance
        
        logger.info(f"Simulation loop ended for {simulation_id}")
    
    def start_simulation_loop(self, simulation_id: str):
        """Start simulation loop"""
        # Debug print removed for performance
        # Debug print removed for performance
        
        if simulation_id in self.simulations:
            task_exists = simulation_id in self.simulation_tasks
            task_done = self.simulation_tasks[simulation_id].done() if task_exists else False
            # Debug print removed for performance
            
            if simulation_id not in self.simulation_tasks or self.simulation_tasks[simulation_id].done():
                logger.info(f"Creating simulation task for {simulation_id}")
                logger.info(f"Creating new simulation task for {simulation_id}")
                self.simulation_tasks[simulation_id] = asyncio.create_task(self.run_simulation_loop(simulation_id))
                logger.info(f"Simulation task created successfully for {simulation_id}")
            else:
                logger.info(f"Simulation task already running for {simulation_id}")
                logger.info(f"Simulation task already running for {simulation_id}")
        else:
            logger.error(f"Simulation {simulation_id} not found!")
            logger.error(f"Cannot start simulation loop: simulation {simulation_id} not found")
    
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
    
    def cleanup_old_simulations(self, max_simulations: int = 5):
        """Clean up old simulations, keeping only the most recent ones"""
        if len(self.simulations) <= max_simulations:
            return
        
        # Sort simulations by ID (which includes timestamp) and keep only the newest ones
        sorted_sims = sorted(self.simulations.keys(), reverse=True)
        sims_to_remove = sorted_sims[max_simulations:]
        
        for sim_id in sims_to_remove:
            logger.info(f"Cleaning up old simulation {sim_id}")
            try:
                # Stop simulation if running
                if sim_id in self.simulations:
                    self.simulations[sim_id].stop()
                
                # Stop broadcast task
                if sim_id in self.broadcast_tasks:
                    task = self.broadcast_tasks[sim_id]
                    if not task.done():
                        task.cancel()
                    self.broadcast_tasks.pop(sim_id, None)
                
                # Stop simulation task
                if sim_id in self.simulation_tasks:
                    task = self.simulation_tasks[sim_id]
                    if not task.done():
                        task.cancel()
                    self.simulation_tasks.pop(sim_id, None)
                
                # Remove from all dictionaries
                self.simulations.pop(sim_id, None)
                self.active_connections.pop(sim_id, None)
                
            except Exception as e:
                logger.warning(f"Error cleaning up simulation {sim_id}: {e}")

# Global server instance
server = Live2Server()
app = server.app

if __name__ == "__main__":
    import asyncio
    import sys
    
    # Fix for Windows asyncio issues
    if sys.platform == "win32":
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
    
    uvicorn.run(app, host="localhost", port=8001)
