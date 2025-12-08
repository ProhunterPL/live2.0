"""
FastAPI server for Live 2.0 simulation
Provides WebSocket streaming and REST API endpoints
"""

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*'api.server' found in sys.modules.*")

import asyncio
import json
import time
import logging
from typing import Dict, List, Optional, Any
from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response
from pydantic import BaseModel
import uvicorn
import numpy as np
import msgpack
import taichi as ti
import multiprocessing

# Setup logging to file with rotation (max 5MB)
from logging.handlers import RotatingFileHandler
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Create logs directory if it doesn't exist
project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
logs_dir = os.path.join(project_root, 'logs')
os.makedirs(logs_dir, exist_ok=True)

# Create rotating file handler with max 5MB per file, keep 3 backup files
log_file_path = os.path.join(logs_dir, 'logs.txt')
file_handler = RotatingFileHandler(
    log_file_path, 
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

# Initialize Taichi: Use CPU for best performance
# Benchmark results show CPU is 728x faster for chemistry operations (Bonds/Clusters)
# GPU: 11659ms, CPU: 16ms for chemistry - CPU wins decisively!
#
# To use GPU instead (not recommended based on benchmarks):
# - Set environment variable: LIVE2_USE_GPU=1
# - Or modify this code to use ti.cuda
import multiprocessing
import os

# Check if GPU mode is explicitly requested
use_gpu = os.environ.get('LIVE2_USE_GPU', '0') == '1'

if use_gpu:
    logger.info("GPU mode requested via LIVE2_USE_GPU environment variable")
    try:
        ti.init(arch=ti.cuda, device_memory_GB=4.0)
        logger.info("Taichi initialized with GPU (CUDA) backend")
        logger.warning("Warning: GPU is slower than CPU for chemistry operations (see benchmarks)")
    except Exception as e:
        logger.error(f"GPU initialization failed: {e}")
        logger.info("Falling back to CPU backend")
        use_gpu = False

if not use_gpu:
    try:
        # Use all CPU cores for maximum performance
        num_threads = multiprocessing.cpu_count()
        ti.init(arch=ti.cpu, cpu_max_num_threads=num_threads)
        logger.info(f"Taichi initialized with CPU backend ({num_threads} threads)")
        logger.info("CPU mode: Optimal for chemistry-heavy workloads (728x faster than GPU)")
    except Exception as e:
        # Fallback to default CPU if thread count fails
        logger.warning(f"Failed to initialize with thread count: {e}")
        ti.init(arch=ti.cpu)
        logger.info("Taichi initialized with CPU backend (default threads)")

from sim.config import SimulationConfig, PresetPrebioticConfig, OpenChemistryConfig
from sim.core.stepper import SimulationStepper
from sim.core.hybrid_stepper import HybridSimulationStepper
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

class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Middleware to add security and cache headers"""
    
    async def dispatch(self, request: Request, call_next):
        response = await call_next(request)
        
        # Add security headers
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        # X-XSS-Protection is deprecated and removed (modern browsers have better protection)
        
        # Ensure Content-Type has charset=utf-8 for text/JSON responses
        content_type = response.headers.get("content-type", "")
        if content_type and ("application/json" in content_type or "text/" in content_type) and "charset=" not in content_type:
            response.headers["content-type"] = f"{content_type}; charset=utf-8"
        
        # Add Cache-Control headers based on endpoint
        if request.url.path.startswith("/simulation/") and request.url.path.endswith("/metrics"):
            # Metrics endpoint - no cache (real-time data)
            response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        elif request.url.path.startswith("/simulation/") and "novel-substances" in request.url.path:
            # Novel substances - short cache (5 minutes)
            response.headers["Cache-Control"] = "public, max-age=300"
        elif request.url.path.startswith("/simulation/") and request.url.path.endswith("/status"):
            # Status endpoint - no cache (real-time data)
            response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        elif request.url.path.startswith("/simulation/create"):
            # Create endpoint - no cache
            response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
        else:
            # Default - short cache for other endpoints
            response.headers["Cache-Control"] = "public, max-age=60"
        
        return response

class Live2Server:
    """Main server class for Live 2.0 simulation"""
    
    def __init__(self):
        self.app = FastAPI(title="Live 2.0 Simulation API", version="1.0.0")
        self.setup_middleware()
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
    
    def setup_middleware(self):
        """Setup CORS and Security headers middleware"""
        # Add Security headers middleware first
        self.app.add_middleware(SecurityHeadersMiddleware)
        
        # Add CORS middleware
        self.app.add_middleware(
            CORSMiddleware,
            allow_origins=["*"],  # In production, specify actual origins
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )
    
    def setup_routes(self):
        """Setup API routes"""
        
        # Mount API v1 (Synthetic Data as a Service)
        try:
            from backend.api.v1.main import app as v1_app
            self.app.mount("/api/v1", v1_app)
            logger.info("Mounted API v1 at /api/v1")
        except Exception as e:
            logger.warning(f"Failed to mount API v1: {e}")
        
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
                    # Performance optimizations for production
                    if not hasattr(config, 'chemistry_snapshot_interval') or config.chemistry_snapshot_interval is None:
                        config.chemistry_snapshot_interval = 200  # Less frequent snapshots = better performance
                    if not hasattr(config, 'metrics_update_interval') or config.metrics_update_interval is None:
                        config.metrics_update_interval = 1000  # Less frequent metrics = better performance
                else:
                    raise ValueError(f"Unknown mode: {request.mode}")
                
                # Create simulation - use HybridSimulationStepper for better performance
                # Hybrid mode runs chemistry in background thread, doesn't block simulation
                try:
                    # Try Hybrid mode first (better performance)
                    simulation = HybridSimulationStepper(config)
                    logger.info(f"Using HybridSimulationStepper for {simulation_id} (GPU physics + CPU chemistry)")
                except Exception as e:
                    logger.warning(f"HybridSimulationStepper failed, falling back to SimulationStepper: {e}")
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
        
        @self.app.get("/simulation/{simulation_id}/substance/{substance_id}/details")
        async def get_substance_details(simulation_id: str, substance_id: str):
            """Get detailed topology data for a specific substance/cluster"""
            logger.info(f"Getting details for substance {substance_id} in simulation {simulation_id}")
            
            if simulation_id not in self.simulations:
                logger.error(f"Simulation {simulation_id} not found")
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            
            try:
                if hasattr(simulation, 'catalog'):
                    # Get substance from catalog
                    substance = None
                    if hasattr(simulation.catalog, 'substances'):
                        for sub in simulation.catalog.substances.values():
                            if str(sub.id) == substance_id:
                                substance = sub
                                break
                    
                    if substance is None:
                        logger.error(f"Substance {substance_id} not found in catalog")
                        raise HTTPException(status_code=404, detail="Substance not found")
                    
                    # Extract graph topology
                    graph = substance.graph
                    
                    # Prepare nodes data
                    nodes = []
                    for particle_idx in graph.particles:
                        attrs = graph.particle_attributes.get(particle_idx, np.array([1.0, 0.0, 0.0, 0.0]))
                        if isinstance(attrs, np.ndarray):
                            attrs = attrs.tolist()
                        
                        nodes.append({
                            "id": int(particle_idx),
                            "label": "A",  # Generic label, can be enhanced later
                            "mass": float(attrs[0]) if len(attrs) > 0 else 1.0,
                            "charge": [float(attrs[1]), float(attrs[2]), float(attrs[3])] if len(attrs) >= 4 else [0.0, 0.0, 0.0],
                            "energy": 0.0  # Can be enhanced with actual energy data
                        })
                    
                    # Prepare bonds data
                    bonds = []
                    for i, j in graph.bonds:
                        bonds.append({
                            "a": int(i),
                            "b": int(j),
                            "order": 1  # Default bond order, can be enhanced
                        })
                    
                    # Prepare metadata
                    properties = substance.properties if hasattr(substance, 'properties') else {}
                    metadata = {
                        "size": int(graph.get_node_count()),
                        "bonds": int(graph.get_edge_count()),
                        "density": float(graph.get_density()),
                        "avg_mass": float(properties.get('avg_mass', 1.0)) if 'avg_mass' in properties else 1.0,
                        "complexity": float(graph.get_complexity()),
                        "timestamp": float(substance.last_seen) if hasattr(substance, 'last_seen') else 0.0
                    }
                    
                    return {
                        "id": substance_id,
                        "nodes": nodes,
                        "bonds": bonds,
                        "metadata": metadata
                    }
                else:
                    logger.error(f"Simulation does not have catalog")
                    raise HTTPException(status_code=500, detail="Simulation catalog not available")
                    
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"Failed to get substance details: {e}", exc_info=True)
                raise HTTPException(status_code=500, detail=f"Failed to get substance details: {str(e)}")
        
        @self.app.post("/simulation/{simulation_id}/substance/{substance_id}/match")
        async def match_substance_to_pubchem(simulation_id: str, substance_id: str):
            """
            Match a substance to PubChem and generate comparison files.
            
            This endpoint:
            1. Gets substance topology
            2. Converts to SMILES
            3. Queries PubChem for similar molecules
            4. Generates comparison panel and chemical format files
            5. Returns match results
            """
            logger.info(f"Matching substance {substance_id} from simulation {simulation_id} to PubChem")
            
            if simulation_id not in self.simulations:
                raise HTTPException(status_code=404, detail="Simulation not found")
            
            simulation = self.simulations[simulation_id]
            
            try:
                # Get substance details (reuse existing logic)
                if not hasattr(simulation, 'catalog'):
                    raise HTTPException(status_code=500, detail="Simulation catalog not available")
                
                # Find substance
                substance = None
                if hasattr(simulation.catalog, 'substances'):
                    for sub in simulation.catalog.substances.values():
                        if str(sub.id) == substance_id:
                            substance = sub
                            break
                
                if substance is None:
                    raise HTTPException(status_code=404, detail="Substance not found")
                
                # Import matcher functions
                import sys
                from pathlib import Path
                matcher_path = Path(__file__).parent.parent.parent / "matcher"
                if str(matcher_path) not in sys.path:
                    sys.path.insert(0, str(matcher_path))
                
                from chem import json_to_mol, mol_to_smiles, pubchem_similar_top, export_all_formats, render_mol_png  # type: ignore
                from compose import compose_panel_with_metadata  # type: ignore
                
                # Prepare cluster data
                graph = substance.graph
                nodes = []
                for particle_idx in graph.particles:
                    attrs = graph.particle_attributes.get(particle_idx, np.array([1.0, 0.0, 0.0, 0.0]))
                    if isinstance(attrs, np.ndarray):
                        attrs = attrs.tolist()
                    
                    nodes.append({
                        "id": int(particle_idx),
                        "label": "A",
                        "mass": float(attrs[0]) if len(attrs) > 0 else 1.0,
                        "charge": [float(attrs[1]), float(attrs[2]), float(attrs[3])] if len(attrs) >= 4 else [0.0, 0.0, 0.0],
                        "energy": 0.0
                    })
                
                bonds = [{"a": int(i), "b": int(j), "order": 1} for i, j in graph.bonds]
                
                cluster_data = {
                    "id": substance_id,
                    "nodes": nodes,
                    "bonds": bonds,
                    "metadata": {
                        "size": int(graph.get_node_count()),
                        "bonds": int(graph.get_edge_count()),
                        "density": float(graph.get_density()),
                        "complexity": float(graph.get_complexity())
                    }
                }
                
                # Convert to RDKit Mol
                mol = json_to_mol(cluster_data)
                smiles = mol_to_smiles(mol)
                logger.info(f"Generated SMILES: {smiles}")
                
                # Create matches directory
                matches_dir = Path(__file__).parent.parent.parent / "matches"
                matches_dir.mkdir(exist_ok=True)
                
                # Generate timestamp
                from datetime import datetime
                timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                base_filename = f"cluster_{timestamp}"
                
                # Export to chemical formats
                base_path = matches_dir / base_filename
                export_all_formats(mol, str(base_path), stamp=substance_id)
                
                # Render LIVE cluster
                live_png = matches_dir / f"{base_filename}_live.png"
                render_mol_png(mol, str(live_png), size=512)
                
                # Query PubChem
                logger.info("Querying PubChem...")
                pubchem_match = pubchem_similar_top(smiles, threshold=80)
                
                pubchem_png = None
                if pubchem_match and pubchem_match.get("smiles"):
                    logger.info(f"Found PubChem match: CID {pubchem_match['cid']}")
                    from rdkit import Chem
                    from chem import render_pubchem_png  # type: ignore
                    
                    pubchem_mol = Chem.MolFromSmiles(pubchem_match["smiles"])
                    if pubchem_mol:
                        # Export PubChem to formats
                        pubchem_base = matches_dir / f"{base_filename}_pubchem"
                        export_all_formats(pubchem_mol, str(pubchem_base), stamp=f"PubChem_CID_{pubchem_match['cid']}")
                        
                        # Render PubChem molecule
                        pubchem_png = matches_dir / f"{base_filename}_pubchem.png"
                        render_pubchem_png(pubchem_match["smiles"], str(pubchem_png), size=512)
                else:
                    logger.info("No PubChem match found")
                    pubchem_png = live_png  # Fallback
                
                # Compose comparison panel
                panel_path = matches_dir / f"{base_filename}_match.png"
                if pubchem_match:
                    compose_panel_with_metadata(
                        str(live_png),
                        str(pubchem_png),
                        cluster_data["metadata"],
                        pubchem_match,
                        str(panel_path)
                    )
                else:
                    from compose import compose_panel  # type: ignore
                    compose_panel(
                        str(live_png),
                        str(live_png),
                        "LIVE 2.0 Cluster",
                        "No PubChem Match",
                        str(panel_path),
                        f"Size: {cluster_data['metadata']['size']}",
                        "No similar molecule found"
                    )
                
                # Prepare result
                result = {
                    "success": True,
                    "cluster_id": substance_id,
                    "smiles": smiles,
                    "timestamp": timestamp,
                    "files": {
                        "panel": f"matches/{base_filename}_match.png",
                        "live_png": f"matches/{base_filename}_live.png",
                        "mol": f"matches/{base_filename}.mol",
                        "xyz": f"matches/{base_filename}.xyz"
                    },
                    "pubchem_match": pubchem_match
                }
                
                if pubchem_match:
                    result["files"]["pubchem_png"] = f"matches/{base_filename}_pubchem.png"
                    result["files"]["pubchem_mol"] = f"matches/{base_filename}_pubchem.mol"
                    result["files"]["pubchem_xyz"] = f"matches/{base_filename}_pubchem.xyz"
                
                logger.info(f"Match completed successfully: {base_filename}")
                return result
                
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"Failed to match substance: {e}", exc_info=True)
                raise HTTPException(status_code=500, detail=f"Failed to match substance: {str(e)}")
        
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
                await asyncio.sleep(0.2)  # 5 FPS broadcast - reduced from 10 FPS for better performance

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
                
                # Wait before next broadcast
                # OPTIMIZATION: 15 FPS broadcasting (was 10 FPS / 0.1s) for smoother visualization
                await asyncio.sleep(0.067)  # ~15 FPS broadcasting
                
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
                    # Log every 100 steps to track progress
                    if simulation.step_count % 100 == 0:
                        logger.info(f"[API] Simulation {simulation_id}: calling step() at step {simulation.step_count}")
                    simulation.step()
                    step_count += 1
                    # Log every 100 steps to avoid spam
                    if simulation.step_count % 100 == 0:
                        logger.info(f"[API] Simulation {simulation_id}: step {simulation.step_count} completed, time {simulation.current_time:.3f}")
                except Exception as e:
                    logger.error(f"[API] Error in simulation step for {simulation_id}: {e}")
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
