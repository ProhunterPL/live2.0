"""
Snapshot management for Live 2.0 simulation
Handles saving and loading simulation states
"""

import json
import time
import os
from typing import Dict, List, Optional, Any
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import taichi as ti

class SnapshotManager:
    """Manages simulation snapshots"""
    
    def __init__(self, snapshot_dir: str = "snapshots"):
        self.snapshot_dir = Path(snapshot_dir)
        self.snapshot_dir.mkdir(exist_ok=True)
        
        # Snapshot metadata
        self.metadata_file = self.snapshot_dir / "metadata.json"
        self.metadata = self.load_metadata()
    
    def load_metadata(self) -> Dict:
        """Load snapshot metadata"""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {"snapshots": {}}
    
    def save_metadata(self):
        """Save snapshot metadata"""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
    
    def create_snapshot(self, simulation_data: Dict, name: str = None, save_images: bool = True) -> str:
        """Create a new snapshot"""
        if name is None:
            name = f"snapshot_{int(time.time())}"
        
        filename = f"{name}.json"
        filepath = self.snapshot_dir / filename
        
        # Add metadata
        snapshot_info = {
            "name": name,
            "filename": filename,
            "created_at": time.time(),
            "created_at_str": time.strftime("%Y-%m-%d %H:%M:%S"),
            "simulation_mode": simulation_data.get("config", {}).get("mode", "unknown"),
            "particle_count": simulation_data.get("particles", {}).get("positions", []).__len__(),
            "current_time": simulation_data.get("current_time", 0.0),
            "step_count": simulation_data.get("step_count", 0),
            "has_images": save_images
        }
        
        # Save snapshot data
        with open(filepath, 'w') as f:
            json.dump(simulation_data, f, indent=2)
        
        # Generate and save visualization images
        if save_images:
            self._generate_visualizations(simulation_data, name)
            snapshot_info["image_files"] = [
                f"{name}_overview.png",
                f"{name}_energy_field.png",
                f"{name}_particles.png"
            ]
        
        # Update metadata
        self.metadata["snapshots"][name] = snapshot_info
        self.save_metadata()
        
        return filename
    
    def load_snapshot(self, name: str) -> Dict:
        """Load a snapshot"""
        if name not in self.metadata["snapshots"]:
            raise ValueError(f"Snapshot '{name}' not found")
        
        snapshot_info = self.metadata["snapshots"][name]
        filename = snapshot_info["filename"]
        filepath = self.snapshot_dir / filename
        
        if not filepath.exists():
            raise FileNotFoundError(f"Snapshot file '{filename}' not found")
        
        with open(filepath, 'r') as f:
            return json.load(f)
    
    def list_snapshots(self) -> List[Dict]:
        """List all available snapshots"""
        snapshots = []
        for name, info in self.metadata["snapshots"].items():
            snapshots.append({
                "name": name,
                "filename": info["filename"],
                "created_at": info["created_at"],
                "created_at_str": info["created_at_str"],
                "simulation_mode": info["simulation_mode"],
                "particle_count": info["particle_count"],
                "current_time": info["current_time"],
                "step_count": info["step_count"]
            })
        
        # Sort by creation time (newest first)
        snapshots.sort(key=lambda x: x["created_at"], reverse=True)
        return snapshots
    
    def delete_snapshot(self, name: str) -> bool:
        """Delete a snapshot"""
        if name not in self.metadata["snapshots"]:
            return False
        
        snapshot_info = self.metadata["snapshots"][name]
        filename = snapshot_info["filename"]
        filepath = self.snapshot_dir / filename
        
        # Delete file
        if filepath.exists():
            filepath.unlink()
        
        # Remove from metadata
        del self.metadata["snapshots"][name]
        self.save_metadata()
        
        return True
    
    def get_snapshot_info(self, name: str) -> Optional[Dict]:
        """Get snapshot information"""
        return self.metadata["snapshots"].get(name)
    
    def cleanup_old_snapshots(self, max_age_days: int = 30):
        """Clean up old snapshots"""
        current_time = time.time()
        max_age_seconds = max_age_days * 24 * 60 * 60
        
        to_delete = []
        for name, info in self.metadata["snapshots"].items():
            if current_time - info["created_at"] > max_age_seconds:
                to_delete.append(name)
        
        for name in to_delete:
            self.delete_snapshot(name)
        
        return len(to_delete)
    
    def _generate_visualizations(self, simulation_data: Dict, name: str):
        """Generate visualization images for snapshot"""
        try:
            # Create overview visualization
            self._create_overview_plot(simulation_data, name)
            
            # Create energy field visualization
            if 'energy_field' in simulation_data:
                self._create_energy_field_plot(simulation_data, name)
            
            # Create particles visualization
            if 'particles' in simulation_data:
                self._create_particles_plot(simulation_data, name)
                
        except Exception as e:
            print(f"Warning: Failed to generate visualizations: {e}")
            import traceback
            traceback.print_exc()
    
    def _create_overview_plot(self, simulation_data: Dict, name: str):
        """Create overview plot showing simulation metrics"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(f'Simulation Overview - {name}', fontsize=14, fontweight='bold')
        
        # Time and step information
        ax = axes[0, 0]
        current_time = simulation_data.get('current_time', 0)
        step_count = simulation_data.get('step_count', 0)
        ax.text(0.1, 0.7, f'Current Time: {current_time:.2f}', transform=ax.transAxes, fontsize=12)
        ax.text(0.1, 0.5, f'Step Count: {step_count:,}', transform=ax.transAxes, fontsize=12)
        ax.text(0.1, 0.3, f'Mode: {simulation_data.get("mode", "unknown")}', transform=ax.transAxes, fontsize=12)
        ax.set_title('Simulation Status')
        ax.axis('off')
        
        # Particle count and metrics
        ax = axes[0, 1]
        particles = simulation_data.get('particles', {})
        particle_count = len(particles.get('positions', []))
        bonds = simulation_data.get('bonds', [])
        bond_count = len(bonds)
        
        ax.text(0.1, 0.7, f'Particles: {particle_count}', transform=ax.transAxes, fontsize=12)
        ax.text(0.1, 0.5, f'Bonds: {bond_count}', transform=ax.transAxes, fontsize=12)
        ax.text(0.1, 0.3, f'Clusters: {len(simulation_data.get("clusters", []))}', transform=ax.transAxes, fontsize=12)
        ax.set_title('System Composition')
        ax.axis('off')
        
        # Metrics visualization
        metrics = simulation_data.get('metrics', {})
        ax = axes[1, 0]
        metric_names = list(metrics.keys())
        metric_values = list(metrics.values())
        if metric_names:
            bars = ax.bar(range(len(metric_names)), metric_values)
            ax.set_xticks(range(len(metric_names)))
            ax.set_xticklabels(metric_names, rotation=45, ha='right')
            ax.set_title('Key Metrics')
            ax.tick_params(axis='x', labelsize=8)
        
        # Energy distribution
        ax = axes[1, 1]
        energy_field = simulation_data.get('energy_field')
        if energy_field:
            energy_array = np.array(energy_field)
            ax.hist(energy_array.flatten(), bins=50, alpha=0.7, color='blue', edgecolor='black')
            ax.set_xlabel('Energy Value')
            ax.set_ylabel('Frequency')
            ax.set_title('Energy Distribution')
        
        plt.tight_layout()
        output_file = self.snapshot_dir / f'{name}_overview.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
    
    def _create_energy_field_plot(self, simulation_data: Dict, name: str):
        """Create energy field heatmap"""
        energy_field = simulation_data.get('energy_field')
        if not energy_field:
            return
            
        energy_array = np.array(energy_field)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        im = ax.imshow(energy_array, cmap='viridis', aspect='auto')
        ax.set_title(f'Energy Field - {name}', fontsize=14, fontweight='bold')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Energy', rotation=270, labelpad=20)
        
        plt.tight_layout()
        output_file = self.snapshot_dir / f'{name}_energy_field.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
    
    def _create_particles_plot(self, simulation_data: Dict, name: str):
        """Create particles and bonds visualization"""
        particles = simulation_data.get('particles', {})
        bonds = simulation_data.get('bonds', [])
        
        positions = particles.get('positions', [])
        if not positions:
            return
            
        positions_array = np.array(positions)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot particles
        if len(positions_array) > 0:
            x_coords = positions_array[:, 0]
            y_coords = positions_array[:, 1]
            ax.scatter(x_coords, y_coords, c='red', s=20, alpha=0.7, label='Particles')
        
        # Plot bonds
        if bonds and len(positions_array) > 0:
            for bond in bonds:
                if len(bond) >= 2:
                    i, j = int(bond[0]), int(bond[1])
                    if i < len(positions_array) and j < len(positions_array):
                        x1, y1 = positions_array[i]
                        x2, y2 = positions_array[j]
                        strength = bond[2] if len(bond) > 2 else 1.0
                        ax.plot([x1, x2], [y1, y2], 'b-', alpha=strength, linewidth=1)
        
        ax.set_title(f'Particles and Bonds - {name}', fontsize=14, fontweight='bold')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = self.snapshot_dir / f'{name}_particles.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()

class SnapshotSerializer:
    """Serializes and deserializes simulation data"""
    
    @staticmethod
    def serialize_simulation(simulation) -> Dict:
        """Serialize simulation state"""
        # Get basic state
        state = simulation.get_simulation_state()
        
        # Get visualization data
        viz_data = simulation.get_visualization_data()
        
        # Get novel substances
        novel_substances = simulation.get_novel_substances(50)
        
        # Get metrics
        metrics = simulation.aggregator.get_aggregated_stats()
        
        # Get catalog data
        catalog_stats = simulation.catalog.get_catalog_stats()
        
        # Combine all data
        snapshot_data = {
            "config": simulation.config.dict(),
            "current_time": state["current_time"],
            "step_count": state["step_count"],
            "is_running": state["is_running"],
            "is_paused": state["is_paused"],
            "mode": state["mode"],
            "particles": viz_data["particles"],
            "energy_field": viz_data["energy_field"],
            "bonds": viz_data["bonds"],
            "clusters": viz_data["clusters"],
            "novel_substances": novel_substances,
            "metrics": metrics,
            "catalog_stats": catalog_stats,
            "timestamp": time.time()
        }
        
        return snapshot_data
    
    @staticmethod
    def deserialize_simulation(snapshot_data: Dict, simulation):
        """Deserialize simulation state"""
        # Load configuration
        config_data = snapshot_data["config"]
        if snapshot_data["mode"] == "preset_prebiotic":
            from sim.config import PresetPrebioticConfig
            config = PresetPrebioticConfig(**config_data)
        else:
            from sim.config import SimulationConfig
            config = SimulationConfig(**config_data)
        
        # Update simulation configuration
        simulation.config = config
        
        # Load simulation state
        simulation.current_time = snapshot_data["current_time"]
        simulation.step_count = snapshot_data["step_count"]
        
        # Load particles
        particles_data = snapshot_data["particles"]
        simulation.particles.reset()
        
        # Add particles back
        positions = np.array(particles_data["positions"])
        attributes = np.array(particles_data["attributes"])
        active_mask = np.array(particles_data["active_mask"])
        
        for i, (pos, attr, active) in enumerate(zip(positions, attributes, active_mask)):
            if active:
                pos_ti = ti.Vector(pos)
                vel_ti = ti.Vector([0.0, 0.0])  # Reset velocity
                attr_ti = ti.Vector(attr)
                
                # Register particle type
                type_id = simulation.particles.register_particle_type(
                    name=f"restored_particle_{i}",
                    mass=attr[0],
                    charge=(attr[1], attr[2], attr[3])
                )
                
                # Add particle
                simulation.particles.add_particle_py(pos_ti, vel_ti, attr_ti, type_id, 2, 1.0)
        
        # Load bonds
        bonds = snapshot_data["bonds"]
        simulation.binding.reset()
        
        for bond in bonds:
            i, j, strength = bond
            if i < simulation.config.max_particles and j < simulation.config.max_particles:
                simulation.binding.form_bond(i, j)
        
        # Load energy field
        energy_field = np.array(snapshot_data["energy_field"])
        simulation.energy_manager.energy_system.energy_field.from_numpy(energy_field)
        
        # Load catalog
        simulation.catalog.clear()
        # Note: Full catalog restoration would require more complex serialization
        
        # Update metrics
        simulation.update_metrics()
    
    @staticmethod
    def validate_snapshot(snapshot_data: Dict) -> bool:
        """Validate snapshot data"""
        required_fields = [
            "config", "current_time", "step_count", "mode",
            "particles", "energy_field", "bonds", "clusters"
        ]
        
        for field in required_fields:
            if field not in snapshot_data:
                return False
        
        # Validate particle data
        particles = snapshot_data["particles"]
        if not all(key in particles for key in ["positions", "attributes", "active_mask"]):
            return False
        
        # Validate data consistency
        positions = particles["positions"]
        attributes = particles["attributes"]
        active_mask = particles["active_mask"]
        
        if len(positions) != len(attributes) or len(positions) != len(active_mask):
            return False
        
        return True

class SnapshotAPI:
    """API for snapshot management"""
    
    def __init__(self, snapshot_manager: SnapshotManager):
        self.snapshot_manager = snapshot_manager
    
    def create_snapshot(self, simulation, name: str = None) -> str:
        """Create snapshot from simulation"""
        snapshot_data = SnapshotSerializer.serialize_simulation(simulation)
        return self.snapshot_manager.create_snapshot(snapshot_data, name)
    
    def load_snapshot(self, simulation, name: str):
        """Load snapshot into simulation"""
        snapshot_data = self.snapshot_manager.load_snapshot(name)
        
        if not SnapshotSerializer.validate_snapshot(snapshot_data):
            raise ValueError("Invalid snapshot data")
        
        SnapshotSerializer.deserialize_simulation(snapshot_data, simulation)
    
    def list_snapshots(self) -> List[Dict]:
        """List available snapshots"""
        return self.snapshot_manager.list_snapshots()
    
    def delete_snapshot(self, name: str) -> bool:
        """Delete snapshot"""
        return self.snapshot_manager.delete_snapshot(name)
    
    def get_snapshot_info(self, name: str) -> Optional[Dict]:
        """Get snapshot information"""
        return self.snapshot_manager.get_snapshot_info(name)
    
    def cleanup_old_snapshots(self, max_age_days: int = 30) -> int:
        """Clean up old snapshots"""
        return self.snapshot_manager.cleanup_old_snapshots(max_age_days)
