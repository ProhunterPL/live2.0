"""
Energy management for Live 2.0 simulation
Handles energy sources, distribution, and dissipation
"""

import taichi as ti
import numpy as np
from typing import List, Tuple, Dict, Optional
from .config import SimulationConfig

@ti.data_oriented
class EnergySystem:
    """Manages energy distribution and sources in the simulation"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.width = config.grid_width
        self.height = config.grid_height
        
        # Energy field
        self.energy_field = ti.field(dtype=ti.f32, shape=(self.width, self.height))
        
        # Energy sources
        self.max_sources = 10
        self.source_positions = ti.Vector.field(2, dtype=ti.f32, shape=(self.max_sources,))
        self.source_intensities = ti.field(dtype=ti.f32, shape=(self.max_sources,))
        self.source_radii = ti.field(dtype=ti.f32, shape=(self.max_sources,))
        self.source_active = ti.field(dtype=ti.i32, shape=(self.max_sources,))
        self.source_count = ti.field(dtype=ti.i32, shape=())
        
        # Energy parameters
        self.energy_decay_rate = ti.field(dtype=ti.f32, shape=())
        self.energy_diffusion_rate = ti.field(dtype=ti.f32, shape=())
        self.energy_threshold = ti.field(dtype=ti.f32, shape=())
        
        # Initialize parameters
        self.energy_decay_rate[None] = config.energy_decay
        self.energy_diffusion_rate[None] = 0.1
        self.energy_threshold[None] = config.energy_threshold
        
        # Initialize
        self.reset()
    
    @ti.kernel
    def reset(self):
        """Reset energy system"""
        # Clear energy field
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] = 0.0
        
        # Clear energy sources
        for i in range(self.max_sources):
            self.source_active[i] = 0
            self.source_intensities[i] = 0.0
            self.source_radii[i] = 0.0
        
        self.source_count[None] = 0
    
    @ti.kernel
    def add_energy_source(self, pos: ti.template(), intensity: ti.f32, 
                         radius: ti.f32) -> ti.i32:
        """Add an energy source. Returns source index or -1 if failed"""
        if self.source_count[None] >= self.max_sources:
            return -1
        
        idx = self.source_count[None]
        self.source_positions[idx] = pos
        self.source_intensities[idx] = intensity
        self.source_radii[idx] = radius
        self.source_active[idx] = 1
        self.source_count[None] += 1
        
        return idx
    
    @ti.kernel
    def remove_energy_source(self, idx: ti.i32):
        """Remove energy source by index"""
        if idx < 0 or idx >= self.max_sources:
            return
        
        self.source_active[idx] = 0
    
    @ti.kernel
    def update_energy_field(self, dt: ti.f32):
        """Update energy field with sources, diffusion, and decay"""
        # Apply energy sources
        for i in range(self.max_sources):
            if self.source_active[i] == 1:
                source_pos = self.source_positions[i]
                intensity = self.source_intensities[i]
                radius = self.source_radii[i]
                
                # Add energy in circular area around source
                for x, y in ti.ndrange(self.width, self.height):
                    dx = x - source_pos[0]
                    dy = y - source_pos[1]
                    
                    # Handle periodic boundary conditions
                    if dx > self.width / 2:
                        dx -= self.width
                    elif dx < -self.width / 2:
                        dx += self.width
                    
                    if dy > self.height / 2:
                        dy -= self.height
                    elif dy < -self.height / 2:
                        dy += self.height
                    
                    distance = ti.sqrt(dx * dx + dy * dy)
                    
                    if distance <= radius:
                        # Energy falls off with distance
                        energy_factor = 1.0 - (distance / radius)
                        energy_factor = ti.max(energy_factor, 0.0)
                        self.energy_field[x, y] += intensity * energy_factor * dt
        
        # Apply diffusion
        self.apply_diffusion(dt)
        
        # Apply decay
        self.apply_decay(dt)
    
    @ti.kernel
    def apply_diffusion(self, dt: ti.f32):
        """Apply energy diffusion using discrete Laplacian"""
        # Create temporary field for diffusion
        temp_field = ti.field(dtype=ti.f32, shape=(self.width, self.height))
        
        # Copy current field
        for i, j in ti.ndrange(self.width, self.height):
            temp_field[i, j] = self.energy_field[i, j]
        
        # Apply diffusion kernel
        diffusion_rate = self.energy_diffusion_rate[None]
        
        for i, j in ti.ndrange(self.width, self.height):
            # Get neighbors with periodic boundary conditions
            left = (i - 1) % self.width
            right = (i + 1) % self.width
            up = (j - 1) % self.height
            down = (j + 1) % self.height
            
            # Discrete Laplacian
            laplacian = (temp_field[left, j] + temp_field[right, j] + 
                        temp_field[i, up] + temp_field[i, down] - 
                        4.0 * temp_field[i, j])
            
            # Update energy field
            self.energy_field[i, j] += diffusion_rate * laplacian * dt
    
    @ti.kernel
    def apply_decay(self, dt: ti.f32):
        """Apply energy decay"""
        decay_rate = self.energy_decay_rate[None]
        
        for i, j in ti.ndrange(self.width, self.height):
            self.energy_field[i, j] *= decay_rate
    
    @ti.kernel
    def add_energy_impulse(self, pos: ti.template(), intensity: ti.f32, 
                          radius: ti.f32):
        """Add a one-time energy impulse"""
        for x, y in ti.ndrange(self.width, self.height):
            dx = x - pos[0]
            dy = y - pos[1]
            
            # Handle periodic boundary conditions
            if dx > self.width / 2:
                dx -= self.width
            elif dx < -self.width / 2:
                dx += self.width
            
            if dy > self.height / 2:
                dy -= self.height
            elif dy < -self.height / 2:
                dy += self.height
            
            distance = ti.sqrt(dx * dx + dy * dy)
            
            if distance <= radius:
                energy_factor = 1.0 - (distance / radius)
                energy_factor = ti.max(energy_factor, 0.0)
                self.energy_field[x, y] += intensity * energy_factor
    
    @ti.kernel
    def get_energy_at_position(self, pos: ti.template()) -> ti.f32:
        """Get energy at a specific position"""
        x = int(pos[0]) % self.width
        y = int(pos[1]) % self.height
        return self.energy_field[x, y]
    
    @ti.kernel
    def get_energy_gradient(self, pos: ti.template()) -> ti.Vector([2], ti.f32):
        """Get energy gradient at a position"""
        x = int(pos[0]) % self.width
        y = int(pos[1]) % self.height
        
        # Get neighbors with periodic boundary conditions
        left = (x - 1) % self.width
        right = (x + 1) % self.width
        up = (y - 1) % self.height
        down = (y + 1) % self.height
        
        # Calculate gradient
        grad_x = (self.energy_field[right, y] - self.energy_field[left, y]) / 2.0
        grad_y = (self.energy_field[x, down] - self.energy_field[x, up]) / 2.0
        
        return ti.Vector([grad_x, grad_y])
    
    @ti.kernel
    def get_high_energy_regions(self, threshold: ti.f32, 
                               regions: ti.template()) -> ti.i32:
        """Get positions of high energy regions"""
        count = 0
        
        for i, j in ti.ndrange(self.width, self.height):
            if self.energy_field[i, j] > threshold and count < regions.shape[0]:
                regions[count] = ti.Vector([float(i), float(j)])
                count += 1
        
        return count
    
    def get_energy_field(self) -> np.ndarray:
        """Get energy field as numpy array"""
        return self.energy_field.to_numpy()
    
    def get_energy_sources(self) -> List[Dict]:
        """Get list of active energy sources"""
        sources = []
        
        for i in range(self.max_sources):
            if self.source_active[i] == 1:
                sources.append({
                    'index': i,
                    'position': self.source_positions[i].to_numpy(),
                    'intensity': self.source_intensities[i],
                    'radius': self.source_radii[i]
                })
        
        return sources
    
    def get_energy_stats(self) -> Dict:
        """Get energy system statistics"""
        energy_field = self.energy_field.to_numpy()
        
        return {
            'total_energy': float(np.sum(energy_field)),
            'max_energy': float(np.max(energy_field)),
            'mean_energy': float(np.mean(energy_field)),
            'energy_std': float(np.std(energy_field)),
            'active_sources': self.source_count[None],
            'energy_decay_rate': self.energy_decay_rate[None],
            'energy_diffusion_rate': self.energy_diffusion_rate[None],
            'energy_threshold': self.energy_threshold[None]
        }
    
    def set_energy_parameters(self, decay_rate: float = None, 
                            diffusion_rate: float = None,
                            threshold: float = None):
        """Set energy system parameters"""
        if decay_rate is not None:
            self.energy_decay_rate[None] = decay_rate
        if diffusion_rate is not None:
            self.energy_diffusion_rate[None] = diffusion_rate
        if threshold is not None:
            self.energy_threshold[None] = threshold

class EnergySource:
    """Individual energy source with temporal behavior"""
    
    def __init__(self, position: Tuple[float, float], intensity: float, 
                 radius: float, duration: float = -1.0):
        self.position = position
        self.intensity = intensity
        self.radius = radius
        self.duration = duration  # -1 for permanent
        self.start_time = time.time()
        self.active = True
    
    def is_active(self) -> bool:
        """Check if source is still active"""
        if not self.active:
            return False
        
        if self.duration > 0:
            return time.time() - self.start_time < self.duration
        
        return True
    
    def get_current_intensity(self) -> float:
        """Get current intensity (may decay over time)"""
        if not self.is_active():
            return 0.0
        
        # Simple linear decay
        if self.duration > 0:
            elapsed = time.time() - self.start_time
            decay_factor = max(0.0, 1.0 - elapsed / self.duration)
            return self.intensity * decay_factor
        
        return self.intensity
    
    def deactivate(self):
        """Deactivate the source"""
        self.active = False

class EnergyManager:
    """High-level energy management"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.energy_system = EnergySystem(config)
        self.energy_sources: List[EnergySource] = []
        
        # Energy events
        self.energy_events = []
        self.event_count = 0
    
    def add_energy_source(self, position: Tuple[float, float], intensity: float,
                         radius: float, duration: float = -1.0) -> int:
        """Add an energy source"""
        source = EnergySource(position, intensity, radius, duration)
        self.energy_sources.append(source)
        
        # Add to energy system
        source_idx = self.energy_system.add_energy_source(
            ti.Vector(position), intensity, radius
        )
        
        return source_idx
    
    def add_energy_impulse(self, position: Tuple[float, float], intensity: float,
                          radius: float):
        """Add a one-time energy impulse"""
        self.energy_system.add_energy_impulse(
            ti.Vector(position), intensity, radius
        )
        
        # Record event
        self.energy_events.append({
            'type': 'impulse',
            'position': position,
            'intensity': intensity,
            'radius': radius,
            'timestamp': time.time()
        })
        self.event_count += 1
    
    def update(self, dt: float):
        """Update energy system"""
        # Update energy sources
        active_sources = []
        for source in self.energy_sources:
            if source.is_active():
                active_sources.append(source)
            else:
                # Remove inactive source from energy system
                # This would require tracking source indices
                pass
        
        self.energy_sources = active_sources
        
        # Update energy system
        self.energy_system.update_energy_field(dt)
    
    def get_energy_at_position(self, position: Tuple[float, float]) -> float:
        """Get energy at a specific position"""
        return self.energy_system.get_energy_at_position(ti.Vector(position))
    
    def get_energy_gradient(self, position: Tuple[float, float]) -> Tuple[float, float]:
        """Get energy gradient at a position"""
        gradient = self.energy_system.get_energy_gradient(ti.Vector(position))
        return gradient.to_numpy()
    
    def get_high_energy_regions(self, threshold: float = None) -> List[Tuple[float, float]]:
        """Get positions of high energy regions"""
        if threshold is None:
            threshold = self.config.energy_threshold
        
        regions = ti.Vector.field(2, dtype=ti.f32, shape=(1000,))
        count = self.energy_system.get_high_energy_regions(threshold, regions)
        
        return [regions[i].to_numpy() for i in range(count)]
    
    def get_stats(self) -> Dict:
        """Get energy management statistics"""
        stats = self.energy_system.get_energy_stats()
        stats.update({
            'active_sources': len(self.energy_sources),
            'total_events': self.event_count,
            'recent_events': len(self.energy_events)
        })
        
        return stats
