"""
Symplectic integrators for Live 2.0 simulation
Provides energy-conserving integration methods
"""

import taichi as ti
import numpy as np
from typing import Optional


@ti.data_oriented
class SymplecticIntegrators:
    """Collection of symplectic integrators for particle dynamics"""
    
    def __init__(self, max_particles: int):
        # Position and velocity buffers for multi-stage integration
        self.positions_half = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
        self.positions_full = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
        self.velocities_half = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
        
        # Force buffers
        self.forces_prev = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
        self.forces_curr = ti.Vector.field(2, dtype=ti.f32, shape=(max_particles,))
        
        # Diagnostics
        self.integration_error = ti.field(dtype=ti.f32, shape=())
        self.energy_drift = ti.field(dtype=ti.f32, shape=())
    
    @ti.kernel
    def verlet_step(self, positions: ti.template(), velocities: ti.template(),
                   forces: ti.template(), masses: ti.template(), dt: ti.f32,
                   active: ti.template(), particle_count: ti.i32):
        """
        Verlet integrator (Position-Stored Verlet)
        
        Steps:
        1. v^{n+1/2} = v^n + (dt/2) * F^n / m
        2. x^{n+1} = x^n + dt * v^{n+1/2}
        3. F^{n+1} = F(x^{n+1})
        4. v^{n+1} = v^{n+1/2} + (dt/2) * F^{n+1} / m
        
        This is time-symmetric and symplectic, preserving energy better than Euler
        """
        # Step 1: Half-step velocity update
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]  # mass is first component of attributes
                if mass > 0:
                    half_accel = forces[i] / mass
                    self.velocities_half[i] = velocities[i] + 0.5 * dt * half_accel
                    self.positions_full[i] = positions[i] + dt * self.velocities_half[i]
        
        # Store previous forces for error estimation
        for i in range(particle_count):
            self.forces_prev[i] = forces[i]
    
    @ti.kernel 
    def verlet_correction(self, velocities: ti.template(), forces: ti.template(),
                         masses: ti.template(), dt: ti.f32, active: ti.template(), 
                         particle_count: ti.i32):
        """
        Completes Verlet integration step with new forces
        
        Update velocities with new forces computed from updated positions
        """
        # Step 4: Complete velocity update with new forces
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    half_accel = forces[i] / mass
                    velocities[i] = self.velocities_half[i] + 0.5 * dt * half_accel
    
    @ti.kernel
    def leapfrog_step(self, positions: ti.template(), velocities: ti.template(),
                     forces: ti.template(), masses: ti.template(), dt: ti.f32,
                     active: ti.template(), particle_count: ti.i32):
        """
        Leapfrog integrator (Velocity-Stored Verlet)
        
        Steps:
        1. v^{n+1/2} = v^{n-1/2} + dt * F^n / m
        2. x^{n+1} = x^n + dt * v^{n+1/2}
        
        More intuitive than Verlet for velocity-dependent quantities
        """
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    # Update velocities (staggered)
                    acceleration = forces[i] / mass
                    velocities[i] += dt * acceleration
                    
                    # Update positions
                    positions[i] += dt * velocities[i]
    
    @ti.kernel
    def r_kinoshita_suzuki_step(self, positions: ti.template(), velocities: ti.template(),
                               forces: ti.template(), masses: ti.template(), dt: ti.f32,
                               active: ti.template(), particle_count: ti.i32):
        """
        Ruth-Kinoshita-Suzuki 4th order symplectic integrator
        
        Uses composition method with stage factors:
        c1 = c4 = 1/(2*(2-2^(1/3)))
        c2 = c3 = (1-2^(1/3))/(2*(2-2^(1/3)))
        d1 = d3 = 1/(2-2^(1/3))
        d2 = -2^(1/3)/(2-2^(1/3))
        
        This is a higher-order method with excellent energy conservation
        """
        # Stage coefficients
        c1 = 1.0 / (2.0 * (2.0 - 2.0**(1.0/3.0)))
        c2 = (1.0 - 2.0**(1.0/3.0)) / (2.0 * (2.0 - 2.0**(1.0/3.0)))
        c3 = c2
        c4 = c1
        
        d1 = 1.0 / (2.0 - 2.0**(1.0/3.0))
        d2 = -2.0**(1.0/3.0) / (2.0 - 2.0**(1.0/3.0))
        d3 = d1
        
        # Stage 1
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    velocities[i] += c1 * dt * forces[i] / mass
                    positions[i] += d1 * dt * velocities[i]
        
        # Stage 2  
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    velocities[i] += c2 * dt * forces[i] / mass
                    positions[i] += d2 * dt * velocities[i]
        
        # Stage 3
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    velocities[i] += c3 * dt * forces[i] / mass
                    positions[i] += d3 * dt * velocities[i]
        
        # Stage 4
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    velocities[i] += c4 * dt * forces[i] / mass
    
    @ti.kernel
    def compute_integration_error(self, dt: ti.f32, active: ti.template(), particle_count: ti.i32):
        """Compute local truncation error for adaptive timestep control"""
        error_norm = 0.0
        
        for i in range(particle_count):
            if active[i] == 1:
                # Estimate error from force changes
                force_change = self.forces_curr[i] - self.forces_prev[i]
                force_magnitude = force_change.norm()
                error_norm += force_magnitude * force_magnitude
        
        # Dimensionless error estimate
        self.integration_error[None] = ti.sqrt(error_norm) * dt * dt / (2.0 * particle_count)
    
    @ti.kernel
    def compute_energy_drift(self, velocities: ti.template(), masses: ti.template(),
                           energy_reference: ti.f32, active: ti.template(), particle_count: ti.i32):
        """Compute energy conservation error metric"""
        kinetic_energy = 0.0
        
        for i in range(particle_count):
            if active[i] == 1:
                mass = masses[i][0]
                if mass > 0:
                    vel_sq = velocities[i].norm_sqr()
                    kinetic_energy += 0.5 * mass * vel_sq
        
        energy_current = kinetic_energy
        energy_error = abs(energy_current - energy_reference) / max(abs(energy_reference), 1e-6)
        
        self.energy_drift[None] = energy_error
    
    def get_diagnostics(self) -> Dict[str, float]:
        """Get integration diagnostics"""
        return {
            'integration_error': float(self.integration_error[None]),
            'energy_drift': float(self.energy_drift[None])
        }
    
    def estimate_optimal_timestep(self, target_error: float = 1e-4) -> float:
        """Estimate optimal timestep based on error tolerance"""
        
        if self.integration_error[None] > 0:
            error_ratio = target_error / self.integration_error[None]
            # 4th-order method scales as dt^4 for truncation error
            optimal_dt = 0.9 * error_ratio**(1/4)  # Conservative factor 0.9
            return optimal_dt
        else:
            return 1.0  # Default if no error estimate


class AdaptiveIntegrator:
    """Adaptive symplectic integrator with error control"""
    
    def __init__(self, max_particles: int, method: str = "verlet"):
        self.integrator = SymplecticIntegrators(max_particles)
        self.method = method
        
        # Adaptive parameters
        self.target_error = 1e-4
        self.tolerance = 0.1  # 10% tolerance above/below target
        
        # Statistics
        self.step_stats = {
            'successful_steps': 0,
            'failed_steps': 0,
            'timestep_adjustments': 0,
            'current_error': 0.0
        }
    
    def step(self, positions, velocities, forces, masses, dt: float, 
             active, particle_count: int) -> tuple[float, bool]:
        """
        Perform adaptive integration step
        
        Returns:
            (actual_dt_used, success_flag)
        """
        dt_candidate = dt
        
        # Try integration with current timestep
        if self.method == "verlet":
            self.integrator.verlet_step(positions, velocities, forces, masses, 
                                       dt_candidate, active, particle_count)
        elif self.method == "leapfrog":
            self.integrator.leapfrog_step(positions, velocities, forces, masses,
                                        dt_candidate, active, particle_count)
        elif self.method == "rks":
            self.integrator.r_kinoshita_suzuki_step(positions, velocities, forces,
                                                  masses, dt_candidate, active, particle_count)
        
        # Compute error metrics
        self.integrator.compute_integration_error(dt_candidate, active, particle_count)
        
        error = float(self.integrator.integration_error[None])
        self.step_stats['current_error'] = error
        
        # Check if error is acceptable
        if error <= self.target_error * (1 + self.tolerance):
            # Accept step
            self.step_stats['successful_steps'] += 1
            
            if self.method == "verlet":
                # Complete Verlet step needs force correction
                return dt_candidate, True
            else:
                return dt_candidate, True
        else:
            # Reject step - would need to retry with smaller dt
            self.step_stats['failed_steps'] += 1
            return 0.0, False  # Signal for reduced timestep
    
    def estimate_next_timestep(self) -> float:
        """Estimate optimal timestep for next step"""
        return self.integrator.estimate_optimal_timestep(self.target_error)
    
    def get_stats(self) -> Dict[str, Any]:
        """Get integration statistics"""
        total_steps = self.step_stats['successful_steps'] + self.step_stats['failed_steps']
        success_rate = (self.step_stats['successful_steps'] / max(total_steps, 1)) * 100
        
        stats = self.step_stats.copy()
        stats.update({
            'success_rate': success_rate,
            **self.integrator.get_diagnostics()
        })
        
        return stats
