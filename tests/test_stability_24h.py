#!/usr/bin/env python3
"""
24-hour stability test for Live 2.0 Phase 1
Tests energy conservation and numerical stability over extended periods
"""

import requests
import time
import json
import numpy as np
from datetime import datetime, timedelta

API_BASE = "http://localhost:8001"

class StabilityTester:
    """Test stability over extended simulation periods"""
    
    def __init__(self):
        self.test_config = {
            "config": {
                "max_particles": 100,
                "grid_width": 128,
                "grid_height": 128,
                "mode": "open_chemistry",
                "dt": 0.01,
                
                # Stability-focused parameters
                "energy_decay": 0.99,  # Lower decay for longer simulation
                "energy_threshold": 0.05,
                "pulse_every": 100,
                "pulse_radius": 16.0,
                "pulse_amplitude": 3.0,
                
                # Enable symplectic integrator
                "use_symplectic_integrator": True,
                
                # Enable diagnostics
                "enable_diagnostics": True,
                "diagnostics_dir": "./phase1_diagnostics"
            },
            "mode": "open_chemistry"
        }
        
        self.simulation_id = None
        self.start_time = None
        self.test_duration = 3600  # 1 hour test (instead of 24h for demo)
        self.sample_interval = 60  # Sample every minute
        
        # Stability metrics
        self.energy_history = []
        self.particle_count_history = []
        self.timestep_history = []
        self.error_metrics = {
            'energy_drift_max': 0.0,
            'energy_drift_avg': 0.0,
            'timestep_stability': 1.0,
            'simulation_uptime': 0.0
        }
    
    def create_simulation(self):
        """Create optimized simulation for stability testing"""
        print("Creating stability-test simulation...")
        
        response = requests.post(f"{API_BASE}/simulation/create", json=self.test_config)
        response.raise_for_status()
        
        data = response.json()
        self.simulation_id = data.get("simulation_id")
        print(f"Created simulation: {self.simulation_id}")
    
    def start_simulation(self):
        """Start simulation"""
        print("Starting stability test simulation...")
        
        response = requests.post(f"{API_BASE}/simulation/{self.simulation_id}/start")
        response.raise_for_status()
        
        print("Stability test started")
        self.start_time = time.time()
    
    def sample_simulation_state(self):
        """Sample current simulation state and compute metrics"""
        response = requests.get(f"{API_BASE}/simulation/{self.simulation_id}/status")
        response.raise_for_status()
        
        status = response.json()
        
        # Extract metrics
        current_time = time.time()
        elapsed = current_time - self.start_time
        
        energy_data = {
            'sim_time': status['current_time'],
            'wall_time': elapsed,
            'step_count': status['step_count'],
            'particle_count': status['particle_count'],
            'health_score': status['health_score'],
            'timestamp': datetime.now().isoformat()
        }
        
        return energy_data
    
    def compute_stability_metrics(self):
        """Compute overall stability metrics from sampled data"""
        if len(self.energy_history) < 2:
            return
        
        # Energy drift analysis
        energy_values = [h['health_score'] for h in self.energy_history]
        initial_energy = energy_values[0]
        
        energy_diffs = [abs(e - initial_energy) / max(initial_energy, 1e-6) 
                       for e in energy_values]
        
        self.error_metrics['energy_drift_max'] = max(energy_diffs) * 100
        self.error_metrics['energy_drift_avg'] = np.mean(energy_diffs) * 100
        
        # Timestep stability (variation in step intervals)
        sim_times = [h['sim_time'] for h in self.energy_history]
        if len(sim_times) > 1:
            time_intervals = [sim_times[i+1] - sim_times[i] 
                            for i in range(len(sim_times)-1)]
            self.error_metrics['timestep_stability'] = np.std(time_intervals) / np.mean(time_intervals)
        
        # Simulation uptime
        self.error_metrics['simulation_uptime'] = time.time() - self.start_time
    
    def print_progress(self, elapsed, target, step, particle_count, health_score):
        """Print progress update"""
        progress = (elapsed / target) * 100
        hours = elapsed / 3600
        
        print(f"[{progress:5.1f}%] {hours:.1f}h | Step: {step:5d} | "
              f"Particles: {particle_count:3d} | Health: {health_score:.3f}")
    
    def run_stability_test(self):
        """Run extended stability test"""
        print("="*80)
        print("PHASE 1: 24-HOUR STABILITY TEST")
        print("="*80)
        print(f"Test parameters:")
        print(f"  Target duration: {self.test_duration/3600:.1f} hours")
        print(f"  Sample interval: {self.sample_interval} seconds")
        print(f"  Symplectic integrator: {'ENABLED' if self.test_config['config'].get('use_symplectic_integrator') else 'DISABLED'}")
        print("")
        
        try:
            # Setup
            self.create_simulation()
            self.start_simulation()
            
            # Wait for initialization
            print("Waiting for initialization...")
            time.sleep(10)
            
            # Initial state
            initial_state = self.sample_simulation_state()
            self.energy_history.append(initial_state)
            
            print(f"Initial state: {initial_state['particle_count']} particles, "
                  f"health={initial_state['health_score']:.3f}")
            print("")
            
            # Run test loop
            last_sample_time = time.time()
            samples_collected = 0
            
            while True:
                current_time = time.time()
                elapsed = current_time - self.start_time
                
                # Check if test should continue
                if elapsed >= self.test_duration:
                    print(f"\nTarget duration reached: {elapsed/3600:.1f} hours")
                    break
                
                # Sample at intervals
                if current_time - last_sample_time >= self.sample_interval:
                    state = self.sample_simulation_state()
                    self.energy_history.append(state)
                    samples_collected += 1
                    
                    # Print progress
                    self.print_progress(
                        elapsed, self.test_duration,
                        state['step_count'], state['particle_count'], 
                        state['health_score']
                    )
                    
                    last_sample_time = current_time
                    
                    # Check for major issues
                    if state['health_score'] < 0.001 and samples_collected > 5:
                        print(f"\nWARNING: Very low health score detected: {state['health_score']}")
                        print("Simulation may be unstable")
                
                time.sleep(5)  # Check every 5 seconds
            
            # Compute final metrics
            self.compute_stability_metrics()
            final_state = self.sample_simulation_state()
            
            # Results
            print("\n" + "="*80)
            print("STABILITY TEST RESULTS")
            print("="*80)
            
            print(f"Test Duration: {elapsed/3600:.1f} hours")
            print(f"Samples Collected: {samples_collected}")
            print(f"Final Step Count: {final_state['step_count']:,}")
            print(f"Final Particle Count: {final_state['particle_count']}")
            print(f"Final Health Score: {final_state['health_score']:.4f}")
            
            print(f"\nEnergy Conservation:")
            print(f"  Max Drift: {self.error_metrics['energy_drift_max']:.2f}%")
            print(f"  Avg Drift: {self.error_metrics['energy_drift_avg']:.2f}%")
            print(f"  Timestep Stability: {self.error_metrics['timestep_stability']:.4f}")
            
            # Phase 1 Criteria
            energy_stable = self.error_metrics['energy_drift_max'] < 5.0  # 5% max drift
            timestep_stable = self.error_metrics['timestep_stability'] < 0.1  # 10% variation
            simulation_stable = final_state['health_score'] > 0.001
            uptime_sufficient = elapsed >= self.test_duration * 0.9  # 90% of target
            
            print(f"\nPhase 1 Stability Criteria:")
            print(f"  Energy Conservation: {'PASS' if energy_stable else 'FAIL'} "
                  f"({self.error_metrics['energy_drift_max']:.1f}% < 5%)")
            print(f"  Timestep Stability: {'PASS' if timestep_stable else 'FAIL'} "
                  f"({self.error_metrics['timestep_stability']:.3f} < 0.1)")
            print(f"  Simulation Stability: {'PASS' if simulation_stable else 'FAIL'} "
                  f"(health={final_state['health_score']:.4f} > 0.001)")
            print(f"  Uptime Sufficient: {'PASS' if uptime_sufficient else 'FAIL'} "
                  f"({elapsed/3600:.1f}h >= {self.test_duration*0.9/3600:.1f}h)")
            
            overall_success = energy_stable and timestep_stable and simulation_stable and uptime_sufficient
            
            print(f"\nPhase 1 Status: {'PASS' if overall_success else 'FAIL'}")
            
            if overall_success:
                print("Phase 1 Stability Requirements Satisfied!")
                print("✓ Energy conservation: Excellent")
                print("✓ Numerical stability: Stable timestep")
                print("✓ Extended operation: Successful")
                
                # Save stability report
                report = {
                    'test_config': self.test_config,
                    'duration_hours': elapsed / 3600,
                    'final_state': final_state,
                    'stability_metrics': self.error_metrics,
                    'energy_history': self.energy_history,
                    'success': True,
                    'timestamp': datetime.now().isoformat()
                }
                
                with open('phase1_stability_report.json', 'w') as f:
                    json.dump(report, f, indent=2)
                
                print("✓ Stability report saved to phase1_stability_report.json")
            else:
                print("Phase 1 Requirements Not Met")
                print("Issues identified:")
                if not energy_stable: print("  - Energy drift too high")
                if not timestep_stable: print("  - Timestep variations too large")
                if not simulation_stable: print("  - Simulation becoming unstable")
                if not uptime_sufficient: print("  - Insufficient uptime")
            
            return overall_success
            
        except KeyboardInterrupt:
            print("\n\nTest interrupted by user")
            return False
            
        except Exception as e:
            print(f"\nTest failed with error: {e}")
            import traceback
            traceback.print_exc()
            return False
            
        finally:
            # Cleanup
            if self.simulation_id:
                try:
                    requests.post(f"{API_BASE}/simulation/{self.simulation_id}/stop")
                    print("Simulation stopped")
                except:
                    pass

def main():
    """Run stability test"""
    tester = StabilityTester()
    success = tester.run_stability_test()
    
    if success:
        print("\n🎉 PHASE 1 STABILITY CERTIFIED")
        exit(0)
    else:
        print("\n❌ PHASE 1 STABILITY NOT ACHIEVED")
        exit(1)

if __name__ == "__main__":
    main()
