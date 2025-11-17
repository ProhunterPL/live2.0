"""
Diagnose why chemistry is still unrealistic despite code fixes
"""
import requests
import json

def diagnose_chemistry(simulation_id='default'):
    """Check what's causing unrealistic molecules"""
    
    base_url = "http://localhost:8000"
    
    print("=" * 70)
    print("CHEMISTRY DIAGNOSTICS")
    print("=" * 70)
    
    try:
        response = requests.get(f"{base_url}/simulation/{simulation_id}/metrics")
        if response.status_code != 200:
            print("\n[ERROR] Cannot connect to simulation")
            print("Make sure backend is running")
            return
        
        data = response.json()
        metrics = data.get('metrics', {})
        
        step = metrics.get('step_count', 0)
        particles = metrics.get('particle_count', 0)
        bonds = metrics.get('bond_count', 0)
        clusters_approx = metrics.get('cluster_count', 0)
        
        print(f"\n[CURRENT STATE]")
        print(f"  Step: {step}")
        print(f"  Particles: {particles}")
        print(f"  Bonds: {bonds}")
        print(f"  Clusters (approx): {clusters_approx}")
        
        if bonds > 0 and particles > 0:
            bonds_per_particle = bonds / particles
            print(f"  Bonds per particle: {bonds_per_particle:.2f}")
            
            if bonds_per_particle < 0.3:
                print(f"\n  ⚠️ LOW BONDING RATIO!")
                print(f"     Possible causes:")
                print(f"     1. binding_threshold too high (check config)")
                print(f"     2. Particles too far apart")
                print(f"     3. Low binding probability")
            elif bonds_per_particle > 2.0:
                print(f"\n  ⚠️ HIGH BONDING RATIO!")
                print(f"     Possible causes:")
                print(f"     1. binding_threshold too low")
                print(f"     2. Bonds forming at long distances")
            else:
                print(f"\n  ✅ Bonding ratio looks OK")
        
        print(f"\n" + "=" * 70)
        print("[PROBLEM ANALYSIS]")
        print("=" * 70)
        
        print(f"\nYour code has the fixes:")
        print(f"  ✅ binding_probability > 0.15 (not 0.005)")
        print(f"  ✅ max_formation_dist: 1.2-1.8 units (not 3.4)")
        print(f"  ✅ max_strain: 0.3-0.8 (not 3.0)")
        
        print(f"\nBut molecules are still unrealistic!")
        print(f"\nMost likely causes:")
        print(f"  1. BACKEND NOT RESTARTED after code changes")
        print(f"  2. Simulation loaded from OLD SNAPSHOT (with old code)")
        print(f"  3. Configuration overrides code parameters")
        
        print(f"\n" + "=" * 70)
        print("[SOLUTIONS]")
        print("=" * 70)
        
        print(f"\n1. RESTART BACKEND (most important!):")
        print(f"   .\\kill_backend.ps1")
        print(f"   .\\start_backend.ps1")
        print(f"\n   This will load the new code with fixed parameters")
        
        print(f"\n2. START NEW SIMULATION (don't load snapshot):")
        print(f"   - In frontend, click 'Create Simulation'")
        print(f"   - Start fresh (not from snapshot)")
        print(f"   - New simulation will use fixed code")
        
        print(f"\n3. MONITOR BOND LENGTHS:")
        print(f"   - Look at 'Largest Connected Cluster' panel")
        print(f"   - Density should be > 0.3 (better > 0.4)")
        print(f"   - If still low density:")
        print(f"     a) Backend wasn't restarted OR")
        print(f"     b) Old snapshot was loaded OR")
        print(f"     c) Config overrides are wrong")
        
        print(f"\n4. CHECK LOGS after restart:")
        print(f"   Get-Content logs\\logs.txt -Tail 30")
        print(f"   Should show simulation initialization")
        
        print(f"\n" + "=" * 70)
        print("[EXPECTED RESULTS after fix]")
        print("=" * 70)
        
        print(f"\nBEFORE (unrealistic):")
        print(f"  - Density: 0.15-0.25 (very low)")
        print(f"  - Long stretched bonds visible")
        print(f"  - 'Spider web' structures")
        print(f"  - Energy often 0.00")
        
        print(f"\nAFTER (realistic):")
        print(f"  - Density: 0.35-0.55 (good)")
        print(f"  - Short compact bonds")
        print(f"  - Tight molecular structures")
        print(f"  - Realistic energy values")
        
        print(f"\n" + "=" * 70)
        print("[VERIFICATION CHECKLIST]")
        print("=" * 70)
        
        print(f"\n□ Backend restarted after code changes?")
        print(f"□ New simulation created (not from snapshot)?")
        print(f"□ Waited 1000+ steps for chemistry to develop?")
        print(f"□ Checked 'Largest Connected Cluster' density > 0.3?")
        print(f"□ Visually inspected bonds (short, not stretched)?")
        
        print(f"\nIf all checked ✓ and still unrealistic:")
        print(f"  → Post screenshot in chat")
        print(f"  → Include metrics (density, size, bonds)")
        print(f"  → Include step count")
        
        return True
        
    except Exception as e:
        print(f"\n[ERROR] {e}")
        print("\nMake sure backend is running:")
        print("  .\\start_backend.ps1")
        return False

if __name__ == "__main__":
    import sys
    sim_id = sys.argv[1] if len(sys.argv) > 1 else 'default'
    diagnose_chemistry(sim_id)

