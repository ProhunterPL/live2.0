"""
Force immediate cluster detection and novel substance catalog update
"""
import requests
import json

def force_detection(simulation_id='default'):
    """
    This script helps debug cluster detection by checking current state
    
    NOTE: There's no direct API to force detection, but we can check:
    1. Current step count
    2. When was last detection
    3. When will be next detection
    4. Current novel substances
    """
    
    base_url = "http://localhost:8000"
    
    print("=" * 60)
    print("FORCE CLUSTER DETECTION - DEBUG HELPER")
    print("=" * 60)
    
    try:
        # 1. Get current metrics
        response = requests.get(f"{base_url}/simulation/{simulation_id}/metrics")
        if response.status_code != 200:
            print(f"[ERROR] Cannot connect to simulation: {response.status_code}")
            print("\nIs the backend running?")
            print("Start it with: .\\start_backend.ps1")
            return
        
        data = response.json()
        metrics = data.get('metrics', {})
        
        step_count = metrics.get('step_count', 0)
        particle_count = metrics.get('particle_count', 0)
        bond_count = metrics.get('bond_count', 0)
        cluster_count_approx = metrics.get('cluster_count', 0)
        
        print(f"\n[CURRENT STATE] Step {step_count}")
        print(f"  Particles: {particle_count}")
        print(f"  Bonds: {bond_count}")
        print(f"  Clusters (approx): {cluster_count_approx}")
        
        # Calculate timing
        novelty_check_interval = 700  # From config
        last_novelty_check = (step_count // novelty_check_interval) * novelty_check_interval
        next_novelty_check = last_novelty_check + novelty_check_interval
        steps_until_detection = next_novelty_check - step_count
        
        print(f"\n[DETECTION TIMING]")
        print(f"  Novelty check interval: every {novelty_check_interval} steps")
        print(f"  Last detection: step {last_novelty_check}")
        print(f"  Next detection: step {next_novelty_check}")
        print(f"  Steps until next: {steps_until_detection}")
        
        if steps_until_detection > novelty_check_interval / 2:
            print(f"\n  >> Far from next detection ({steps_until_detection} steps)")
            print(f"     Consider waiting or manually checking at step {next_novelty_check}")
        else:
            print(f"\n  >> Close to next detection! ({steps_until_detection} steps)")
            print(f"     Novel substances should appear soon")
        
        # 2. Check novel substances
        response = requests.get(f"{base_url}/simulation/{simulation_id}/novel-substances?count=50")
        if response.status_code == 200:
            data = response.json()
            substances = data.get('substances', [])
            
            print(f"\n[NOVEL SUBSTANCES]")
            print(f"  Count: {len(substances)}")
            
            if len(substances) > 0:
                print(f"\n  Top 5 recent:")
                for i, sub in enumerate(substances[:5]):
                    print(f"    {i+1}. Size: {sub['size']} atoms, Bonds: {sub['properties']['bonds']}, Density: {sub['properties']['graph_density']:.3f}")
            else:
                print(f"\n  [NO SUBSTANCES YET]")
                print(f"  Possible reasons:")
                print(f"    1. Haven't reached step {next_novelty_check} (next detection)")
                print(f"    2. No clusters with >=3 particles (min_cluster_size)")
                print(f"    3. Clusters are too sparse (low density)")
                print(f"    4. Bonds are unstable (breaking before detection)")
        
        # 3. Analyze bonding
        if particle_count > 0 and bond_count > 0:
            bonds_per_particle = bond_count / particle_count
            print(f"\n[BONDING ANALYSIS]")
            print(f"  Bonds per particle: {bonds_per_particle:.2f}")
            
            if bonds_per_particle < 0.5:
                print(f"  [WARNING] Low bonding ratio!")
                print(f"    Clusters are sparse or bonds are breaking quickly")
                print(f"    Consider:")
                print(f"      - Reducing unbinding_threshold")
                print(f"      - Increasing bond lifetime (max_age)")
                print(f"      - Reducing max_bond_length")
            elif bonds_per_particle > 2.0:
                print(f"  [INFO] High bonding ratio")
                print(f"    Many bonds forming - good for cluster formation")
            else:
                print(f"  [OK] Moderate bonding ratio")
        
        # 4. Recommendations
        print(f"\n" + "=" * 60)
        print(f"[RECOMMENDATIONS]")
        print(f"=" * 60)
        
        if len(substances) == 0:
            if steps_until_detection <= 100:
                print("1. WAIT - Detection will happen in ~100 steps")
                print(f"   Check again at step {next_novelty_check}")
            else:
                print(f"1. WAIT for step {next_novelty_check} (novelty detection)")
                print("2. Or restart simulation with:")
                print("     novelty_check_interval: 200  # More frequent")
                
            if bonds_per_particle < 0.5:
                print("\n3. LOW BONDING DETECTED:")
                print("   Bonds may be breaking before detection")
                print("   Consider adjusting bond parameters:")
                print("     - Increase bond lifetime")
                print("     - Reduce unbinding threshold")
                print("     - Reduce max bond length (currently 6.8 * radius)")
        else:
            print(f"[OK] {len(substances)} novel substances detected!")
            print("PubChem Matcher should show them in frontend")
            print("\nIf not visible:")
            print("1. Refresh frontend (Ctrl+R)")
            print("2. Check browser console for errors")
            print("3. Wait for next frontend update (15 seconds)")
        
        print("\n" + "=" * 60)
        
    except Exception as e:
        print(f"[ERROR] {e}")
        print("\nMake sure backend is running:")
        print("  .\\start_backend.ps1")

if __name__ == "__main__":
    import sys
    sim_id = sys.argv[1] if len(sys.argv) > 1 else 'default'
    force_detection(sim_id)

