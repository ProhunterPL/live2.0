"""
Diagnostic script to check real cluster status in running simulation
"""
import requests
import json
import sys

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')

def check_clusters(simulation_id='default'):
    """Check real cluster and novel substance status"""
    
    base_url = "http://localhost:8000"
    
    print("=" * 60)
    print("LIVE 2.0 CLUSTER DIAGNOSTICS")
    print("=" * 60)
    
    # 1. Check metrics
    try:
        response = requests.get(f"{base_url}/simulation/{simulation_id}/metrics")
        if response.status_code == 200:
            data = response.json()
            metrics = data.get('metrics', {})
            
            print("\n[METRICS] (approximated):")
            print(f"  Particles: {metrics.get('particle_count', 0)}")
            print(f"  Bonds: {metrics.get('bond_count', 0)}")
            print(f"  Clusters (approx): {metrics.get('cluster_count', 0)}")
            print(f"  Step: {metrics.get('step_count', 0)}")
            
            print(f"\n[CATALOG] Stats:")
            print(f"  Total Novel: {metrics.get('total_novel', 0)}")
            print(f"  Total Discovered: {metrics.get('total_discovered', 0)}")
            print(f"  Novelty Rate: {metrics.get('novelty_rate', 0):.4f}")
        else:
            print(f"[ERROR] Failed to get metrics: {response.status_code}")
            return
    except Exception as e:
        print(f"[ERROR] Error getting metrics: {e}")
        return
    
    # 2. Check novel substances
    try:
        response = requests.get(f"{base_url}/simulation/{simulation_id}/novel-substances?count=50")
        if response.status_code == 200:
            data = response.json()
            substances = data.get('substances', [])
            
            print(f"\n[NOVEL SUBSTANCES] (detected by catalog):")
            print(f"  Count: {len(substances)}")
            
            if len(substances) > 0:
                print(f"\n  Recent discoveries:")
                for i, sub in enumerate(substances[:5]):
                    print(f"    {i+1}. ID: {sub['id'][-8:]}")
                    print(f"       Size: {sub['size']} atoms, Bonds: {sub['properties']['bonds']}")
                    print(f"       Density: {sub['properties']['graph_density']:.3f}")
                    print(f"       Complexity: {sub['complexity']:.3f}")
            else:
                print("\n  [WARNING] No novel substances detected yet!")
                print(f"  Possible reasons:")
                print(f"    1. Novelty detection runs every 700 steps (last run: step {(metrics.get('step_count', 0) // 700) * 700})")
                print(f"    2. Clusters must have >= 3 particles (min_cluster_size)")
                print(f"    3. Clusters are cached (updated every 500 steps)")
        else:
            print(f"[ERROR] Failed to get novel substances: {response.status_code}")
    except Exception as e:
        print(f"[ERROR] Error getting novel substances: {e}")
    
    # 3. Check if data caching is causing issues
    step_count = metrics.get('step_count', 0)
    last_cluster_update = (step_count // 500) * 500
    last_novelty_check = (step_count // 700) * 700
    
    print(f"\n[TIMING] Analysis:")
    print(f"  Current step: {step_count}")
    print(f"  Last cluster update: {last_cluster_update} (every 500 steps)")
    print(f"  Last novelty check: {last_novelty_check} (every 700 steps)")
    print(f"  Next cluster update: {last_cluster_update + 500}")
    print(f"  Next novelty check: {last_novelty_check + 700}")
    
    if step_count < 700:
        print(f"\n[WARNING] Simulation hasn't reached first novelty check yet!")
        print(f"   Wait until step 700 for first novel substance detection")
    
    print("\n" + "=" * 60)
    print("[RECOMMENDATIONS]")
    print("=" * 60)
    
    if len(substances) == 0:
        print("1. Wait for novelty check at step", last_novelty_check + 700)
        print("2. Check if clusters have >= 3 particles")
        print("3. Clusters with long bonds may be unstable")
        print("4. Consider reducing bond_distance threshold")
    else:
        print("[OK] Novel substances are being detected!")
        print(f"   {len(substances)} substances found")
    
    # 4. Check bond statistics
    bond_count = metrics.get('bond_count', 0)
    particle_count = metrics.get('particle_count', 0)
    
    if bond_count > 0 and particle_count > 0:
        avg_bonds_per_particle = bond_count / particle_count
        print(f"\n[BOND STATS]:")
        print(f"  Average bonds per particle: {avg_bonds_per_particle:.2f}")
        
        if avg_bonds_per_particle < 0.5:
            print(f"  [WARNING] Low bonding ratio - clusters may be sparse or unstable")
        elif avg_bonds_per_particle > 2.0:
            print(f"  [WARNING] High bonding ratio - many bonds forming")

if __name__ == "__main__":
    import sys
    sim_id = sys.argv[1] if len(sys.argv) > 1 else 'default'
    check_clusters(sim_id)

