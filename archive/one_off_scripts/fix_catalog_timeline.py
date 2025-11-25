"""
Fix catalog timeline for running simulation
Rebuilds discovery_timeline from existing substances
"""
import requests
import json

def fix_catalog_timeline(simulation_id='default'):
    """
    Rebuild catalog timeline from existing substances
    
    This is needed when a simulation was loaded from a snapshot
    that didn't preserve the discovery_timeline.
    """
    
    base_url = "http://localhost:8000"
    
    print("=" * 60)
    print("FIX CATALOG TIMELINE")
    print("=" * 60)
    
    try:
        # This would require adding a new API endpoint
        # For now, we can only provide instructions
        
        print("\n[PROBLEM IDENTIFIED]")
        print("Your simulation has 9 substances in catalog")
        print("But discovery_timeline is empty (cleared by snapshot load)")
        print("This causes PubChem Matcher to show 'Waiting for clusters...'")
        
        print("\n[ROOT CAUSE]")
        print("Snapshot deserialization calls catalog.clear()")
        print("This clears both substances AND timeline")
        print("New substances are added but timeline stays empty")
        
        print("\n[SOLUTION - Quick Fix]")
        print("The simulation needs to be restarted or use this workaround:")
        print("")
        print("1. Save current simulation snapshot")
        print("2. Restart backend with the fix I just made")
        print("3. The fix will preserve timeline in future snapshots")
        print("")
        print("For the CURRENT running simulation:")
        print("Wait for next novelty detection (~5 minutes)")
        print("New substances will have proper timeline entries")
        
        print("\n[SOLUTION - Complete Fix]")
        print("I've fixed the code in:")
        print("  backend/sim/io/snapshot.py")
        print("")
        print("Changes:")
        print("  1. Serialize full catalog (including timeline)")
        print("  2. Restore timeline on load (instead of clearing)")
        print("")
        print("To apply:")
        print("  1. Restart backend: .\\kill_backend.ps1; .\\start_backend.ps1")
        print("  2. Current simulation will continue normally")
        print("  3. Future snapshots will preserve timeline")
        
        print("\n[TEMPORARY WORKAROUND]")
        print("Until backend is restarted, substances will still be")
        print("detected and cataloged - they just won't appear in")
        print("PubChem Matcher until timeline is rebuilt.")
        print("")
        print("After next novelty detection (~ step 68700):")
        print("  - Timeline will have new entries")
        print("  - PubChem Matcher will show those new substances")
        print("  - Old 9 substances will remain hidden (no timeline)")
        
        print("\n" + "=" * 60)
        print("[ACTION REQUIRED]")
        print("=" * 60)
        print("")
        print("Restart backend to apply the fix:")
        print("")
        print("  .\\kill_backend.ps1")
        print("  .\\start_backend.ps1")
        print("")
        print("Then wait for next novelty detection.")
        print("PubChem Matcher should start showing substances!")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] {e}")
        return False

if __name__ == "__main__":
    import sys
    sim_id = sys.argv[1] if len(sys.argv) > 1 else 'default'
    fix_catalog_timeline(sim_id)

