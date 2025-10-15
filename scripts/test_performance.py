#!/usr/bin/env python3
"""
Performance Test Script for Phase 2 Optimizations
=================================================

Tests simulation performance with optimized settings.
Target: >100 steps/second for 10,000 step run

Usage:
    python scripts/test_performance.py
"""

import sys
import time
import logging
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import subprocess

def run_performance_test():
    """Run optimized performance test"""
    
    print("=" * 70)
    print("PHASE 2 PERFORMANCE TEST")
    print("=" * 70)
    print()
    print("Configuration: configs/phase2_performance_test.yaml")
    print("Target: >100 steps/second")
    print("Duration: 10,000 steps")
    print()
    print("=" * 70)
    print()
    
    # Run simulation
    start_time = time.time()
    
    cmd = [
        'python',
        'scripts/run_phase2_full.py',
        '--config', 'configs/phase2_performance_test.yaml',
        '--output', 'results/performance_test',
        '--steps', '10000',
        '--seed', '42'
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print()
    print("Running test...")
    print("-" * 70)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=False,  # Show output in real-time
            text=True
        )
        
        elapsed_time = time.time() - start_time
        
        print()
        print("=" * 70)
        print("PERFORMANCE TEST RESULTS")
        print("=" * 70)
        
        if result.returncode == 0:
            print("✅ Test completed successfully!")
            print()
            print(f"Total time: {elapsed_time:.2f} seconds")
            print(f"Steps completed: 10,000")
            print(f"Performance: {10000 / elapsed_time:.2f} steps/second")
            print()
            
            if 10000 / elapsed_time >= 100:
                print("🎉 SUCCESS: Performance target ACHIEVED (>100 steps/s)!")
                multiplier = (10000 / elapsed_time) / 2  # Current vs before (2 steps/s)
                print(f"   Speed improvement: {multiplier:.1f}x faster than before")
                print()
                print(f"📊 Projected time for 10M steps: {(10000000 * elapsed_time / 10000) / 3600:.1f} hours")
                print()
                return True
            else:
                print("⚠️  WARNING: Performance below target")
                print(f"   Current: {10000 / elapsed_time:.2f} steps/s")
                print(f"   Target: 100 steps/s")
                print(f"   Shortfall: {100 / (10000 / elapsed_time):.1f}x too slow")
                print()
                return False
        else:
            print("❌ Test FAILED")
            print(f"   Exit code: {result.returncode}")
            print()
            return False
            
    except Exception as e:
        print()
        print("❌ ERROR during test:")
        print(f"   {e}")
        print()
        return False

def main():
    success = run_performance_test()
    
    print("=" * 70)
    if success:
        print("✅ PERFORMANCE TEST: PASSED")
        print()
        print("Next steps:")
        print("  1. Update Phase 2 configs with optimized settings")
        print("  2. Run full 30 simulations (10M steps each)")
        print("  3. Expected completion: ~9 days (with 4 parallel runs)")
    else:
        print("❌ PERFORMANCE TEST: FAILED")
        print()
        print("Next steps:")
        print("  1. Review bottlenecks in logs")
        print("  2. Apply additional optimizations")
        print("  3. Re-run performance test")
    print("=" * 70)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

