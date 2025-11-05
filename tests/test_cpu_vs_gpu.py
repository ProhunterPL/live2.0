#!/usr/bin/env python3
"""
Quick CPU vs GPU Performance Test
==================================

Runs a short simulation (10K steps) on both CPU and GPU to compare performance.
"""

import sys
import time
from pathlib import Path
import subprocess

sys.path.insert(0, str(Path(__file__).parent))

def run_test(arch, output_dir, steps=10000):
    """Run a quick test simulation"""
    
    print(f"\n{'='*70}")
    print(f"Testing {arch.upper()}...")
    print(f"{'='*70}")
    
    # Create modified script that uses specified arch
    test_script = f"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti
import time
from scripts.run_phase2_full import Phase2FullRunner

# Force architecture
ti.init(arch=ti.{arch})
print(f"Initialized: {{arch}}")

# Run short test
runner = Phase2FullRunner(
    config_path="aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml",
    output_dir="{output_dir}",
    max_steps={steps},
    seed=42
)

start = time.time()
runner.run()
elapsed = time.time() - start

print(f"\\n{'='*70}")
print(f"{arch.upper()} Test Complete!")
print(f"Time: {{elapsed/60:.1f}} minutes")
print(f"Speed: {{runner.max_steps/elapsed:.1f}} steps/second")
print(f"{'='*70}")
"""
    
    # Write and run
    test_file = Path("temp_test_arch.py")
    test_file.write_text(test_script)
    
    try:
        start_time = time.time()
        result = subprocess.run(
            [sys.executable, str(test_file)],
            capture_output=True,
            text=True,
            timeout=1800  # 30 min timeout
        )
        elapsed = time.time() - start_time
        
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        return elapsed, result.returncode == 0
        
    finally:
        if test_file.exists():
            test_file.unlink()

def main():
    print("=" * 70)
    print("CPU vs GPU Performance Test")
    print("=" * 70)
    print("\nThis will run 10,000 steps on both CPU and GPU to compare.")
    print("Estimated time: 10-30 minutes per test\n")
    
    results = {}
    
    # Test CPU
    try:
        cpu_time, cpu_success = run_test("cpu", "results/test_cpu_perf", steps=10000)
        results['CPU'] = {'time': cpu_time, 'success': cpu_success}
    except Exception as e:
        print(f"CPU test failed: {e}")
        results['CPU'] = {'time': None, 'success': False}
    
    # Test GPU
    try:
        gpu_time, gpu_success = run_test("cuda", "results/test_gpu_perf", steps=10000)
        results['GPU'] = {'time': gpu_time, 'success': gpu_success}
    except Exception as e:
        print(f"GPU test failed: {e}")
        results['GPU'] = {'time': None, 'success': False}
    
    # Summary
    print("\n" + "=" * 70)
    print("PERFORMANCE COMPARISON")
    print("=" * 70)
    
    for arch, data in results.items():
        if data['success'] and data['time']:
            speed = 10000 / data['time']
            eta_500k = (500000 / speed) / 3600  # hours
            print(f"\n{arch}:")
            print(f"  Time (10K steps): {data['time']/60:.1f} minutes")
            print(f"  Speed: {speed:.1f} steps/second")
            print(f"  ETA for 500K steps: {eta_500k:.1f} hours ({eta_500k/24:.1f} days)")
    
    if results['CPU']['success'] and results['GPU']['success']:
        if results['CPU']['time'] and results['GPU']['time']:
            ratio = results['GPU']['time'] / results['CPU']['time']
            print(f"\n{'='*70}")
            if ratio > 1:
                print(f"CPU is {ratio:.1f}x FASTER than GPU!")
                print("Recommendation: Use CPU mode")
            else:
                print(f"GPU is {1/ratio:.1f}x FASTER than CPU!")
                print("Recommendation: Use GPU mode")
            print("=" * 70)

if __name__ == "__main__":
    main()

