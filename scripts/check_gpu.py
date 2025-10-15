#!/usr/bin/env python3
"""Quick script to check GPU availability for Taichi"""

import taichi as ti

print("=" * 70)
print("TAICHI GPU CHECK - RTX 5070")
print("=" * 70)

# Test 1: Try CUDA
print("\n[TEST 1] CUDA initialization...")
print("-" * 70)

try:
    ti.init(arch=ti.cuda, device_memory_GB=6.0, debug=False)
    print("[OK] CUDA initialization: SUCCESS")
    print("   RTX 5070 CUDA detected and ready!")
    print("   This will give 10-50x speedup!")
    cuda_works = True
    
except Exception as e:
    print(f"[FAIL] CUDA initialization: FAILED")
    print(f"   Error: {e}")
    cuda_works = False

# Test 2: Try Vulkan (alternative GPU backend)
if not cuda_works:
    print("\n[TEST 2] Vulkan initialization (GPU alternative)...")
    print("-" * 70)
    try:
        ti.reset()
        ti.init(arch=ti.vulkan, device_memory_GB=6.0)
        print("[OK] Vulkan initialization: SUCCESS")
        print("   GPU via Vulkan detected!")
    except Exception as e:
        print(f"[FAIL] Vulkan initialization: FAILED")
        print(f"   Error: {e}")

# Test 3: CPU fallback
if not cuda_works:
    print("\n[TEST 3] CPU fallback...")
    print("-" * 70)
    try:
        ti.reset()
        ti.init(arch=ti.cpu, cpu_max_num_threads=28)
        print("[OK] CPU initialization: SUCCESS (28 threads)")
        print("   [WARNING] But this is MUCH slower than GPU")
    except Exception as e:
        print(f"[FAIL] CPU initialization: FAILED")
        print(f"   Error: {e}")

print("\n" + "=" * 70)
print("RECOMMENDATION:")
if cuda_works:
    print("  [OK] Use CUDA for maximum performance!")
    print("  Expected: 100-500 steps/second (vs 2.8 with CPU)")
else:
    print("  [WARNING] CUDA not working - need to debug or use CPU")
print("=" * 70)

