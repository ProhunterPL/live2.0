#!/usr/bin/env python3
"""
Quick script to check which Taichi backend is currently being used
"""

import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import taichi as ti

print("="*70)
print("TAICHI BACKEND DETECTION")
print("="*70)

# Try to detect what's available
print("\n1️⃣ Checking available backends...")

backends_available = []

# Check CUDA
print("\n  Testing CUDA (NVIDIA GPU)...")
try:
    ti.reset()
    ti.init(arch=ti.cuda, device_memory_GB=2.0)
    print("    ✅ CUDA is AVAILABLE")
    backends_available.append("cuda")
    ti.reset()
except Exception as e:
    print(f"    ❌ CUDA not available: {str(e)[:60]}")

# Check Vulkan
print("\n  Testing Vulkan (Generic GPU)...")
try:
    ti.reset()
    ti.init(arch=ti.vulkan, device_memory_GB=2.0)
    print("    ✅ Vulkan is AVAILABLE")
    backends_available.append("vulkan")
    ti.reset()
except Exception as e:
    print(f"    ❌ Vulkan not available: {str(e)[:60]}")

# CPU is always available
print("\n  Testing CPU...")
try:
    ti.reset()
    import multiprocessing
    max_threads = multiprocessing.cpu_count()
    ti.init(arch=ti.cpu, cpu_max_num_threads=max_threads)
    print(f"    ✅ CPU is AVAILABLE ({max_threads} threads)")
    backends_available.append("cpu")
    ti.reset()
except Exception as e:
    print(f"    ❌ CPU error: {str(e)[:60]}")

# Check what the server.py uses
print("\n2️⃣ Checking what backend/api/server.py will use...")
print("   (based on server.py lines 56-65)")

if "cuda" in backends_available:
    print("    → Will use: CUDA (GPU)")
    current_backend = "CUDA"
elif "vulkan" in backends_available:
    print("    → Will use: Vulkan (GPU)")
    current_backend = "Vulkan"
elif "cpu" in backends_available:
    print("    → Will use: CPU")
    current_backend = "CPU"
else:
    print("    → ERROR: No backend available!")
    current_backend = "NONE"

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"\n✅ Available backends: {', '.join(backends_available) if backends_available else 'NONE'}")
print(f"✅ Current backend: {current_backend}")

# Recommendations
print("\n" + "="*70)
print("RECOMMENDATIONS")
print("="*70)

if "cuda" in backends_available:
    print("\n🚀 You have CUDA (NVIDIA GPU)!")
    print("   For best performance, GPU is usually 10-50x faster than CPU")
    print("   Run the benchmark to confirm: python tests/benchmark_gpu_vs_cpu.py")
elif "vulkan" in backends_available:
    print("\n🚀 You have Vulkan (Generic GPU)!")
    print("   GPU should be faster than CPU for most operations")
    print("   Run the benchmark to confirm: python tests/benchmark_gpu_vs_cpu.py")
else:
    print("\n⚠️  No GPU detected, using CPU")
    print("   CPU can be fast with many threads, but typically slower than GPU")
    print("   Consider running on a machine with NVIDIA GPU for better performance")

print("\n💡 To run full benchmark:")
print("   PowerShell: .\\run_benchmark.ps1")
print("   Python: python tests/benchmark_gpu_vs_cpu.py")

print("\n" + "="*70)

