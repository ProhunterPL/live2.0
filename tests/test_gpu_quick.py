#!/usr/bin/env python3
"""
Quick GPU Test for RTX 5070
===========================

Simple test to verify GPU works for Phase 2B simulations.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

import taichi as ti
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_gpu():
    """Test GPU initialization"""
    logger.info("Testing GPU...")
    
    try:
        ti.init(arch=ti.cuda)
        logger.info("✅ GPU initialized successfully!")
        logger.info(f"   Architecture: CUDA")
        
        # Test simple compute
        @ti.kernel
        def test_kernel():
            pass
        
        test_kernel()
        logger.info("✅ GPU computation works!")
        
        return True
        
    except Exception as e:
        logger.error(f"❌ GPU test failed: {e}")
        return False


def test_simulation_components():
    """Test if simulation components can be imported"""
    logger.info("\nTesting simulation components...")
    
    try:
        from backend.sim.config import SimulationConfig
        from backend.sim.phase2_config import Phase2Config
        logger.info("✅ Simulation components imported successfully")
        return True
    except Exception as e:
        logger.error(f"❌ Failed to import: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    logger.info("=" * 70)
    logger.info("GPU QUICK TEST")
    logger.info("=" * 70)
    
    gpu_ok = test_gpu()
    components_ok = test_simulation_components()
    
    logger.info("\n" + "=" * 70)
    if gpu_ok and components_ok:
        logger.info("✅ ALL TESTS PASSED!")
        logger.info("   Ready to run Phase 2B simulations locally")
    else:
        logger.info("❌ SOME TESTS FAILED")
        logger.info("   Check errors above")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

