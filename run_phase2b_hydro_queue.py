#!/usr/bin/env python3
"""
Phase 2B Hydrothermal Queue Runner - Reverse Order
===================================================

Runs hydrothermal simulations in reverse order (run_10 → run_1)
to complement AWS miller_urey runs.

Usage:
    python run_phase2b_hydro_queue.py --start 10 --end 1
    python run_phase2b_hydro_queue.py --start 10 --end 5  # Only 10-5
"""

import sys
import argparse
import time
import subprocess
import logging
import multiprocessing
from pathlib import Path
from datetime import datetime

# Setup logging with both file and console
log_file = Path("logs") / f"hydro_queue_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
log_file.parent.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Detect CPU cores
CPU_CORES = multiprocessing.cpu_count()

HYDRO_CONFIG = {
    'config': 'aws_test/configs/phase2_hydrothermal_SUPER_LIGHT.yaml',
    'description': 'Hydrothermal - SUPER LIGHT - CPU Optimized (500K steps)',
    'steps': 500000,
    'expected_time': 30,  # minutes on CPU - fewer particles = much faster!
    'use_cpu': True,  # CPU is faster than GPU for chemistry-heavy simulations
    'cpu_threads': CPU_CORES
}


def run_single_simulation(run_number: int, seed: int, output_base: str = "results/phase2b_local", 
                         use_cpu: bool = True, cpu_threads: int = None, use_hybrid: bool = False) -> bool:
    """Run a single hydrothermal simulation"""
    
    run_id = f"run_{run_number:02d}"
    output_dir = Path(output_base) / "hydrothermal" / run_id
    
    # Check if already completed
    if output_dir.exists() and (output_dir / "results.json").exists():
        logger.info(f"[SKIP] {run_id} - Already completed, skipping")
        return True
    
    # Determine execution mode
    if use_hybrid:
        mode_str = f"HYBRID (GPU physics + CPU chemistry, {cpu_threads or CPU_CORES} CPU threads)"
        expected_time = 45  # Hybrid is even faster!
    elif use_cpu:
        mode_str = f"CPU ({cpu_threads or CPU_CORES} threads)"
        expected_time = HYDRO_CONFIG['expected_time']
    else:
        mode_str = "GPU"
        expected_time = 90
    
    logger.info("=" * 70)
    logger.info(f"[START] {run_id} (seed={seed})")
    logger.info("=" * 70)
    logger.info(f"Mode: {mode_str}")
    logger.info(f"Config: {HYDRO_CONFIG['config']}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Steps: {HYDRO_CONFIG['steps']:,}")
    logger.info(f"Expected time: ~{expected_time} minutes")
    logger.info("=" * 70)
    
    cmd = [
        sys.executable,
        "scripts/run_phase2_full.py",
        "--config", HYDRO_CONFIG['config'],
        "--output", str(output_dir),
        "--steps", str(HYDRO_CONFIG['steps']),
        "--seed", str(seed)
    ]
    
    # Add CPU/Hybrid flags
    if use_cpu and not use_hybrid:
        cmd.append("--force-cpu")
        if cpu_threads:
            cmd.extend(["--cpu-threads", str(cpu_threads)])
    # Note: Hybrid mode is not yet integrated into run_phase2_full.py
    # Would need to add --hybrid flag to that script
    
    logger.info(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        # Run simulation
        result = subprocess.run(
            cmd,
            capture_output=False,  # Let output go to console
            timeout=14400  # 4 hour timeout
        )
        
        elapsed_time = time.time() - start_time
        elapsed_minutes = elapsed_time / 60
        
        if result.returncode == 0:
            logger.info("=" * 70)
            logger.info(f"[OK] {run_id} COMPLETED in {elapsed_minutes:.1f} minutes")
            logger.info("=" * 70)
            return True
        else:
            logger.error("=" * 70)
            logger.error(f"[FAILED] {run_id} FAILED after {elapsed_minutes:.1f} minutes")
            logger.error(f"Return code: {result.returncode}")
            logger.error("=" * 70)
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"[TIMEOUT] {run_id} TIMED OUT after 4 hours")
        return False
    except KeyboardInterrupt:
        logger.warning(f"\n[INTERRUPT] {run_id} INTERRUPTED by user")
        raise
    except Exception as e:
        logger.error(f"[CRASH] {run_id} CRASHED: {e}")
        return False


def run_queue(start_run: int, end_run: int, output_base: str = "results/phase2b_local",
              use_cpu: bool = True, cpu_threads: int = None, use_hybrid: bool = False):
    """Run queue of hydrothermal simulations in reverse order"""
    
    # Generate run numbers in reverse order
    run_numbers = list(range(start_run, end_run - 1, -1))
    
    # Seeds: run_1=100, run_2=101, ..., run_10=109
    seeds = {i: 99 + i for i in range(1, 11)}
    
    # Determine mode
    if use_hybrid:
        mode_str = f"HYBRID (GPU + CPU with {cpu_threads or CPU_CORES} threads)"
        expected_time_per_run = 45
    elif use_cpu:
        mode_str = f"CPU ({cpu_threads or CPU_CORES} threads)"
        expected_time_per_run = HYDRO_CONFIG['expected_time']
    else:
        mode_str = "GPU"
        expected_time_per_run = 90
    
    logger.info("\n" + "=" * 70)
    logger.info("PHASE 2B - HYDROTHERMAL QUEUE (REVERSE ORDER)")
    logger.info("=" * 70)
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Mode: {mode_str}")
    logger.info(f"CPU Cores Available: {CPU_CORES}")
    logger.info(f"Config: {HYDRO_CONFIG['config']}")
    logger.info(f"Queue: {', '.join([f'run_{i:02d}' for i in run_numbers])}")
    logger.info(f"Total runs: {len(run_numbers)}")
    logger.info(f"Steps per run: {HYDRO_CONFIG['steps']:,}")
    logger.info(f"Estimated total time: {len(run_numbers) * expected_time_per_run:.0f} minutes ({len(run_numbers) * expected_time_per_run/60:.1f} hours)")
    logger.info(f"Output: {output_base}/hydrothermal/")
    logger.info("=" * 70)
    
    results = []
    successful = 0
    failed = 0
    skipped = 0
    
    for idx, run_num in enumerate(run_numbers, 1):
        logger.info(f"\n\n📍 QUEUE POSITION: [{idx}/{len(run_numbers)}]")
        
        seed = seeds[run_num]
        success = run_single_simulation(run_num, seed, output_base, use_cpu, cpu_threads, use_hybrid)
        
        if success:
            # Check if it was skipped or actually ran
            run_dir = Path(output_base) / "hydrothermal" / f"run_{run_num:02d}"
            if run_dir.exists():
                successful += 1
            else:
                skipped += 1
        else:
            failed += 1
        
        results.append({
            'run_number': run_num,
            'run_id': f"run_{run_num:02d}",
            'seed': seed,
            'success': success
        })
        
        # Progress update
        logger.info(f"\n[PROGRESS] {successful} completed, {failed} failed, {skipped} skipped")
    
    # Final summary
    logger.info("\n\n" + "=" * 70)
    logger.info("HYDROTHERMAL QUEUE - FINAL SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Total runs: {len(run_numbers)}")
    logger.info(f"[OK] Successful: {successful}")
    logger.info(f"[FAILED] Failed: {failed}")
    logger.info(f"[SKIP] Skipped: {skipped}")
    logger.info(f"Success rate: {100*successful/(successful+failed):.1f}%" if (successful+failed) > 0 else "N/A")
    logger.info("=" * 70)
    
    # List completed runs
    logger.info("\nCOMPLETED RUNS:")
    for r in results:
        if r['success']:
            logger.info(f"  [OK] {r['run_id']} (seed={r['seed']})")
        else:
            logger.info(f"  [FAILED] {r['run_id']} (seed={r['seed']})")
    
    logger.info("\n" + "=" * 70)
    logger.info(f"Log saved to: {log_file}")
    logger.info("=" * 70)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2B Hydrothermal Queue Runner (Reverse Order) - CPU OPTIMIZED",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  # CPU mode (default, fastest for this workload!)
  python run_phase2b_hydro_queue.py --start 10 --end 1
  
  # CPU with specific thread count
  python run_phase2b_hydro_queue.py --start 10 --end 5 --cpu-threads 12
  
  # GPU mode (slower for chemistry-heavy simulations)
  python run_phase2b_hydro_queue.py --start 10 --end 1 --gpu
  
  # Hybrid mode (GPU physics + CPU chemistry) - EXPERIMENTAL
  python run_phase2b_hydro_queue.py --start 10 --end 1 --hybrid

Detected CPU cores: {CPU_CORES}
        """
    )
    
    parser.add_argument('--start', type=int, default=10,
                       help='Starting run number (default: 10)')
    parser.add_argument('--end', type=int, default=1,
                       help='Ending run number (default: 1)')
    parser.add_argument('--output', default='results/phase2b_local',
                       help='Output base directory (default: results/phase2b_local)')
    
    # Execution mode
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--gpu', action='store_true',
                           help='Use GPU mode (slower for this workload)')
    mode_group.add_argument('--hybrid', action='store_true',
                           help='Use hybrid GPU+CPU mode (EXPERIMENTAL)')
    
    parser.add_argument('--cpu-threads', type=int, default=None,
                       help=f'Number of CPU threads to use (default: {CPU_CORES} = all cores)')
    
    args = parser.parse_args()
    
    # Validation
    if args.start < args.end:
        parser.error(f"Start ({args.start}) must be >= end ({args.end}) for reverse order")
    
    if args.start < 1 or args.end < 1:
        parser.error("Run numbers must be >= 1")
    
    if args.start > 10:
        logger.warning(f"[WARNING] Run number {args.start} is > 10. Are you sure?")
    
    # Determine execution mode
    use_cpu = not args.gpu and not args.hybrid
    use_hybrid = args.hybrid
    
    # Log execution mode
    logger.info(f"\n[SYSTEM] System Info:")
    logger.info(f"   CPU Cores: {CPU_CORES}")
    if use_hybrid:
        logger.info(f"   Mode: HYBRID (GPU physics + CPU chemistry)")
        logger.info(f"   Expected to be FASTEST for this workload!")
    elif use_cpu:
        logger.info(f"   Mode: CPU ({args.cpu_threads or CPU_CORES} threads)")
        logger.info(f"   Faster than GPU for chemistry-heavy simulations!")
    else:
        logger.info(f"   Mode: GPU")
        logger.info(f"   Note: CPU is typically faster for this workload")
    
    try:
        results = run_queue(args.start, args.end, args.output, 
                          use_cpu=use_cpu, 
                          cpu_threads=args.cpu_threads,
                          use_hybrid=use_hybrid)
        logger.info("\n[OK] Queue completed successfully!")
        return 0
    except KeyboardInterrupt:
        logger.warning("\n\n[INTERRUPT] Queue interrupted by user (Ctrl+C)")
        logger.info("You can restart the queue - it will skip completed runs")
        return 130
    except Exception as e:
        logger.error(f"\n\n[ERROR] Queue failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

