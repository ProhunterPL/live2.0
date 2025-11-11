# Phase 2B Cluster Detection Deadlock - Fix Instructions

## Problem Summary

Simulations are getting stuck in an infinite loop inside the cluster detection kernel (`update_clusters_kernel` in `backend/sim/core/binding.py`).

**Evidence:**
- ✅ Processes show high CPU usage (373-2239% each)
- ❌ Progress stopped hours ago (logs frozen at specific steps)
- 😴 Process state "Sl" (sleeping/waiting on thread synchronization)
- 🧵 196 threads per process (Taichi parallel execution)
- ⚠️ No new snapshots/checkpoints created for 6-41 hours

**Root Cause:**
The cluster detection uses an O(N²) nested loop that checks all particle pairs. When bond networks become complex (high connectivity), the Union-Find algorithm triggers pathological cases causing exponential slowdown or deadlock.

## Solution: Two-Step Fix

### Step 1: Kill Stuck Simulations

Run on AWS:

```bash
cd ~/live2.0

# Make script executable
chmod +x aws_test/scripts/kill_stuck_simulations.sh

# Kill stuck simulations (keeps run_9 which is still progressing)
bash aws_test/scripts/kill_stuck_simulations.sh

# Verify only run_9 remains
ps aux | grep "run_phase2_full.py" | grep -v grep
```

### Step 2: Apply Code Hotfix

The hotfix makes cluster detection optional via config:

```bash
cd ~/live2.0

# Make script executable
chmod +x aws_test/scripts/apply_cluster_fix.sh

# Apply the fix (creates backup automatically)
bash aws_test/scripts/apply_cluster_fix.sh

# Verify the fix
grep -A 10 "Update clusters - now configurable" backend/sim/core/stepper.py
```

**What the fix does:**
- Reads `cluster_check_interval` from config (defaults to 1200 if not specified)
- If `cluster_check_interval >= 999999999`, cluster detection is completely disabled
- This prevents the infinite loop while maintaining all other simulation features

### Step 3: Monitor run_9

Let run_9 finish to see if it also hits the deadlock:

```bash
# Monitor progress every 5 minutes
watch -n 300 "python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"

# Or check manually
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

**Expected behavior:**
- If run_9 completes successfully: cluster detection works for simple networks
- If run_9 also gets stuck around 300,000+ steps: need to disable cluster detection

### Step 4: Restart with Safer Configuration

After run_9 finishes (or gets stuck), restart with cluster detection disabled:

```bash
cd ~/live2.0

# Make script executable
chmod +x aws_test/scripts/restart_phase2b_safe.sh

# Start 9 new runs (to complete 30 total with run_1 + run_9)
# Seeds 110-118 (avoiding conflicts with previous runs)
bash aws_test/scripts/restart_phase2b_safe.sh 9 110
```

This uses the SAFER config which:
- Sets `cluster_check_interval: 999999999` (effectively disabled)
- Maintains all other simulation parameters
- Uses `python3 -u` for unbuffered logging (real-time progress visibility)

## Monitoring New Simulations

```bash
# Check overall progress
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py

# Check file-based progress (bypasses log buffering)
python3 ~/live2.0/aws_test/scripts/monitor_by_filesize.py

# Check specific run logs
tail -f ~/live2.0/logs/phase2b_safe/run_10.log

# Quick status
python3 ~/live2.0/aws_test/scripts/quick_diagnose.py
```

## Recovery Strategy

### Current Status (Before Fix)
- ✅ **1 completed**: run_1
- ⏳ **1 in progress**: run_9 (at 15%, ~6-8 hours remaining)
- 🚫 **7 stuck**: runs 2, 3, 4, 5, 6, 7, 8 (deadlocked for hours)

### After Fix
- ✅ **1 completed**: run_1
- ⏳ **1 in progress**: run_9 (monitoring)
- 🔄 **9 new runs**: runs 10-18 (with cluster detection disabled)

### Final Result
- **Target**: 30 simulations total for statistical significance
- **After recovery**: 10-11 completed (1 existing + 1 monitored + 9 new)
- **Still need**: ~19-20 more runs in subsequent batches

## Expected Timeline

With cluster detection disabled:
- **Each simulation**: 500K steps at ~140 steps/sec/core = ~60-90 minutes per run
- **9 parallel runs** on 96-core instance: ~60-90 minutes total
- **Total recovery time**: ~2 hours (including monitoring run_9)

## Verification Checklist

- [ ] Stuck simulations killed (only run_9 should remain)
- [ ] Code hotfix applied and verified
- [ ] run_9 monitored until completion or gets stuck
- [ ] New simulations started with SAFER config
- [ ] Progress monitoring showing recent file updates
- [ ] No processes stuck in 'Sl' state for >1 hour
- [ ] Snapshots being created every 50K steps
- [ ] Logs updating in real-time (with python3 -u)

## Troubleshooting

### If simulations still get stuck:

1. **Check process state:**
   ```bash
   ps aux | grep "run_phase2_full.py" | grep -v grep
   ```
   If state is 'D' → I/O issue
   If state is 'Sl' for hours → Still deadlocked (shouldn't happen with fix)

2. **Verify fix was applied:**
   ```bash
   grep "cluster_interval = getattr" backend/sim/core/stepper.py
   ```
   Should show the new configurable code

3. **Check config is correct:**
   ```bash
   grep "cluster_check_interval" aws_test/configs/phase2_miller_urey_extended_SAFER.yaml
   ```
   Should show: `cluster_check_interval: 999999999`

4. **Last resort - disable in code directly:**
   Edit `backend/sim/core/stepper.py` line 509 and comment out the entire cluster update block

## Files Created by This Fix

- `aws_test/patches/disable_cluster_detection.patch` - Patch file (for reference)
- `aws_test/scripts/apply_cluster_fix.sh` - Automated fix application
- `aws_test/scripts/kill_stuck_simulations.sh` - Kill deadlocked processes
- `aws_test/scripts/restart_phase2b_safe.sh` - Restart with safe config
- `aws_test/configs/phase2_miller_urey_extended_SAFER.yaml` - Config without cluster detection
- `docs/aws_test/CLUSTER_FIX_INSTRUCTIONS.md` - This file

## Contact

If issues persist, check:
1. System resources: `htop`, `free -h`, `df -h`
2. Taichi compilation: Look for LLVM errors in logs
3. Memory pressure: Might need to reduce `n_particles` further

Good luck! 🚀

