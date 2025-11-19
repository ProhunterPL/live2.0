---
date: 2025-11-19
label: guide
---

# Hydrothermal Deployment - Summary 🚀

**Date**: 2025-11-19  
**Status**: ✅ READY FOR AWS DEPLOYMENT  
**Next Action**: Start when miller_urey run_4 completes

---

## 📦 What Was Prepared

### 1. AWS Production Config
**File**: `aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml`

**Key features**:
- ✅ Safe particle count: 1000 (tested locally)
- ✅ Cluster detection DISABLED (prevents deadlock)
- ✅ Mutations DISABLED (prevents LLVM crash)
- ✅ Conservative timestep: dt=0.002 (stability)
- ✅ Output to: `results/phase2b_additional/hydrothermal_extended/`
- ✅ Based on local SUPER_LIGHT config (proven stable)

### 2. Auto-Restart Queue System
**File**: `aws_test/scripts/auto_queue_restart_hydro.sh`

**Features**:
- 🔄 Auto-restart completed runs
- 📊 Max 4 parallel simulations
- 🔍 Progress monitoring every 5 minutes
- 🚨 Stuck process detection
- 📝 Detailed logging
- 🎯 Runs 1-17 (seeds 100-116)

### 3. Quick Launcher
**File**: `aws_test/scripts/start_hydro_queue_aws.sh`

**Features**:
- ✅ Pre-flight checks
- 📊 System info display
- ❓ Interactive confirmation
- 🚀 One-command start
- 📋 Post-start monitoring commands

### 4. Documentation
- **`aws_test/HYDRO_QUEUE_AWS_READY.md`**: Complete deployment guide
- **`aws_test/DEPLOY_HYDRO_COMMANDS.txt`**: Copy-paste commands

---

## 🎯 Local vs AWS Strategy

### Local (Your Machine)
- **Runs**: 10 → 1 (reverse order)
- **Config**: `phase2_hydrothermal_SUPER_LIGHT.yaml`
- **Script**: `run_phase2b_hydro_queue.py`
- **Started**: Already running
- **Purpose**: Quick results while AWS finishes miller_urey

### AWS (Cloud)
- **Runs**: 1 → 17 (forward order)
- **Config**: `phase2_hydrothermal_AWS_OPTIMIZED.yaml` (same as SUPER_LIGHT)
- **Script**: `auto_queue_restart_hydro.sh`
- **Start when**: miller_urey run_4 completes (~6 hours)
- **Purpose**: Complete Phase 2B dataset

**Combined total**: 27 hydrothermal runs (17 AWS + 10 local)

---

## ⏱️ Timeline

| Time | Event | Location | Action |
|------|-------|----------|--------|
| **T+0h** (now) | Local hydro started | Local | Running (10→1) |
| **T+6h** | Miller run_4 done | AWS | Start hydro queue (1→17) |
| **T+10h** | Local hydro done | Local | Download results |
| **T+41h** | AWS hydro done | AWS | Download results |
| **T+42h** | **PHASE 2B COMPLETE** | 🎉 | 17 miller + 27 hydro |

---

## 🚀 Deployment Steps (AWS)

### When to Start
✅ **After** miller_urey run_4 completes (currently at 84%)

### Commands (Copy-Paste)

```bash
# 1. Connect to AWS
ssh ubuntu@<YOUR_AWS_IP>

# 2. Navigate and update
cd ~/live2.0
git pull origin main

# 3. Start queue (interactive)
bash aws_test/scripts/start_hydro_queue_aws.sh

# 4. Monitor
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

**Full command list**: See `aws_test/DEPLOY_HYDRO_COMMANDS.txt`

---

## 📊 Expected Results

### Per-Run Expectations
- **Runtime**: 6-8 hours
- **Snapshots**: 10 (every 50K steps)
- **Checkpoints**: 5 (every 100K steps)
- **File size**: 1-2 MB
- **Molecules**: 20-100 unique species
- **Bonds**: 50-200 at completion

### Total Expectations (17 runs)
- **Total time**: 35-45 hours
- **Total size**: 20-35 MB
- **Success rate**: 17/17 (100%)
- **Completion**: ~T+41h from start

---

## 🔍 Monitoring

### Quick Check (Every Few Hours)
```bash
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

### Detailed Logs
```bash
# Main queue log
tail -f logs/auto_restart_hydro_main.log

# Individual run
tail -f logs/hydro_run_1.log
```

### Progress Indicators
✅ **Good signs**:
- 4 processes running
- Log files growing
- Snapshots appearing every ~1 hour
- No errors in logs

⚠️ **Warning signs**:
- Processes stuck (no log updates >1 hour)
- "Too many particles" errors
- LLVM crash messages
- Out of memory errors

---

## ✅ Success Criteria

Phase 2B Hydrothermal is successful when:

1. ✅ All 17 AWS runs complete (results.json exists)
2. ✅ All 10 local runs complete
3. ✅ Each run has 10 snapshots
4. ✅ Molecules detected (>0 bonds in most runs)
5. ✅ No deadlocks or crashes
6. ✅ Total: 27 hydrothermal runs for statistical analysis

---

## 🆘 If Something Goes Wrong

### Queue Not Starting
```bash
# Check if already running
pgrep -f auto_queue_restart_hydro

# Check logs
cat logs/auto_restart_hydro_launcher.log

# Try manual start
bash aws_test/scripts/auto_queue_restart_hydro.sh
```

### Run Getting Stuck
```bash
# Find stuck process
ps aux | grep run_phase2_full | grep hydrothermal

# Kill it (queue will restart)
kill -9 <PID>
```

### Config Error
```bash
# Test with 1000 steps
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml \
    --output /tmp/test \
    --steps 1000 \
    --seed 999 \
    --force-cpu
```

---

## 📝 Files Changed/Created

### New Files
1. `aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml` - Production config
2. `aws_test/scripts/auto_queue_restart_hydro.sh` - Queue manager
3. `aws_test/scripts/start_hydro_queue_aws.sh` - Launcher
4. `aws_test/HYDRO_QUEUE_AWS_READY.md` - Deployment guide
5. `aws_test/DEPLOY_HYDRO_COMMANDS.txt` - Quick commands
6. `HYDRO_DEPLOYMENT_SUMMARY.md` - This file

### Existing Files (Unchanged)
- `run_phase2b_hydro_queue.py` - Local runner (already working)
- `start_hydro_queue.ps1` - Local launcher (already working)
- `aws_test/configs/phase2_hydrothermal_SUPER_LIGHT.yaml` - Local config

---

## 🎓 Key Differences: Local vs AWS Config

Both configs are **almost identical**, just different output paths:

| Setting | SUPER_LIGHT (local) | AWS_OPTIMIZED (AWS) |
|---------|---------------------|---------------------|
| particles | 1000 | 1000 |
| dt | 0.002 | 0.002 |
| grid | 96×96 | 96×96 |
| cluster_check | disabled | disabled |
| mutations | disabled | disabled |
| output | `phase2b_local` | `phase2b_additional` |

**Why two configs?** Different output directories (local vs AWS results structure)

---

## 🎯 After Phase 2B Completes

### Data You'll Have
- ✅ 17 miller_urey runs (500K steps each)
- ✅ 27 hydrothermal runs (17 AWS + 10 local, 500K steps each)
- ✅ Total: 44 runs × 500K = 22 million steps simulated
- ✅ ~2-4 GB total data
- ✅ Statistical validation complete

### Next Steps
1. Download all results from AWS
2. Run comprehensive analysis
3. Compare miller_urey vs hydrothermal chemistry
4. Write Phase 3 paper sections
5. Generate publication figures

---

## 📚 Reference Documentation

- **AWS Setup**: `aws_test/AUTO_RESTART_GUIDE.md`
- **Cluster Fix**: `CLUSTER_DEADLOCK_FIX.md`
- **Project Rules**: `.cursorrules` (main project doc)
- **Deployment Guide**: `aws_test/HYDRO_QUEUE_AWS_READY.md`
- **Quick Commands**: `aws_test/DEPLOY_HYDRO_COMMANDS.txt`

---

## ✨ Summary

**What's Ready**:
- ✅ Config tested and optimized
- ✅ Queue system implemented
- ✅ Documentation complete
- ✅ Safety checks in place
- ✅ Monitoring tools ready

**What to Do**:
1. ⏰ Wait for miller_urey run_4 (~6h)
2. 🔌 Connect to AWS
3. 🚀 Run: `bash aws_test/scripts/start_hydro_queue_aws.sh`
4. 📊 Monitor occasionally
5. ☕ Let it run (~35-45h)

**You're all set! Just wait for miller_urey to finish.** 🎉

---

*Good luck with Phase 2B! 🧪🔬*

