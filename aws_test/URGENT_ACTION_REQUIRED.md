# 🚨 URGENT ACTION REQUIRED - Phase 2B Deadlock

**Date**: 2025-11-13 21:28:02  
**Status**: 4 out of 5 processes are DEADLOCKED

---

## 📊 Current Situation

### ✅ Active Processes (1/5)
- **PID 79914**: State `Rsl` (Running) - ✅ **ACTIVE**
  - Likely one of: run_3, run_5, run_7, or run_8
  - Recent log activity (1-44 min ago)

### 🚨 Deadlocked Processes (4/5)
- **PID 73085**: State `Ssl`, CPU 995%, Running 8150:59 - **DEADLOCKED**
- **PID 76917**: State `Ssl`, CPU 1367%, Running 8023:24 - **DEADLOCKED**
- **PID 79713**: State `Ssl`, CPU 1241%, Running 5201:50 - **DEADLOCKED**
- **PID 80149**: State `Ssl`, CPU 1139%, Running 4739:48 - **DEADLOCKED**

**Pattern**: High CPU (>1000%) but sleeping state = **Cluster detection deadlock**

---

## 🔍 Analysis

### Active Runs (Recent Log Activity)
- **run_3**: Step 117,000 (22%) - Log 44 min ago ✅
- **run_5**: Step 215,000 (43%) - Log 1.2 min ago ✅
- **run_7**: Step 134,000 (27%) - Log 2.1 min ago ✅
- **run_8**: Step 134,000 (27%) - Log 1.8 min ago ✅

### Stuck Runs (Old Logs)
- **run_4**: Step 198,000 (40%) - Log 2+ days ago ❌
- **run_6**: Step 128,000 (26%) - Log 8+ hours ago ❌
- **run_10-18**: Step 160,000-165,000 (32-33%) - Logs 2+ days ago ❌ **DEADLOCKED AT 160K**

---

## ⚡ IMMEDIATE ACTIONS

### Step 1: Identify Which Runs Are Deadlocked

Run this on AWS:
```bash
cd ~/live2.0
bash aws_test/scripts/identify_process_runs.sh
```

This will show which PID corresponds to which run.

### Step 2: Kill Deadlocked Processes

**Option A: Smart Kill (Recommended)**
```bash
cd ~/live2.0
bash aws_test/scripts/kill_stuck_smart.sh
```

This will:
- ✅ Identify processes in `Ssl` state with high CPU
- ✅ Check log files to confirm they're stuck
- ✅ Kill only confirmed deadlocked processes
- ✅ Preserve actively running processes

**Option B: Manual Kill (If smart script doesn't work)**
```bash
# Kill specific PIDs (replace with actual stuck PIDs)
kill -9 73085 76917 79713 80149

# Or kill by run pattern (if you know which runs are stuck)
pkill -f "run_phase2_full.py.*run_1[0-8]"
pkill -f "run_phase2_full.py.*run_[46]"
```

### Step 3: Verify Fix Applied

```bash
cd ~/live2.0
python3 aws_test/scripts/diagnose_phase2b_issues.py
```

This checks:
- ✅ If cluster detection fix is in code
- ✅ If config files have `cluster_check_interval: 999999999`
- ✅ Current process count vs limit (5)

### Step 4: Restart Killed Runs

**IMPORTANT**: Only restart if you have <5 processes running!

```bash
cd ~/live2.0

# Check current count
CURRENT=$(ps aux | grep run_phase2_full.py | grep -v grep | wc -l)
echo "Currently running: $CURRENT processes"

# If < 5, restart using auto queue
if [ "$CURRENT" -lt 5 ]; then
    bash aws_test/scripts/auto_queue_restart.sh
else
    echo "⚠️  Already at limit (5 processes). Wait for some to complete."
fi
```

---

## 📈 Expected Outcome

### After Killing Deadlocked Processes:
- **Remaining**: 1-4 active processes (run_3, run_5, run_7, run_8)
- **Available slots**: 1-4 slots for new runs
- **Action**: Restart killed runs one by one as slots open

### Timeline:
- **T+0h**: Kill deadlocked processes → 1-4 processes remain
- **T+0h**: Restart 1-4 killed runs (up to limit of 5)
- **T+5h**: run_7 completes → slot opens → restart next
- **T+9h**: run_8 completes → slot opens → restart next
- **T+15h**: run_5 completes → slot opens → restart next
- **T+23h**: run_3 completes → slot opens → restart next

---

## ⚠️ Critical Warnings

1. **DO NOT kill PID 79914** - This is the only actively running process
2. **DO NOT restart run_3, run_5, run_7, run_8** - They're making progress
3. **Maintain ≤5 total processes** - Check count before restarting
4. **Verify fix before restarting** - Otherwise they'll deadlock again

---

## 🔍 Verification Commands

### Check Process States
```bash
ps aux | grep run_phase2_full.py | grep -v grep | awk '{print "PID", $2, "State:", $8, "CPU:", $3"%", $NF}'
```

### Check Which Runs Are Active
```bash
for run in 3 5 7 8; do
    echo "run_$run:"
    tail -5 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_$run/simulation.log | grep "Step"
done
```

### Check Which Runs Are Stuck
```bash
for run in 4 6 10 11 12 13 14 15 16 17 18; do
    echo "run_$run:"
    stat -c "%y" ~/live2.0/results/phase2b_additional/miller_urey_extended/run_$run/simulation.log 2>/dev/null || echo "No log"
done
```

---

## 📝 Notes

- **Process State "Ssl"**: Sleeping but with high CPU (>1000%) = deadlock
- **160K Steps**: Classic deadlock point (cluster detection runs every 1200 steps)
- **Only 1 process running**: This is critical - 4 processes are wasting CPU cycles
- **Limit 5 processes**: After killing, restart maintaining this limit

---

**Last Updated**: 2025-11-13 21:28:02  
**Next Action**: Run `kill_stuck_smart.sh` on AWS instance

