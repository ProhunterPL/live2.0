# Phase 2B Fix - Quick Reference Card

## 🚀 ONE-LINE DEPLOY

### On Windows (local):
```powershell
.\copy_fix_to_aws.ps1
```

### On AWS:
```bash
ssh ubuntu@ip-172-31-0-42
cd ~/live2.0
bash aws_test/DEPLOY_FIX_NOW.sh
```

---

## 📊 MONITORING COMMANDS

```bash
# Full status check
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py

# File-based monitoring (bypasses log buffering)
python3 ~/live2.0/aws_test/scripts/monitor_by_filesize.py

# Watch real-time (updates every 5 min)
watch -n 300 "python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"

# Check running processes
ps aux | grep "run_phase2_full.py" | grep -v grep

# Quick diagnose
python3 ~/live2.0/aws_test/scripts/quick_diagnose.py
```

---

## 📝 LOG COMMANDS

```bash
# Check new simulation logs (real-time)
tail -f ~/live2.0/logs/phase2b_safe/run_10.log

# Check all new logs
ls -lh ~/live2.0/logs/phase2b_safe/

# Check old simulation (run_9)
tail -f ~/live2.0/results/phase2b_additional/miller_urey_extended/run_9/simulation.log

# Find errors in all logs
grep -i error ~/live2.0/logs/phase2b_safe/*.log
```

---

## 🔧 MANUAL FIX STEPS

If automatic deployment fails:

```bash
# 1. Kill stuck processes
bash aws_test/scripts/kill_stuck_simulations.sh

# 2. Apply code fix
bash aws_test/scripts/apply_cluster_fix.sh

# 3. Verify fix applied
grep "cluster_interval = getattr" backend/sim/core/stepper.py

# 4. Start new simulations
bash aws_test/scripts/restart_phase2b_safe.sh 9 110
```

---

## 🩺 HEALTH CHECK

### Healthy simulation shows:
- ✅ Process state: `R` or `Sl`
- ✅ Logs updating every 1-2 minutes
- ✅ Snapshots created every 50K steps
- ✅ CPU usage 100-300% per process
- ✅ Files modified within last hour

### Stuck simulation shows:
- ❌ Process state: `Sl` for >1 hour
- ❌ Logs frozen for >1 hour
- ❌ High CPU but no new files
- ❌ No snapshots for >2 hours
- ❌ Last file modified >6 hours ago

**Check with:**
```bash
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

---

## 🛑 EMERGENCY COMMANDS

```bash
# Kill ALL simulations
pkill -9 -f "run_phase2_full.py"

# Kill specific run
kill -9 <PID>

# Check if anything is stuck in D state (I/O wait)
ps aux | grep "run_phase2_full.py" | grep " D"

# Check system resources
htop
free -h
df -h

# Check Taichi version
python3 -c "import taichi; print(taichi.__version__)"
```

---

## 📁 KEY FILES

### Scripts
- `aws_test/DEPLOY_FIX_NOW.sh` - Automatic deployment
- `aws_test/scripts/kill_stuck_simulations.sh` - Kill deadlocked processes
- `aws_test/scripts/apply_cluster_fix.sh` - Apply code hotfix
- `aws_test/scripts/restart_phase2b_safe.sh` - Restart with safe config
- `aws_test/scripts/check_actual_progress.py` - Real progress monitor
- `copy_fix_to_aws.ps1` - Deploy from Windows

### Configs
- `aws_test/configs/phase2_miller_urey_extended_SAFER.yaml` - Safe config (no cluster detection)

### Docs
- `aws_test/QUICK_REFERENCE.md` - This file
- `aws_test/EMERGENCY_SUMMARY.md` - Full problem description
- `aws_test/PODSUMOWANIE_NAPRAWY_PL.md` - Polish version
- `aws_test/CLUSTER_FIX_INSTRUCTIONS.md` - Detailed instructions

---

## ⏱️ EXPECTED TIMELINE

| Phase | Duration | Result |
|-------|----------|--------|
| Copy files to AWS | 1 min | Files deployed |
| Kill stuck processes | 1 min | 7 processes killed |
| Apply code fix | 1 min | Code patched |
| Monitor run_9 | 1-2 hours | 1 more completed |
| Run 9 new simulations | 60-90 min | 9 completed |
| **TOTAL** | **2-3 hours** | **10-11 completed** |

---

## ✅ SUCCESS CHECKLIST

- [ ] Files copied to AWS
- [ ] Stuck processes killed
- [ ] Code fix applied and verified
- [ ] run_9 completed or monitored
- [ ] 9 new simulations started
- [ ] Logs showing real-time progress
- [ ] Snapshots being created
- [ ] No processes stuck >1 hour
- [ ] All processes in healthy state

---

## 🆘 TROUBLESHOOTING

### Problem: Scripts fail to execute
**Solution:**
```bash
chmod +x aws_test/DEPLOY_FIX_NOW.sh
chmod +x aws_test/scripts/*.sh
```

### Problem: Python imports fail
**Solution:**
```bash
cd ~/live2.0
export PYTHONPATH=$PWD:$PYTHONPATH
```

### Problem: Simulations still getting stuck
**Solution:**
1. Verify fix: `grep "cluster_interval" backend/sim/core/stepper.py`
2. Check config: `grep "cluster_check_interval" aws_test/configs/phase2_miller_urey_extended_SAFER.yaml`
3. Nuclear option: Comment out lines 509-514 in `backend/sim/core/stepper.py`

### Problem: Out of memory
**Solution:**
```bash
free -h
# If low memory, reduce parallel runs:
bash aws_test/scripts/restart_phase2b_safe.sh 4 110  # Only 4 at a time
```

---

## 📞 SUPPORT

**If everything fails:**
1. Save current state: `tar -czf ~/phase2b_backup_$(date +%Y%m%d).tar.gz ~/live2.0/results/`
2. Reboot instance: `sudo reboot`
3. Re-deploy fix after reboot

**Check system health:**
```bash
# System load
uptime

# CPU usage
top -bn1 | head -20

# Memory usage
free -h

# Disk usage
df -h

# Recent errors
dmesg | tail -50
```

---

## 🎯 FINAL GOAL

**Target**: 30 completed simulations for Phase 2B statistical analysis

**Current**: 1 completed (run_1) + 1 in progress (run_9) + 7 stuck (to be killed)

**After fix**: 10-11 completed → Need 19-20 more in next batch

**Total time to 30**: ~5-6 hours from now

---

**Remember**: The fix is solid. Cluster detection is just for metrics, doesn't affect chemistry. Trust the process! 🚀

