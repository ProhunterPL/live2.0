# Hydrothermal Queue - AWS Deployment Ready 🚀

**Status**: READY TO DEPLOY  
**Date**: 2025-11-19  
**Config**: `phase2_hydrothermal_AWS_OPTIMIZED.yaml`

---

## 📋 Pre-Flight Checklist

Before starting hydrothermal queue on AWS:

- [ ] ✅ Miller-Urey run_4 completed (17/17 runs done)
- [ ] ✅ Old hydrothermal processes killed (cleanup done)
- [ ] ✅ Local hydrothermal test validated (runs 10→1 working)
- [ ] ✅ Git pushed to main (latest code on AWS)
- [ ] ✅ Disk space available (~5-10 GB for 17 runs)

---

## 🚀 Quick Start Commands

### 1. Connect to AWS

```bash
ssh ubuntu@<YOUR_AWS_IP>
cd ~/live2.0
```

### 2. Pull Latest Code

```bash
git pull origin main
```

### 3. Verify Config

```bash
# Check config exists
ls -lh aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml

# Quick peek at critical settings
grep -E "n_particles|dt|cluster_check_interval|mutation" aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml
```

Expected output:
```
n_particles: 1000
dt: 0.002
cluster_check_interval: 999999999  # DISABLED
enable_mutations: false
```

### 4. Start Queue (Simple Method)

```bash
# Interactive launcher (recommended for first time)
bash aws_test/scripts/start_hydro_queue_aws.sh
```

**OR** Manual start:

```bash
# Make executable (if needed)
chmod +x aws_test/scripts/auto_queue_restart_hydro.sh
chmod +x aws_test/scripts/start_hydro_queue_aws.sh

# Start in background
nohup bash aws_test/scripts/auto_queue_restart_hydro.sh > logs/auto_restart_hydro.log 2>&1 &

echo "Queue started! PID: $!"
```

### 5. Verify Started

```bash
# Check queue process
ps aux | grep auto_queue_restart_hydro

# Check if simulations starting
sleep 30
ps aux | grep run_phase2_full | grep hydrothermal

# Check logs
tail -f logs/auto_restart_hydro_main.log
```

---

## 📊 Monitoring

### Quick Status Check

```bash
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

### Detailed Monitoring

```bash
# Main queue log (system status, starts/stops)
tail -f logs/auto_restart_hydro_main.log

# Progress log (high-level progress)
tail -f logs/auto_restart_hydro_progress.log

# Individual run logs
tail -f logs/hydro_run_1.log  # Replace 1 with run number

# Watch all running processes
watch -n 60 'ps aux | grep run_phase2_full | grep hydrothermal | grep -v grep'
```

### Progress Table

```bash
# Check completed runs
ls -lh results/phase2b_additional/hydrothermal_extended/*/results.json | wc -l

# Check snapshots
du -sh results/phase2b_additional/hydrothermal_extended/run_*
```

---

## ⏱️ Timeline Expectations

| Time | Event | Status |
|------|-------|--------|
| T+0h | Queue starts | 4 runs begin (1-4) |
| T+7h | First batch done | Start runs 5-8 |
| T+14h | Second batch done | Start runs 9-12 |
| T+21h | Third batch done | Start runs 13-16 |
| T+28h | Fourth batch done | Start run 17 |
| T+35h | **ALL COMPLETE** | 17/17 runs done ✅ |

**Total expected time**: ~35-45 hours (depends on chemistry complexity)

---

## 🛑 How to Stop

### Graceful Stop

```bash
# Stop queue (no new runs will start)
pkill -f auto_queue_restart_hydro.sh

# Wait for current runs to finish (recommended)
watch -n 60 'ps aux | grep run_phase2_full | grep hydrothermal | grep -v grep'
```

### Force Stop (Emergency)

```bash
# Kill queue
pkill -f auto_queue_restart_hydro.sh

# Kill all hydrothermal simulations
pkill -f 'run_phase2_full.py.*hydrothermal'

# Verify all stopped
ps aux | grep hydrothermal
```

---

## 🔧 Troubleshooting

### Problem: Queue not starting

```bash
# Check if already running
pgrep -f auto_queue_restart_hydro

# Check launcher log
cat logs/auto_restart_hydro_launcher.log

# Try manual start
bash aws_test/scripts/auto_queue_restart_hydro.sh
```

### Problem: Runs getting stuck

```bash
# Check which runs are stuck
python3 aws_test/scripts/check_actual_progress.py

# Find stuck process
ps aux | grep run_phase2_full | grep hydrothermal

# Kill stuck process (example: PID 12345)
kill -9 12345

# Queue will auto-restart the run
```

### Problem: Out of disk space

```bash
# Check space
df -h

# Clean old logs (if safe)
rm logs/hydro_run_*.log.old

# Archive completed runs
tar -czf hydro_runs_backup.tar.gz results/phase2b_additional/hydrothermal_extended/
```

### Problem: Config error

```bash
# Test config with dry-run (1000 steps)
python3 scripts/run_phase2_full.py \
    --config aws_test/configs/phase2_hydrothermal_AWS_OPTIMIZED.yaml \
    --output /tmp/test_hydro \
    --steps 1000 \
    --seed 999 \
    --force-cpu

# Check log
cat /tmp/test_hydro/simulation.log
```

---

## 📁 Output Structure

After completion, you'll have:

```
results/phase2b_additional/hydrothermal_extended/
├── run_1/
│   ├── results.json          # Summary statistics
│   ├── molecules.json        # Detected molecules
│   ├── simulation.log        # Detailed log
│   ├── snapshots/            # 10 snapshots (50K intervals)
│   │   ├── snapshot_50000.json
│   │   ├── snapshot_100000.json
│   │   └── ...
│   └── checkpoints/          # 4 checkpoints (100K intervals)
├── run_2/
├── ...
└── run_17/
```

---

## 🎯 Post-Completion

After all 17 runs complete:

### 1. Verify Results

```bash
# Check all runs completed
python3 aws_test/scripts/check_actual_progress.py

# Count results files
ls results/phase2b_additional/hydrothermal_extended/*/results.json | wc -l
# Should show: 17
```

### 2. Download to Local Machine

```bash
# On your local machine (Windows PowerShell):
scp -r ubuntu@<AWS_IP>:~/live2.0/results/phase2b_additional/hydrothermal_extended ./aws_results/
```

### 3. Analyze Results

```bash
# On AWS or locally
python3 aws_test/scripts/analyze_results.py --scenario hydrothermal
```

### 4. Clean Up AWS (Optional)

```bash
# Stop instance (to save costs)
# From AWS Console or CLI:
aws ec2 stop-instances --instance-ids i-XXXXXXXXX
```

---

## 📊 Expected Results

Based on local tests and miller_urey data:

- **Completion rate**: 17/17 runs (100%)
- **Average runtime**: 6-8 hours per run
- **Molecules detected**: 20-100 unique species per run
- **Bonds formed**: 50-200 bonds at final snapshot
- **File size**: ~1-2 MB per run (total ~20-35 MB)

---

## ✅ Success Criteria

Queue is successful if:

1. ✅ All 17 runs complete (results.json exists)
2. ✅ Each run has 10 snapshots
3. ✅ No stuck/zombie processes
4. ✅ Molecules detected in most runs (>0 bonds)
5. ✅ No LLVM crashes or deadlocks

---

## 🆘 Emergency Contacts

If something goes wrong:

1. **Check logs first**: `logs/auto_restart_hydro_main.log`
2. **Check progress**: `python3 aws_test/scripts/check_actual_progress.py`
3. **Review this doc**: Troubleshooting section above
4. **Last resort**: Kill all, cleanup, restart queue

---

## 📝 Notes

- **CPU mode**: Proven faster than GPU for hydrothermal chemistry
- **Particle count**: 1000 safe limit (tested locally)
- **Cluster detection**: DISABLED to prevent deadlock
- **Mutations**: DISABLED to prevent LLVM crash
- **Config source**: Based on `SUPER_LIGHT` with local test validation

---

**Good luck! 🚀**

*The queue will handle everything automatically - just monitor occasionally and let it run.*

