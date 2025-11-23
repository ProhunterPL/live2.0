---
date: 2025-11-19
label: guide
---

# Formamide Queue - AWS Deployment Ready 🚀

**Status**: READY TO DEPLOY  
**Date**: 2025-11-19  
**Config**: `phase2_formamide_AWS_OPTIMIZED.yaml`  
**Runs**: 8 (run_1 to run_8)

---

## 📋 Pre-Flight Checklist

Before starting formamide queue on AWS:

- [ ] ✅ Hydrothermal runs completed (all 17 runs done)
- [ ] ✅ Old formamide processes killed (cleanup done)
- [ ] ✅ Git pushed to main (latest code on AWS)
- [ ] ✅ Disk space available (~2-3 GB for 8 runs)
- [ ] ✅ Config verified (cluster_check_interval disabled, mutations disabled)

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
ls -lh aws_test/configs/phase2_formamide_AWS_OPTIMIZED.yaml

# Quick peek at critical settings
grep -E "n_particles|dt|cluster_check_interval|enable_mutations" aws_test/configs/phase2_formamide_AWS_OPTIMIZED.yaml
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
bash aws_test/scripts/start_formamide_queue_aws.sh
```

**OR** Manual start:

```bash
# Make executable (if needed)
chmod +x aws_test/scripts/auto_queue_restart_formamide.sh
chmod +x aws_test/scripts/start_formamide_queue_aws.sh

# Start in background
nohup bash aws_test/scripts/auto_queue_restart_formamide.sh > logs/auto_restart_formamide.log 2>&1 &

echo "Queue started! PID: $!"
```

### 5. Verify Started

```bash
# Check queue process
ps aux | grep auto_queue_restart_formamide

# Check if simulations starting
sleep 30
ps aux | grep run_phase2_full | grep formamide

# Check logs
tail -f logs/auto_restart_formamide_main.log
```

---

## 📊 Monitoring

### Quick Status Check

```bash
python3 ~/live2.0/aws_test/scripts/check_actual_progress.py
```

### Detailed Monitoring

```bash
# Main queue log
tail -f logs/auto_restart_formamide_main.log

# Progress log
tail -f logs/auto_restart_formamide_progress.log

# Specific run log
tail -f logs/formamide_run_1.log

# Running processes
ps aux | grep run_phase2_full | grep formamide
```

### Check Completion Status

```bash
# Count completed runs
ls -d results/phase2b_additional/formamide_extended/run_*/results.json 2>/dev/null | wc -l

# List all runs
ls -lh results/phase2b_additional/formamide_extended/
```

---

## ⚙️ Configuration Details

### Key Settings (Applied Fixes)

1. **Cluster Detection DISABLED**
   - `cluster_check_interval: 999999999`
   - **Why**: Prevents deadlock in complex bond networks (learned from Miller-Urey issues)

2. **Mutations DISABLED**
   - `enable_mutations: false`
   - **Why**: Prevents LLVM crash on CPU (learned from hydrothermal issues)

3. **Conservative Timestep**
   - `dt: 0.002`
   - **Why**: Stability for long runs (same as hydrothermal)

4. **Safe Particle Count**
   - `n_particles: 1000`
   - **Why**: Proven stable in previous tests

5. **CPU Optimized**
   - `grid_width: 96`, `grid_height: 96`
   - **Why**: Optimal for 28-core AWS instance

### Run Configuration

- **Total runs**: 8 (run_1 to run_8)
- **Seeds**: 100-107 (run_1=100, run_2=101, etc.)
- **Max parallel**: 4 simulations
- **Steps per run**: 500,000
- **Expected time**: 6-8 hours per run
- **Total time**: ~12-16 hours (with 4 parallel)

---

## 🛠️ Troubleshooting

### Queue Not Starting

```bash
# Check if already running
ps aux | grep auto_queue_restart_formamide

# Check launcher log
cat logs/auto_restart_formamide_launcher.log

# Verify script permissions
ls -lh aws_test/scripts/auto_queue_restart_formamide.sh
```

### Simulations Not Starting

```bash
# Check main log for errors
tail -50 logs/auto_restart_formamide_main.log

# Verify config exists
ls -lh aws_test/configs/phase2_formamide_AWS_OPTIMIZED.yaml

# Check Python script
python3 scripts/run_phase2_full.py --help
```

### Stuck Processes

```bash
# Check for stuck processes (no log update for 1+ hour)
# The queue script will warn about these automatically

# Manual check
ps aux | grep run_phase2_full | grep formamide

# Kill stuck process (if needed)
pkill -f "run_phase2_full.py.*formamide.*run_X"
```

### Out of Memory

```bash
# Check memory usage
free -h

# Reduce parallel runs (edit script)
# Change MAX_PARALLEL=4 to MAX_PARALLEL=2 in auto_queue_restart_formamide.sh
```

---

## 🛑 Stopping the Queue

### Graceful Stop

```bash
# Stop queue manager (won't kill running simulations)
pkill -f auto_queue_restart_formamide.sh

# Wait for current runs to finish, or kill them manually:
pkill -f "run_phase2_full.py.*formamide"
```

### Emergency Stop

```bash
# Kill everything
pkill -f formamide
pkill -f auto_queue_restart_formamide
```

---

## 📁 Output Structure

After completion, results will be in:

```
results/phase2b_additional/formamide_extended/
├── run_1/
│   ├── results.json
│   ├── molecules.json
│   ├── snapshots/
│   │   ├── snapshot_000000.json
│   │   ├── snapshot_050000.json
│   │   └── ... (10 snapshots total)
│   └── checkpoints/
├── run_2/
│   └── ...
└── ...
```

---

## ✅ Post-Completion Steps

1. **Verify All Runs Completed**
   ```bash
   ls -d results/phase2b_additional/formamide_extended/run_*/results.json | wc -l
   # Should show: 8
   ```

2. **Extract Molecules** (if needed)
   ```python
   from backend.sim.molecule_extractor import extract_molecules_from_results
   for run_num in range(1, 9):
       results = extract_molecules_from_results(
           f"results/phase2b_additional/formamide_extended/run_{run_num}"
       )
   ```

3. **Analyze Results**
   ```bash
   python3 aws_test/scripts/analyze_results.py
   ```

4. **Backup Results** (before terminating AWS instance)
   ```bash
   # SCP to local machine
   scp -r ubuntu@<AWS_IP>:~/live2.0/results/phase2b_additional/formamide_extended ./
   ```

---

## 📊 Expected Timeline

With 4 parallel runs:

- **T+0h**: Runs 1-4 start
- **T+6-8h**: Runs 1-4 complete → Runs 5-8 start
- **T+12-16h**: All 8 runs complete ✅

---

## 🎯 Key Differences from Hydrothermal

1. **Fewer runs**: 8 vs 17 (faster completion)
2. **Same stability**: Uses same conservative settings
3. **Same fixes**: Cluster detection disabled, mutations disabled
4. **Different chemistry**: Formamide-rich environment vs hydrothermal

---

## 📚 Related Documentation

- **Hydrothermal Guide**: `docs/aws_test/HYDRO_QUEUE_AWS_READY.md`
- **Cluster Deadlock Fix**: `docs/CLUSTER_DEADLOCK_FIX.md`
- **Auto-Restart Guide**: `docs/aws_test/AUTO_RESTART_GUIDE.md`

---

**Last Updated**: 2025-11-19  
**Status**: Ready for deployment after hydrothermal completes

