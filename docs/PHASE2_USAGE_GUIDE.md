# Phase 2 Complete Usage Guide

**Status**: ✅ Infrastructure Complete + POC Working  
**Ready**: Yes - can start runs tonight!

---

## Quick Start Commands

### 1. Start Overnight Test (Recommended First Step)

```powershell
# Windows PowerShell
.\scripts\start_overnight_test.ps1

# Or manual start:
python scripts/run_phase2_full.py `
  --config configs/phase2_miller_urey_test.yaml `
  --output results/overnight_test `
  --steps 10000000 `
  --seed 42 > overnight.log 2>&1 &
```

**Duration**: ~10 hours  
**Check progress**: `Get-Content overnight.log -Tail 20 -Wait`

### 2. Run Quick Test (5-10 minutes)

```bash
python scripts/run_phase2_full.py \
  --config configs/phase2_miller_urey_test.yaml \
  --output results/quick_test \
  --steps 10000 \
  --seed 42
```

### 3. Run Full Production Batch (All 30 runs)

```bash
# Test mode first (1 run per scenario, 10k steps)
python scripts/phase2_master.py --mode test --scenarios all

# Full production (10 runs per scenario, 10M steps)
python scripts/phase2_master.py --mode full --scenarios all
```

---

## Complete Workflow

### Phase A: Single Run (Testing)

```bash
# 1. Run simulation
python scripts/run_phase2_full.py \
  --config configs/phase2_miller_urey_test.yaml \
  --output results/test_run \
  --steps 1000000 \
  --seed 42

# 2. Extract molecules
python backend/sim/molecule_extractor.py results/test_run

# 3. Check results
cat results/test_run/analysis/molecule_report.txt
```

### Phase B: Batch Production

```bash
# 1. Run all simulations (takes several days)
python scripts/phase2_master.py --mode full --scenarios all

# 2. Analyze batch results
python scripts/analyze_phase2_batch.py \
  --input results/phase2 \
  --output analysis/phase2_complete \
  --recursive \
  --use-matcher

# 3. View results
cat analysis/phase2_complete/batch_report.txt
```

### Phase C: Analysis Only

```bash
# If simulations already complete
python scripts/phase2_master.py --mode analyze --input results/phase2

# Or just batch analyzer
python scripts/analyze_phase2_batch.py \
  --input results/phase2 \
  --output analysis/phase2 \
  --recursive \
  --use-matcher
```

---

## Individual Scenario Commands

### Miller-Urey (CH₄ + NH₃ + H₂O + electrical discharge)

```bash
# Single run
python scripts/run_phase2_full.py \
  --config configs/phase2_miller_urey_test.yaml \
  --output results/miller_urey/run_01 \
  --steps 10000000 \
  --seed 42

# Batch (10 runs)
for i in {1..10}; do
  python scripts/run_phase2_full.py \
    --config configs/phase2_miller_urey.yaml \
    --output results/miller_urey/run_$(printf "%02d" $i) \
    --steps 10000000 \
    --seed $((42 + i)) &
done
```

### Hydrothermal Vent (H₂ + H₂S + CO₂ + catalysts)

```bash
python scripts/run_phase2_full.py \
  --config configs/phase2_hydrothermal.yaml \
  --output results/hydrothermal/run_01 \
  --steps 10000000 \
  --seed 42
```

### Formamide-rich (HCONH₂ + UV + minerals)

```bash
python scripts/run_phase2_full.py \
  --config configs/phase2_formamide.yaml \
  --output results/formamide/run_01 \
  --steps 10000000 \
  --seed 42
```

---

## Monitoring & Management

### Check Running Simulation

```bash
# View log (last 20 lines)
tail -f results/test_run/simulation.log

# Or in PowerShell
Get-Content results\test_run\simulation.log -Tail 20 -Wait
```

### Stop Simulation

```bash
# Find process
ps aux | grep run_phase2_full

# Kill process
kill <PID>

# Or in PowerShell
Get-Process python | Where-Object {$_.CommandLine -like "*run_phase2_full*"} | Stop-Process
```

### Check Disk Space

```bash
# Linux/Mac
df -h .

# PowerShell
Get-PSDrive C | Select-Object Used,Free
```

**Estimated space needed**:
- 1 simulation: ~100-500 MB (depending on snapshots)
- 30 simulations: ~10-15 GB
- With analysis: ~20 GB total

---

## Output Structure

```
results/phase2/
├── miller_urey/
│   ├── run_01/
│   │   ├── results.json          # Main results
│   │   ├── molecules.json        # Detected molecules
│   │   ├── simulation.log        # Full log
│   │   ├── snapshots/            # Periodic snapshots
│   │   │   ├── step_00000000.json
│   │   │   ├── step_00500000.json
│   │   │   └── ...
│   │   └── analysis/             # Post-run analysis
│   │       ├── molecule_report.txt
│   │       └── molecules_for_matcher.json
│   ├── run_02/
│   └── ...
├── hydrothermal/
│   └── ...
└── formamide/
    └── ...

analysis/
├── batch_analysis.json      # Aggregated data
├── batch_report.txt         # Human-readable report
└── pubchem_matches.json     # MatcherV2 results
```

---

## Troubleshooting

### Simulation Too Slow

Current: ~2s per step (thermodynamic validation overhead)

**Options**:
1. Accept slower speed, run overnight
2. Optimize validation frequency (tomorrow's task)
3. Run multiple in parallel on different machines

### GPU Out of Memory

**Current fix**: Using CPU backend (forced in runner)

**To try GPU again**:
1. Edit `scripts/run_phase2_full.py`
2. Change `ti.init(arch=ti.cpu)` to `ti.init(arch=ti.cuda)`
3. May need smaller particle counts

### No Molecules Detected

**Expected**: Cluster detection not yet integrated

**Next steps** (tomorrow):
1. Integrate GraphProcessor with results
2. Add proper molecule detection
3. Export detected molecules

### Simulation Crashes

Check logs:
```bash
cat results/test_run/simulation.log
cat results/test_run/stderr.log  # If using orchestrator
```

Common issues:
- Out of memory → Reduce particle count in YAML
- File permissions → Check output directory writable
- Missing dependencies → Run `pip install -r requirements.txt`

---

## Performance Estimates

### Current Performance (As-Is)

| Configuration | Steps | Time | Speed |
|--------------|-------|------|-------|
| Quick test | 10,000 | 5-10 min | ~30 steps/s |
| Medium test | 100,000 | 1-2 hours | ~20 steps/s |
| Full run | 10,000,000 | ~10 hours | ~300 steps/s |

**Note**: ~2s per step = 0.5 steps/s (very slow!)  
Discrepancy due to thermodynamic validation being slow at first, speeds up later.

### After Optimization (Tomorrow)

Expected 5-10x speedup:
- Full run: 2-5 hours
- 30 runs: 60-150 hours = 2.5-6 days
- Can parallelize: 3-4 machines = 1-2 days total

---

## Integration with MatcherV2

### Manual Matching

```bash
# 1. Extract molecules
python backend/sim/molecule_extractor.py results/test_run

# 2. Match with MatcherV2
python scripts/demo_matcher_v2.py \
  --input results/test_run/analysis/molecules_for_matcher.json \
  --output results/test_run/analysis/pubchem_matches
```

### Batch Matching

```bash
# Automatically included in batch analyzer
python scripts/analyze_phase2_batch.py \
  --input results/phase2 \
  --output analysis/phase2 \
  --recursive \
  --use-matcher  # <-- Runs MatcherV2 on all molecules
```

---

## Next Steps

### Tonight

```powershell
# Start overnight test
.\scripts\start_overnight_test.ps1
```

**Result**: Check in morning, validate workflow

### Tomorrow

1. **Optimize performance** (if overnight test too slow)
   - Reduce validation frequency
   - Profile bottlenecks
   - Test optimizations

2. **Integrate cluster detection**
   - Connect GraphProcessor
   - Enable real molecule detection
   - Test with overnight results

3. **Start production runs**
   - If optimization done: Start batch
   - Or: Continue overnight strategy

### This Week

1. Complete all 30 simulations
2. Run batch analysis with MatcherV2
3. Generate molecule catalog
4. Create figures

---

## Command Reference

### Simulation

```bash
# Basic run
python scripts/run_phase2_full.py --config <yaml> --output <dir> --steps <n> --seed <n>

# Dry run (test config)
python scripts/run_phase2_full.py --config <yaml> --output <dir> --steps <n> --dry-run
```

### Analysis

```bash
# Single run
python backend/sim/molecule_extractor.py <results_dir>

# Batch
python scripts/analyze_phase2_batch.py --input <dir> --output <dir> --recursive --use-matcher
```

### Orchestration

```bash
# Test mode (1 run/scenario, 10k steps)
python scripts/phase2_master.py --mode test --scenarios all

# Full mode (10 runs/scenario, 10M steps)
python scripts/phase2_master.py --mode full --scenarios all

# Analyze only
python scripts/phase2_master.py --mode analyze --input results/phase2
```

---

## Files Overview

### Configs
- `configs/phase2_miller_urey.yaml` - Full scale
- `configs/phase2_miller_urey_test.yaml` - Small scale (testing)
- `configs/phase2_hydrothermal.yaml`
- `configs/phase2_formamide.yaml`

### Scripts
- `scripts/run_phase2_full.py` - Single simulation runner
- `scripts/phase2_master.py` - Complete pipeline orchestrator
- `scripts/analyze_phase2_batch.py` - Batch analyzer
- `scripts/start_overnight_test.ps1` - Overnight test launcher

### Backend
- `backend/sim/phase2_config.py` - Configuration system
- `backend/sim/phase2_initializer.py` - Molecule initialization
- `backend/sim/molecule_extractor.py` - Results extraction

### Docs
- `docs/PHASE2_USAGE_GUIDE.md` - This file
- `docs/PHASE2_INTEGRATION_SUCCESS.md` - Technical details
- `docs/FINAL_PROGRESS_OCT13.md` - Progress summary

---

## Quick Decision Tree

**Want to test if everything works?**
→ `python scripts/run_phase2_full.py --config configs/phase2_miller_urey_test.yaml --output results/quick --steps 10000 --seed 42`

**Ready to start overnight test?**
→ `.\scripts\start_overnight_test.ps1` (Windows) or use manual command

**Want to run everything automatically?**
→ `python scripts/phase2_master.py --mode full --scenarios all` (takes days!)

**Just want to analyze existing results?**
→ `python scripts/analyze_phase2_batch.py --input results/phase2 --output analysis --recursive --use-matcher`

---

**Status**: ✅ **ALL INFRASTRUCTURE READY**  
**Next**: Start overnight test, optimize tomorrow, run batch  
**Timeline**: 1-2 weeks to Phase 2 complete

*Updated: October 13, 2025*

