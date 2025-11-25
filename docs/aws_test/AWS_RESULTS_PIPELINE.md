# AWS Results Pipeline Guide
**Complete workflow for downloading and analyzing AWS simulation results**

*Last updated: October 16, 2025*

---

## 📋 Overview

When running large-scale simulations on AWS, you need an automated pipeline to:
1. **Download** results from AWS instance
2. **Analyze** all simulations in batch
3. **Generate** publication-ready figures and reports

This guide covers the complete workflow.

---

## 🚀 Quick Start (One Command)

```bash
# Complete pipeline: Download → Analyze → Report
bash scripts/aws_pipeline.sh <aws-host> <ssh-key>
```

**Example**:
```bash
bash scripts/aws_pipeline.sh 54.123.45.67 ~/.ssh/aws_key.pem
```

This will:
- ✅ Download all completed simulations
- ✅ Extract molecules from each run
- ✅ Build reaction networks per scenario
- ✅ Detect autocatalytic cycles
- ✅ Generate comparison statistics
- ✅ Create publication figures

---

## 📥 Step 1: Download Results

### Using the downloader script:

```bash
python scripts/aws_results_downloader.py \
    --host 54.123.45.67 \
    --key ~/.ssh/aws_key.pem \
    --local-base ./results/aws_batch
```

### Options:

- `--host`: AWS instance IP or hostname (required)
- `--key`: Path to SSH private key (required)
- `--remote-base`: Remote results directory (default: `~/live2.0/results`)
- `--local-base`: Local download directory (default: `./results/aws_batch`)
- `--force`: Re-download existing files
- `--status-only`: Check status without downloading

### Check download status:

```bash
python scripts/aws_results_downloader.py \
    --host 54.123.45.67 \
    --key ~/.ssh/aws_key.pem \
    --status-only
```

**Output**:
```
📊 DOWNLOAD STATUS
====================================
Last update: 2025-10-16T10:30:00
Total downloaded: 24
Valid simulations: 24
Invalid simulations: 0

By scenario:
  miller_urey: 8
  hydrothermal: 8
  formamide: 8
```

---

## 📊 Step 2: Analyze Results

### Run batch analysis:

```bash
python scripts/aws_results_analyzer.py \
    --input ./results/aws_batch
```

### What it does:

1. **Discovers** all completed simulations
2. **Extracts** molecules from each run using `quick_analyze.py`
3. **Builds** reaction networks per scenario using `reaction_network_analyzer.py`
4. **Detects** autocatalytic cycles using `autocatalytic_detector.py`
5. **Compares** scenarios using `compare_scenarios.py`
6. **Generates** figures using plot scripts
7. **Creates** comprehensive JSON report

### Output structure:

```
results/aws_batch/
├── miller_urey/
│   ├── run_1/
│   ├── run_2/
│   └── ...
├── hydrothermal/
│   ├── run_1/
│   └── ...
├── formamide/
│   └── ...
└── analysis/
    ├── batch_analysis_report.json  ← Main report
    ├── miller_urey/
    │   ├── run_1_molecules.json
    │   ├── network/
    │   │   └── reaction_network.json
    │   └── cycles/
    │       └── cycles_report.json
    ├── hydrothermal/
    ├── formamide/
    ├── comparison/
    │   └── scenario_comparison.json
    └── figures/
        ├── molecular_diversity.png
        ├── reaction_networks.png
        ├── autocatalytic_cycles.png
        └── emergence_timeline.png
```

---

## 📈 Step 3: Review Results

### Main analysis report:

```bash
cat results/aws_batch/analysis/batch_analysis_report.json | python -m json.tool
```

**Example output**:
```json
{
  "timestamp": "2025-10-16T10:45:00",
  "scenarios": {
    "miller_urey": {
      "network_file": "results/aws_batch/analysis/miller_urey/network/reaction_network.json",
      "cycles": 5
    },
    "hydrothermal": {
      "network_file": "results/aws_batch/analysis/hydrothermal/network/reaction_network.json",
      "cycles": 8
    },
    "formamide": {
      "network_file": "results/aws_batch/analysis/formamide/network/reaction_network.json",
      "cycles": 12
    }
  },
  "summary": {
    "total_scenarios": 3,
    "total_cycles": 25,
    "figures_generated": 4,
    "errors_count": 0
  }
}
```

### View generated figures:

```bash
# Open figures directory
cd results/aws_batch/analysis/figures
ls -lh
```

### Scenario comparison:

```bash
cat results/aws_batch/analysis/comparison/scenario_comparison.json | python -m json.tool
```

---

## 🔍 Monitoring AWS Progress (Real-time)

While simulations are running on AWS:

```bash
# Check how many are completed
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67 \
    'cd live2.0/results && find . -name "summary.txt" | wc -l'

# Monitor system resources
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67 \
    'free -h && df -h && ps aux | grep python | grep -v grep | wc -l'

# Watch logs
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67 \
    'tail -f live2.0/logs/*.log'
```

Or use the monitoring script:

```bash
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67 \
    'bash live2.0/monitor_aws_runs.sh'
```

---

## 🛠️ Troubleshooting

### Problem: Connection refused

**Solution**: Check security group allows SSH (port 22) from your IP:
```bash
# Test connection
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67 echo "OK"
```

### Problem: Download timeout

**Solution**: Large files may timeout. Increase timeout or use `--force` to resume:
```bash
python scripts/aws_results_downloader.py \
    --host 54.123.45.67 \
    --key ~/.ssh/aws_key.pem \
    --force
```

### Problem: Invalid simulations

**Solution**: Re-run incomplete simulations:
```bash
# Check which are invalid
python scripts/aws_results_downloader.py \
    --host 54.123.45.67 \
    --key ~/.ssh/aws_key.pem \
    --status-only

# On AWS, re-run failed ones
ssh -i ~/.ssh/aws_key.pem ubuntu@54.123.45.67
cd live2.0
# Re-run specific simulation...
```

### Problem: Analysis fails

**Solution**: Check error log in `batch_analysis_report.json`:
```bash
cat results/aws_batch/analysis/batch_analysis_report.json | \
    python -m json.tool | grep -A5 "errors"
```

---

## 📊 Integration with Paper

After pipeline completes, use results for paper:

### 1. Generate all publication figures:

```bash
python scripts/generate_all_figures.py \
    --input results/aws_batch/analysis \
    --output paper/figures
```

### 2. Extract statistics for text:

```python
import json

# Load report
with open('results/aws_batch/analysis/batch_analysis_report.json') as f:
    report = json.load(f)

# For Results section:
print(f"Total cycles detected: {report['summary']['total_cycles']}")
print(f"Miller-Urey cycles: {report['scenarios']['miller_urey']['cycles']}")
# etc.
```

### 3. Copy figures to paper directory:

```bash
cp results/aws_batch/analysis/figures/*.png paper/figures/
```

---

## 🔄 Incremental Updates

If AWS is still running and you want to analyze partial results:

```bash
# Download only new results
python scripts/aws_results_downloader.py \
    --host 54.123.45.67 \
    --key ~/.ssh/aws_key.pem

# Analyze what's available
python scripts/aws_results_analyzer.py \
    --input ./results/aws_batch
```

The pipeline is **incremental** - it won't re-process already analyzed runs.

---

## ⏱️ Expected Timings

| Step | Duration | Notes |
|------|----------|-------|
| Download 24 runs | 5-15 min | Depends on file size and bandwidth |
| Extract molecules | 10-20 min | ~30 sec per run |
| Build networks | 2-5 min | Per scenario |
| Detect cycles | 1-3 min | Per scenario |
| Generate figures | 2-5 min | All figures |
| **Total** | **~20-45 min** | For 24 simulations |

For 72 simulations (3 rounds): ~60-90 minutes total.

---

## 📋 Checklist: Ready for Paper?

- [ ] All simulations downloaded (check count)
- [ ] No invalid simulations
- [ ] All molecules extracted
- [ ] Reaction networks built (3 scenarios)
- [ ] Autocatalytic cycles detected (>10 total)
- [ ] Scenario comparison complete
- [ ] All 7 figures generated
- [ ] Analysis report reviewed
- [ ] Statistics extracted for text
- [ ] Figures copied to paper directory

---

## 🎯 Next Steps

Once pipeline is complete:

1. **Review figures** - Check quality and clarity
2. **Extract statistics** - For Results section
3. **Select top molecules** - For detailed analysis
4. **Write Results section** - Using generated data
5. **DFT validation** - For top 5 molecules (optional)

---

## 📞 Quick Reference Commands

```bash
# Complete pipeline
bash scripts/aws_pipeline.sh 54.123.45.67 ~/.ssh/aws_key.pem

# Just download
python scripts/aws_results_downloader.py --host <ip> --key <key>

# Just analyze
python scripts/aws_results_analyzer.py --input results/aws_batch

# Check status
python scripts/aws_results_downloader.py --host <ip> --key <key> --status-only

# Monitor AWS
ssh -i <key> ubuntu@<ip> 'bash live2.0/monitor_aws_runs.sh'
```

---

*For questions or issues, refer to the main documentation or open an issue on GitHub.*

