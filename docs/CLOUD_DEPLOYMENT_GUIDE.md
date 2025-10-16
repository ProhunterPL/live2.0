# CLOUD DEPLOYMENT GUIDE
# =======================
**Phase 2 Prebiotic Chemistry Simulation**  
**Optimized for AWS/Azure/GCP**

## 📊 PERFORMANCE SUMMARY

### **Local Performance (Baseline):**
- **CPU:** Intel/AMD 28 threads
- **Speed:** 4.8 steps/s (baseline), 11.4 steps/s (optimized)
- **1M steps:** 1.0 day
- **150 simulations @ 4 parallel:** 38 days

### **Target Cloud Performance:**
- **CPU:** 64+ cores (16-32 parallel)
- **Expected:** 2-8 days for full Phase 2 (150 sims)
- **Cost:** $130-500 depending on configuration

---

## 🌩️ RECOMMENDED CLOUD INSTANCES

### **AWS Options:**

#### **1. c6i.16xlarge (RECOMMENDED)**
- **vCPUs:** 64 (32 physical cores)
- **RAM:** 128 GB
- **Cost:** ~$2.72/hour (~$65/day)
- **Performance:** 16 parallel jobs
- **Timeline:** 2.4 days for 150 sims
- **Total cost:** ~$156

#### **2. c6i.32xlarge (FASTEST)**
- **vCPUs:** 128 (64 physical cores)
- **RAM:** 256 GB
- **Cost:** ~$5.44/hour (~$131/day)
- **Performance:** 32 parallel jobs
- **Timeline:** 1.2 days for 150 sims
- **Total cost:** ~$157

#### **3. c6i.8xlarge (BUDGET)**
- **vCPUs:** 32 (16 physical cores)
- **RAM:** 64 GB
- **Cost:** ~$1.36/hour (~$33/day)
- **Performance:** 8 parallel jobs
- **Timeline:** 4.8 days for 150 sims
- **Total cost:** ~$158

### **Azure Options:**

#### **1. F64s v2 (RECOMMENDED)**
- **vCPUs:** 64
- **RAM:** 128 GB
- **Cost:** ~$2.72/hour
- **Equivalent to AWS c6i.16xlarge**

### **GCP Options:**

#### **1. c2-standard-60 (RECOMMENDED)**
- **vCPUs:** 60
- **RAM:** 240 GB
- **Cost:** ~$2.52/hour
- **High single-thread performance**

---

## 🚀 QUICK START

### **1. Launch Instance:**

```bash
# AWS CLI example (c6i.16xlarge)
aws ec2 run-instances \
  --image-id ami-0c55b159cbfafe1f0 \
  --instance-type c6i.16xlarge \
  --key-name your-key \
  --security-group-ids sg-xxxxx \
  --subnet-id subnet-xxxxx \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":200}}]'
```

### **2. Setup Environment:**

```bash
# SSH into instance
ssh -i your-key.pem ubuntu@<instance-ip>

# Update system
sudo apt update && sudo apt upgrade -y

# Install dependencies
sudo apt install -y python3.11 python3-pip git

# Clone repository
git clone https://github.com/ProhunterPL/live2.0.git
cd live2.0

# Install Python packages
pip3 install -r requirements.txt
```

### **3. Verify Setup:**

```bash
# Test run (1000 steps)
python3 scripts/run_phase2_full.py \
  --config configs/phase2_final_test.yaml \
  --output results/test \
  --steps 1000 \
  --seed 42

# Check performance
# Expected: 4-6 steps/s (depends on instance)
```

### **4. Launch Production Runs:**

```bash
# Option A: Manual parallel (4 jobs)
for i in {1..4}; do
  python3 scripts/run_phase2_full.py \
    --config configs/phase2_miller_urey_1M.yaml \
    --output results/run_$i \
    --steps 200000 \
    --seed $((42 + i)) &
done

# Option B: Use orchestrator
python3 scripts/phase2_master_1M.py \
  --mode full \
  --scenarios all \
  --max-parallel 16
```

---

## 📋 REQUIREMENTS

### **Python Packages:**
```
taichi>=1.7.0
numpy>=1.24.0
pyyaml>=6.0
scipy>=1.10.0
networkx>=3.0
pydantic>=2.0
```

### **System:**
- **Python:** 3.11+
- **RAM:** 4 GB per parallel job (64 GB for 16 parallel)
- **Disk:** 100 GB (for results)
- **OS:** Ubuntu 22.04 LTS (recommended)

---

## 📊 COST ESTIMATES

### **For 150 Simulations (200k steps each):**

| Instance | Parallel | Days | $/hour | Total Cost |
|----------|----------|------|--------|------------|
| c6i.8xlarge | 8 | 4.8 | $1.36 | **$157** |
| c6i.16xlarge | 16 | 2.4 | $2.72 | **$157** |
| c6i.32xlarge | 32 | 1.2 | $5.44 | **$157** |

**Note:** Cost is similar across instance types! Choose based on timeline.

### **For 300 Simulations (better statistics):**

| Instance | Days | Total Cost |
|----------|------|------------|
| c6i.16xlarge | 4.8 | **$314** |
| c6i.32xlarge | 2.4 | **$314** |

---

## 💰 COST OPTIMIZATION

### **1. Use Spot Instances (50-90% savings):**
```bash
# AWS Spot request
aws ec2 request-spot-instances \
  --instance-count 1 \
  --type "one-time" \
  --launch-specification file://spot-config.json

# Savings: $157 -> $15-80 (depends on availability)
```

### **2. Use Reserved Instances (30-70% savings):**
- Good if running multiple experiments
- 1-year commitment

### **3. Use Preemptible VMs (GCP):**
- 60-80% cheaper
- May be interrupted (need checkpointing)

---

## 🔧 PERFORMANCE TUNING

### **1. Maximize CPU Usage:**
```yaml
# In config file
simulation:
  dt: 0.003  # Larger timestep
  spatial_hash_cell_size: 15.0  # Larger cells
```

### **2. Parallel Execution:**
```bash
# Match parallel jobs to CPU count
# c6i.16xlarge (64 cores) -> 16-32 parallel jobs
PARALLEL_JOBS=$(nproc)
echo "Running $PARALLEL_JOBS jobs"
```

### **3. Monitor Performance:**
```bash
# Install htop
sudo apt install htop

# Monitor CPU usage
htop

# Should see 90-100% usage across all cores
```

---

## 📈 EXPECTED TIMELINES

### **200k Steps per Simulation:**

| Scenario | Sims | Parallel | Days | Cost (c6i.16xlarge) |
|----------|------|----------|------|---------------------|
| Quick Test | 30 | 16 | 0.5 | **$33** |
| Standard | 150 | 16 | 2.4 | **$157** |
| Extended | 300 | 16 | 4.8 | **$314** |

### **1M Steps per Simulation:**

| Scenario | Sims | Parallel | Days | Cost |
|----------|------|----------|------|------|
| Standard | 150 | 16 | 9.5 | **$622** |
| Extended | 300 | 16 | 19 | **$1,244** |

**Recommendation:** Use 200k steps (still scientifically valid, 5x cheaper!)

---

## ⚠️ IMPORTANT NOTES

### **1. Data Transfer:**
- Results will be ~10-50 GB
- Download to local: `scp -r ubuntu@<ip>:~/live2.0/results ./`
- Or upload to S3: `aws s3 sync results/ s3://your-bucket/`

### **2. Monitoring:**
```bash
# Check progress
tail -f results/*/simulation.log

# Count completed runs
ls results/*/results.json | wc -l
```

### **3. Checkpointing:**
- Current code doesn't support resume
- Use shorter runs (200k instead of 1M)
- Or implement checkpointing if needed

### **4. GPU:**
- **NOT RECOMMENDED** (caused system crash)
- Spatial hash works well on CPU
- GPU may work after further optimization

---

## 🎯 RECOMMENDED WORKFLOW

### **Day 1: Setup & Test**
1. Launch instance
2. Setup environment (30 min)
3. Run test (1000 steps, 5 min)
4. Validate performance (>4 steps/s)

### **Day 2-3: Production Runs**
5. Launch 16 parallel jobs
6. Monitor progress (htop, logs)
7. Let run overnight

### **Day 4: Download & Cleanup**
8. Download results to local/S3
9. Terminate instance
10. Total cost: ~$157

---

## 📞 SUPPORT

### **If Performance Lower Than Expected:**
- Check CPU usage (should be 90-100%)
- Check for thermal throttling: `sensors`
- Try larger dt (0.005 instead of 0.003)
- Use c-series instances (compute-optimized)

### **If Out of Memory:**
- Reduce parallel jobs
- Use larger instance (more RAM)
- Set `gc_interval` lower in config

### **If Crashes:**
- Check logs: `tail -100 results/*/simulation.log`
- Reduce parallel jobs
- GPU should be disabled (CPU only)

---

## ✅ CHECKLIST

Before production run:
- [ ] Instance launched and SSH working
- [ ] Repository cloned
- [ ] Dependencies installed
- [ ] Test run completed (1000 steps)
- [ ] Performance validated (>4 steps/s)
- [ ] Config files prepared
- [ ] Parallel jobs configured
- [ ] Monitoring setup (htop, logs)
- [ ] Data backup plan (S3/local)

---

## 🎉 SUCCESS METRICS

**You'll know it's working when:**
- ✅ CPU usage: 90-100% across all cores
- ✅ Speed: 4-6 steps/s per job
- ✅ Memory: <4 GB per job
- ✅ No errors in logs
- ✅ Results files being created

**Timeline check:**
- 1 simulation (200k steps): ~5-7 hours
- 16 parallel: 5-7 hours for 16 sims
- 150 total: ~2.4 days

**Good luck! 🚀**


