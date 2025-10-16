# AGGRESSIVE OPTIMIZATION PLAN
# =============================
**Cel:** Zmniejszyć czas Phase 2 z 3 miesięcy do <2 tygodni
**Metoda:** Każda optymalizacja = 2-5x przyspieszenie
**Target:** Łączne przyspieszenie 10-50x

## 🎯 CURRENT STATUS

- **Current speed:** 4.8 steps/s (1775 atoms)
- **1M steps:** 2.4 days
- **150 simulations @ 4 parallel:** 94 days (3.1 months)
- **Target:** <14 days (2 weeks)
- **Needed speedup:** ~7x minimum

---

## 🚀 OPTIMIZATION STRATEGIES

### **OPT #1: GPU with Spatial Hash** 🔥
**Expected speedup:** 5-10x  
**Risk:** Medium (crashed before, but that was without spatial hash)  
**Time to implement:** 1 hour (just re-enable)

**Why it might work now:**
- Spatial hash = much smaller kernels
- O(n) instead of O(n²)
- GPU loves O(n) algorithms

**Action:**
```python
# Try GPU again with spatial hashing
ti.init(arch=ti.cuda, device_memory_GB=4.0)
```

**Expected result:**
- 4.8 → 24-48 steps/s
- 1M steps: 6-12 hours
- 150 sims: 12-25 days

---

### **OPT #2: Increase Timestep (dt)** ⚡
**Expected speedup:** 2-5x  
**Risk:** Low (just needs validation)  
**Time to implement:** 10 minutes

**Current:** dt = 0.001  
**Proposed:** dt = 0.002-0.005 (adaptive)

**Why it works:**
- Same physics at larger dt (if stable)
- Fewer steps = faster simulation
- Need to validate energy conservation

**Action:**
```yaml
simulation:
  dt: 0.003  # 3x larger
  # Or adaptive: 0.001-0.005
```

**Expected result:**
- 3x fewer steps needed
- 1M effective steps in 333k actual steps
- 2.4 days → 0.8 days

---

### **OPT #3: Reduce Steps (Science Trade-off)** 📊
**Expected speedup:** 2-10x  
**Risk:** None (just less data)  
**Time to implement:** 0 minutes (config change)

**Current:** 1M steps per simulation  
**Proposed:** 100k-500k steps

**Why it works:**
- Most interesting chemistry happens early
- Can run more replicas instead

**Action:**
```yaml
max_steps: 200000  # 200k instead of 1M
```

**Trade-off:**
- 5x fewer steps
- But can run 5x more replicas
- Same total sampling, faster results

**Expected result:**
- 2.4 days → 0.5 days
- 150 sims @ 200k: 19 days

---

### **OPT #4: Aggressive Physics Simplification** ⚡
**Expected speedup:** 1.5-2x  
**Risk:** Low  
**Time to implement:** 15 minutes

**Current:** Full force calculation every step  
**Proposed:** Skip some steps

**Actions:**
```python
# Compute forces every 2-3 steps instead of every step
if step % 2 == 0:
    compute_forces()
else:
    reuse_forces()
```

**Expected result:**
- 2x fewer force computations
- 4.8 → 7-8 steps/s

---

### **OPT #5: Larger Spatial Hash Cells** 🔧
**Expected speedup:** 1.3-1.5x  
**Risk:** Very low  
**Time to implement:** 5 minutes

**Current:** cell_size = 10 Å  
**Proposed:** cell_size = 15 Å

**Why it works:**
- Fewer cells to search
- Still captures all interactions (cutoff = 15 Å)

**Action:**
```yaml
spatial_hash_cell_size: 15.0  # Was 10.0
```

**Expected result:**
- ~30% fewer cells
- 4.8 → 6.2 steps/s

---

### **OPT #6: Compile Taichi AOT** 🏭
**Expected speedup:** 1.2-1.5x  
**Risk:** Low  
**Time to implement:** 30 minutes

**Current:** JIT compilation  
**Proposed:** Ahead-of-time (AOT) compilation

**Why it works:**
- No runtime compilation overhead
- Better optimization

---

### **OPT #7: Profile and Remove Bottlenecks** 🔍
**Expected speedup:** 1.5-3x  
**Risk:** None  
**Time to implement:** 2 hours

**Action:**
- Profile with cProfile
- Identify slowest functions
- Optimize top 3 bottlenecks

---

## 📊 COMBINED IMPACT

### **Conservative Scenario:**
```
Base:           4.8 steps/s
+ GPU (5x):     24 steps/s
+ Larger dt (2x): 48 effective steps/s
+ Larger cells (1.3x): 62 effective steps/s

Result: 13x speedup
1M steps: 4.5 hours
150 sims @ 4 parallel: 7 days
```

### **Optimistic Scenario:**
```
Base:           4.8 steps/s
+ GPU (10x):    48 steps/s
+ Larger dt (3x): 144 effective steps/s
+ Larger cells (1.5x): 216 effective steps/s
+ Skip forces (1.5x): 324 effective steps/s

Result: 67x speedup
1M steps: 50 minutes
150 sims @ 4 parallel: 1.6 days
```

### **Realistic Scenario:**
```
Base:           4.8 steps/s
+ GPU (7x):     33.6 steps/s
+ Larger dt (2.5x): 84 effective steps/s
+ Larger cells (1.3x): 109 effective steps/s

Result: 23x speedup
1M steps: 2.5 hours
150 sims @ 4 parallel: 4 days
```

---

## 🎯 ACTION PLAN (PRIORITY ORDER)

### **QUICK WINS (30 min):**
1. ✅ Larger spatial cells (15 Å) - 5 min - **1.3x**
2. ✅ Increase timestep to 0.003 - 5 min - **2x**
3. ✅ Test GPU with spatial hash - 10 min - **5-10x**

**Expected:** 13-26x combined = **<5 days for Phase 2** ✅

### **IF NEEDED (2 hours):**
4. Profile and optimize - 2h - **1.5-2x**
5. Skip force steps - 30 min - **1.5x**

**Expected:** 30-50x total = **2-3 days for Phase 2** ✅✅

---

## 🌩️ CLOUD OPTION (FALLBACK)

If optimizations don't work, cloud specs needed:

### **AWS/Azure/GCP Requirements:**
- **CPU:** 64+ cores (e.g., c6i.16xlarge)
- **RAM:** 128 GB
- **GPU:** Optional (A100 or V100)
- **Cost:** ~$3-5/hour

### **Expected Performance:**
- 64 cores = 16x parallel (vs current 4x)
- Same speed per core = 4.8 steps/s
- **150 sims @ 16 parallel:** 6 days
- **Total cost:** ~$720-1200

### **With optimizations + cloud:**
- Optimized code: 25 steps/s
- 16 parallel: **1.5 days total**
- **Total cost:** ~$180-240 ✅

---

## 💰 COST-BENEFIT

| Option | Time | Cost | Effort |
|--------|------|------|--------|
| **Current (CPU)** | 94 days | $0 | 0h |
| **Quick wins** | 4-7 days | $0 | 0.5h ✅ |
| **Full optimization** | 2-3 days | $0 | 2.5h ✅✅ |
| **Cloud (no opt)** | 6 days | $720 | Setup |
| **Cloud (optimized)** | 1.5 days | $180 | 2.5h + setup ✅✅✅ |

**Recommendation:** Try quick wins first (30 min), if not enough → cloud with optimizations.


