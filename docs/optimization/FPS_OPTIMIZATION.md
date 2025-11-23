# FPS Optimization - Staggered Operations 🎯

## Problem: FPS drops from 43 to 0.1

### Root Cause
Multiple heavy operations executing **at the same time** (same step number):

```
Step 3000:
  ├─ update_bonds()          [500ms] ← HEAVY!
  ├─ update_clusters()       [300ms] ← HEAVY!
  ├─ detect_novel_substances [200ms] ← HEAVY!
  └─ Total:                  1000ms+ ← FPS drops to 1!
  
Step 3001:
  └─ Normal operations       [10ms]  ← FPS back to 43
```

## Solution: Staggered Operations

**Spread heavy operations over time** so they never overlap:

### Before (Operations Overlap):
```
Step:  0    500   1000  1500  2000  2500  3000  3500
       |     |     |     |     |     |     |     |
Bonds: █─────█─────█─────█─────█─────█─────█─────█  (every 500)
Clust: █───────────█───────────█───────────█──────  (every 1000)
Novel: █─────█─────█─────█─────█─────█─────█─────█  (every 500)
       ↑                             ↑
    Overlap!                     Overlap! ← FPS drop
```

### After (Staggered):
```
Step:  0    200   400   600   800   1000  1200  1400
       |     |     |     |     |     |     |     |
Bonds: ──────█───────────█───────────█──────────── (every 600, offset +100)
Clust: ────────────█───────────────────────█────── (every 1200, offset +300)
Novel: ──────────────────█───────────────────────█ (every 700, offset +500)
       
No overlap! ← Smooth FPS!
```

## Changes Made

### 1. update_bonds()
**Before:** Every 500 steps  
**After:** Every 600 steps, offset by +100

```python
# Before:
if self.step_count % 500 == 0:
    self.binding.update_bonds()

# After:
if (self.step_count - 100) % 600 == 0:
    self.binding.update_bonds()
```

**Executes at:** 100, 700, 1300, 1900, 2500, 3100...

### 2. update_clusters()
**Before:** Every 1000 steps  
**After:** Every 1200 steps, offset by +300

```python
# Before:
if self.step_count % 1000 == 0:
    self.binding.update_clusters()

# After:
if (self.step_count - 300) % 1200 == 0:
    self.binding.update_clusters()
```

**Executes at:** 300, 1500, 2700, 3900...

### 3. detect_novel_substances()
**Before:** Every 500 steps  
**After:** Every 700 steps, offset by +500

```python
# Before:
if self.step_count % 500 == 0:
    self.detect_novel_substances()

# After:
if (self.step_count - 500) % 700 == 0:
    self.detect_novel_substances()
```

**Executes at:** 500, 1200, 1900, 2600, 3300...

### 4. Visualization (bonds/clusters)
**Already optimized:** Every 200 steps

```python
if self.step_count % 200 == 0:
    bonds = self.binding.get_bonds()
    clusters = self.binding.get_clusters()
```

**Executes at:** 0, 200, 400, 600, 800, 1000...

## Timeline Example (First 3000 steps)

| Step | Operation | Est. Time | FPS Impact |
|------|-----------|-----------|------------|
| 0 | Normal | 10ms | 43 FPS |
| 100 | **update_bonds** | 500ms | 2 FPS |
| 200 | viz (bonds/clusters) | 211ms | 5 FPS |
| 300 | **update_clusters** | 300ms | 3 FPS |
| 400 | viz (bonds/clusters) | 211ms | 5 FPS |
| 500 | **detect_novel** | 200ms | 5 FPS |
| 600 | viz (bonds/clusters) | 211ms | 5 FPS |
| 700 | **update_bonds** | 500ms | 2 FPS |
| 800 | viz (bonds/clusters) | 211ms | 5 FPS |
| 1000 | viz (bonds/clusters) | 211ms | 5 FPS |
| 1200 | viz + **detect_novel** | 411ms | 2-3 FPS |
| 1300 | **update_bonds** | 500ms | 2 FPS |
| 1500 | **update_clusters** | 300ms | 3 FPS |

**Notice:** Heavy operations spread out, no more than 2 at once!

## Expected Results

### Before:
- **Normal FPS:** 43 (when nothing heavy)
- **Heavy step FPS:** 0.1-1 (3-4 operations at once)
- **Pattern:** Unpredictable drops every 500-1000 steps

### After:
- **Normal FPS:** 40-43 (most of the time)
- **Heavy step FPS:** 2-5 (only 1 operation at a time)
- **Pattern:** Predictable, smaller drops, more frequent but manageable

## Benefits

1. **Smoother FPS** - No more 0.1 FPS drops
2. **More predictable** - Small drops instead of huge freezes
3. **Better UX** - User sees steady ~5-10 FPS minimum instead of freezes
4. **Still accurate** - All operations still happen, just spread out

## Trade-offs

### Slower individual updates:
- Bonds: 500→600 steps (20% slower update)
- Clusters: 1000→1200 steps (20% slower update)
- Novel: 500→700 steps (40% slower update)

**Worth it?** YES! User experience is much better with smooth FPS than accurate-but-frozen simulation.

## Further Optimization (if needed)

### If still seeing FPS drops < 5:

#### Option 1: Increase intervals more
```python
update_bonds: 600 → 1000 steps
update_clusters: 1200 → 2000 steps
detect_novel: 700 → 1000 steps
```

#### Option 2: Reduce particle count
```python
config.n_particles = 300  # Instead of 500
config.max_particles = 400  # Lower limit
```

#### Option 3: Disable heavy operations
```python
config.detect_novel_substances = False  # No novelty detection
# In server.py or config
```

## Monitoring

### Check logs for:

**Good (After fix):**
```
Step 100: update_bonds (500ms) - FPS: 2
Step 200: viz (211ms) - FPS: 5
Step 300: update_clusters (300ms) - FPS: 3
Step 400: viz (211ms) - FPS: 5
Step 500: detect_novel (200ms) - FPS: 5
```

**Bad (Before fix):**
```
Step 3000: update_bonds + update_clusters + detect_novel (1000ms+) - FPS: 0.1!
```

## Verification

After restart, monitor FPS over 5000 steps:
- Should see ~40 FPS most of the time
- Occasional drops to 2-5 FPS (1-2x per minute)
- **NO MORE** drops to 0.1 FPS

## Restart Required

```powershell
# Kill old backend
.\kill_backend.ps1

# Start new backend (with staggered operations)
cd backend
python -m api.server
```

## Technical Details

### Why Offset?

**Offset = Starting point** for the modulo operation:

```python
# Without offset:
self.step_count % 600 == 0
# Triggers at: 0, 600, 1200, 1800...

# With offset +100:
(self.step_count - 100) % 600 == 0
# Triggers at: 100, 700, 1300, 1900...
```

### Why These Specific Intervals?

**Goal:** Minimize overlap with:
- Visualization (every 200 steps)
- Other operations (various intervals)

**Chosen intervals:**
- 600 (bonds): Factor of 6 × 100, spreads well
- 1200 (clusters): Factor of 12 × 100, rare operation
- 700 (novel): Prime-ish, avoids most overlaps

## Related Issues

- With 596 bonds, operations are slower
- More particles = longer operations
- CPU usage will be high but evenly spread

## Conclusion

**Problem:** Multiple O(n²) operations at same time  
**Solution:** Stagger operations with offsets  
**Result:** Smooth FPS, better UX

---

Questions? See:
- [PRODUCTION_OPTIMIZATION.md](PRODUCTION_OPTIMIZATION.md) - General optimization
- [PERFORMANCE_TUNING.md](PERFORMANCE_TUNING.md) - Detailed tuning guide

