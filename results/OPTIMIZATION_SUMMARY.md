# Quick Optimization Summary
# ==========================

## Baseline (with spatial hash):
- dt = 0.001
- cell_size = 10.0
- Speed: 4.8 steps/s
- 1M steps: 2.4 days

## OPT #1: Larger Timestep (dt=0.003)
- Speed: 3.8 steps/s (raw)
- But: 3x more physics per step
- Effective: 11.4 steps/s equivalent
- **Speedup: 2.4x**
- 1M steps: 1.0 day ✅

## Target with all optimizations:
- dt = 0.003: 2.4x
- cell_size = 15.0: 1.3x (pending test)
- GPU: 5-10x (risky)

**Conservative estimate: 3-4x total = 0.6-0.8 days per simulation**
**150 sims @ 4 parallel = 23-30 days (1 month)** ✅

**With GPU: 15-30x total = 2-4 hours per simulation**
**150 sims @ 4 parallel = 3-6 days** ✅✅✅

## Recommendation:
1. Test with cell_size=15.0 (expect 1.3x more)
2. If total <1 month acceptable: USE THIS
3. If need faster: Try GPU carefully or use cloud


