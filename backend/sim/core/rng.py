"""
Random Number Generator for Live 2.0 simulation
Provides deterministic, seeded random number generation
"""

import taichi as ti
import numpy as np
from typing import Optional, Tuple

@ti.data_oriented
class RNG:
    """Deterministic random number generator using Taichi"""
    
    def __init__(self, seed: Optional[int] = None):
        if seed is None:
            seed = np.random.randint(0, 2**31)
        
        self.seed = seed
        self.state = ti.field(dtype=ti.u32, shape=())
        self.state[None] = seed
    
    @ti.kernel
    def reset(self, new_seed: ti.i32):
        """Reset RNG with new seed"""
        self.state[None] = new_seed
    
    @ti.func
    def next(self) -> ti.f32:
        """Generate next random float in [0, 1)"""
        # Linear congruential generator
        self.state[None] = (self.state[None] * 1664525 + 1013904223) % (2**32)
        return self.state[None] / (2**32)
    
    @ti.func
    def next_int(self, max_val: ti.i32) -> ti.i32:
        """Generate next random integer in [0, max_val)"""
        return int(self.next() * max_val)
    
    @ti.func
    def next_range(self, min_val: ti.f32, max_val: ti.f32) -> ti.f32:
        """Generate next random float in [min_val, max_val)"""
        return min_val + self.next() * (max_val - min_val)
    
    @ti.func
    def next_gaussian(self, mean: ti.f32 = 0.0, std: ti.f32 = 1.0) -> ti.f32:
        """Generate next random float from Gaussian distribution (Box-Muller)"""
        # Box-Muller transform
        u1 = self.next()
        u2 = self.next()
        
        if u1 == 0.0:
            u1 = 1e-10  # Avoid log(0)
        
        z0 = ti.sqrt(-2.0 * ti.log(u1)) * ti.cos(2.0 * 3.14159265359 * u2)
        return mean + std * z0
    
    @ti.func
    def next_vector2(self, min_val: ti.f32 = 0.0, max_val: ti.f32 = 1.0) -> ti.Vector([2], ti.f32):
        """Generate next random 2D vector"""
        return ti.Vector([self.next_range(min_val, max_val), 
                         self.next_range(min_val, max_val)])
    
    @ti.func
    def next_vector3(self, min_val: ti.f32 = 0.0, max_val: ti.f32 = 1.0) -> ti.Vector([3], ti.f32):
        """Generate next random 3D vector"""
        return ti.Vector([self.next_range(min_val, max_val),
                         self.next_range(min_val, max_val),
                         self.next_range(min_val, max_val)])
    
    @ti.func
    def next_vector4(self, min_val: ti.f32 = 0.0, max_val: ti.f32 = 1.0) -> ti.Vector([4], ti.f32):
        """Generate next random 4D vector"""
        return ti.Vector([self.next_range(min_val, max_val),
                         self.next_range(min_val, max_val),
                         self.next_range(min_val, max_val),
                         self.next_range(min_val, max_val)])
    
    @ti.func
    def bernoulli(self, p: ti.f32) -> ti.i32:
        """Generate Bernoulli random variable (0 or 1) with probability p"""
        return 1 if self.next() < p else 0
    
    @ti.func
    def choice(self, weights: ti.template()) -> ti.i32:
        """Choose index based on weights (roulette wheel selection)"""
        total_weight = 0.0
        for i in range(weights.shape[0]):
            total_weight += weights[i]
        
        if total_weight == 0.0:
            return 0
        
        r = self.next() * total_weight
        cumulative = 0.0
        
        for i in range(weights.shape[0]):
            cumulative += weights[i]
            if cumulative >= r:
                return i
        
        return weights.shape[0] - 1
    
    def get_state(self) -> int:
        """Get current RNG state"""
        return self.state[None]
    
    def set_state(self, state: int):
        """Set RNG state"""
        self.state[None] = state

    # -------- Python-scope helpers (safe to call from normal Python) --------
    def py_next(self) -> float:
        rnd = np.random.default_rng(self.seed)
        val = float(rnd.random())
        # advance seed deterministically
        self.seed = int((self.seed * 1664525 + 1013904223) % (2**31 - 1))
        return val

    def py_next_range(self, min_val: float, max_val: float) -> float:
        return min_val + (max_val - min_val) * self.py_next()

    def py_next_vector2(self, min_val: float, max_val: float) -> Tuple[float, float]:
        return (
            self.py_next_range(min_val, max_val),
            self.py_next_range(min_val, max_val),
        )

# Global RNG instance
_global_rng = None

def get_global_rng() -> RNG:
    """Get global RNG instance"""
    global _global_rng
    if _global_rng is None:
        _global_rng = RNG()
    return _global_rng

def set_global_rng_seed(seed: int):
    """Set seed for global RNG"""
    global _global_rng
    if _global_rng is None:
        _global_rng = RNG(seed)
    else:
        _global_rng.reset(seed)
