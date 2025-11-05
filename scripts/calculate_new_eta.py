#!/usr/bin/env python3
"""Calculate new ETA based on CPU performance"""

steps_total = 500000
steps_completed = 1000
steps_remaining = steps_total - steps_completed

# CPU performance from logs
time_per_step = 0.1881  # seconds (188.1ms)

print("=" * 70)
print("CPU PERFORMANCE - NEW ETA")
print("=" * 70)
print(f"\nPerformance:")
print(f"  Time per step: {time_per_step*1000:.1f} ms")
print(f"  Steps per second: {1/time_per_step:.1f}")
print(f"  Steps per minute: {60/time_per_step:.0f}")

print(f"\nProgress:")
print(f"  Steps completed: {steps_completed:,} / {steps_total:,}")
print(f"  Steps remaining: {steps_remaining:,}")

print(f"\nTime remaining:")
remaining_seconds = steps_remaining * time_per_step
remaining_minutes = remaining_seconds / 60
remaining_hours = remaining_minutes / 60
remaining_days = remaining_hours / 24

print(f"  {remaining_seconds:,.0f} seconds")
print(f"  {remaining_minutes:,.0f} minutes")
print(f"  {remaining_hours:,.1f} hours")
print(f"  {remaining_days:.2f} days")

print(f"\n{'='*70}")
print(f"ETA: ~{int(remaining_hours)} hours ({remaining_days:.1f} days)")
print("=" * 70)

print(f"\nComparison:")
print(f"  Previous (GPU): 4.5 days")
print(f"  Current (CPU):  {remaining_days:.1f} days")
print(f"  Improvement: {4.5/remaining_days:.1f}x FASTER!")

