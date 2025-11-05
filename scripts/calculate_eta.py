#!/usr/bin/env python3
"""Calculate ETA for simulation"""

# From logs:
# Step 1000: 776.0ms per step
# Step 2000: 790.8ms per step
# Between step 1000-2000: ~13 minutes for 1000 steps

steps_completed = 2000
steps_total = 500000
steps_remaining = steps_total - steps_completed

# Average time per step
avg_time_per_step = 0.79  # seconds (average of 776ms and 790.8ms)

# Or use time between checkpoints
time_between_steps = 13 * 60  # seconds (13 minutes for 1000 steps)
steps_between_checkpoints = 1000
time_per_step = time_between_steps / steps_between_checkpoints

print("=" * 70)
print("SIMULATION ETA CALCULATION")
print("=" * 70)
print(f"\nSteps completed: {steps_completed:,} / {steps_total:,}")
print(f"Steps remaining: {steps_remaining:,}")
print(f"\nPerformance metrics:")
print(f"  Time per step (checkpoint): {time_per_step:.2f} seconds")
print(f"  Steps per second: {1/time_per_step:.2f}")
print(f"  Steps per minute: {60/time_per_step:.1f}")
print(f"\nTime remaining:")
remaining_seconds = steps_remaining * time_per_step
remaining_minutes = remaining_seconds / 60
remaining_hours = remaining_minutes / 60
remaining_days = remaining_hours / 24

print(f"  {remaining_seconds:,.0f} seconds")
print(f"  {remaining_minutes:,.1f} minutes")
print(f"  {remaining_hours:,.1f} hours")
print(f"  {remaining_days:,.2f} days")
print(f"\nETA: ~{int(remaining_hours)} hours ({remaining_days:.1f} days)")
print("=" * 70)

