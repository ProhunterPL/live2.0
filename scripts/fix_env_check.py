#!/usr/bin/env python3
"""
Check and fix .env file for Stripe webhook secret.
"""

import os
from pathlib import Path

env_file = Path(".env")

if not env_file.exists():
    print("ERROR: .env file not found")
    exit(1)

print("=" * 60)
print("Checking .env file for STRIPE_WEBHOOK_SECRET")
print("=" * 60)

with open(env_file, "r", encoding="utf-8") as f:
    lines = f.readlines()

webhook_lines = []
for i, line in enumerate(lines, 1):
    if "STRIPE_WEBHOOK_SECRET" in line.upper():
        webhook_lines.append((i, line))

print(f"\nFound {len(webhook_lines)} line(s) with STRIPE_WEBHOOK_SECRET:\n")

for line_num, line in webhook_lines:
    stripped = line.strip()
    is_comment = stripped.startswith("#")
    
    print(f"Line {line_num}:")
    print(f"  Is comment: {is_comment}")
    
    if not is_comment and "=" in stripped:
        key, value = stripped.split("=", 1)
        value = value.strip()
        
        print(f"  Key: {key}")
        print(f"  Value length: {len(value)}")
        print(f"  Value preview: {value[:20]}...{value[-10:] if len(value) > 30 else ''}")
        
        if value == "***" or len(value) < 10:
            print(f"  STATUS: [ERROR] Value is masked or too short")
            print(f"  ACTION: Replace with actual webhook secret from Stripe")
        elif value.startswith("whsec_"):
            print(f"  STATUS: [OK] Format looks correct")
        else:
            print(f"  STATUS: [WARNING] Format might be incorrect (should start with 'whsec_')")
    else:
        print(f"  STATUS: {'Commented out' if is_comment else 'No value'}")
    print()

print("=" * 60)
print("If STRIPE_WEBHOOK_SECRET is masked as '***', you need to:")
print("1. Get webhook secret from Stripe (see docs/guides/STRIPE_WEBHOOK_SECRET_FIX.md)")
print("2. Replace '***' with actual value in .env file")
print("=" * 60)

