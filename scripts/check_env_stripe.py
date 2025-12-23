#!/usr/bin/env python3
"""
Check Stripe environment variables from .env file.
"""

import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env
project_root = Path(__file__).parent.parent
env_file = project_root / ".env"
load_dotenv(env_file)

print("=" * 60)
print("Stripe Environment Variables Check")
print("=" * 60)
print(f"\n.env file location: {env_file}")
print(f".env exists: {env_file.exists()}")
print()

# Check each variable
vars_to_check = [
    "STRIPE_SECRET_KEY",
    "STRIPE_WEBHOOK_SECRET",
    "STRIPE_PRICE_ID_HOBBY",
    "STRIPE_PRICE_ID_RESEARCH",
    "STRIPE_PRICE_ID_PRO",
]

for var_name in vars_to_check:
    value = os.getenv(var_name, "")
    if value:
        # Show first/last chars for security
        if len(value) > 20:
            display = f"{value[:10]}...{value[-10:]}"
        else:
            display = value
        print(f"{var_name}:")
        print(f"  Status: SET")
        print(f"  Length: {len(value)}")
        print(f"  Preview: {display}")
        print(f"  Starts with expected prefix: ", end="")
        if var_name == "STRIPE_WEBHOOK_SECRET":
            print("YES" if value.startswith("whsec_") else f"NO (expected 'whsec_', got '{value[:6]}')")
        elif var_name == "STRIPE_SECRET_KEY":
            print("YES" if value.startswith("sk_") else f"NO (expected 'sk_', got '{value[:3]}')")
        else:
            print("YES" if value.startswith("price_") else f"NO (expected 'price_', got '{value[:6]}')")
    else:
        print(f"{var_name}:")
        print(f"  Status: NOT SET")
    print()

# Also check raw .env file content (if exists)
if env_file.exists():
    print("=" * 60)
    print("Raw .env file content (STRIPE_WEBHOOK_SECRET line):")
    print("=" * 60)
    with open(env_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            if "STRIPE_WEBHOOK_SECRET" in line:
                # Mask the secret value
                if "=" in line:
                    key, value = line.split("=", 1)
                    masked_value = value[:10] + "..." + value[-10:] if len(value) > 20 else "***"
                    print(f"Line {line_num}: {key}={masked_value}")
                    print(f"  Full line length: {len(line)}")
                    print(f"  Has trailing whitespace: {line.rstrip() != line}")
                    print(f"  Value length: {len(value.strip())}")

