#!/usr/bin/env python3
"""
Check system logs for OOM killer or process kills
"""
import subprocess
import os

print("=" * 80)
print("🔍 CHECKING SYSTEM LOGS FOR PROCESS KILLS")
print("=" * 80)

# Check dmesg for OOM kills
print("\n📋 Checking dmesg for OOM kills:")
print("-" * 80)
try:
    result = subprocess.run(['dmesg', '-T'], capture_output=True, text=True, timeout=5)
    oom_lines = [line for line in result.stdout.split('\n') if 'oom' in line.lower() or 'killed process' in line.lower()]
    if oom_lines:
        print("⚠️  Found OOM/kill messages:")
        for line in oom_lines[-10:]:  # Last 10
            print(f"  {line}")
    else:
        print("✅ No OOM kills found in dmesg")
except Exception as e:
    print(f"⚠️  Could not check dmesg: {e}")

# Check system logs
print("\n📋 Checking system logs:")
print("-" * 80)
log_files = [
    '/var/log/syslog',
    '/var/log/kern.log',
    '/var/log/messages'
]

for log_file in log_files:
    if os.path.exists(log_file):
        try:
            result = subprocess.run(['tail', '-50', log_file], capture_output=True, text=True, timeout=5)
            kill_lines = [line for line in result.stdout.split('\n') 
                         if 'python' in line.lower() and ('killed' in line.lower() or 'oom' in line.lower() or 'signal' in line.lower())]
            if kill_lines:
                print(f"\n⚠️  Found kill messages in {log_file}:")
                for line in kill_lines[-5:]:
                    print(f"  {line}")
        except Exception as e:
            pass

# Check if there are any Python processes that were killed recently
print("\n📋 Checking process history:")
print("-" * 80)
print("Run these commands manually to check:")
print("  sudo journalctl -u systemd --since '5 hours ago' | grep -i 'killed\\|oom\\|python'")
print("  dmesg | grep -i 'killed process' | tail -20")

print("\n" + "=" * 80)
print("💡 Most likely causes:")
print("  1. SSH connection dropped and processes were killed")
print("  2. Out of Memory (OOM) killer terminated processes")
print("  3. System reboot or instance restart")
print("  4. Process timeout or resource limit")
print("=" * 80)

