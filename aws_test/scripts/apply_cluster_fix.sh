#!/bin/bash
# Apply hotfix to disable cluster detection deadlock
# ====================================================

set -e

echo "================================================================================"
echo "🔧 APPLYING CLUSTER DETECTION HOTFIX"
echo "================================================================================"
echo ""

# Backup original file
STEPPER_FILE="$HOME/live2.0/backend/sim/core/stepper.py"
BACKUP_FILE="$HOME/live2.0/backend/sim/core/stepper.py.backup_$(date +%Y%m%d_%H%M%S)"

if [ ! -f "$STEPPER_FILE" ]; then
    echo "❌ ERROR: $STEPPER_FILE not found!"
    exit 1
fi

echo "📦 Creating backup: $BACKUP_FILE"
cp "$STEPPER_FILE" "$BACKUP_FILE"

echo "🔧 Applying inline patch..."

# Use Python to apply the fix (more reliable than sed on different systems)
python3 - <<'EOF'
import re

stepper_file = "backend/sim/core/stepper.py"

with open(stepper_file, 'r') as f:
    content = f.read()

# Find and replace the cluster update section
old_pattern = r'''(\s+)# OPTIMIZATION: Update clusters every 1200 steps - STAGGERED
(\s+)# Offset by 300 to spread load over time
(\s+)if \(self\.step_count - 300\) % 1200 == 0:
(\s+)self\.binding\.update_clusters\(
(\s+)self\.particles\.positions,
(\s+)self\.particles\.active,
(\s+)self\.particles\.particle_count\[None\]
(\s+)\)'''

new_code = r'''\1# OPTIMIZATION: Update clusters - now configurable to prevent infinite loops
\2# Check if cluster detection is enabled (default interval = 1200)
\3cluster_interval = getattr(self.config, 'cluster_check_interval', 1200)
\3if cluster_interval < 999999999:  # Only if not disabled
\3    # Offset by 300 to spread load over time
\3    if (self.step_count - 300) % cluster_interval == 0:
\3        self.binding.update_clusters(
\3            self.particles.positions,
\3            self.particles.active,
\3            self.particles.particle_count[None]
\3        )'''

if re.search(old_pattern, content, re.MULTILINE):
    content = re.sub(old_pattern, new_code, content, flags=re.MULTILINE)
    print("✅ Pattern found and replaced")
else:
    print("⚠️  Pattern not found - checking if already patched...")
    if 'cluster_interval = getattr(self.config' in content:
        print("✅ Already patched!")
    else:
        print("❌ Cannot find expected pattern. Manual intervention needed.")
        exit(1)

with open(stepper_file, 'w') as f:
    f.write(content)

print("✅ File updated successfully")
EOF

echo ""
echo "================================================================================"
echo "✅ HOTFIX APPLIED SUCCESSFULLY"
echo "================================================================================"
echo ""
echo "Backup saved to: $BACKUP_FILE"
echo ""
echo "Now you can use cluster_check_interval: 999999999 in config to disable"
echo "cluster detection completely!"
echo ""
echo "To verify the fix:"
echo "  grep -A 10 'Update clusters - now configurable' $STEPPER_FILE"
echo ""

