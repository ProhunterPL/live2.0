#!/bin/bash
# Backup completed simulation results to compressed archive
# =========================================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/live2.0}"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional"
BACKUP_DIR="$PROJECT_ROOT/backups"

cd "$PROJECT_ROOT" || exit 1

echo "================================================================================"
echo "📦 BACKUP COMPLETED RESULTS"
echo "================================================================================"
echo ""

# Create backup directory
mkdir -p "$BACKUP_DIR"

# Find all completed runs
COMPLETED_RUNS=()
for scenario in miller_urey_extended hydrothermal_extended; do
    SCENARIO_DIR="$RESULTS_DIR/$scenario"
    if [ -d "$SCENARIO_DIR" ]; then
        for run_dir in "$SCENARIO_DIR"/run_*; do
            if [ -f "$run_dir/results.json" ]; then
                RUN_NAME=$(basename "$run_dir")
                COMPLETED_RUNS+=("$scenario/$RUN_NAME")
            fi
        done
    fi
done

if [ ${#COMPLETED_RUNS[@]} -eq 0 ]; then
    echo "❌ No completed runs found!"
    exit 1
fi

echo "Found ${#COMPLETED_RUNS[@]} completed runs:"
for run in "${COMPLETED_RUNS[@]}"; do
    echo "  ✅ $run"
done
echo ""

# Generate timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
ARCHIVE_NAME="phase2b_completed_${TIMESTAMP}.tar.gz"
ARCHIVE_PATH="$BACKUP_DIR/$ARCHIVE_NAME"

echo "Creating archive: $ARCHIVE_NAME"
echo ""

# Create tar command with all completed runs
cd "$RESULTS_DIR" || exit 1
tar -czf "$ARCHIVE_PATH" "${COMPLETED_RUNS[@]}"

if [ $? -eq 0 ]; then
    ARCHIVE_SIZE=$(du -h "$ARCHIVE_PATH" | awk '{print $1}')
    echo "✅ Archive created successfully!"
    echo ""
    echo "📊 Details:"
    echo "  Path: $ARCHIVE_PATH"
    echo "  Size: $ARCHIVE_SIZE"
    echo "  Runs: ${#COMPLETED_RUNS[@]}"
    echo ""
    echo "📥 To download to local machine:"
    echo "  scp -i your-key.pem ubuntu@<AWS-IP>:$ARCHIVE_PATH ."
    echo ""
    echo "📦 To extract locally:"
    echo "  tar -xzf $(basename "$ARCHIVE_PATH") -C results/phase2b_additional/"
    echo ""
else
    echo "❌ Archive creation failed!"
    exit 1
fi

# Optional: Create a manifest file
MANIFEST_PATH="$BACKUP_DIR/manifest_${TIMESTAMP}.txt"
echo "Backup created: $(date)" > "$MANIFEST_PATH"
echo "Archive: $ARCHIVE_NAME" >> "$MANIFEST_PATH"
echo "Runs included:" >> "$MANIFEST_PATH"
for run in "${COMPLETED_RUNS[@]}"; do
    echo "  - $run" >> "$MANIFEST_PATH"
done

echo "📝 Manifest created: $MANIFEST_PATH"
echo ""
echo "================================================================================"
echo "✅ DONE"
echo "================================================================================"

