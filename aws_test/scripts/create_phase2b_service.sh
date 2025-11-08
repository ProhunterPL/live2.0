#!/bin/bash
# Script to create systemd service for Phase 2B simulations

SERVICE_FILE="/etc/systemd/system/phase2b.service"
WORK_DIR="/home/ubuntu/live2.0/aws_test"
PYTHON_CMD="/usr/bin/python3"
SCRIPT_PATH="$WORK_DIR/run_phase2b_master.py"
LOG_FILE="$WORK_DIR/phase2b_service.log"

echo "Creating systemd service file..."

sudo tee $SERVICE_FILE > /dev/null <<EOF
[Unit]
Description=Phase 2B Simulations
After=network.target

[Service]
Type=simple
User=ubuntu
WorkingDirectory=$WORK_DIR
ExecStart=$PYTHON_CMD $SCRIPT_PATH --mode run
Restart=on-failure
RestartSec=10
StandardOutput=append:$LOG_FILE
StandardError=append:$LOG_FILE

[Install]
WantedBy=multi-user.target
EOF

echo "Service file created: $SERVICE_FILE"
echo ""
echo "Now run:"
echo "  sudo systemctl daemon-reload"
echo "  sudo systemctl enable phase2b"
echo "  sudo systemctl start phase2b"
echo "  sudo systemctl status phase2b"

