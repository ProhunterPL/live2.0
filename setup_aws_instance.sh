#!/bin/bash
# AWS Instance Setup Script for Phase 2 Simulations
# Usage: bash setup_aws_instance.sh

set -e  # Exit on error

echo "================================"
echo "AWS Instance Setup - Phase 2"
echo "================================"

# Update system
echo "1. Updating system..."
sudo apt update && sudo apt upgrade -y

# Install dependencies
echo "2. Installing dependencies..."
sudo apt install -y python3.11 python3-pip git htop

# Clone repository
echo "3. Cloning repository..."
if [ ! -d "live2.0" ]; then
    git clone https://github.com/ProhunterPL/live2.0.git
else
    echo "Repository already exists, pulling latest changes..."
    cd live2.0
    git pull
    cd ..
fi

cd live2.0

# Install Python packages
echo "4. Installing Python packages..."
pip3 install --upgrade pip
pip3 install -r requirements.txt

# Verify installation
echo "5. Verifying installation..."
python3 -c "import taichi as ti; print(f'Taichi version: {ti.__version__}')"
python3 -c "import numpy as np; print(f'NumPy version: {np.__version__}')"

echo ""
echo "================================"
echo "✅ Setup complete!"
echo "================================"
echo ""
echo "Next steps:"
echo "1. Test run:    bash test_aws_instance.sh"
echo "2. Production:  bash run_aws_production.sh"
echo ""

