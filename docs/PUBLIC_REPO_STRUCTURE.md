---
date: 2025-12-04
label: guide
---

# Public Repository Structure Guide

**Purpose**: This document explains what is included in the public repository and what is intentionally excluded.

---

## 📋 Overview

This repository was made public as part of the paper submission process and is linked via DOI. The repository contains:

- ✅ **Core simulation engine** (`backend/sim/`)
- ✅ **Configuration files** (`aws_test/configs/`, `configs/`)
- ✅ **Analysis scripts** (`scripts/`)
- ✅ **Documentation** (`docs/`)
- ✅ **Paper materials** (`paper/`)
- ✅ **Archive** (`archive/`) - historical reference only

---

## 🚫 What is NOT in the Public Repository

The following types of files are intentionally excluded from the public repository to protect strategic information, future research directions, and competitive advantages:

### 1. Work-in-Progress Files
- `AWS_*.txt` - AWS deployment notes (may contain IP addresses, commands)
- `*DEPLOY*.txt`, `*DEPLOY*.sh` - Deployment command files
- `HYDRO_*.txt` - Hydrothermal setup notes
- `structure.txt`, `docs_structure.txt` - Temporary structure files

### 2. Local PowerShell Scripts
- `*.ps1` files in root directory and `scripts/` (except those in `archive/`)
- These are local development tools, not needed for public use

### 3. AWS Test Results (Duplicates)
- `aws_test/results/` - Duplicate results (main results are in `results/`)
- `aws_test/results_*/` - Backup/duplicate result directories

### 4. Backup Files
- `*.backup` - Backup copies of scripts
- `*_backup.py` - Backup Python files

### 5. Strategic Plans and Future Publications
- `docs/plans/` - All strategic development plans and roadmaps
- `docs/wniosek/` - Grant proposals and funding applications
- `docs/LIVE2_QUANTUM_AI_EXPANSION.md` - Future expansion plans
- `docs/VALIDATION_ROADMAP.md` - Validation roadmap
- `docs/phase3/PAPER2_OUTLINE.md` - Future publication outlines
- `paper/POST_SUBMISSION_PLAN.md` - Post-submission strategy
- `paper/QUANTUM_AI_EXPANSION_ANALYSIS.md` - Expansion analysis
- `paper/WORK_PLAN.md` - Work plans

### 6. Editor Configuration (Project-Specific)
- `.cursor/` - Cursor editor rules and project-specific configurations
- `.cursorrules` - Cursor editor rules file
- These contain project strategies, architecture decisions, and development workflows

### 7. Sensitive Information
- SSH keys (`.pem`, `.key` files)
- Environment files (`.env`)
- Hardcoded IP addresses or credentials in scripts

---

## 📁 Repository Structure

```
live2.0/
├── backend/              # Core simulation engine
│   └── sim/             # Simulation code (Taichi kernels, physics, chemistry)
├── aws_test/            # AWS deployment scripts and configs
│   ├── configs/        # Simulation configuration files (YAML)
│   └── scripts/        # AWS deployment and monitoring scripts
├── scripts/             # Analysis and utility scripts
├── docs/                # Documentation
│   ├── aws_test/       # AWS deployment documentation
│   ├── phase2b/        # Phase 2B specific documentation
│   └── technical/      # Technical documentation
├── paper/               # Paper materials (manuscript, figures, data)
├── archive/             # Historical reference (old scripts, deprecated code)
└── results/             # Simulation results (if included)
```

---

## 🔒 Security Considerations

### What We Checked
- ✅ No AWS credentials (access keys, secret keys)
- ✅ No API keys or tokens
- ✅ No SSH private keys
- ✅ No hardcoded IP addresses in active scripts
- ✅ No passwords or sensitive data

### Scripts Using Credentials
Scripts that require AWS credentials use:
- **Environment variables**: `AWS_IP`, `AWS_SSH_KEY_PATH`
- **Command-line arguments**: `--host`, `--key`
- **Never hardcoded values** in active scripts

### Example Usage
```bash
# Set environment variables
export AWS_IP="your-aws-ip"
export AWS_SSH_KEY_PATH="/path/to/key.pem"

# Run script
python3 aws_test/scripts/final_status_check.py
```

---

## 📊 Data Availability

### What Data is Available
- ✅ **Configuration files** - All simulation configs (YAML)
- ✅ **Analysis scripts** - Complete analysis pipeline
- ✅ **Documentation** - Full documentation of methods
- ⚠️ **Results** - May be excluded due to size (see `.gitignore`)

### How to Reproduce Results
1. Use configuration files from `aws_test/configs/`
2. Run simulation using `scripts/run_phase2_full.py`
3. Analyze results using scripts in `scripts/`
4. See `docs/` for detailed instructions

---

## 🛠️ For Reviewers

### If You Need Full Data
- Contact the authors for access to full simulation results
- Results are large (500K+ steps × 43 runs × multiple scenarios)
- Main results are available in `paper/results_data/` for figures

### If You Want to Run Simulations
1. See `README.md` for setup instructions
2. Use configs from `aws_test/configs/*_SUPER_FAST.yaml` for quick tests
3. See `docs/aws_test/` for AWS deployment guide
4. See `docs/local/` for local development guide

---

## 📝 Notes

### Archive Directory
The `archive/` directory contains:
- Historical scripts (may contain old hardcoded paths - these are archived for reference only)
- Deprecated code
- One-off scripts used during development

**Note**: Files in `archive/` are kept for historical reference but are not actively used.

### Documentation
All documentation in `docs/` is public and may contain:
- Example commands with placeholders (`<AWS_IP>`, `<your-key.pem>`)
- Historical notes about deployment
- Troubleshooting guides

These are safe to share as they use placeholders, not actual credentials.

---

## ✅ Verification

To verify what is tracked in git:
```bash
git ls-files | grep -E "(AWS_|DEPLOY|HYDRO_|\.ps1$|\.backup)"
```

If this returns any files (except in `archive/`), they should be removed.

---

**Last Updated**: 2025-12-04  
**Status**: Repository cleaned and ready for public access

