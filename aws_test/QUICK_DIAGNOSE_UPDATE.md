# 🔄 quick_diagnose.py Update - Nov 8, 2025

## 🎯 Problem Solved

**Old behavior**: Script showed "STOPPED" for simulations that were actually running, because it only checked log file timestamps. Logs are buffered, so timestamps can be 60-90 minutes old even while process is actively computing.

**New behavior**: Script now checks CPU usage via `ps aux` to determine real status, not just log timestamps.

---

## ✨ New Features

### **1. Process Detection**:
```python
def get_running_processes():
    """Get CPU usage and info for running simulation processes"""
```

- Runs `ps aux` to get all running processes
- Identifies simulation processes by `run_phase2_full.py`
- Extracts CPU%, memory%, CPU time, PID
- Determines if process is active (CPU > 100%)

### **2. Smart Status Detection**:

**Priority order**:
1. ✅ **Process actively computing** (CPU > 100%) → `🔄 RUNNING`
2. ⏸️ **Process idle** (CPU < 100%) → `⏸️ PAUSED`
3. 🔄 **Recent log update** (< 10 min, no process) → `🔄 RUNNING`
4. ⏸️ **Old log, no process** → `⏸️ STOPPED`

### **3. Enhanced Output**:

**Before**:
```
⏸️ STOPPED miller_urey_extended/run_2: Step 24,000/500,000 (4.8%)
    Last update: 78.4 minutes ago
```

**After**:
```
🔄 RUNNING miller_urey_extended/run_2: Step 24,000/500,000 (4.8%)
    Status: CPU: 1106%, TIME: 1829:18
    ℹ️  Log buffered (78min old) but process actively computing
```

### **4. Process Summary**:

Now shows at the top:
```
💻 ACTIVE PROCESSES:
   Total processes: 4
   Actively computing: 4 (CPU > 100%)
   Total CPU usage: 5387%
```

### **5. Informative Footer**:

```
ℹ️  NOTE: Status based on CPU usage (>100% = actively computing)
   Log timestamps may be delayed due to buffering - this is normal!
```

---

## 🔧 Technical Changes

### **Added**:
- `import subprocess` - to run `ps aux`
- `get_running_processes()` function - parse process list
- Process correlation logic - match processes to simulation runs
- CPU-based status determination
- Log buffering notes

### **Improved**:
- Status detection: Process-first, log-second
- Error handling: Graceful fallback if `ps` fails
- User feedback: Clear explanation of what's happening

---

## 📊 Example Output

```
================================================================================
🔍 PHASE 2B DIAGNOSTICS (Enhanced)
================================================================================
Time: 2025-11-08 21:33:19

💻 ACTIVE PROCESSES:
   Total processes: 4
   Actively computing: 4 (CPU > 100%)
   Total CPU usage: 5387%

📄 results.json files: 0

📊 Status from phase2b_results.json:
   Completed: 0/30
   Failed: 4/30

📝 SIMULATION STATUS:
  🔄 RUNNING miller_urey_extended/run_1: Step 60,000/500,000 (12.0%)
      Status: CPU: 1487%, TIME: 2460:19
  🔄 RUNNING miller_urey_extended/run_2: Step 24,000/500,000 (4.8%)
      Status: CPU: 1106%, TIME: 1829:18
      ℹ️  Log buffered (78min old) but process actively computing
  🔄 RUNNING miller_urey_extended/run_3: Step 60,000/500,000 (12.0%)
      Status: CPU: 1489%, TIME: 2463:04
  🔄 RUNNING miller_urey_extended/run_4: Step 26,000/500,000 (5.2%)
      Status: CPU: 1073%, TIME: 1775:09
      ℹ️  Log buffered (73min old) but process actively computing

================================================================================
ℹ️  NOTE: Status based on CPU usage (>100% = actively computing)
   Log timestamps may be delayed due to buffering - this is normal!
================================================================================
```

---

## 🚀 Usage

### **On AWS**:
```bash
cd ~/live2.0
python3 aws_test/scripts/quick_diagnose.py
```

### **Benefits**:
- ✅ Accurate status (no more false "STOPPED" alerts)
- ✅ Shows CPU usage and time (proof of activity)
- ✅ Explains log buffering (reduces confusion)
- ✅ Works even with delayed log writes

---

## 🔄 Backward Compatibility

- Still shows all old information (steps, progress, errors)
- Just adds CPU info and improves status logic
- Same command line usage
- No breaking changes

---

## 💡 Why This Matters

**Before**: Users would see "STOPPED" and panic, restart service, lose progress  
**After**: Users see "RUNNING" with CPU proof, can relax and let it work

**Key insight**: On Linux with Python logging, file writes are buffered. A process can compute for hours without writing to log. CPU usage is the **only reliable indicator** of activity.

---

**Created**: November 8, 2025  
**Version**: 2.0 (Enhanced with CPU detection)  
**Status**: Ready to deploy

