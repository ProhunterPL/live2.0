#!/usr/bin/env python3
"""
Quick AWS Status Check
======================

Commands to check Phase 2B simulation status on AWS.
Run these commands on AWS instance via SSH.
"""

print("""
🔍 SPRAWDZANIE STATUSU PHASE 2B NA AWS
=======================================

1. Sprawdź procesy Python (szczegółowo):
   ps aux | grep python | grep -v grep

2. Sprawdź tylko procesy związane z symulacjami:
   ps aux | grep -E "run_phase2|run_phase2b" | grep -v grep

3. Sprawdź postęp symulacji:
   cd ~/live2.0
   python3 aws_test/scripts/check_phase2b_progress.py --results-dir results/phase2b_additional

4. Sprawdź ostatnie kroki w logach:
   tail -5 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_*/simulation.log | grep "Step"

5. Sprawdź kiedy ostatnio były aktualizowane logi:
   find ~/live2.0/results/phase2b_additional -name "simulation.log" -exec ls -lh {} \\; | awk '{print $6, $7, $8, $9}'

6. Sprawdź zużycie CPU przez procesy Python:
   top -b -n 1 | grep python

7. Sprawdź ile procesów Python działa:
   ps aux | grep python | grep -v grep | wc -l

8. Sprawdź czy master runner działa:
   ps aux | grep run_phase2b_additional | grep -v grep

=======================================
""")

