#!/bin/bash
# Diagnostyka Round 1 - sprawdź co poszło nie tak

echo "=========================================="
echo "DIAGNOSTYKA ROUND 1"
echo "=========================================="
echo ""

echo "1. KTÓRE LOGI ISTNIEJĄ:"
echo "------------------------"
ls -lh logs/*.log 2>/dev/null | wc -l
echo ""

echo "2. KTÓRE PROCESY DZIAŁAJĄ:"
echo "------------------------"
ps aux | grep run_phase2_full | grep -v grep
echo ""

echo "3. KTÓRE SUMMARY ISTNIEJĄ:"
echo "------------------------"
find results -name "summary.txt" -exec echo {} \;
echo ""

echo "4. SPRAWDŹ OSTATNIE 3 LINIE KAŻDEGO LOGA (szukaj błędów):"
echo "------------------------"
for log in logs/miller_*.log logs/hydro_*.log logs/form_*.log; do
    if [ -f "$log" ]; then
        echo "=== $log ==="
        tail -3 "$log"
        echo ""
    fi
done

echo ""
echo "5. SPRAWDŹ CZY BYŁY OUT OF MEMORY:"
echo "------------------------"
dmesg | grep -i "out of memory" | tail -10

echo ""
echo "6. SPRAWDŹ CZY BYŁY KILL:"
echo "------------------------"
dmesg | grep -i "killed process" | tail -10

echo ""
echo "7. AKTUALNE ZUŻYCIE PAMIĘCI:"
echo "------------------------"
free -h

echo ""
echo "8. PRĘDKOŚĆ SYMULACJI (steps/s):"
echo "------------------------"
for log in logs/miller_*.log logs/hydro_*.log logs/form_*.log; do
    if [ -f "$log" ]; then
        speed=$(tail -20 "$log" | grep "steps/s" | tail -1 | grep -oP '\d+\.\d+ steps/s')
        if [ ! -z "$speed" ]; then
            echo "$log: $speed"
        fi
    fi
done

echo ""
echo "=========================================="
echo "DIAGNOZA ZAKOŃCZONA"
echo "=========================================="

