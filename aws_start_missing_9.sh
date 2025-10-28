#!/bin/bash
# Quick script to start missing 9 simulations with 50k steps

echo "=========================================="
echo "URUCHAMIANIE 9 BRAKUJĄCYCH SYMULACJI"
echo "=========================================="
echo ""

echo "1. Sprawdzam aktualne procesy..."
CURRENT=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)
echo "   Aktualnie: $CURRENT procesów"
echo ""

echo "2. Uruchamiam Miller Urey 9-16 (8 symulacji, 50k steps)..."
for i in {9..16}; do 
    nohup python3 scripts/run_phase2_full.py \
        --config configs/phase2_miller_urey.yaml \
        --output results/miller_urey/run_$i \
        --steps 50000 \
        --seed $((42 + i)) \
        > logs/miller_new_$i.log 2>&1 &
    echo "   ✓ Miller Urey run_$i uruchomione"
done
echo ""

echo "3. Uruchamiam Formamide 2 (1 symulacja, 50k steps)..."
nohup python3 scripts/run_phase2_full.py \
    --config configs/phase2_formamide.yaml \
    --output results/formamide/run_2 \
    --steps 50000 \
    --seed 60 \
    > logs/form_new_2.log 2>&1 &
echo "   ✓ Formamide run_2 uruchomione"
echo ""

sleep 3

echo "4. Sprawdzam nową liczbę procesów..."
NEW=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)
echo "   Teraz: $NEW procesów"
echo ""

if [ $NEW -eq 24 ]; then
    echo "✅ SUKCES! Wszystkie 24 procesy działają"
else
    echo "⚠️  Oczekiwano 24, jest $NEW"
fi
echo ""

echo "5. Pamięć:"
free -h | grep Mem
echo ""

echo "=========================================="
echo "GOTOWE!"
echo "=========================================="
echo ""
echo "Monitoruj nowe symulacje:"
echo "  tail -f logs/miller_new_9.log"
echo ""
echo "Sprawdź prędkość (po 10 min):"
echo "  tail -20 logs/miller_new_9.log | grep 'steps/s'"
echo ""
echo "Sprawdź wszystkie procesy:"
echo "  ps aux | grep run_phase2_full | grep -v grep"
echo ""

