# Naprawa Skryptów Phase 2B na AWS

## Problem
Skrypty w `aws_test/scripts/` używają `"python"` zamiast `"python3"`.

## Szybka Naprawa - Na AWS:

```bash
# Napraw wszystkie skrypty na raz
cd ~/live2.0/aws_test/scripts

# Napraw debug_formamide.py
sed -i 's/"python"/"python3"/g' debug_formamide.py

# Napraw run_phase2b_additional.py  
sed -i 's/"python"/"python3"/g' run_phase2b_additional.py

# Napraw download_phase2b_results.py
sed -i 's/"python"/"python3"/g' download_phase2b_results.py

# Sprawdź czy są jeszcze wywołania python (nie python3)
grep -n '"python"' *.py | grep -v python3
```

## Następnie uruchom ponownie:
```bash
cd ~/live2.0/aws_test
python3 run_phase2b_master.py --mode all
```

## Sprawdź czy faktycznie były uruchomione:

```bash
# Sprawdź foldery wyników
ls -la results/phase2b_additional/
ls -la results/phase2b_additional/miller_urey_extended/
ls -la results/phase2b_additional/formamide_debug/

# Sprawdź logi
tail -50 results/phase2b_additional/logs/phase2b_runner.log

# Sprawdź czy są pliki molecules.json
find results/phase2b_additional -name "molecules.json" | wc -l
```

