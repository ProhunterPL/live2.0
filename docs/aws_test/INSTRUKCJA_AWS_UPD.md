# Instrukcja - Naprawa na AWS

## Szybka Naprawa (Zalecana)

### Krok 1: Połącz się z AWS
```bash
ssh -i ~/.ssh/your-key.pem ubuntu@<twoj-ip>
```

### Krok 2: Napraw plik bezpośrednio na AWS
```bash
cd ~/live2.0/aws_test

# Napraw wszystkie wywołania python -> python3
sed -i 's/"python"/"python3"/g' run_phase2b_master.py

# Sprawdź poprawki
grep -n '"python' run_phase2b_master.py
# Powinno pokazać 0 wyników (lub tylko w komentarzach)
```

### Krok 3: Uruchom Phase 2B
```bash
python3 run_phase2b_master.py --mode all
```

---

## Alternatywnie: Wgraj przez Git

### Na lokalnym komputerze:
```bash
cd D:\live2.0

# Commit zmian
git add aws_test/run_phase2b_master.py
git commit -m "Fix: python -> python3 dla AWS Ubuntu"
git push origin main
```

### Na AWS:
```bash
cd ~/live2.0
git pull origin main
cd aws_test
python3 run_phase2b_master.py --mode all
```

---

## Co Wybrać?

- **Opcja 1 (sed)**: Szybka, doraźna naprawa ✅
- **Opcja 2 (git)**: Jeśli chcesz zaktualizować repo

**Rekomendacja**: Opcja 1 (sed) - szybsze i wystarczające.

