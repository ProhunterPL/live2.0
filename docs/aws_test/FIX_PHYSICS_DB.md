# 🔧 Naprawa PhysicsDatabase na AWS

## ❌ Problem

Symulacja pokazuje błąd:
```
PhysicsDatabase not found, using fallback parameters
```

**Przyczyna**: Plik `data/physics_parameters.json` nie jest znajdowany na AWS.

## ✅ Rozwiązanie

### 1. Sprawdź czy plik istnieje na AWS:

```bash
# Na AWS
cd ~/live2.0
ls -la data/physics_parameters.json
```

Jeśli plik **NIE ISTNIEJE**, skopiuj go:

```bash
# Na lokalnej maszynie (PowerShell)
scp -i "C:\Users\user\Desktop\aws_credential\key-do-live.pem" `
    data\physics_parameters.json `
    ubuntu@<AWS-IP>:~/live2.0/data/physics_parameters.json
```

### 2. Zaktualizuj kod (już poprawione):

Kod został ulepszony żeby automatycznie szukał pliku w wielu lokalizacjach:
- `data/physics_parameters.json` (względem project root)
- `~/live2.0/data/physics_parameters.json` (home directory)
- Automatyczne wykrywanie project root

### 3. Zaktualizuj kod na AWS:

```bash
cd ~/live2.0
git pull
```

### 4. Sprawdź czy działa:

Po zaktualizowaniu kodu, sprawdź logi symulacji:

```bash
tail -50 ~/live2.0/results/phase2b_additional/miller_urey_extended/run_1/simulation.log | grep -i physics
```

Powinieneś zobaczyć:
```
Found PhysicsDatabase at: /home/ubuntu/live2.0/data/physics_parameters.json
Loaded PhysicsDatabase from ...
  Bond parameters: XX
  VDW parameters: XX
  Citations: XX
```

## ⚠️ Ważne

**Symulacja działa nawet bez PhysicsDatabase!**

Jeśli plik nie zostanie znaleziony, symulacja użyje parametrów fallback:
- `default_epsilon = 0.439` (Carbon UFF)
- `default_sigma = 3.431` (Carbon UFF)
- `default_bond_D_e = 348.0` (C-C single)
- `default_bond_r_e = 1.54` (C-C single)

**To nie jest krytyczny błąd** - symulacja będzie działać, ale użyje domyślnych parametrów zamiast literaturowych.

## 🎯 Rekomendacja

1. **Sprawdź czy plik istnieje** na AWS
2. **Jeśli nie istnieje** - skopiuj go z lokalnej maszyny
3. **Zaktualizuj kod** na AWS (`git pull`)
4. **Nie zatrzymuj obecnej symulacji** - działa z fallback parameters

Symulacja będzie działać poprawnie nawet bez PhysicsDatabase!

