# Naprawa physics_parameters.json

## Problem
Symulacja nie może znaleźć `physics_parameters.json`.

## Sprawdź co jest w data/:

```bash
find ~/live2.0/data -name "*.json"
ls -la ~/live2.0/data/
```

## Jeśli brakuje pliku, stwórz prosty fallback

```bash
cd ~/live2.0/data

# Stwórz minimalny physics_parameters.json
cat > physics_parameters.json << 'EOF'
{
    "lj_epsilon": {
        "H-H": 0.12,
        "C-C": 0.15,
        "N-N": 0.13,
        "O-O": 0.15
    },
    "lj_sigma": {
        "H-H": 2.4,
        "C-C": 3.5,
        "N-N": 3.0,
        "O-O": 3.0
    },
    "bond_params": {
        "C-H": {"k": 350.0, "r0": 1.1},
        "C-C": {"k": 350.0, "r0": 1.5},
        "N-H": {"k": 350.0, "r0": 1.0},
        "O-H": {"k": 350.0, "r0": 0.95}
    }
}
EOF
```

Ale najpierw - niech symulacja się zakończy z fallback parameters!

