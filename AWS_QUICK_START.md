# AWS Quick Start - Instrukcje krok po kroku

## 🚀 Krok 1: Połącz się z instancją

```bash
ssh -i twoj-klucz.pem ubuntu@<TWOJ-IP-AWS>
```

## 📦 Krok 2: Skopiuj skrypty instalacyjne

**Opcja A: Pobierz z repozytorium**
```bash
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/setup_aws_instance.sh
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/test_aws_instance.sh
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/run_aws_production.sh
wget https://raw.githubusercontent.com/ProhunterPL/live2.0/main/monitor_aws_runs.sh
chmod +x *.sh
```

**Opcja B: Ręczne skopiowanie**
Z lokalnego komputera:
```bash
scp -i twoj-klucz.pem setup_aws_instance.sh ubuntu@<IP>:~/
scp -i twoj-klucz.pem test_aws_instance.sh ubuntu@<IP>:~/
scp -i twoj-klucz.pem run_aws_production.sh ubuntu@<IP>:~/
scp -i twoj-klucz.pem monitor_aws_runs.sh ubuntu@<IP>:~/
```

Na instancji AWS:
```bash
chmod +x *.sh
```

## ⚙️ Krok 3: Instalacja środowiska

```bash
bash setup_aws_instance.sh
```

To zainstaluje:
- Python 3.11
- Git
- Wszystkie zależności z requirements.txt
- Sklonuje repozytorium

**Czas: ~5-10 minut**

## ✅ Krok 4: Test wydajności

```bash
bash test_aws_instance.sh
```

To uruchomi:
- Test symulacji (1000 kroków)
- Sprawdzi czy wszystko działa
- Pokaże wydajność (powinna być 4-6 steps/s)

**Czas: ~5 minut**

Sprawdź wyniki:
```bash
cat ~/live2.0/results/aws_test/summary.txt
```

## 🎯 Krok 5: Produkcyjne uruchomienie

### Opcja A: Automatyczny (ZALECANE)
```bash
bash run_aws_production.sh
```

Skrypt automatycznie:
- Wykryje liczbę CPU
- Ustawi optymalną liczbę równoległych zadań
- Uruchomi wszystkie scenariusze

### Opcja B: Ręczne ustawienie liczby zadań równoległych
```bash
# Dla instancji z 64 CPU -> 16 zadań równoległych
bash run_aws_production.sh 16

# Dla instancji z 32 CPU -> 8 zadań równoległych
bash run_aws_production.sh 8
```

## 📊 Krok 6: Monitorowanie

### Terminal 1: Monitor zasobów
```bash
htop
```
Powinieneś zobaczyć:
- ✅ CPU: 90-100% na wszystkich rdzeniach
- ✅ Memory: <80%

### Terminal 2: Monitor postępu
```bash
bash monitor_aws_runs.sh
```

Lub:
```bash
# Logi na żywo
tail -f ~/live2.0/results/phase2_aws/master.log
```

## 💾 Krok 7: Pobieranie wyników

### Podczas działania symulacji (w innym terminalu lokalnie):
```bash
# Sprawdź ile miejsca zajmują wyniki
ssh -i twoj-klucz.pem ubuntu@<IP> "du -sh ~/live2.0/results/"
```

### Po zakończeniu:
```bash
# Pobierz wszystkie wyniki
scp -r -i twoj-klucz.pem ubuntu@<IP>:~/live2.0/results/ ./aws_results/

# Lub prześlij do S3
ssh -i twoj-klucz.pem ubuntu@<IP>
cd live2.0
aws s3 sync results/ s3://twoj-bucket/phase2-results/
```

## 🧹 Krok 8: Cleanup

Po pobraniu wyników:
```bash
# Zakończ instancję w AWS Console
# LUB przez CLI:
aws ec2 terminate-instances --instance-ids <INSTANCE-ID>
```

---

## 📋 Szybki Checklist

- [ ] Połączono z instancją przez SSH
- [ ] Skopiowano skrypty
- [ ] Uruchomiono `setup_aws_instance.sh`
- [ ] Uruchomiono `test_aws_instance.sh`
- [ ] Sprawdzono wydajność (>4 steps/s)
- [ ] Uruchomiono `run_aws_production.sh`
- [ ] Monitoruje się postęp (htop, logs)
- [ ] Pobrano wyniki
- [ ] Zakończono instancję

---

## 🆘 Troubleshooting

### Problem: Niska wydajność (<3 steps/s)
```bash
# Sprawdź typ instancji
curl http://169.254.169.254/latest/meta-data/instance-type

# Powinna być c6i.16xlarge lub podobna (compute-optimized)
```

### Problem: Brak pamięci
```bash
# Zmniejsz liczbę równoległych zadań
bash run_aws_production.sh 8  # zamiast 16
```

### Problem: Skrypty nie działają
```bash
# Upewnij się że są wykonywalne
chmod +x *.sh

# Uruchom z bash
bash setup_aws_instance.sh
```

---

## 📞 Szybkie komendy

```bash
# Status symulacji
bash monitor_aws_runs.sh

# Logi na żywo
tail -f ~/live2.0/results/phase2_aws/master.log

# Zasoby systemu
htop

# Ile zajmują wyniki
du -sh ~/live2.0/results/

# Sprawdź postęp
find ~/live2.0/results -name "summary.txt" | wc -l
```

---

## ⏱️ Oczekiwane czasy (dla 200k kroków)

| Instancja | Równoległe | Czas dla 150 symulacji |
|-----------|------------|------------------------|
| c6i.8xlarge (32 CPU) | 8 | ~4.8 dni |
| c6i.16xlarge (64 CPU) | 16 | **~2.4 dni** |
| c6i.32xlarge (128 CPU) | 32 | ~1.2 dni |

---

## 💰 Szacowany koszt

- **c6i.16xlarge**: ~$2.72/godz = ~$65/dzień
- **2.4 dni**: ~**$157 total**
- **Spot instances**: ~$15-80 (oszczędność 50-90%)

---

**Powodzenia! 🚀**

