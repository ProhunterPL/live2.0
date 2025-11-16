#!/bin/bash
# Skrypt do uruchamiania testów lokalnie tak samo jak w CI
# Użycie: bash .github/scripts/run_local_ci.sh

set -e  # Zatrzymaj przy pierwszym błędzie

echo "================================"
echo "🧪 Live 2.0 Local CI Tests"
echo "================================"
echo ""

# Kolory dla outputu
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Sprawdź czy jesteśmy w głównym katalogu projektu
if [ ! -f "requirements.txt" ]; then
    echo -e "${RED}❌ Error: Must be run from project root${NC}"
    exit 1
fi

# Ustaw Taichi na CPU mode
export TI_ARCH=cpu
export PYTHONPATH=$(pwd)

echo -e "${YELLOW}📦 Step 1: Checking dependencies...${NC}"
if ! python -c "import pytest" 2>/dev/null; then
    echo -e "${YELLOW}Installing dependencies...${NC}"
    pip install -r requirements.txt
fi

echo ""
echo -e "${YELLOW}🎨 Step 2: Code Quality Checks${NC}"
echo "--------------------------------"

# Black
echo -e "\n${YELLOW}→ Checking code formatting (black)...${NC}"
if black --check backend/ scripts/ matcher/ 2>/dev/null; then
    echo -e "${GREEN}✓ Code formatting OK${NC}"
else
    echo -e "${RED}✗ Code formatting issues found${NC}"
    echo -e "  Run: ${YELLOW}black backend/ scripts/ matcher/${NC} to fix"
fi

# isort
echo -e "\n${YELLOW}→ Checking import sorting (isort)...${NC}"
if isort --check-only backend/ scripts/ matcher/ 2>/dev/null; then
    echo -e "${GREEN}✓ Import sorting OK${NC}"
else
    echo -e "${RED}✗ Import sorting issues found${NC}"
    echo -e "  Run: ${YELLOW}isort backend/ scripts/ matcher/${NC} to fix"
fi

# mypy (non-blocking)
echo -e "\n${YELLOW}→ Type checking (mypy)...${NC}"
if mypy backend/sim/ --ignore-missing-imports 2>/dev/null; then
    echo -e "${GREEN}✓ Type checking OK${NC}"
else
    echo -e "${YELLOW}⚠ Type checking found issues (non-blocking)${NC}"
fi

echo ""
echo -e "${YELLOW}🧪 Step 3: Unit Tests${NC}"
echo "--------------------------------"

# Backend tests
echo -e "\n${YELLOW}→ Running backend tests (excluding slow tests)...${NC}"
cd backend
if pytest tests/ -v -m "not slow" --tb=short --color=yes --maxfail=5; then
    echo -e "${GREEN}✓ Backend tests passed${NC}"
else
    echo -e "${RED}✗ Backend tests failed${NC}"
    cd ..
    exit 1
fi
cd ..

# Root tests
echo -e "\n${YELLOW}→ Running root tests (excluding stability tests)...${NC}"
if pytest tests/ -v -k "not stability and not 24h" --tb=short --color=yes --maxfail=5; then
    echo -e "${GREEN}✓ Root tests passed${NC}"
else
    echo -e "${RED}✗ Root tests failed${NC}"
    exit 1
fi

echo ""
echo "================================"
echo -e "${GREEN}✅ All CI checks passed!${NC}"
echo "================================"
echo ""
echo -e "Your code is ready to push to main! 🚀"

