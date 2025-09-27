# Live 2.0 Makefile

.PHONY: help install dev build test clean docker-up docker-down

# Default target
help:
	@echo "Live 2.0 - Available commands:"
	@echo ""
	@echo "Development:"
	@echo "  install     - Install all dependencies"
	@echo "  dev         - Start development servers"
	@echo "  build       - Build production versions"
	@echo "  test        - Run all tests"
	@echo "  clean       - Clean build artifacts"
	@echo ""
	@echo "Docker:"
	@echo "  docker-up   - Start Docker containers"
	@echo "  docker-down - Stop Docker containers"
	@echo ""
	@echo "Backend:"
	@echo "  backend-install - Install backend dependencies"
	@echo "  backend-dev     - Start backend server"
	@echo "  backend-test    - Run backend tests"
	@echo "  backend-lint    - Lint backend code"
	@echo ""
	@echo "Frontend:"
	@echo "  frontend-install - Install frontend dependencies"
	@echo "  frontend-dev     - Start frontend server"
	@echo "  frontend-test    - Run frontend tests"
	@echo "  frontend-lint    - Lint frontend code"

# Install all dependencies
install: backend-install frontend-install

# Development
dev:
	@echo "Starting development servers..."
	@echo "Backend: http://localhost:8000"
	@echo "Frontend: http://localhost:3000"
	@echo "API Docs: http://localhost:8000/docs"
	@echo ""
	@echo "Press Ctrl+C to stop all servers"
	@trap 'kill %1 %2' INT; \
	cd backend && python -m api.server & \
	cd frontend && npm run dev & \
	wait

# Build production versions
build: backend-build frontend-build

# Run all tests
test: backend-test frontend-test

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf backend/__pycache__ backend/*/__pycache__ backend/*/*/__pycache__
	rm -rf backend/.pytest_cache backend/.coverage backend/htmlcov
	rm -rf frontend/dist frontend/node_modules/.vite
	rm -rf snapshots/*.json
	@echo "Clean complete!"

# Docker commands
docker-up:
	@echo "Starting Docker containers..."
	docker-compose up -d
	@echo "Services started!"
	@echo "Backend: http://localhost:8000"
	@echo "Frontend: http://localhost:3000"

docker-down:
	@echo "Stopping Docker containers..."
	docker-compose down
	@echo "Services stopped!"

# Backend commands
backend-install:
	@echo "Installing backend dependencies..."
	cd backend && pip install -r requirements.txt
	@echo "Backend dependencies installed!"

backend-dev:
	@echo "Starting backend server..."
	cd backend && python -m api.server

backend-test:
	@echo "Running backend tests..."
	cd backend && python -m pytest tests/ -v

backend-lint:
	@echo "Linting backend code..."
	cd backend && python -m black sim/ api/ tests/
	cd backend && python -m isort sim/ api/ tests/
	cd backend && python -m mypy sim/ api/

backend-build:
	@echo "Building backend..."
	cd backend && python -m pip install -r requirements.txt
	@echo "Backend build complete!"

# Frontend commands
frontend-install:
	@echo "Installing frontend dependencies..."
	cd frontend && npm install
	@echo "Frontend dependencies installed!"

frontend-dev:
	@echo "Starting frontend server..."
	cd frontend && npm run dev

frontend-test:
	@echo "Running frontend tests..."
	cd frontend && npm run test

frontend-lint:
	@echo "Linting frontend code..."
	cd frontend && npm run lint

frontend-build:
	@echo "Building frontend..."
	cd frontend && npm run build
	@echo "Frontend build complete!"

# Database commands (if needed)
db-migrate:
	@echo "Running database migrations..."
	# Add database migration commands here

db-seed:
	@echo "Seeding database..."
	# Add database seeding commands here

# Deployment commands
deploy-staging:
	@echo "Deploying to staging..."
	# Add staging deployment commands here

deploy-production:
	@echo "Deploying to production..."
	# Add production deployment commands here

# Monitoring commands
logs:
	@echo "Showing application logs..."
	docker-compose logs -f

logs-backend:
	@echo "Showing backend logs..."
	docker-compose logs -f backend

logs-frontend:
	@echo "Showing frontend logs..."
	docker-compose logs -f frontend

# Health checks
health:
	@echo "Checking service health..."
	@curl -f http://localhost:8000/ && echo "Backend: OK" || echo "Backend: FAIL"
	@curl -f http://localhost:3000/ && echo "Frontend: OK" || echo "Frontend: FAIL"

# Backup commands
backup-snapshots:
	@echo "Backing up snapshots..."
	tar -czf snapshots-backup-$(shell date +%Y%m%d-%H%M%S).tar.gz snapshots/
	@echo "Snapshots backed up!"

# Performance testing
perf-test:
	@echo "Running performance tests..."
	cd backend && python -m pytest tests/test_performance.py -v

# Security scanning
security-scan:
	@echo "Running security scan..."
	cd backend && python -m bandit -r sim/ api/
	cd frontend && npm audit

# Documentation
docs:
	@echo "Generating documentation..."
	cd backend && python -m sphinx-build -b html docs/ docs/_build/html
	@echo "Documentation generated in backend/docs/_build/html/"

# Release commands
release-patch:
	@echo "Creating patch release..."
	# Add release commands here

release-minor:
	@echo "Creating minor release..."
	# Add release commands here

release-major:
	@echo "Creating major release..."
	# Add release commands here
