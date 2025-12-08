# Billing Module Migrations

Alembic migrations for the billing module.

## Setup

1. Ensure PostgreSQL is running (via Docker Compose or local installation)
2. Set `DATABASE_URL` environment variable:
   ```bash
   export DATABASE_URL="postgresql://live2:password@localhost:5432/live2_billing"
   ```

## Running Migrations

### From project root:

```bash
# Run migrations
cd backend/billing
python -m alembic upgrade head

# Create new migration
python -m alembic revision --autogenerate -m "description"

# Rollback one migration
python -m alembic downgrade -1
```

### Using Alembic directly:

```bash
# From backend/billing directory
alembic upgrade head
alembic revision --autogenerate -m "description"
alembic downgrade -1
```

## Migration Files

- `001_initial.py` - Initial schema (users, subscriptions, usage tables)

## Notes

- All migrations use PostgreSQL UUID type for primary keys
- Foreign keys have CASCADE delete
- Indexes are created for frequently queried columns

