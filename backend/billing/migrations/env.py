"""
Alembic environment configuration for billing module.
"""

from logging.config import fileConfig
from sqlalchemy import create_engine, engine_from_config
from sqlalchemy import pool
from alembic import context
import os
import sys
import importlib.util

# Add parent directory to path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Load .env file first to ensure DATABASE_URL is available
from dotenv import load_dotenv
env_path = os.path.join(project_root, '.env')
if os.path.exists(env_path):
    load_dotenv(env_path)

# Import DATABASE_URL from config directly (avoid __init__.py)
# But ensure we use the value from environment (which may have been loaded from .env)
config_path = os.path.join(project_root, 'backend', 'billing', 'config.py')
config_spec = importlib.util.spec_from_file_location("billing_config", config_path)
config_module = importlib.util.module_from_spec(config_spec)
config_spec.loader.exec_module(config_module)
DATABASE_URL = config_module.DATABASE_URL

# If DATABASE_URL is still default, try to get it directly from environment
if DATABASE_URL == "postgresql://live2:password@localhost:5432/live2_billing":
    DATABASE_URL = os.getenv("DATABASE_URL", DATABASE_URL)

# Import models directly (avoid circular imports from __init__.py)
# Import models.py directly without going through __init__.py
models_path = os.path.join(project_root, 'backend', 'billing', 'models.py')
models_spec = importlib.util.spec_from_file_location("billing_models", models_path)
models_module = importlib.util.module_from_spec(models_spec)
models_spec.loader.exec_module(models_module)
Base = models_module.Base

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
# This line sets up loggers basically.
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Set database URL from config
# Ensure proper encoding for connection string
if isinstance(DATABASE_URL, bytes):
    DATABASE_URL = DATABASE_URL.decode('utf-8')
# Don't set in config - we'll use create_engine directly to avoid encoding issues
# This prevents UnicodeDecodeError when Alembic tries to parse the URL

# add your model's MetaData object here
# for 'autogenerate' support
target_metadata = Base.metadata

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    # Use DATABASE_URL directly instead of config.get_main_option to avoid encoding issues
    from urllib.parse import urlparse, unquote
    parsed = urlparse(DATABASE_URL)
    password = unquote(parsed.password) if parsed.password else None
    # Reconstruct URL without password for offline mode (password not needed)
    url = f"{parsed.scheme}://{parsed.username}@{parsed.hostname}:{parsed.port}{parsed.path}"
    
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """
    # Create engine directly from DATABASE_URL to avoid encoding issues
    # Parse URL and use connect_args to avoid UnicodeDecodeError
    from urllib.parse import urlparse, unquote
    import psycopg2
    
    # Parse connection string
    parsed = urlparse(DATABASE_URL)
    # Decode password if URL-encoded
    password = unquote(parsed.password) if parsed.password else None
    
    # Create connection creator function that uses psycopg2.connect directly
    # This avoids SQLAlchemy URL parsing which can cause UnicodeDecodeError
    # Use DSN string instead of keyword arguments to avoid encoding issues
    def create_connection():
        # Build DSN string manually to avoid encoding issues
        dsn_parts = []
        dsn_parts.append(f"host={parsed.hostname}")
        dsn_parts.append(f"port={parsed.port}")
        dsn_parts.append(f"user={parsed.username}")
        if password:
            # Password is already decoded, use as-is
            dsn_parts.append(f"password={password}")
        db_name = parsed.path[1:] if parsed.path.startswith('/') else parsed.path
        dsn_parts.append(f"dbname={db_name}")
        dsn_parts.append("connect_timeout=10")
        dsn = " ".join(dsn_parts)
        return psycopg2.connect(dsn)
    
    # Use create_engine with connection creator to bypass URL parsing
    connectable = create_engine(
        "postgresql://",  # Dummy URL, we'll use connection creator
        poolclass=pool.NullPool,
        echo=False,
        creator=create_connection
    )

    with connectable.connect() as connection:
        context.configure(
            connection=connection, target_metadata=target_metadata
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
