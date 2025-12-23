#!/usr/bin/env python3
"""
Test database connection for billing module (without importing full billing module).
"""

import sys
import os
from dotenv import load_dotenv

load_dotenv()

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

def test_db_connection():
    """Test database connection directly."""
    print("=" * 60)
    print("Database Connection Test (Billing)")
    print("=" * 60)
    
    try:
        from sqlalchemy import create_engine, text
        from sqlalchemy.pool import NullPool
        
        # Get DATABASE_URL from env
        database_url = os.getenv("DATABASE_URL")
        if not database_url:
            print("\nFAILED: DATABASE_URL not set in .env")
            return False
        
        print(f"\n[1/3] Creating engine...")
        engine = create_engine(
            database_url,
            poolclass=NullPool,
            echo=False
        )
        
        print(f"[2/3] Testing connection...")
        with engine.connect() as conn:
            result = conn.execute(text("SELECT 1"))
            row = result.fetchone()
            if row and row[0] == 1:
                print("  OK - Database is responding")
            else:
                print("  FAILED - Unexpected response")
                return False
        
        print(f"[3/3] Checking if tables exist...")
        with engine.connect() as conn:
            # Check if users table exists
            result = conn.execute(text("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' 
                    AND table_name = 'users'
                );
            """))
            exists = result.fetchone()[0]
            if exists:
                print("  OK - Tables exist (users table found)")
                
                # Count users
                result = conn.execute(text("SELECT COUNT(*) FROM users"))
                count = result.fetchone()[0]
                print(f"  Users in DB: {count}")
            else:
                print("  WARNING - Tables don't exist (run migrations!)")
                print("\n  To create tables, run:")
                print("    alembic -c backend/billing/migrations/alembic.ini upgrade head")
        
        print("\n" + "=" * 60)
        print("SUCCESS: Database connection working!")
        print("=" * 60)
        return True
        
    except Exception as e:
        print(f"\nFAILED: Database connection error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_db_connection()
    sys.exit(0 if success else 1)

