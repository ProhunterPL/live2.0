"""initial

Revision ID: 001_initial
Revises: 
Create Date: 2025-12-04 12:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = '001_initial'
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Create users table
    op.create_table(
        'users',
        sa.Column('id', postgresql.UUID(as_uuid=True), primary_key=True),
        sa.Column('email', sa.String(255), nullable=False, unique=True),
        sa.Column('password_hash', sa.String(255), nullable=False),
        sa.Column('api_key', sa.String(64), nullable=False, unique=True),
        sa.Column('tier', sa.String(20), nullable=False, server_default='hobby'),
        sa.Column('subscription_status', sa.String(20), nullable=False, server_default='trial'),
        sa.Column('stripe_customer_id', sa.String(255), nullable=True),
        sa.Column('stripe_subscription_id', sa.String(255), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('updated_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
    )
    op.create_index('ix_users_email', 'users', ['email'])
    op.create_index('ix_users_api_key', 'users', ['api_key'])
    op.create_index('ix_users_stripe_customer_id', 'users', ['stripe_customer_id'])
    op.create_index('ix_users_stripe_subscription_id', 'users', ['stripe_subscription_id'])

    # Create subscriptions table
    op.create_table(
        'subscriptions',
        sa.Column('id', postgresql.UUID(as_uuid=True), primary_key=True),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('tier', sa.String(20), nullable=False),
        sa.Column('status', sa.String(20), nullable=False),
        sa.Column('current_period_start', sa.Date(), nullable=False),
        sa.Column('current_period_end', sa.Date(), nullable=False),
        sa.Column('cancel_at_period_end', sa.Boolean(), nullable=False, server_default='false'),
        sa.Column('stripe_subscription_id', sa.String(255), nullable=True, unique=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('updated_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ondelete='CASCADE'),
    )
    op.create_index('ix_subscriptions_user_id', 'subscriptions', ['user_id'])
    op.create_index('ix_subscriptions_stripe_subscription_id', 'subscriptions', ['stripe_subscription_id'])

    # Create usage table
    op.create_table(
        'usage',
        sa.Column('id', postgresql.UUID(as_uuid=True), primary_key=True),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('date', sa.Date(), nullable=False),
        sa.Column('reactions_count', sa.Integer(), nullable=False, server_default='0'),
        sa.Column('api_calls_count', sa.Integer(), nullable=False, server_default='0'),
        sa.Column('tier', sa.String(20), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('updated_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ondelete='CASCADE'),
        sa.UniqueConstraint('user_id', 'date', name='uq_user_date'),
    )
    op.create_index('ix_usage_user_id', 'usage', ['user_id'])
    op.create_index('ix_usage_date', 'usage', ['date'])


def downgrade() -> None:
    op.drop_table('usage')
    op.drop_table('subscriptions')
    op.drop_table('users')

