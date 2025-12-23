#!/bin/bash
# Setup script dla DigitalOcean Droplet
# Uruchom jako root: bash setup_do_droplet.sh

set -e

echo "🚀 Live 2.0 - DigitalOcean Setup"
echo "=================================="

# Update system
echo "📦 Updating system..."
apt update && apt upgrade -y

# Install dependencies
echo "📦 Installing dependencies..."
apt install -y \
    python3.11 \
    python3.11-venv \
    python3-pip \
    nginx \
    git \
    certbot \
    python3-certbot-nginx \
    curl \
    wget \
    build-essential

# Create user (if not exists)
if ! id "live2" &>/dev/null; then
    echo "👤 Creating user 'live2'..."
    adduser --disabled-password --gecos "" live2
    usermod -aG sudo live2
fi

# Setup project directory
echo "📁 Setting up project directory..."
mkdir -p /opt
cd /opt

if [ ! -d "live2.0" ]; then
    echo "⚠️  Repository not found. Please clone it manually:"
    echo "   cd /opt && git clone https://github.com/YOUR_REPO/live2.0.git"
else
    echo "✅ Repository found at /opt/live2.0"
fi

# Setup Python venv
echo "🐍 Setting up Python environment..."
cd /opt/live2.0
python3.11 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# Create .env if not exists
if [ ! -f ".env" ]; then
    echo "📝 Creating .env template..."
    cat > .env <<EOF
# DATABASE (Supabase)
DATABASE_URL=

# REDIS
REDIS_HOST=
REDIS_PORT=6379
REDIS_USERNAME=
REDIS_PASSWORD=

# JWT
JWT_SECRET_KEY=

# STRIPE
STRIPE_SECRET_KEY=
STRIPE_PUBLISHABLE_KEY=
STRIPE_WEBHOOK_SECRET=
STRIPE_PRICE_ID_HOBBY=
STRIPE_PRICE_ID_RESEARCH=
STRIPE_PRICE_ID_PRO=

# AWS
AWS_ACCESS_KEY_ID=
AWS_SECRET_ACCESS_KEY=
AWS_REGION=us-east-1
AWS_BATCH_JOB_QUEUE=
AWS_BATCH_JOB_DEFINITION=

# SUPABASE
SUPABASE_URL=
SUPABASE_SERVICE_KEY=

# APP
ENV=prod
API_BASE_URL=
EOF
    echo "⚠️  Please edit /opt/live2.0/.env with your production values"
    chmod 600 .env
fi

# Create systemd service
echo "⚙️  Creating systemd service..."
cat > /etc/systemd/system/live2-backend.service <<EOF
[Unit]
Description=Live 2.0 Backend API
After=network.target

[Service]
Type=simple
User=live2
WorkingDirectory=/opt/live2.0
Environment="PATH=/opt/live2.0/venv/bin"
ExecStart=/opt/live2.0/venv/bin/gunicorn \
  -w 4 \
  -k uvicorn.workers.UvicornWorker \
  -b 127.0.0.1:8000 \
  backend.api.server:app
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

# Setup Nginx (basic config)
echo "🌐 Setting up Nginx..."
cat > /etc/nginx/sites-available/live2 <<EOF
server {
    listen 80;
    server_name _;

    location /api {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }

    location /status {
        proxy_pass http://127.0.0.1:8000;
    }
}
EOF

# Enable site
ln -sf /etc/nginx/sites-available/live2 /etc/nginx/sites-enabled/
rm -f /etc/nginx/sites-enabled/default
nginx -t

echo ""
echo "✅ Setup complete!"
echo ""
echo "📋 Next steps:"
echo "1. Edit /opt/live2.0/.env with your production values"
echo "2. Run migrations: cd /opt/live2.0 && source venv/bin/activate && alembic -c backend/billing/migrations/alembic.ini upgrade head"
echo "3. Start service: systemctl enable live2-backend && systemctl start live2-backend"
echo "4. Setup SSL: certbot --nginx -d your-domain.com"
echo ""

