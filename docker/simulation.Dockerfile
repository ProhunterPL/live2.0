FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install AWS CLI (for S3 uploads)
RUN pip install --no-cache-dir boto3 awscli

# Install Supabase client (for status updates)
RUN pip install --no-cache-dir supabase

# Copy backend
COPY backend/ ./backend/

# Set Python path
ENV PYTHONPATH=/app

# Default command (może być nadpisane przez Batch)
CMD ["python", "-m", "backend.sim.run_simulation"]

