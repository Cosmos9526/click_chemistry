# Use slim Python 3.12 as base image for efficiency
FROM python:3.12-slim

# Set working directory
WORKDIR /app

# Install system dependencies (for RDKit Draw and matplotlib)
RUN apt-get update && apt-get install -y \
    build-essential \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    libgomp1 \
    libexpat1 \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip first, then install Python dependencies
RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    rdkit \
    pandas \
    matplotlib \
    numpy \
    && pip cache purge

# Copy all Python files and JSON to /app
COPY *.py .
COPY merged_starred_linkers.json .

# Create outputs dir for mounted volumes
RUN mkdir -p /app/outputs

# Default command: run main.py (can override with args)
CMD ["python", "main.py"]