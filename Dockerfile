FROM ubuntu:24.04

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.10 \
    python3-pip \
    python3-venv \
    git \
    build-essential \
    cmake \
    gcc \
    g++ \
    libstdc++6 \
    wget \
    curl \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /workspace

# Create Python virtual environment
RUN python3 -m venv /opt/spinnaker_env

# Activate virtual environment for subsequent commands
ENV PATH="/opt/spinnaker_env/bin:$PATH"

# Upgrade pip
RUN pip install --upgrade pip setuptools wheel

# Install common Python packages (add your specific requirements here)
RUN pip install numpy scipy matplotlib

# Copy your local files (uncomment and adjust as needed)
# COPY ./SpiNNMan2 /workspace/SpiNNMan2
# COPY ./py-spinnaker2 /workspace/py-spinnaker2

# Set up network access for SpiNNaker board
# You'll need to use --network host when running the container

# Set entrypoint to bash with venv activated
CMD ["/bin/bash"]

