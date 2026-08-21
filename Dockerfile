# 1. Base Image: NVIDIA Official CUDA 11.7.1 Runtime on Ubuntu 22.04
FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

# Prevent interactive prompts during package installation and unbuffer Python logs
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/conda/bin:${PATH}"

WORKDIR /workspace/PHOCI

# 2. Install essential system dependencies and build libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    git \
    build-essential \
    libgl1-mesa-glx \
    libglib2.0-0 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 3. Install lightweight Miniconda to manage Python 3.8 environment
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_23.11.0-2-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    conda clean -afy

# 4. Install Python dependencies step-by-step

# Step 2: Upgrade basic build tools
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Step 3: Install PyTorch and CUDA 11.7 binaries
RUN pip install --no-cache-dir \
    torch==2.0.1 \
    torchvision==0.15.2 \
    torchaudio==2.0.2

# Step 4: Install core data science and computational packages
RUN pip install --no-cache-dir \
    pandas==1.5.3 \
    numpy==1.24.3 \
    scipy==1.10.1 \
    scikit-learn==1.3.2 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    h5py==3.11.0 \
    notebook==6.5.7

# Step 5: Install PyTorch Geometric (PyG) and pre-compiled C++ extensions
RUN pip install --no-cache-dir torch-geometric==2.5.2 && \
    pip install --no-cache-dir \
    torch-scatter==2.1.2+pt20cu117 \
    torch-sparse==0.6.18+pt20cu117 \
    torch-cluster==1.6.3+pt20cu117 \
    torch-spline-conv==1.2.2+pt20cu117 \
    -f https://data.pyg.org/whl/torch-2.0.0+cu117.html

# Step 6: Install 3D genomics and Hi-C data analysis tools
RUN pip install --no-cache-dir \
    cooler==0.9.3 \
    bioframe==0.3.3 \
    pybigwig==0.3.22 \
    pyfaidx==0.8.1.2 \
    trackc==0.0.18

# Copy local third-party wheel into the build stage for installation
COPY third_party_packages/hic_straw-1.3.1-cp38-cp38-linux_x86_64.whl /tmp/

# Install hicstraw directly from the local pre-compiled wheel
RUN pip install --no-cache-dir /tmp/hic_straw-1.3.1-cp38-cp38-linux_x86_64.whl && \
    rm /tmp/hic_straw-1.3.1-cp38-cp38-linux_x86_64.whl

# Step 7: Install utility and helper packages
RUN pip install --no-cache-dir \
    torchmetrics==1.4.1 \
    lightning-utilities==0.11.6 \
    aiohttp==3.10.3 \
    pyyaml==6.0.2 \
    tqdm==4.66.5 \
    efficient-apriori==2.0.5

# 5. Copy current project repository into the container workspace
COPY . /workspace/PHOCI

# Default entry command
CMD ["/bin/bash"]
