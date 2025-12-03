#  RAPTOR v2.1.0 Installation Guide

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

Complete installation instructions for all platforms and use cases.

---

##  Table of Contents

1. [System Requirements](#system-requirements)
2. [Quick Installation](#quick-installation)
3. [Detailed Installation](#detailed-installation)
4. [Pipeline Tools Installation](#pipeline-tools-installation)
5. [Platform-Specific Instructions](#platform-specific-instructions)
6. [Docker Installation](#docker-installation)
7. [Cloud Installation](#cloud-installation)
8. [Verification](#verification)
9. [Troubleshooting](#troubleshooting)
10. [Updating](#updating)

---

##  System Requirements

### Minimum Requirements

**Hardware:**
- CPU: 4 cores
- RAM: 8 GB
- Storage: 50 GB free space
- Network: Internet connection (for installation)

**Software:**
- Python: 3.8 or higher
- pip: 20.0 or higher
- OS: Linux, macOS, or Windows (WSL2)

### Recommended Requirements

**Hardware:**
- CPU: 16+ cores
- RAM: 32 GB or more
- Storage: 500 GB free space (SSD preferred)
- Network: High-speed internet

**Software:**
- Python: 3.9 or 3.10
- conda: 4.10 or higher (optional)
- Docker: 20.10 or higher (optional)

### Platform Support

| Platform | Support Level | Notes |
|----------|--------------|-------|
| Linux (Ubuntu 20.04+) | ✅ Full | Recommended |
| Linux (CentOS 7+) | ✅ Full | Tested |
| macOS (11+) | ✅ Full | Intel and Apple Silicon |
| Windows (WSL2) | ✅ Good | Ubuntu 20.04+ on WSL2 |
| Windows (Native) | ⚠️ Limited | Use WSL2 instead |
| HPC Clusters | ✅ Full | SLURM/PBS supported |
| Cloud (AWS/GCP/Azure) | ✅ Full | Native support |

---

##  Quick Installation

### Option 1: pip (Fastest)

```bash
# Install RAPTOR
pip install raptor-rnaseq

# Verify installation
raptor --version
```

**Time:** 5-10 minutes

### Option 2: conda (Recommended)

```bash
# Create environment
conda create -n raptor python=3.9
conda activate raptor

# Install RAPTOR
pip install raptor-rnaseq

# Verify installation
raptor --version
```

**Time:** 10-15 minutes

### Option 3: Docker (Easiest)

```bash
# Pull image
docker pull ayehblk/raptor:2.1.0

# Run RAPTOR
docker run -it ayehblk/raptor:2.1.0 raptor --version
```

**Time:** 15-20 minutes (depending on download speed)

---

##  Detailed Installation

### Step 1: Prepare Your System

#### Update System Packages

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install -y \
    build-essential \
    python3-dev \
    python3-pip \
    git \
    wget \
    curl \
    libgomp1
```

**CentOS/RHEL:**
```bash
sudo yum update
sudo yum groupinstall "Development Tools"
sudo yum install -y \
    python3-devel \
    python3-pip \
    git \
    wget \
    curl \
    libgomp
```

**macOS:**
```bash
# Install Homebrew if not installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install python@3.9 git wget
```

### Step 2: Install Python (if needed)

#### Check Python Version

```bash
python3 --version
# Should be 3.8 or higher
```

#### Install Python 3.9

**Using conda (Recommended):**
```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# Create Python environment
conda create -n raptor python=3.9
conda activate raptor
```

**Using apt (Ubuntu):**
```bash
sudo apt-get install python3.9 python3.9-venv python3.9-dev
```

### Step 3: Create Virtual Environment

#### Option A: venv (Built-in)

```bash
# Create environment
python3 -m venv raptor_env

# Activate environment
source raptor_env/bin/activate  # Linux/macOS
# OR
raptor_env\Scripts\activate.bat  # Windows

# Upgrade pip
pip install --upgrade pip
```

#### Option B: conda (Recommended)

```bash
# Create environment
conda create -n raptor python=3.9 -y

# Activate environment
conda activate raptor

# Install pip
conda install pip
```

### Step 4: Install RAPTOR

#### From PyPI (Stable Release)

```bash
pip install raptor-rnaseq
```

#### From GitHub (Latest Development)

```bash
pip install git+https://github.com/AyehBlk/RAPTOR.git@v2.1.0
```

#### From Source (For Development)

```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install in development mode
pip install -e .

# Install development dependencies (optional)
pip install -r requirements-dev.txt
```

### Step 5: Verify Installation

```bash
# Check version
raptor --version

# Run self-test
raptor test

# Check system info
raptor sysinfo
```

---

##  Pipeline Tools Installation

RAPTOR requires external RNA-seq tools. Install the ones you need:

### Option 1: conda (Easiest)

```bash
# Activate your environment
conda activate raptor

# Install all tools at once
conda install -c bioconda \
    salmon \
    kallisto \
    star \
    rsem \
    hisat2 \
    samtools \
    -y

# Verify installations
salmon --version
kallisto version
STAR --version
```

**Time:** 20-30 minutes

### Option 2: Individual Tool Installation

#### Salmon

```bash
# conda
conda install -c bioconda salmon

# OR from binary
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar xzf salmon-1.10.0_linux_x86_64.tar.gz
export PATH=$PATH:$(pwd)/salmon-latest_linux_x86_64/bin
```

#### Kallisto

```bash
# conda
conda install -c bioconda kallisto

# OR from binary
wget https://github.com/pachterlab/kallisto/releases/download/v0.50.0/kallisto_linux-v0.50.0.tar.gz
tar xzf kallisto_linux-v0.50.0.tar.gz
export PATH=$PATH:$(pwd)/kallisto
```

#### STAR

```bash
# conda
conda install -c bioconda star

# OR from source
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make STAR
export PATH=$PATH:$(pwd)
```

#### RSEM

```bash
# conda
conda install -c bioconda rsem
```

#### HISAT2

```bash
# conda
conda install -c bioconda hisat2
```

### Option 3: Containerized Tools

```bash
# Use Docker/Singularity for tools
# RAPTOR can work with containerized tools
raptor config --use-containers
```

---

##  Platform-Specific Instructions

### Ubuntu/Debian

```bash
# 1. Update system
sudo apt-get update && sudo apt-get upgrade

# 2. Install dependencies
sudo apt-get install -y python3-pip python3-venv build-essential

# 3. Create virtual environment
python3 -m venv raptor_env
source raptor_env/bin/activate

# 4. Install RAPTOR
pip install raptor-rnaseq

# 5. Install pipeline tools
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create -n raptor python=3.9
conda activate raptor
conda install -c bioconda salmon kallisto star rsem
```

### CentOS/RHEL

```bash
# 1. Enable EPEL
sudo yum install epel-release

# 2. Install dependencies
sudo yum groupinstall "Development Tools"
sudo yum install python39 python39-devel

# 3. Create virtual environment
python3.9 -m venv raptor_env
source raptor_env/bin/activate

# 4. Install RAPTOR
pip install raptor-rnaseq

# 5. Install tools via conda (recommended for CentOS)
# Follow conda installation above
```

### macOS

```bash
# 1. Install Xcode Command Line Tools
xcode-select --install

# 2. Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 3. Install Python
brew install python@3.9

# 4. Create virtual environment
python3 -m venv raptor_env
source raptor_env/bin/activate

# 5. Install RAPTOR
pip install raptor-rnaseq

# 6. Install conda for tools
brew install miniconda
conda init "$(basename "${SHELL}")"
# Restart terminal
conda create -n raptor python=3.9
conda activate raptor
conda install -c bioconda salmon kallisto star rsem
```

**Apple Silicon (M1/M2) Note:**
```bash
# Some tools may need Rosetta 2
softwareupdate --install-rosetta

# Use x86_64 conda environment
CONDA_SUBDIR=osx-64 conda create -n raptor python=3.9
conda activate raptor
conda config --env --set subdir osx-64
conda install -c bioconda salmon kallisto
```

### Windows (WSL2)

```bash
# 1. Install WSL2 (PowerShell as Administrator)
wsl --install -d Ubuntu-22.04

# 2. Launch Ubuntu and update
sudo apt-get update && sudo apt-get upgrade

# 3. Install Python
sudo apt-get install python3-pip python3-venv

# 4. Follow Ubuntu instructions above
```

**Windows Terminal Configuration:**
```json
{
  "profiles": {
    "list": [
      {
        "name": "Ubuntu RAPTOR",
        "commandline": "wsl.exe -d Ubuntu-22.04",
        "startingDirectory": "//wsl$/Ubuntu-22.04/home/YOUR_USERNAME"
      }
    ]
  }
}
```

### HPC Clusters

```bash
# 1. Load required modules
module load python/3.9
module load gcc/10.2
module load git

# 2. Install in user space
pip install --user raptor-rnaseq

# 3. Install conda in user space
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate

# 4. Install tools
conda create -n raptor python=3.9
conda activate raptor
conda install -c bioconda salmon kallisto star rsem

# 5. Add to .bashrc
echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
echo 'source $HOME/miniconda3/bin/activate' >> ~/.bashrc
echo 'conda activate raptor' >> ~/.bashrc
```

**SLURM Job Script:**
```bash
#!/bin/bash
#SBATCH --job-name=raptor
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load python/3.9
source $HOME/miniconda3/bin/activate
conda activate raptor

raptor analyze --config config.yaml
```

---

##  Docker Installation

### Basic Docker Setup

```bash
# 1. Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# 2. Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# 3. Pull RAPTOR image
docker pull ayehblk/raptor:2.1.0

# 4. Test installation
docker run --rm ayehblk/raptor:2.1.0 raptor --version
```

### Using Docker

```bash
# Run interactive session
docker run -it \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  ayehblk/raptor:2.1.0 bash

# Run analysis directly
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/config.yaml:/config.yaml \
  ayehblk/raptor:2.1.0 \
  raptor analyze --config /config.yaml

# Start dashboard
docker run -p 8501:8501 \
  -v $(pwd)/results:/results \
  ayehblk/raptor:2.1.0 \
  raptor dashboard --host 0.0.0.0
```

### Building Custom Image

```dockerfile
# Dockerfile
FROM ayehblk/raptor:2.1.0

# Add custom tools
RUN conda install -c bioconda additional-tool

# Add custom scripts
COPY my_script.py /usr/local/bin/

# Set working directory
WORKDIR /workspace
```

```bash
# Build image
docker build -t my-raptor:latest .

# Run custom image
docker run -it my-raptor:latest bash
```

---

##  Cloud Installation

### AWS

```bash
# Launch EC2 instance (Amazon Linux 2)
# t3.2xlarge or larger recommended

# 1. Update system
sudo yum update -y

# 2. Install Python 3.9
sudo yum install python39 -y

# 3. Install RAPTOR
pip3.9 install raptor-rnaseq

# 4. Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/miniconda3/bin/activate

# 5. Install tools
conda create -n raptor python=3.9
conda activate raptor
conda install -c bioconda salmon kallisto star rsem

# 6. Configure AWS credentials
aws configure
```

### Google Cloud Platform

```bash
# Launch Compute Engine instance (Ubuntu 22.04)
# n2-standard-8 or larger recommended

# 1. Update system
sudo apt-get update && sudo apt-get upgrade -y

# 2. Install RAPTOR
sudo apt-get install python3-pip -y
pip3 install raptor-rnaseq

# 3. Install conda and tools
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/miniconda3/bin/activate
conda create -n raptor python=3.9
conda activate raptor
conda install -c bioconda salmon kallisto star rsem

# 4. Configure GCP credentials
gcloud auth login
```

### Microsoft Azure

```bash
# Launch VM (Ubuntu 22.04)
# Standard_D8s_v3 or larger recommended

# Follow Ubuntu installation instructions above

# Configure Azure credentials
az login
```

---

##  Verification

### Quick Check

```bash
# 1. Check RAPTOR version
raptor --version
# Expected: RAPTOR v2.1.0

# 2. Check Python packages
pip list | grep raptor

# 3. Check pipeline tools
which salmon
which kallisto
which STAR

# 4. Run self-test
raptor test
```

### Comprehensive Verification

```bash
# 1. System information
raptor sysinfo

# 2. Test each pipeline
raptor test --pipeline salmon
raptor test --pipeline kallisto
raptor test --pipeline star

# 3. Test ML features
raptor ml test

# 4. Test dashboard
raptor dashboard --test-mode

# 5. Run example analysis
raptor analyze --demo
```

### Expected Output

```
✓ RAPTOR v2.1.0 installed successfully
✓ Python 3.9.7 found
✓ All required packages installed
✓ Salmon v1.10.0 detected
✓ Kallisto v0.50.0 detected
✓ STAR v2.7.11a detected
✓ ML model loaded successfully
✓ Dashboard components ready
✓ System tests passed (23/23)

Installation verified!   
```

---

##  Troubleshooting

### Common Issues

#### Issue: pip install fails

**Error:**
```
ERROR: Could not find a version that satisfies the requirement raptor-rnaseq
```

**Solution:**
```bash
# Update pip
pip install --upgrade pip

# Try with --user flag
pip install --user raptor-rnaseq

# Check Python version
python --version  # Must be 3.8+
```

#### Issue: Import errors

**Error:**
```python
ImportError: No module named 'raptor'
```

**Solution:**
```bash
# Verify installation
pip show raptor-rnaseq

# Check Python path
python -c "import sys; print(sys.path)"

# Reinstall
pip uninstall raptor-rnaseq
pip install raptor-rnaseq
```

#### Issue: Tool not found

**Error:**
```
FileNotFoundError: salmon not found in PATH
```

**Solution:**
```bash
# Check if tool is installed
which salmon

# Add to PATH
export PATH=$PATH:/path/to/salmon/bin

# Or reinstall with conda
conda install -c bioconda salmon
```

#### Issue: Permission denied

**Error:**
```
PermissionError: [Errno 13] Permission denied
```

**Solution:**
```bash
# Install in user space
pip install --user raptor-rnaseq

# Or use virtual environment
python3 -m venv raptor_env
source raptor_env/bin/activate
pip install raptor-rnaseq
```

#### Issue: Memory error during installation

**Error:**
```
MemoryError: Unable to allocate array
```

**Solution:**
```bash
# Increase swap space (Linux)
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Install with --no-cache-dir
pip install --no-cache-dir raptor-rnaseq
```

### Getting Help

**If issues persist:**

1. Check [FAQ](FAQ.md)
2. Search [GitHub issues](https://github.com/AyehBlk/RAPTOR/issues)
3. Create new issue with:
   - Error message
   - System information (`raptor sysinfo`)
   - Installation method
   - Python version

---

##  Updating

### Update RAPTOR

```bash
# Update from PyPI
pip install --upgrade raptor-rnaseq

# Update from GitHub
pip install --upgrade git+https://github.com/AyehBlk/RAPTOR.git

# Check new version
raptor --version
```

### Update Pipeline Tools

```bash
# Update all conda packages
conda update --all

# Update specific tool
conda update salmon

# Check versions
salmon --version
kallisto version
STAR --version
```

### Migrate Configuration

```bash
# Migrate from v2.0.0 to v2.1.0
raptor migrate --from v2.0.0 --to v2.1.0

# Backup old config
cp config.yaml config.yaml.bak

# Validate new config
raptor config --validate
```

---

##  Next Steps

After installation:

1. **Read the User Guide**: [USER_GUIDE.md](USER_GUIDE.md)
2. **Configure RAPTOR**: [CONFIGURATION.md](CONFIGURATION.md)
3. **Try tutorials**: [examples/](examples/)
4. **Join community**: [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)

---

##  Support

**Need help with installation?**

- **Documentation**: [GitHub Wiki](https://github.com/AyehBlk/RAPTOR/wiki)
- **Issues**: [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- **Email**: ayehbolouki1988@gmail.com

---

**Author:** Ayeh Bolouki   
**Version:** 2.1.0  
**License:** MIT

---

*"Installation is just the beginning of the journey!"* 
