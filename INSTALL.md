# Installation Instructions

## Quick Start

### Option 1: Full Environment (Recommended)
```bash
# Clone repository
git clone https://github.com/Nfumia/moringa_pangenome.git
cd moringa_pangenome

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate moringa_pangenome

# Verify installation
bash scripts/verify_installation.sh
```

**Size:** ~5 GB  
**Install time:** 10-20 minutes

---

### Option 2: Minimal Environment

For systems with limited space or if you only need core assembly tools:
```bash
# Create minimal environment
conda env create -f environment.minimal.yml

# Activate environment
conda activate moringa_pangenome_minimal
```

**Size:** ~2 GB  
**Install time:** 5-10 minutes

---

## Manual Installation

If you prefer to install tools individually:
```bash
# Create new environment
conda create -n moringa_pangenome python=3.10
conda activate moringa_pangenome

# Install assembly tools
conda install -c bioconda hifiasm=0.19.8
conda install -c bioconda yahs=1.2a.2
conda install -c bioconda bwa-mem2=2.2.1
conda install -c bioconda samtools=1.17

# Install QC tools
conda install -c bioconda seqkit=2.3.1
conda install -c bioconda busco=5.4.7
conda install -c bioconda quast=5.2.0

# Install analysis tools
conda install -c bioconda blast=2.14.0
conda install -c bioconda minimap2=2.24
conda install -c bioconda meryl=1.3
conda install -c bioconda merqury=1.3

# Install Python packages
conda install -c conda-forge pandas numpy matplotlib seaborn
```

---

## Verification

After installation, verify all tools are working:
```bash
# Check versions
hifiasm --version
yahs --version
samtools --version
seqkit version
blastn -version
minimap2 --version

# Should all return version information without errors
```

---

## Troubleshooting

### Common Issues

**1. Conda is slow**
```bash
# Use mamba instead (much faster)
conda install mamba -n base -c conda-forge
mamba env create -f environment.yml
```

**2. Channel conflicts**
```bash
# Set channel priority
conda config --set channel_priority strict
```

**3. Out of disk space**
```bash
# Clean conda cache
conda clean --all

# Use minimal environment instead
conda env create -f environment.minimal.yml
```

**4. Specific package fails to install**
```bash
# Try installing without version constraints
conda install -c bioconda package_name

# Or use pip for Python packages
pip install package_name
```

---

## System Requirements

**Minimum:**
- 64 GB RAM
- 24 CPU cores  
- 2 TB storage
- Linux OS (Ubuntu 22.04+)

**Recommended:**
- 128 GB RAM
- 48+ CPU cores
- 5 TB storage
- Fast SSD for temporary files

---

## Updating

To update the environment:
```bash
# Update all packages
conda env update -f environment.yml --prune

# Or recreate from scratch
conda env remove -n moringa_pangenome
conda env create -f environment.yml
```

---

## Uninstallation
```bash
# Remove conda environment
conda env remove -n moringa_pangenome

# Remove cloned repository
rm -rf moringa_pangenome/
```

---

## Alternative: Docker (Coming Soon)

For reproducibility, we will provide a Docker image:
```bash
# Pull Docker image (future)
docker pull username/moringa_pangenome:latest

# Run container
docker run -it -v $(pwd):/data username/moringa_pangenome:latest
```

---

**Last Updated:** January 22, 2026
