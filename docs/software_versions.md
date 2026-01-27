# Software Versions

Complete list of software used in the Moringa pangenome project.

---

## Assembly Pipeline

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **hifiasm** | v0.19.8 | HiFi genome assembly & haplotype phasing | `conda install -c bioconda hifiasm=0.19.8` |
| **YaHS** | v1.2a.2 | Hi-C scaffolding | `conda install -c bioconda yahs=1.2a.2` |
| **bwa-mem2** | v2.2.1 | Hi-C read alignment | `conda install -c bioconda bwa-mem2=2.2.1` |
| **samtools** | v1.17 | BAM file processing | `conda install -c bioconda samtools=1.17` |

---

## Quality Control & Statistics

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **seqkit** | v2.3.1 | Sequence statistics and manipulation | `conda install -c bioconda seqkit=2.3.1` |
| **BUSCO** | v5.4.7 | Gene completeness assessment | `conda install -c bioconda busco=5.4.7` |
| **Quast** | v5.2.0 | Assembly quality evaluation | `conda install -c bioconda quast=5.2.0` |
| **GenomeScope2** | v2.0 | Heterozygosity estimation | `conda install -c bioconda genomescope2=2.0` |

---

## Haplotype Validation

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **BLAST+** | v2.14.0 | Homology analysis & contamination screening | `conda install -c bioconda blast=2.14.0` |
| **minimap2** | v2.24 | Synteny alignment | `conda install -c bioconda minimap2=2.24` |
| **meryl** | v1.3 | K-mer counting | `conda install -c bioconda meryl=1.3` |
| **merqury** | v1.3 | K-mer based assembly QC | `conda install -c bioconda merqury=1.3` |

---

## Annotation (Planned)

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **BRAKER3** | v3.0.3 | Gene prediction pipeline | `conda install -c bioconda braker3=3.0.3` |
| **AUGUSTUS** | v3.5.0 | Ab initio gene prediction | `conda install -c bioconda augustus=3.5.0` |
| **GeneMark-ES** | v4.69 | Gene modeling | License required |

---

## Pangenome Construction (Planned)

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **Minigraph-Cactus** | v2.6.0 | Pangenome graph construction | `conda install -c bioconda cactus=2.6.0` |

---

## Utilities

| Software | Version | Purpose | Installation |
|----------|---------|---------|--------------|
| **Python** | v3.10 | Scripting and analysis | `conda install python=3.10` |
| **pandas** | v2.0.0 | Data manipulation | `conda install pandas=2.0.0` |
| **numpy** | v1.24.0 | Numerical computing | `conda install numpy=1.24.0` |
| **matplotlib** | v3.7.0 | Plotting | `conda install matplotlib=3.7.0` |
| **seaborn** | v0.12.0 | Statistical visualization | `conda install seaborn=0.12.0` |

---

## Complete Conda Environment

### Create environment from scratch:
```bash
# Create new environment
conda create -n moringa_pangenome python=3.10

# Activate environment
conda activate moringa_pangenome

# Install assembly tools
conda install -c bioconda -c conda-forge \
    hifiasm=0.19.8 \
    yahs=1.2a.2 \
    bwa-mem2=2.2.1 \
    samtools=1.17

# Install QC tools
conda install -c bioconda -c conda-forge \
    seqkit=2.3.1 \
    busco=5.4.7 \
    quast=5.2.0 \
    genomescope2=2.0

# Install analysis tools
conda install -c bioconda -c conda-forge \
    blast=2.14.0 \
    minimap2=2.24 \
    meryl=1.3 \
    merqury=1.3

# Install Python packages
conda install -c conda-forge \
    pandas=2.0.0 \
    numpy=1.24.0 \
    matplotlib=3.7.0 \
    seaborn=0.12.0
```

### Or use environment.yml:
```bash
conda env create -f environment.yml
conda activate moringa_pangenome
```

---

## System Requirements

**Minimum:**
- 64 GB RAM (for assembly)
- 24 CPU cores
- 2 TB storage

**Recommended:**
- 128 GB RAM
- 48+ CPU cores
- 5 TB storage (for raw data + assemblies)

**Operating System:**
- Linux (Ubuntu 22.04 or later recommended)
- MacOS (with modifications)

---

## Verification

Check installed versions:
```bash
# Assembly tools
hifiasm --version
yahs --version
bwa-mem2 version
samtools --version

# QC tools
seqkit version
busco --version
quast --version

# Analysis tools
blastn -version
minimap2 --version
meryl --version

# Python
python --version
conda list | grep -E "pandas|numpy|matplotlib|seaborn"
```

---

## Citation Information

If using this software environment, please cite:

**Hifiasm:**
> Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*, 18, 170-175.

**YaHS:**
> Zhou, C., McCarthy, S.A., Durbin, R. (2023) YaHS: yet another Hi-C scaffolding tool. *Bioinformatics*, 39(1), btac808.

**BUSCO:**
> Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A., Zdobnov, E.M. (2021) BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and Evolution*, 38(10), 4647-4654.

**Merqury:**
> Rhie, A., Walenz, B.P., Koren, S., Phillippy, A.M. (2020) Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biology*, 21, 245.

---

**Last Updated:** January 22, 2026  
**Environment Name:** `moringa_pangenome`  
**Conda Version:** 24.1.2
