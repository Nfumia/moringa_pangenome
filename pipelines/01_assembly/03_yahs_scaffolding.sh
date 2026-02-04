#!/bin/bash

# YaHS Hi-C Scaffolding Pipeline
# Usage: bash 03_yahs_scaffolding.sh <SAMPLE_NAME>

# set -e  # Disabled to prevent exit on minor errors

SAMPLE=$1

if [ -z "$SAMPLE" ]; then
    echo "Usage: bash 03_yahs_scaffolding.sh <SAMPLE_NAME>"
    echo "Example: bash 03_yahs_scaffolding.sh Mo-TH-30"
    exit 1
fi

# Directories
BASE_DIR="$HOME/moringa_pangenome"
ASSEMBLY_DIR="$BASE_DIR/assemblies/production/$SAMPLE"
HIC_DIR="$BASE_DIR/raw_data/hic/$SAMPLE"
SCAFFOLD_DIR="$BASE_DIR/assemblies/scaffolded/$SAMPLE"
LOG_DIR="$BASE_DIR/logs"

# Create output directory
mkdir -p "$SCAFFOLD_DIR"
mkdir -p "$LOG_DIR"

# Files
ASSEMBLY="$ASSEMBLY_DIR/${SAMPLE}.primary.fasta"
HIC_R1="$HIC_DIR/${SAMPLE}_R1.fq.gz"
HIC_R2="$HIC_DIR/${SAMPLE}_R2.fq.gz"
LOG_FILE="$LOG_DIR/yahs_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"

echo "========================================" | tee "$LOG_FILE"
echo "YaHS Hi-C Scaffolding" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Sample: $SAMPLE" | tee -a "$LOG_FILE"
echo "Assembly: $ASSEMBLY" | tee -a "$LOG_FILE"
echo "Hi-C R1: $HIC_R1" | tee -a "$LOG_FILE"
echo "Hi-C R2: $HIC_R2" | tee -a "$LOG_FILE"
echo "Output: $SCAFFOLD_DIR" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Check files exist
if [ ! -f "$ASSEMBLY" ]; then
    echo "ERROR: Assembly not found: $ASSEMBLY" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$HIC_R1" ] || [ ! -f "$HIC_R2" ]; then
    echo "ERROR: Hi-C reads not found" | tee -a "$LOG_FILE"
    echo "Expected: $HIC_R1" | tee -a "$LOG_FILE"
    echo "Expected: $HIC_R2" | tee -a "$LOG_FILE"
    exit 1
fi

cd "$SCAFFOLD_DIR"

# Step 1: Index assembly with BWA
echo "Step 1: Indexing assembly with BWA..." | tee -a "$LOG_FILE"
bwa index "$ASSEMBLY" 2>&1 | tee -a "$LOG_FILE"

# Step 1b: Create FASTA index for YaHS (CRITICAL - YaHS requires .fai file)
echo "" | tee -a "$LOG_FILE"
echo "Step 1b: Creating FASTA index for YaHS..." | tee -a "$LOG_FILE"
samtools faidx "$ASSEMBLY" 2>&1 | tee -a "$LOG_FILE"

# Verify .fai file was created
if [ ! -f "${ASSEMBLY}.fai" ]; then
    echo "ERROR: Failed to create FASTA index (.fai file)" | tee -a "$LOG_FILE"
    exit 1
fi
echo "✓ FASTA index created: ${ASSEMBLY}.fai" | tee -a "$LOG_FILE"

# Step 2: Align Hi-C reads to assembly
echo "" | tee -a "$LOG_FILE"
echo "Step 2: Aligning Hi-C reads..." | tee -a "$LOG_FILE"
bwa mem -5SP -t 32 "$ASSEMBLY" "$HIC_R1" "$HIC_R2" | \
    samtools view -@ 8 -bS - | \
    samtools sort -@ 8 -o ${SAMPLE}_hic.bam - 2>&1 | tee -a "$LOG_FILE"

# Index BAM
echo "Indexing BAM file..." | tee -a "$LOG_FILE"
samtools index -@ 8 ${SAMPLE}_hic.bam 2>&1 | tee -a "$LOG_FILE"

# Step 3: Run YaHS scaffolding
echo "" | tee -a "$LOG_FILE"
echo "Step 3: Running YaHS scaffolding..." | tee -a "$LOG_FILE"
yahs "$ASSEMBLY" ${SAMPLE}_hic.bam -o ${SAMPLE}_scaffolds 2>&1 | tee -a "$LOG_FILE"

# Step 4: Generate statistics
echo "" | tee -a "$LOG_FILE"
echo "Step 4: Assembly statistics..." | tee -a "$LOG_FILE"
echo "=== BEFORE SCAFFOLDING ===" | tee -a "$LOG_FILE"
seqkit stats -a "$ASSEMBLY" 2>&1 | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "=== AFTER SCAFFOLDING ===" | tee -a "$LOG_FILE"
if [ -f "${SAMPLE}_scaffolds_scaffolds_final.fa" ]; then
    seqkit stats -a ${SAMPLE}_scaffolds_scaffolds_final.fa 2>&1 | tee -a "$LOG_FILE"
else
    echo "WARNING: Scaffolded assembly not found!" | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Scaffolding complete!" | tee -a "$LOG_FILE"
echo "Output directory: $SCAFFOLD_DIR" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
