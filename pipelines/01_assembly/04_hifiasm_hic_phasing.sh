#!/bin/bash

# Hifiasm with Hi-C Integration for Haplotype Phasing
# Usage: bash 04_hifiasm_hic_phasing.sh <SAMPLE_NAME>

# set -e  # Disabled to prevent exit on errors

SAMPLE=$1

if [ -z "$SAMPLE" ]; then
    echo "Usage: bash 04_hifiasm_hic_phasing.sh <SAMPLE_NAME>"
    echo "Example: bash 04_hifiasm_hic_phasing.sh Mo-TH-55"
    exit 1
fi

# Directories
BASE_DIR="$HOME/moringa_pangenome"
HIFI_DIR="$BASE_DIR/raw_data/hifi"
HIC_DIR="$BASE_DIR/raw_data/hic/$SAMPLE"
OUTPUT_DIR="$BASE_DIR/assemblies/phased/$SAMPLE"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Files
HIFI="$HIFI_DIR/${SAMPLE}_hifi.fastq.gz"
HIC_R1="$HIC_DIR/${SAMPLE}_R1.fq.gz"
HIC_R2="$HIC_DIR/${SAMPLE}_R2.fq.gz"
LOG_FILE="$LOG_DIR/hifiasm_hic_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"

echo "========================================" | tee "$LOG_FILE"
echo "Hifiasm Hi-C Integrated Assembly" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Sample: $SAMPLE" | tee -a "$LOG_FILE"
echo "HiFi reads: $HIFI" | tee -a "$LOG_FILE"
echo "Hi-C R1: $HIC_R1" | tee -a "$LOG_FILE"
echo "Hi-C R2: $HIC_R2" | tee -a "$LOG_FILE"
echo "Output: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Check files exist
if [ ! -f "$HIFI" ]; then
    echo "ERROR: HiFi reads not found: $HIFI" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$HIC_R1" ] || [ ! -f "$HIC_R2" ]; then
    echo "ERROR: Hi-C reads not found" | tee -a "$LOG_FILE"
    echo "Expected: $HIC_R1" | tee -a "$LOG_FILE"
    echo "Expected: $HIC_R2" | tee -a "$LOG_FILE"
    exit 1
fi

cd "$OUTPUT_DIR"

# Run hifiasm with Hi-C integration
echo "Running hifiasm with Hi-C integration..." | tee -a "$LOG_FILE"
echo "This will take 3-6 hours..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

hifiasm \
    -o ${SAMPLE} \
    --h1 "$HIC_R1" \
    --h2 "$HIC_R2" \
    -t 48 \
    -z 240m \
    -l 2 \
    "$HIFI" 2>&1 | tee -a "$LOG_FILE"

# Convert GFA to FASTA
echo "" | tee -a "$LOG_FILE"
echo "Converting GFA to FASTA..." | tee -a "$LOG_FILE"

# Check which GFA files were created
if [ -f "${SAMPLE}.hic.hap1.p_ctg.gfa" ]; then
    echo "Converting haplotype 1..." | tee -a "$LOG_FILE"
    awk '/^S/{print ">"$2"\n"$3}' ${SAMPLE}.hic.hap1.p_ctg.gfa > ${SAMPLE}.hap1.fasta
    
    echo "Converting haplotype 2..." | tee -a "$LOG_FILE"
    awk '/^S/{print ">"$2"\n"$3}' ${SAMPLE}.hic.hap2.p_ctg.gfa > ${SAMPLE}.hap2.fasta
    
    # Statistics
    echo "" | tee -a "$LOG_FILE"
    echo "=== HAPLOTYPE 1 ===" | tee -a "$LOG_FILE"
    seqkit stats -a ${SAMPLE}.hap1.fasta 2>&1 | tee -a "$LOG_FILE"
    
    echo "" | tee -a "$LOG_FILE"
    echo "=== HAPLOTYPE 2 ===" | tee -a "$LOG_FILE"
    seqkit stats -a ${SAMPLE}.hap2.fasta 2>&1 | tee -a "$LOG_FILE"
else
    echo "WARNING: Expected haplotype GFA files not found!" | tee -a "$LOG_FILE"
    echo "Checking for alternative output files..." | tee -a "$LOG_FILE"
    ls -lh ${SAMPLE}*.gfa 2>&1 | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Haplotype phasing complete!" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
