#!/bin/bash

##########################################################################
# Script: 01_setup_and_convert_bam.sh
# Purpose: Setup test sample and convert HiFi BAM to FASTQ for assembly
# Author: Nathan Fumia - Moringa Pangenome Project
# Date: November 2025
##########################################################################

# Exit on any error
# set -e # Temporally disabled for debugging
set -u
set -o pipefail

# Print commands as they execute
set -x

#########################################################################
# Configuration
#########################################################################

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate hifi_assembly

# Base directories
PROJECT_DIR=~/moringa_pangenome
HIFI_BASE=~/10193_78f4e76953b74a288de0f6e0bbd93029/central_memory2/pm_storage/caolinzhi/shujuxiazai/zhikong/163_20251121/SG251015HK01S67N1_福建农林大学Maringa_HiFi、HiC及转录组建库测序/BC2025100554-PB-hifi-9samples

# Test sample
SAMPLE="Mo-TH-30"

# Output Directories
RAW_DATA_DIR="${PROJECT_DIR}/raw_data/hifi"
TEST_DIR="${PROJECT_DIR}/assemblies/test/${SAMPLE}"
LOG_DIR="${PROJECT_DIR}/logs"

# Create directories
mkdir -p "${RAW_DATA_DIR}"
mkdir -p "${TEST_DIR}"
mkdir -p "${LOG_DIR}"

echo "=================================="
echo "HiFi Data Setup and Conversion"
echo "=================================="
echo "Sample: ${SAMPLE}"
echo "Start time: $(date)"
echo ""

#######################################################################
# STEP 1: Create symlink to HiFi BAM file
#######################################################################

echo "Step 1: Creating Symlink to HiFi BAM file..."

BAM_SOURCE="${HIFI_BASE}/${SAMPLE}/${SAMPLE}.bam"
BAM_LINK="${RAW_DATA_DIR}/${SAMPLE}.bam"

# Check if source BAM exists
if [ ! -f "${BAM_SOURCE}" ]; then
    echo "ERROR: Source BAM file not found: ${BAM_SOURCE}"
    exit 1
fi

# Create symlink (remove if already exists)
if [ -L "${BAM_LINK}" ]; then
    echo "Symlink already exists, removing old link..."
    rm "${BAM_LINK}"
fi

ln -s "${BAM_SOURCE}" "${BAM_LINK}"
echo "Symlink created: ${BAM_LINK} -> ${BAM_SOURCE}"
echo ""

################################################################################
# STEP 2: Check BAM file statistics
################################################################################

echo "Step 2: Checking BAM file statistics..."

# Get basic BAM info
echo "BAM file size: $(ls -lh ${BAM_LINK} | awk '{print $5}')"
echo ""

# Count reads in BAM
echo "Counting reads in BAM file (this may take a minute)..."
NUM_READS=$(samtools view -c "${BAM_LINK}")
echo "Total number of reads: ${NUM_READS}"
echo ""

# Get read length statistics (sample first 10000 reads)
echo "Sampling read lengths from first 10,000 reads..."
samtools view "${BAM_LINK}" | head -10000 | awk '{print length($10)}' | \
    awk '{sum+=$1; sumsq+=$1*$1; if(NR==1){min=max=$1}} 
         {if($1<min) min=$1; if($1>max) max=$1} 
         END {print "Mean length: "sum/NR; 
              print "Min length: "min; 
              print "Max length: "max}' | \
    tee "${LOG_DIR}/${SAMPLE}_read_stats.txt"

echo ""

################################################################################
# STEP 3: Convert BAM to FASTQ
################################################################################

echo "Step 3: Converting BAM to FASTQ..."
echo "This will take several minutes for a 7GB BAM file..."
echo ""

FASTQ_OUTPUT="${TEST_DIR}/${SAMPLE}.fastq.gz"

# Check if FASTQ already exists
if [ -f "${FASTQ_OUTPUT}" ]; then
    echo "FASTQ file already exists: ${FASTQ_OUTPUT}"
    echo "Skipping conversion. Delete file to reconvert."
else
    echo "Converting ${SAMPLE}.bam to FASTQ format..."
    echo "Output: ${FASTQ_OUTPUT}"
    
    # Convert BAM to FASTQ with bam2fq
    # -T '*' means include all tags (quality scores, etc)
    samtools bam2fq -T '*' "${BAM_LINK}" > "${TEST_DIR}/${SAMPLE}.fastq" 2>&1
    
    # Compress the FASTQ
    gzip "${TEST_DIR}/${SAMPLE}.fastq"
    
    echo "Conversion complete!"
fi

echo ""

################################################################################
# STEP 4: Verify FASTQ file
################################################################################

echo "Step 4: Verifying FASTQ file..."

# Check file size
echo "FASTQ file size: $(ls -lh ${FASTQ_OUTPUT} | awk '{print $5}')"

# Count reads in FASTQ (quick check)
echo "Counting reads in FASTQ..."
FASTQ_READS=$(zcat "${FASTQ_OUTPUT}" | wc -l | awk '{print $1/4}')
echo "Number of reads in FASTQ: ${FASTQ_READS}"

# Compare with BAM read count
if [ "${NUM_READS}" -eq "${FASTQ_READS}" ]; then
    echo "✓ Read counts match! Conversion successful."
else
    echo "⚠ WARNING: Read counts don't match!"
    echo "  BAM reads: ${NUM_READS}"
    echo "  FASTQ reads: ${FASTQ_READS}"
fi

echo ""

################################################################################
# STEP 5: Get detailed FASTQ statistics with seqkit
################################################################################

echo "Step 5: Generating detailed FASTQ statistics with seqkit..."

STATS_FILE="${LOG_DIR}/${SAMPLE}_fastq_stats.txt"

seqkit stats "${FASTQ_OUTPUT}" > "${STATS_FILE}"
cat "${STATS_FILE}"

echo ""
echo "Statistics saved to: ${STATS_FILE}"
echo ""

################################################################################
# STEP 6: Calculate estimated coverage
################################################################################

echo "Step 6: Calculating estimated genome coverage..."

# Moringa genome size
GENOME_SIZE=240000000  # 315 Mb

# Extract total bases from seqkit stats
TOTAL_BASES=$(grep -v "^file" "${STATS_FILE}" | awk '{print $5}' | sed 's/,//g')

# Calculate coverage
COVERAGE=$(awk "BEGIN {printf \"%.2f\", $TOTAL_BASES / $GENOME_SIZE}")

echo "Genome size: ${GENOME_SIZE} bp (315 Mb)"
echo "Total bases sequenced: ${TOTAL_BASES} bp"
echo "Estimated coverage: ${COVERAGE}x"
echo ""

# Check if coverage is adequate
if (( $(echo "$COVERAGE >= 20" | bc -l) )); then
    echo "✓ Coverage is adequate for high-quality assembly (>20x)"
elif (( $(echo "$COVERAGE >= 15" | bc -l) )); then
    echo "⚠ Coverage is marginal (15-20x). Assembly possible but may have gaps."
else
    echo "✗ WARNING: Coverage is low (<15x). Assembly quality may be poor."
fi

echo ""

################################################################################
# COMPLETION SUMMARY
################################################################################

echo "========================================"
echo "SETUP AND CONVERSION COMPLETE!"
echo "========================================"
echo "End time: $(date)"
echo ""
echo "Files created:"
echo "  - Symlink: ${BAM_LINK}"
echo "  - FASTQ: ${FASTQ_OUTPUT}"
echo "  - Stats: ${STATS_FILE}"
echo "  - Read stats: ${LOG_DIR}/${SAMPLE}_read_stats.txt"
echo ""
echo "Sample: ${SAMPLE}"
echo "Total reads: ${NUM_READS}"
echo "Coverage: ${COVERAGE}x"
echo ""
echo "Next step: Run assembly with hifiasm"
echo "  Script: 02_run_hifiasm_assembly.sh"
echo ""
