#!/bin/bash

################################################################################
# Script: 02_run_hifiasm_assembly.sh
# Purpose: Run hifiasm assembly on HiFi BAM data
# Author: Moringa Pangenome Project
# Date: December 2025
################################################################################

# Exit on error
set -e
set -u
set -o pipefail

# Print commands (for debugging)
set -x

################################################################################
# CONFIGURATION
################################################################################

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate hifi_assembly

# Base directories
PROJECT_DIR=~/moringa_pangenome
HIFI_BASE=~/10193_78f4e76953b74a288de0f6e0bbd93029/central_memory2/pm_storage/caolinzhi/shujuxiazai/zhikong/163_20251121/SG251015HK01S67N1_福建农林大学Maringa_HiFi、HiC及转录组建库测序/BC2025100554-PB-hifi-9samples

# Test sample
SAMPLE="Mo-TH-30"

# Input/Output directories
TEST_DIR="${PROJECT_DIR}/assemblies/test/${SAMPLE}"
LOG_DIR="${PROJECT_DIR}/logs"

# HiFi BAM file
HIFI_BAM="${HIFI_BASE}/${SAMPLE}/${SAMPLE}.bam"

# Output prefix for hifiasm
OUTPUT_PREFIX="${TEST_DIR}/${SAMPLE}_hifiasm"

# Number of threads
THREADS=32

# Create directories
mkdir -p "${TEST_DIR}"
mkdir -p "${LOG_DIR}"

echo "========================================"
echo "hifiasm Assembly for Moringa"
echo "========================================"
echo "Sample: ${SAMPLE}"
echo "Threads: ${THREADS}"
echo "Start time: $(date)"
echo ""

################################################################################
# STEP 1: Verify input BAM file exists
################################################################################

echo "Step 1: Verifying HiFi BAM file..."

if [ ! -f "${HIFI_BAM}" ]; then
    echo "ERROR: HiFi BAM file not found: ${HIFI_BAM}"
    exit 1
fi

echo "HiFi BAM found: ${HIFI_BAM}"
BAM_SIZE=$(ls -lh "${HIFI_BAM}" | awk '{print $5}')
echo "BAM file size: ${BAM_SIZE}"
echo ""

################################################################################
# STEP 2: Check if assembly already exists
################################################################################

echo "Step 2: Checking for existing assembly..."

# hifiasm creates several output files:
# .bp.p_ctg.gfa - primary contigs (phased)
# .bp.hap1.p_ctg.gfa - haplotype 1
# .bp.hap2.p_ctg.gfa - haplotype 2

if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
    echo "WARNING: Assembly output already exists!"
    echo "Files found:"
    ls -lh "${OUTPUT_PREFIX}"*.gfa 2>/dev/null || echo "  (some files exist)"
    echo ""
    echo "Delete these files to re-run assembly, or use a different output prefix."
    echo "Exiting to prevent overwriting existing assembly."
    exit 0
fi

echo "No existing assembly found. Proceeding with assembly..."
echo ""

################################################################################
# STEP 3: Run hifiasm assembly
################################################################################

echo "========================================"
echo "Step 3: Running hifiasm assembly"
echo "========================================"
echo ""
echo "Input: ${HIFI_BAM}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Threads: ${THREADS}"
echo ""
echo "This will take several hours (2-6 hours typical for 315 Mb genome)"
echo "Assembly started at: $(date)"
echo ""

# Run hifiasm
# -o: output prefix
# -t: number of threads
# --primary: output primary contigs in addition to haplotype-resolved contigs
# The BAM file is the last argument

hifiasm \
    -o "${OUTPUT_PREFIX}" \
    -t ${THREADS} \
    --primary \
    "${HIFI_BAM}" \
    2>&1 | tee "${LOG_DIR}/${SAMPLE}_hifiasm.log"

echo ""
echo "Assembly completed at: $(date)"
echo ""

################################################################################
# STEP 4: Check assembly outputs and get basic statistics
################################################################################

echo "========================================"
echo "Step 4: Assembly output summary"
echo "========================================"
echo ""

# List all output files
echo "Output files created:"
ls -lh "${OUTPUT_PREFIX}"* | awk '{print "  " $9 " - " $5}'
echo ""

# Check for expected outputs
if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
    echo "✓ Primary contigs: ${OUTPUT_PREFIX}.bp.p_ctg.gfa"
else
    echo "✗ WARNING: Primary contigs not found!"
fi

if [ -f "${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa" ]; then
    echo "✓ Haplotype 1: ${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa"
else
    echo "✗ WARNING: Haplotype 1 not found!"
fi

if [ -f "${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa" ]; then
    echo "✓ Haplotype 2: ${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa"
else
    echo "✗ WARNING: Haplotype 2 not found!"
fi

echo ""

################################################################################
# STEP 5: Convert GFA to FASTA for easier analysis
################################################################################

echo "========================================"
echo "Step 5: Converting GFA to FASTA format"
echo "========================================"
echo ""

# Function to convert GFA to FASTA
convert_gfa_to_fasta() {
    local gfa_file=$1
    local fasta_file=$2
    
    if [ -f "${gfa_file}" ]; then
        echo "Converting: ${gfa_file}"
        awk '/^S/{print ">"$2; print $3}' "${gfa_file}" > "${fasta_file}"
        echo "  Created: ${fasta_file}"
        
        # Get basic stats
        NUM_CONTIGS=$(grep -c "^>" "${fasta_file}")
        echo "  Number of contigs: ${NUM_CONTIGS}"
    else
        echo "File not found: ${gfa_file}"
    fi
}

# Convert each assembly
convert_gfa_to_fasta "${OUTPUT_PREFIX}.bp.p_ctg.gfa" "${OUTPUT_PREFIX}.primary.fasta"
echo ""
convert_gfa_to_fasta "${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa" "${OUTPUT_PREFIX}.hap1.fasta"
echo ""
convert_gfa_to_fasta "${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa" "${OUTPUT_PREFIX}.hap2.fasta"
echo ""

################################################################################
# STEP 6: Get assembly statistics with seqkit
################################################################################

echo "========================================"
echo "Step 6: Detailed assembly statistics"
echo "========================================"
echo ""

STATS_FILE="${LOG_DIR}/${SAMPLE}_assembly_stats.txt"

echo "Assembly Statistics for ${SAMPLE}" > "${STATS_FILE}"
echo "Generated: $(date)" >> "${STATS_FILE}"
echo "========================================" >> "${STATS_FILE}"
echo "" >> "${STATS_FILE}"

# Get stats for each assembly
for FASTA in "${OUTPUT_PREFIX}.primary.fasta" \
             "${OUTPUT_PREFIX}.hap1.fasta" \
             "${OUTPUT_PREFIX}.hap2.fasta"; do
    
    if [ -f "${FASTA}" ]; then
        BASENAME=$(basename "${FASTA}")
        echo "Statistics for ${BASENAME}:" | tee -a "${STATS_FILE}"
        echo "----------------------------------------" | tee -a "${STATS_FILE}"
        seqkit stats "${FASTA}" | tee -a "${STATS_FILE}"
        echo "" | tee -a "${STATS_FILE}"
        
        # Calculate N50
        echo "N50 statistics:" | tee -a "${STATS_FILE}"
        seqkit stats -a "${FASTA}" | grep -v "^file" | \
            awk '{printf "  N50: %s\n  Longest contig: %s\n", $11, $7}' | \
            tee -a "${STATS_FILE}"
        echo "" | tee -a "${STATS_FILE}"
    fi
done

echo "Detailed statistics saved to: ${STATS_FILE}"
echo ""

################################################################################
# COMPLETION SUMMARY
################################################################################

echo "========================================"
echo "ASSEMBLY COMPLETE!"
echo "========================================"
echo "End time: $(date)"
echo ""
echo "Sample: ${SAMPLE}"
echo "Input: ${HIFI_BAM}"
echo ""
echo "Output files:"
echo "  - GFA assemblies: ${OUTPUT_PREFIX}*.gfa"
echo "  - FASTA assemblies: ${OUTPUT_PREFIX}*.fasta"
echo "  - Assembly log: ${LOG_DIR}/${SAMPLE}_hifiasm.log"
echo "  - Statistics: ${STATS_FILE}"
echo ""
echo "Assembly types generated:"
echo "  1. Primary assembly (collapsed haplotypes)"
echo "  2. Haplotype 1 (phased)"
echo "  3. Haplotype 2 (phased)"
echo ""
echo "Next steps:"
echo "  1. Review assembly statistics"
echo "  2. Run BUSCO quality assessment"
echo "  3. Align to reference genome for comparison"
echo "  4. Repeat for remaining 8 samples"
echo ""
echo "To view statistics:"
echo "  cat ${STATS_FILE}"
echo ""
