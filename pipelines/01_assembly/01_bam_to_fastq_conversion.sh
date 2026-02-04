#!/bin/bash
set -e

# ============================================================================
# Script: 01_bam_to_fastq_conversion.sh
# Purpose: Convert PacBio CCS BAM files to FASTQ format using AWK method
# Author: Moringa Pangenome Project
# Date: December 2, 2025
# ============================================================================

# Configuration
BAM="/home/hawaii-agriculture-research-cent/10193_78f4e76953b74a288de0f6e0bbd93029/central_memory2/pm_storage/caolinzhi/shujuxiazai/zhikong/163_20251121/SG251015HK01S67N1_福建农林大学Maringa_HiFi、HiC及转录组建库测序/BC2025100554-PB-hifi-9samples/Mo-TH-30/Mo-TH-30.bam"
SAMTOOLS="/home/hawaii-agriculture-research-cent/miniconda3/envs/hifi_assembly/bin/samtools"
OUTPUT_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome/assemblies/test/Mo-TH-30"
SAMPLE="Mo-TH-30"
LOG_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome/logs"

# Create directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Start logging
LOG_FILE="${LOG_DIR}/01_bam_to_fastq_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "${LOG_FILE}")
exec 2>&1

echo "=========================================="
echo "BAM to FASTQ Conversion - AWK Method"
echo "=========================================="
echo "Start time: $(date)"
echo "Sample: ${SAMPLE}"
echo "BAM file: ${BAM}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Step 1: Test with 100 reads
echo "[$(date +%H:%M:%S)] Step 1: Testing with first 100 reads..."
${SAMTOOLS} view ${BAM} | head -n 100 | \
  awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${OUTPUT_DIR}/${SAMPLE}_test100.fastq

echo "  Created: ${OUTPUT_DIR}/${SAMPLE}_test100.fastq"
ls -lh ${OUTPUT_DIR}/${SAMPLE}_test100.fastq

# Step 2: Validate test FASTQ
echo ""
echo "[$(date +%H:%M:%S)] Step 2: Validating test FASTQ..."
NUM_READS=$(grep -c "^@" ${OUTPUT_DIR}/${SAMPLE}_test100.fastq)
echo "  Number of reads: ${NUM_READS}"

echo ""
echo "  First read in FASTQ:"
head -n 4 ${OUTPUT_DIR}/${SAMPLE}_test100.fastq

# Validate format
NUM_AT=$(grep -c "^@" ${OUTPUT_DIR}/${SAMPLE}_test100.fastq)
NUM_PLUS=$(grep -c "^+$" ${OUTPUT_DIR}/${SAMPLE}_test100.fastq)
TOTAL_LINES=$(wc -l < ${OUTPUT_DIR}/${SAMPLE}_test100.fastq)
EXPECTED_LINES=$((NUM_AT * 4))

echo ""
echo "  Format check:"
echo "    @ lines: ${NUM_AT}"
echo "    + lines: ${NUM_PLUS}"
echo "    Total lines: ${TOTAL_LINES}"
echo "    Expected lines: ${EXPECTED_LINES}"

if [ ${TOTAL_LINES} -eq ${EXPECTED_LINES} ] && [ ${NUM_AT} -eq ${NUM_PLUS} ]; then
    echo "  ✓ FASTQ format is CORRECT!"
else
    echo "  ✗ FASTQ format has issues!"
    exit 1
fi

# Step 3: Check with seqkit
echo ""
echo "[$(date +%H:%M:%S)] Step 3: Validating with seqkit..."
if command -v seqkit &> /dev/null; then
    seqkit stats ${OUTPUT_DIR}/${SAMPLE}_test100.fastq
    echo "  ✓ seqkit validation passed!"
else
    echo "  Note: seqkit not found, skipping"
fi

# Step 4: User confirmation
echo ""
echo "=========================================="
echo "Test successful! Ready for full conversion."
echo "=========================================="
echo ""
echo "The full conversion will:"
echo "  - Extract all ~974,000 reads from the BAM file"
echo "  - Create a compressed FASTQ file (~7-8 GB)"
echo "  - Take approximately 10-15 minutes"
echo ""
read -p "Proceed with full conversion? (yes/no): " PROCEED

if [ "${PROCEED}" != "yes" ]; then
    echo "Conversion cancelled. Test file kept at: ${OUTPUT_DIR}/${SAMPLE}_test100.fastq"
    exit 0
fi

# Step 5: Full conversion
echo ""
echo "[$(date +%H:%M:%S)] Step 4: Starting full BAM to FASTQ conversion..."
echo "  This may take 10-15 minutes..."

${SAMTOOLS} view ${BAM} | \
  awk '{print "@"$1"\n"$10"\n+\n"$11}' | \
  gzip > ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz

echo "  ✓ Conversion complete!"
ls -lh ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz

# Step 6: Validate full file
echo ""
echo "[$(date +%H:%M:%S)] Step 5: Validating full FASTQ file..."
TOTAL_READS=$(zcat ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz | grep -c "^@")
echo "  Total reads extracted: ${TOTAL_READS}"
echo "  Expected reads: ~973,980"

if [ ${TOTAL_READS} -gt 900000 ] && [ ${TOTAL_READS} -lt 1000000 ]; then
    echo "  ✓ Read count looks correct!"
else
    echo "  ⚠ Read count seems unusual, please verify"
fi

# Step 7: Final seqkit stats
if command -v seqkit &> /dev/null; then
    echo ""
    echo "[$(date +%H:%M:%S)] Step 6: Final statistics with seqkit..."
    seqkit stats ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz
fi

# Clean up test file
rm ${OUTPUT_DIR}/${SAMPLE}_test100.fastq

echo ""
echo "=========================================="
echo "SUCCESS! Conversion complete."
echo "=========================================="
echo "End time: $(date)"
echo "Output file: ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz"
echo "Log file: ${LOG_FILE}"
echo ""
echo "Next steps:"
echo "  1. Validate: seqkit stats ${OUTPUT_DIR}/${SAMPLE}_hifi.fastq.gz"
echo "  2. Run assembly with: bash pipelines/02_run_hifiasm_assembly.sh"
echo ""
