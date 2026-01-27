#!/bin/bash

################################################################################
# Pre-flight Check for Hi-C Alignment Pipeline
# 
# This script verifies that all required files and tools are available
# before running the main alignment pipeline
################################################################################

echo "========================================"
echo "Hi-C Pipeline Pre-flight Check"
echo "========================================"
echo ""

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate hic_align

ERRORS=0

################################################################################
# Check 1: Conda environment and tools
################################################################################

echo "[1/6] Checking conda environment and tools..."

if command -v bwa &> /dev/null; then
    echo "  ✓ BWA found: $(bwa 2>&1 | grep Version)"
else
    echo "  ✗ ERROR: BWA not found"
    ERRORS=$((ERRORS+1))
fi

if command -v samtools &> /dev/null; then
    echo "  ✓ samtools found: $(samtools --version | head -1)"
else
    echo "  ✗ ERROR: samtools not found"
    ERRORS=$((ERRORS+1))
fi

echo ""

################################################################################
# Check 2: Reference genome
################################################################################

echo "[2/6] Checking reference genome..."

REFERENCE_DIR=~/moringa_reference
REFERENCE_GENOME="${REFERENCE_DIR}/MoringaV2.genome.fa"

if [ -f "${REFERENCE_GENOME}" ]; then
    SIZE=$(ls -lh "${REFERENCE_GENOME}" | awk '{print $5}')
    echo "  ✓ Reference genome found: ${REFERENCE_GENOME} (${SIZE})"
else
    echo "  ✗ ERROR: Reference genome not found at ${REFERENCE_GENOME}"
    ERRORS=$((ERRORS+1))
fi

echo ""

################################################################################
# Check 3: Sample files
################################################################################

echo "[3/6] Checking sample files..."

ARCHIVE1=~/9842_04198ee0ed804489b6136a9ec70ffb70
ARCHIVE2=~/9843_502ded7b18ba4ee6ab7982104bc3b562

SAMPLES=(
    "Mo-TH-30"
    "Mo-TH-55"
    "Mo-TH-16"
    "Mo-TH-43"
    "Mo-TH-6"
    "Mo-TH-63"
    "Mo-TH-66"
    "Mo-TW-13"
    "Mo-US-5"
)

FOUND=0

for SAMPLE in "${SAMPLES[@]}"; do
    # Search in both archives
    R1=$(find "${ARCHIVE1}" "${ARCHIVE2}" -name "${SAMPLE}_R1.fq.gz" 2>/dev/null | head -1)
    R2=$(find "${ARCHIVE1}" "${ARCHIVE2}" -name "${SAMPLE}_R2.fq.gz" 2>/dev/null | head -1)
    
    if [ -n "$R1" ] && [ -n "$R2" ]; then
        R1_SIZE=$(ls -lh "$R1" | awk '{print $5}')
        R2_SIZE=$(ls -lh "$R2" | awk '{print $5}')
        echo "  ✓ ${SAMPLE}: R1=${R1_SIZE}, R2=${R2_SIZE}"
        FOUND=$((FOUND+1))
    else
        echo "  ✗ ${SAMPLE}: Files not found"
        ERRORS=$((ERRORS+1))
    fi
done

echo ""
echo "  Summary: Found ${FOUND}/9 sample pairs"

echo ""

################################################################################
# Check 4: Disk space
################################################################################

echo "[4/6] Checking available disk space..."

AVAILABLE=$(df -h ~ | tail -1 | awk '{print $4}')
AVAILABLE_GB=$(df ~ | tail -1 | awk '{print $4}')
AVAILABLE_GB=$((AVAILABLE_GB / 1024 / 1024))

echo "  Available disk space in home directory: ${AVAILABLE} (~${AVAILABLE_GB} GB)"

if [ ${AVAILABLE_GB} -lt 500 ]; then
    echo "  ⚠ WARNING: Less than 500 GB available. Hi-C alignment may need 500-1000 GB"
    echo "             Consider using a different output location with more space"
else
    echo "  ✓ Sufficient disk space available"
fi

echo ""

################################################################################
# Check 5: System resources
################################################################################

echo "[5/6] Checking system resources..."

CORES=$(nproc)
MEM_GB=$(free -g | grep Mem | awk '{print $2}')

echo "  CPU cores: ${CORES}"
echo "  Total RAM: ${MEM_GB} GB"

if [ ${CORES} -ge 8 ] && [ ${MEM_GB} -ge 32 ]; then
    echo "  ✓ System resources look good for Hi-C alignment"
elif [ ${CORES} -lt 4 ] || [ ${MEM_GB} -lt 16 ]; then
    echo "  ⚠ WARNING: Limited system resources. Pipeline may be slow."
    echo "             Recommended: 8+ cores, 32+ GB RAM"
else
    echo "  ⚠ System resources are adequate but not optimal"
fi

echo ""

################################################################################
# Check 6: Estimated runtime and storage
################################################################################

echo "[6/6] Estimating pipeline requirements..."

echo "  Estimated per-sample processing:"
echo "    - Time: 2-6 hours (depending on coverage)"
echo "    - Storage: 50-100 GB per sample for BAM files"
echo ""
echo "  Total for 9 samples:"
echo "    - Time: 18-54 hours (can run in parallel)"
echo "    - Storage: ~450-900 GB"
echo ""

################################################################################
# Summary
################################################################################

echo "========================================"
echo "Pre-flight Check Summary"
echo "========================================"

if [ ${ERRORS} -eq 0 ]; then
    echo "✓ All checks passed! Ready to run the pipeline."
    echo ""
    echo "To run the pipeline:"
    echo "  bash ~/hic_alignment_pipeline.sh"
    echo ""
    echo "Or run in background with logging:"
    echo "  nohup bash ~/hic_alignment_pipeline.sh > ~/hic_pipeline.log 2>&1 &"
    echo ""
else
    echo "✗ ${ERRORS} error(s) found. Please fix before running the pipeline."
    echo ""
fi

echo "========================================"
