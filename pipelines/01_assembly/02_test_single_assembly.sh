#!/bin/bash
set -e

# ============================================================================
# Script: 02_test_single_assembly.sh
# Purpose: Test hifiasm assembly on Mo-TH-30 sample
# Estimated time: 2-4 hours
# ============================================================================

# Configuration
BASE_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome"
SAMPLE="Mo-TH-30"
FASTQ="${BASE_DIR}/assemblies/test/${SAMPLE}/${SAMPLE}_hifi.fastq.gz"
OUTPUT_DIR="${BASE_DIR}/assemblies/test/${SAMPLE}"
OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}"
LOG_DIR="${BASE_DIR}/logs"
HIFIASM="hifiasm"

# Assembly parameters (using more threads for single job)
THREADS=48
GENOME_SIZE="240m"
PURGE_LEVEL=2

# Create directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Log file
LOG_FILE="${LOG_DIR}/test_hifiasm_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "${LOG_FILE}")
exec 2>&1

echo "============================================================================"
echo "Test hifiasm Assembly - ${SAMPLE}"
echo "============================================================================"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  Sample: ${SAMPLE}"
echo "  Input: ${FASTQ}"
echo "  Output: ${OUTPUT_DIR}"
echo "  Threads: ${THREADS} (using more since single job)"
echo "  Genome size: ${GENOME_SIZE}"
echo "  Purge level: ${PURGE_LEVEL}"
echo ""
echo "Estimated time: 2-4 hours"
echo "============================================================================"
echo ""

# Check if FASTQ exists
if [ ! -f "${FASTQ}" ]; then
    echo "ERROR: FASTQ file not found: ${FASTQ}"
    exit 1
fi

# Check if assembly already exists
if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
    echo "Assembly already exists!"
    ls -lh ${OUTPUT_PREFIX}.*
    echo ""
    read -p "Overwrite? (yes/no): " OVERWRITE
    if [ "${OVERWRITE}" != "yes" ]; then
        echo "Exiting."
        exit 0
    fi
    echo "Removing old assembly files..."
    rm -f ${OUTPUT_PREFIX}.*
fi

# Display input statistics
echo "Input FASTQ statistics:"
if command -v seqkit &> /dev/null; then
    seqkit stats ${FASTQ}
else
    echo "  File: ${FASTQ}"
    echo "  Size: $(du -h ${FASTQ} | cut -f1)"
    NUM_READS=$(zcat ${FASTQ} | grep -c "^@")
    echo "  Reads: ${NUM_READS}"
fi

echo ""
echo "============================================================================"
echo "Running hifiasm assembly..."
echo "============================================================================"
echo "Command:"
echo "  hifiasm --primary -t ${THREADS} -z ${GENOME_SIZE} -l ${PURGE_LEVEL} -o ${OUTPUT_PREFIX} ${FASTQ}"
echo ""
echo "This will take approximately 2-4 hours. You can monitor progress with:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "Start time: $(date)"
echo "----------------------------------------------------------------------------"

START_TIME=$(date +%s)

${HIFIASM} --primary \
           -t ${THREADS} \
           -z ${GENOME_SIZE} \
           -l ${PURGE_LEVEL} \
           -o ${OUTPUT_PREFIX} \
           ${FASTQ}

EXIT_CODE=$?
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "----------------------------------------------------------------------------"
echo "End time: $(date)"
echo "Duration: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""

# Check if successful
if [ ${EXIT_CODE} -eq 0 ] && [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
    echo "============================================================================"
    echo "✓ Assembly completed successfully!"
    echo "============================================================================"
    echo ""
    
    # List output files
    echo "Output files:"
    ls -lh ${OUTPUT_PREFIX}.*.gfa
    echo ""
    
    # Convert GFA to FASTA
    echo "Converting GFA to FASTA format..."
    
    if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
        awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.p_ctg.gfa > ${OUTPUT_PREFIX}.primary.fasta
        echo "  ✓ Primary: ${OUTPUT_PREFIX}.primary.fasta"
    fi
    
    if [ -f "${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa" ]; then
        awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa > ${OUTPUT_PREFIX}.hap1.fasta
        echo "  ✓ Haplotype 1: ${OUTPUT_PREFIX}.hap1.fasta"
    fi
    
    if [ -f "${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa" ]; then
        awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa > ${OUTPUT_PREFIX}.hap2.fasta
        echo "  ✓ Haplotype 2: ${OUTPUT_PREFIX}.hap2.fasta"
    fi
    
    echo ""
    echo "Assembly statistics:"
    if command -v seqkit &> /dev/null; then
        echo ""
        echo "Primary assembly:"
        seqkit stats ${OUTPUT_PREFIX}.primary.fasta
        
        if [ -f "${OUTPUT_PREFIX}.hap1.fasta" ]; then
            echo ""
            echo "Haplotype 1:"
            seqkit stats ${OUTPUT_PREFIX}.hap1.fasta
        fi
        
        if [ -f "${OUTPUT_PREFIX}.hap2.fasta" ]; then
            echo ""
            echo "Haplotype 2:"
            seqkit stats ${OUTPUT_PREFIX}.hap2.fasta
        fi
    else
        echo "  (install seqkit for detailed statistics)"
    fi
    
    echo ""
    echo "============================================================================"
    echo "SUCCESS! Test assembly complete."
    echo "============================================================================"
    echo ""
    echo "Key metrics to check:"
    echo "  - Assembly size: Should be ~240-260 Mb"
    echo "  - N50: Should be >1 Mb for HiFi assemblies"
    echo "  - Number of contigs: Typically 50-500"
    echo ""
    echo "If this looks good, proceed with batch assembly:"
    echo "  bash pipelines/02_batch_hifiasm_assembly.sh"
    echo ""
    
else
    echo "============================================================================"
    echo "✗ Assembly FAILED!"
    echo "============================================================================"
    echo "Exit code: ${EXIT_CODE}"
    echo ""
    echo "Check the log for errors: ${LOG_FILE}"
    echo ""
    exit 1
fi

echo "Log file: ${LOG_FILE}"
echo ""
