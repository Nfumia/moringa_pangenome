#!/bin/bash
set -e

# ============================================================================
# Script: 01_batch_bam_to_fastq.sh
# Purpose: Batch convert all 9 PacBio CCS BAM files to FASTQ format
# Author: Moringa Pangenome Project
# Date: December 2, 2025
# ============================================================================

# Configuration
BASE_DIR="/home/hawaii-agriculture-research-cent"
BAM_BASE="${BASE_DIR}/10193_78f4e76953b74a288de0f6e0bbd93029/central_memory2/pm_storage/caolinzhi/shujuxiazai/zhikong/163_20251121/SG251015HK01S67N1_福建农林大学Maringa_HiFi、HiC及转录组建库测序/BC2025100554-PB-hifi-9samples"
SAMTOOLS="${BASE_DIR}/miniconda3/envs/hifi_assembly/bin/samtools"
OUTPUT_BASE="${BASE_DIR}/moringa_pangenome/raw_data/hifi"
LOG_DIR="${BASE_DIR}/moringa_pangenome/logs"

# Sample list (9 samples)
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

# Create directories
mkdir -p ${OUTPUT_BASE}
mkdir -p ${LOG_DIR}

# Main log file
MAIN_LOG="${LOG_DIR}/batch_bam_to_fastq_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "${MAIN_LOG}")
exec 2>&1

echo "============================================================================"
echo "Batch BAM to FASTQ Conversion"
echo "============================================================================"
echo "Start time: $(date)"
echo "Number of samples: ${#SAMPLES[@]}"
echo "Output directory: ${OUTPUT_BASE}"
echo ""

# Summary arrays
declare -a SUCCESS_SAMPLES
declare -a FAILED_SAMPLES

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "========================================================================"
    echo "Processing sample: ${SAMPLE}"
    echo "========================================================================"
    echo "Start time: $(date)"
    
    # Define paths
    BAM="${BAM_BASE}/${SAMPLE}/${SAMPLE}.bam"
    OUTPUT="${OUTPUT_BASE}/${SAMPLE}_hifi.fastq.gz"
    SAMPLE_LOG="${LOG_DIR}/bam_to_fastq_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"
    
    # Check if BAM file exists
    if [ ! -f "${BAM}" ]; then
        echo "ERROR: BAM file not found: ${BAM}"
        FAILED_SAMPLES+=("${SAMPLE}")
        continue
    fi
    
    # Check if output already exists
    if [ -f "${OUTPUT}" ]; then
        echo "Output file already exists: ${OUTPUT}"
        read -p "Overwrite? (yes/no): " OVERWRITE
        if [ "${OVERWRITE}" != "yes" ]; then
            echo "Skipping ${SAMPLE}"
            continue
        fi
        rm "${OUTPUT}"
    fi
    
    # Get BAM stats
    echo "Getting BAM statistics..."
    NUM_READS=$(${SAMTOOLS} view -c ${BAM})
    echo "  Total reads in BAM: ${NUM_READS}"
    
    # Run conversion
    echo "Converting BAM to FASTQ..."
    echo "  Command: samtools view | awk | gzip"
    
    START_TIME=$(date +%s)
    
    ${SAMTOOLS} view ${BAM} | \
      awk '{print "@"$1"\n"$10"\n+\n"$11}' | \
      gzip > ${OUTPUT} 2>&1 | tee -a "${SAMPLE_LOG}"
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    
    echo "  Conversion completed in ${DURATION} seconds"
    
    # Validate output
    echo "Validating FASTQ file..."
    if [ -f "${OUTPUT}" ]; then
        FILE_SIZE=$(du -h ${OUTPUT} | cut -f1)
        echo "  File size: ${FILE_SIZE}"
        
        # Count reads in FASTQ
        echo "  Counting reads (this may take a moment)..."
        FASTQ_READS=$(zcat ${OUTPUT} | grep -c "^@")
        echo "  Reads in FASTQ: ${FASTQ_READS}"
        
        # Compare counts
        if [ ${FASTQ_READS} -eq ${NUM_READS} ]; then
            echo "  ✓ Read count matches! Conversion successful."
            SUCCESS_SAMPLES+=("${SAMPLE}")
            
            # Run seqkit stats if available
            if command -v seqkit &> /dev/null; then
                echo "  seqkit statistics:"
                seqkit stats ${OUTPUT}
            fi
        else
            echo "  ✗ Read count mismatch!"
            echo "    Expected: ${NUM_READS}"
            echo "    Got: ${FASTQ_READS}"
            FAILED_SAMPLES+=("${SAMPLE}")
        fi
    else
        echo "  ✗ Output file not created!"
        FAILED_SAMPLES+=("${SAMPLE}")
    fi
    
    echo "End time: $(date)"
done

# Final summary
echo ""
echo "============================================================================"
echo "Batch Conversion Complete"
echo "============================================================================"
echo "End time: $(date)"
echo ""
echo "Summary:"
echo "  Total samples: ${#SAMPLES[@]}"
echo "  Successful: ${#SUCCESS_SAMPLES[@]}"
echo "  Failed: ${#FAILED_SAMPLES[@]}"
echo ""

if [ ${#SUCCESS_SAMPLES[@]} -gt 0 ]; then
    echo "Successful samples:"
    for SAMPLE in "${SUCCESS_SAMPLES[@]}"; do
        echo "  ✓ ${SAMPLE}"
    done
    echo ""
fi

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "Failed samples:"
    for SAMPLE in "${FAILED_SAMPLES[@]}"; do
        echo "  ✗ ${SAMPLE}"
    done
    echo ""
fi

echo "Log file: ${MAIN_LOG}"
echo "Output directory: ${OUTPUT_BASE}"
echo ""
echo "Next steps:"
echo "  1. Review conversion logs in: ${LOG_DIR}"
echo "  2. Verify all FASTQ files: ls -lh ${OUTPUT_BASE}"
echo "  3. Run assembly with: bash pipelines/02_batch_hifiasm_assembly.sh"
echo ""
