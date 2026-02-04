#!/bin/bash
set -e

# ============================================================================
# Script: 02_batch_hifiasm_assembly.sh
# Purpose: Batch hifiasm assembly for all 9 Moringa samples
# Author: Moringa Pangenome Project
# Date: December 2, 2025
# ============================================================================

# Configuration
BASE_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome"
FASTQ_DIR="${BASE_DIR}/raw_data/hifi"
OUTPUT_DIR="${BASE_DIR}/assemblies/production"
LOG_DIR="${BASE_DIR}/logs"
HIFIASM="hifiasm"  # Assuming it's in PATH from conda env

# Assembly parameters
THREADS=32           # Use 32 threads per job (allows 2 jobs simultaneously)
GENOME_SIZE="240m"   # Moringa genome size (AOCC v2 2022 chromosome-scale assembly)
PURGE_LEVEL=2        # Purge level for heterozygosity (0-3)

# Sample list
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
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Main log file
MAIN_LOG="${LOG_DIR}/batch_hifiasm_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "${MAIN_LOG}")
exec 2>&1

echo "============================================================================"
echo "Batch hifiasm Assembly - Moringa Pangenome Project"
echo "============================================================================"
echo "Start time: $(date)"
echo "Number of samples: ${#SAMPLES[@]}"
echo ""
echo "Assembly parameters:"
echo "  Threads per job: ${THREADS}"
echo "  Genome size: ${GENOME_SIZE}"
echo "  Purge level: ${PURGE_LEVEL}"
echo "  Mode: --primary (creates primary + hap1 + hap2)"
echo ""
echo "Expected outputs per sample:"
echo "  - .bp.p_ctg.gfa (primary assembly)"
echo "  - .bp.hap1.p_ctg.gfa (haplotype 1)"
echo "  - .bp.hap2.p_ctg.gfa (haplotype 2)"
echo ""
echo "Estimated time per sample: 2-6 hours"
echo "Total estimated time: 18-54 hours for all 9 samples"
echo "============================================================================"
echo ""

# Check if running in parallel
echo "Assembly mode options:"
echo "  1. Sequential (one at a time, safest, ~18-54 hours total)"
echo "  2. Parallel 2x (two simultaneously, faster, ~9-27 hours total)"
echo "  3. Parallel 3x (three simultaneously, fastest, uses all cores)"
echo ""
read -p "Select mode (1/2/3): " MODE

if [ "${MODE}" == "2" ]; then
    MAX_PARALLEL=2
    echo "Running 2 assemblies in parallel"
elif [ "${MODE}" == "3" ]; then
    MAX_PARALLEL=3
    THREADS=21  # Adjust threads: 64/3 = ~21 per job
    echo "Running 3 assemblies in parallel (threads adjusted to ${THREADS})"
else
    MAX_PARALLEL=1
    THREADS=48  # Use more threads for sequential
    echo "Running assemblies sequentially (threads increased to ${THREADS})"
fi

echo ""

# Summary arrays
declare -a SUCCESS_SAMPLES
declare -a FAILED_SAMPLES
declare -a SKIPPED_SAMPLES

# Function to run single assembly
run_assembly() {
    local SAMPLE=$1
    local FASTQ="${FASTQ_DIR}/${SAMPLE}_hifi.fastq.gz"
    local SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
    local OUTPUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}"
    local SAMPLE_LOG="${LOG_DIR}/hifiasm_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"
    
    echo ""
    echo "========================================================================"
    echo "Assembling: ${SAMPLE}"
    echo "========================================================================"
    echo "Start time: $(date)"
    
    # Create sample directory
    mkdir -p ${SAMPLE_DIR}
    
    # Check if FASTQ exists
    if [ ! -f "${FASTQ}" ]; then
        echo "ERROR: FASTQ file not found: ${FASTQ}"
        echo "Skipping ${SAMPLE}"
        SKIPPED_SAMPLES+=("${SAMPLE}")
        return 1
    fi
    
    # Check if assembly already exists
    if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
        echo "Assembly already exists: ${OUTPUT_PREFIX}.bp.p_ctg.gfa"
        read -p "Overwrite? (yes/no): " OVERWRITE
        if [ "${OVERWRITE}" != "yes" ]; then
            echo "Skipping ${SAMPLE}"
            SKIPPED_SAMPLES+=("${SAMPLE}")
            return 0
        fi
        echo "Removing old assembly files..."
        rm -f ${OUTPUT_PREFIX}.*
    fi
    
    # Get FASTQ stats
    echo "Input FASTQ statistics:"
    if command -v seqkit &> /dev/null; then
        seqkit stats ${FASTQ}
    else
        echo "  File: ${FASTQ}"
        echo "  Size: $(du -h ${FASTQ} | cut -f1)"
    fi
    
    echo ""
    echo "Running hifiasm..."
    echo "Command: hifiasm --primary -t ${THREADS} -z ${GENOME_SIZE} -l ${PURGE_LEVEL} -o ${OUTPUT_PREFIX} ${FASTQ}"
    
    START_TIME=$(date +%s)
    
    # Run hifiasm
    ${HIFIASM} --primary \
               -t ${THREADS} \
               -z ${GENOME_SIZE} \
               -l ${PURGE_LEVEL} \
               -o ${OUTPUT_PREFIX} \
               ${FASTQ} \
               2>&1 | tee ${SAMPLE_LOG}
    
    EXIT_CODE=${PIPESTATUS[0]}
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    
    echo ""
    echo "Assembly completed in ${HOURS}h ${MINUTES}m"
    
    # Check outputs
    if [ ${EXIT_CODE} -eq 0 ] && [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
        echo "✓ Assembly successful!"
        
        # Convert GFA to FASTA
        echo "Converting GFA to FASTA..."
        
        # Primary assembly
        if [ -f "${OUTPUT_PREFIX}.bp.p_ctg.gfa" ]; then
            awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.p_ctg.gfa > ${OUTPUT_PREFIX}.primary.fasta
            echo "  Primary assembly: ${OUTPUT_PREFIX}.primary.fasta"
        fi
        
        # Haplotype 1
        if [ -f "${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa" ]; then
            awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.hap1.p_ctg.gfa > ${OUTPUT_PREFIX}.hap1.fasta
            echo "  Haplotype 1: ${OUTPUT_PREFIX}.hap1.fasta"
        fi
        
        # Haplotype 2
        if [ -f "${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa" ]; then
            awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.bp.hap2.p_ctg.gfa > ${OUTPUT_PREFIX}.hap2.fasta
            echo "  Haplotype 2: ${OUTPUT_PREFIX}.hap2.fasta"
        fi
        
        # Get assembly statistics
        echo ""
        echo "Assembly statistics:"
        if command -v seqkit &> /dev/null; then
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
        fi
        
        SUCCESS_SAMPLES+=("${SAMPLE}")
    else
        echo "✗ Assembly failed!"
        FAILED_SAMPLES+=("${SAMPLE}")
    fi
    
    echo "End time: $(date)"
}

# Export function for parallel execution
export -f run_assembly
export FASTQ_DIR OUTPUT_DIR LOG_DIR HIFIASM THREADS GENOME_SIZE PURGE_LEVEL
export SUCCESS_SAMPLES FAILED_SAMPLES SKIPPED_SAMPLES

# Run assemblies
if [ ${MAX_PARALLEL} -eq 1 ]; then
    # Sequential
    for SAMPLE in "${SAMPLES[@]}"; do
        run_assembly ${SAMPLE}
    done
else
    # Parallel using GNU parallel or xargs
    if command -v parallel &> /dev/null; then
        printf "%s\n" "${SAMPLES[@]}" | parallel -j ${MAX_PARALLEL} run_assembly {}
    else
        # Fallback: simple background jobs
        RUNNING=0
        for SAMPLE in "${SAMPLES[@]}"; do
            run_assembly ${SAMPLE} &
            RUNNING=$((RUNNING + 1))
            
            if [ ${RUNNING} -ge ${MAX_PARALLEL} ]; then
                wait -n  # Wait for any job to finish
                RUNNING=$((RUNNING - 1))
            fi
        done
        wait  # Wait for all remaining jobs
    fi
fi

# Final summary
echo ""
echo "============================================================================"
echo "Batch Assembly Complete"
echo "============================================================================"
echo "End time: $(date)"
echo ""
echo "Summary:"
echo "  Total samples: ${#SAMPLES[@]}"
echo "  Successful: ${#SUCCESS_SAMPLES[@]}"
echo "  Failed: ${#FAILED_SAMPLES[@]}"
echo "  Skipped: ${#SKIPPED_SAMPLES[@]}"
echo ""

if [ ${#SUCCESS_SAMPLES[@]} -gt 0 ]; then
    echo "Successful assemblies:"
    for SAMPLE in "${SUCCESS_SAMPLES[@]}"; do
        echo "  ✓ ${SAMPLE}"
        echo "    Primary: ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.primary.fasta"
    done
    echo ""
fi

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "Failed assemblies:"
    for SAMPLE in "${FAILED_SAMPLES[@]}"; do
        echo "  ✗ ${SAMPLE}"
    done
    echo ""
fi

if [ ${#SKIPPED_SAMPLES[@]} -gt 0 ]; then
    echo "Skipped samples:"
    for SAMPLE in "${SKIPPED_SAMPLES[@]}"; do
        echo "  - ${SAMPLE}"
    done
    echo ""
fi

echo "Main log: ${MAIN_LOG}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Next steps:"
echo "  1. Review assembly logs in: ${LOG_DIR}"
echo "  2. Check assembly quality: seqkit stats ${OUTPUT_DIR}/*/*.primary.fasta"
echo "  3. Compare assemblies to reference"
echo ""
