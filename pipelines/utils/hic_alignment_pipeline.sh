#!/bin/bash

################################################################################
# Hi-C Alignment Pipeline for Moringa oleifera Samples
# 
# This script aligns Hi-C paired-end reads to the Moringa v2 reference genome
# 
# Author: Analysis pipeline for pangenomics project
# Date: November 2025
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

################################################################################
# CONFIGURATION
################################################################################

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate hic_align

# Set paths
BASE_DIR=~/moringa_pangenome
REFERENCE_DIR=~/moringa_pangenome/moringa_reference
REFERENCE_GENOME="${REFERENCE_DIR}/MoringaV2.genome.fa"
WORK_DIR="${BASE_DIR}/hic_analysis"
OUTPUT_DIR="${WORK_DIR}/alignments"
LOG_DIR="${WORK_DIR}/logs"
STATS_DIR="${WORK_DIR}/stats"

# Archive paths
ARCHIVE1=~/9842_04198ee0ed804489b6136a9ec70ffb70
ARCHIVE2=~/9843_502ded7b18ba4ee6ab7982104bc3b562

# Number of threads (adjust based on your system - you have 64 cores!)
THREADS=16

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

################################################################################
# SETUP
################################################################################

echo "========================================"
echo "Hi-C Alignment Pipeline for Moringa"
echo "========================================"
echo "Start time: $(date)"
echo ""

# Create output directories
mkdir -p "${WORK_DIR}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${STATS_DIR}"

echo "Working directory: ${WORK_DIR}"
echo "Reference genome: ${REFERENCE_GENOME}"
echo "Number of threads: ${THREADS}"
echo "Number of samples: ${#SAMPLES[@]}"
echo ""

################################################################################
# STEP 1: INDEX REFERENCE GENOME
################################################################################

echo "========================================"
echo "STEP 1: Indexing reference genome"
echo "========================================"

if [ ! -f "${REFERENCE_GENOME}.bwt" ]; then
    echo "Creating BWA index for reference genome..."
    bwa index "${REFERENCE_GENOME}" 2>&1 | tee "${LOG_DIR}/bwa_index.log"
    echo "BWA indexing complete!"
else
    echo "BWA index already exists, skipping..."
fi

# Create samtools fai index
if [ ! -f "${REFERENCE_GENOME}.fai" ]; then
    echo "Creating samtools fai index..."
    samtools faidx "${REFERENCE_GENOME}"
    echo "Samtools indexing complete!"
else
    echo "Samtools fai index already exists, skipping..."
fi

echo ""

################################################################################
# FUNCTION: Find sample fastq files
################################################################################

find_sample_files() {
    local sample=$1
    local r1_file=""
    local r2_file=""
    
    # Search in archive 1
    r1_file=$(find "${ARCHIVE1}" -name "${sample}_R1.fq.gz" 2>/dev/null | head -1)
    r2_file=$(find "${ARCHIVE1}" -name "${sample}_R2.fq.gz" 2>/dev/null | head -1)
    
    # If not found, search in archive 2
    if [ -z "$r1_file" ]; then
        r1_file=$(find "${ARCHIVE2}" -name "${sample}_R1.fq.gz" 2>/dev/null | head -1)
        r2_file=$(find "${ARCHIVE2}" -name "${sample}_R2.fq.gz" 2>/dev/null | head -1)
    fi
    
    echo "$r1_file|$r2_file"
}

################################################################################
# STEP 2: ALIGN SAMPLES
################################################################################

echo "========================================"
echo "STEP 2: Aligning Hi-C reads to reference"
echo "========================================"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "Processing sample: ${SAMPLE}"
    echo "----------------------------------------"
    
    # Find input files
    FILES=$(find_sample_files "${SAMPLE}")
    R1=$(echo "$FILES" | cut -d'|' -f1)
    R2=$(echo "$FILES" | cut -d'|' -f2)
    
    if [ -z "$R1" ] || [ -z "$R2" ]; then
        echo "ERROR: Could not find fastq files for ${SAMPLE}"
        echo "Skipping this sample..."
        continue
    fi
    
    echo "R1: $R1"
    echo "R2: $R2"
    
    # Output files
    RAW_BAM="${OUTPUT_DIR}/${SAMPLE}.raw.bam"
    SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}.sorted.bam"
    FINAL_BAM="${OUTPUT_DIR}/${SAMPLE}.final.bam"
    
    # Skip if final BAM already exists
    if [ -f "${FINAL_BAM}" ]; then
        echo "Final BAM already exists for ${SAMPLE}, skipping..."
        continue
    fi
    
    # STEP 2A: Alignment with BWA-MEM
    echo "Running BWA-MEM alignment..."
    echo "Start: $(date)"
    bwa mem \
        -t ${THREADS} \
        -SP \
        "${REFERENCE_GENOME}" \
        "${R1}" \
        "${R2}" \
        2> "${LOG_DIR}/${SAMPLE}_bwa.log" \
        | samtools view -@ ${THREADS} -bS - > "${RAW_BAM}"
    
    echo "Alignment complete for ${SAMPLE}"
    
    # STEP 2B: Sort BAM
    echo "Sorting BAM file..."
    samtools sort \
        -@ ${THREADS} \
        -o "${SORTED_BAM}" \
        "${RAW_BAM}"
    
    # Remove raw BAM to save space
    rm "${RAW_BAM}"
    
    # STEP 2C: Mark duplicates (using samtools)
    echo "Marking duplicates..."
    samtools markdup \
        -@ ${THREADS} \
        -s \
        "${SORTED_BAM}" \
        "${FINAL_BAM}" \
        2> "${LOG_DIR}/${SAMPLE}_markdup.log"
    
    # Index final BAM
    echo "Indexing final BAM..."
    samtools index -@ ${THREADS} "${FINAL_BAM}"
    
    # Remove sorted BAM to save space
    rm "${SORTED_BAM}"
    
    # STEP 2D: Generate alignment statistics
    echo "Generating alignment statistics..."
    samtools flagstat "${FINAL_BAM}" > "${STATS_DIR}/${SAMPLE}_flagstat.txt"
    samtools stats "${FINAL_BAM}" > "${STATS_DIR}/${SAMPLE}_stats.txt"
    samtools idxstats "${FINAL_BAM}" > "${STATS_DIR}/${SAMPLE}_idxstats.txt"
    
    echo "Sample ${SAMPLE} processing complete!"
    echo "End: $(date)"
    echo ""
    
done

################################################################################
# STEP 3: SUMMARY STATISTICS
################################################################################

echo "========================================"
echo "STEP 3: Generating summary statistics"
echo "========================================"
echo ""

# Create summary file
SUMMARY="${STATS_DIR}/alignment_summary.txt"
echo "Moringa Hi-C Alignment Summary" > "${SUMMARY}"
echo "Generated: $(date)" >> "${SUMMARY}"
echo "========================================" >> "${SUMMARY}"
echo "" >> "${SUMMARY}"

for SAMPLE in "${SAMPLES[@]}"; do
    FLAGSTAT="${STATS_DIR}/${SAMPLE}_flagstat.txt"
    
    if [ -f "${FLAGSTAT}" ]; then
        echo "Sample: ${SAMPLE}" >> "${SUMMARY}"
        echo "----------------------------------------" >> "${SUMMARY}"
        cat "${FLAGSTAT}" >> "${SUMMARY}"
        echo "" >> "${SUMMARY}"
    fi
done

echo "Summary statistics saved to: ${SUMMARY}"
echo ""

################################################################################
# STEP 4: QUICK COVERAGE CALCULATION
################################################################################

echo "========================================"
echo "STEP 4: Calculating coverage statistics"
echo "========================================"
echo ""

COVERAGE_SUMMARY="${STATS_DIR}/coverage_summary.txt"
echo "Sample,Total_Reads,Mapped_Reads,Properly_Paired,Mapping_Rate,Estimated_Coverage" > "${COVERAGE_SUMMARY}"

GENOME_SIZE=315000000  # 315 Mb

for SAMPLE in "${SAMPLES[@]}"; do
    FLAGSTAT="${STATS_DIR}/${SAMPLE}_flagstat.txt"
    
    if [ -f "${FLAGSTAT}" ]; then
        # Extract key statistics
        TOTAL=$(grep "in total" "${FLAGSTAT}" | awk '{print $1}')
        MAPPED=$(grep "mapped (" "${FLAGSTAT}" | head -1 | awk '{print $1}')
        PAIRED=$(grep "properly paired" "${FLAGSTAT}" | awk '{print $1}')
        
        # Calculate mapping rate
        MAPPING_RATE=$(awk "BEGIN {printf \"%.2f\", ($MAPPED/$TOTAL)*100}")
        
        # Estimate coverage (assuming 150bp reads)
        READ_LENGTH=150
        COVERAGE=$(awk "BEGIN {printf \"%.2f\", ($MAPPED*$READ_LENGTH)/$GENOME_SIZE}")
        
        echo "${SAMPLE},${TOTAL},${MAPPED},${PAIRED},${MAPPING_RATE}%,${COVERAGE}x" >> "${COVERAGE_SUMMARY}"
        
        echo "Sample ${SAMPLE}: ${MAPPING_RATE}% mapped, ~${COVERAGE}x coverage"
    fi
done

echo ""
echo "Coverage summary saved to: ${COVERAGE_SUMMARY}"
echo ""

################################################################################
# COMPLETION
################################################################################

echo "========================================"
echo "PIPELINE COMPLETE!"
echo "========================================"
echo "End time: $(date)"
echo ""
echo "Results locations:"
echo "  - BAM files: ${OUTPUT_DIR}"
echo "  - Statistics: ${STATS_DIR}"
echo "  - Logs: ${LOG_DIR}"
echo ""
echo "Summary files:"
echo "  - Alignment summary: ${SUMMARY}"
echo "  - Coverage summary: ${COVERAGE_SUMMARY}"
echo ""
echo "Next steps:"
echo "  1. Review alignment statistics in ${STATS_DIR}"
echo "  2. Check coverage levels"
echo "  3. These BAM files will be useful for variant calling later"
echo "  4. Wait for HiFi data to arrive for de novo assembly"
echo ""
