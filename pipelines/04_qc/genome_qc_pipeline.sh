#!/bin/bash

#####################################################################
# GENOME QC PIPELINE - Comprehensive Quality Assessment
# Runs all assemblies through 5 QC methods in optimal order:
#   1. Assembly Statistics (seqkit)
#   2. Q30 Scores (seqkit)  
#   3. Coverage Analysis (mosdepth)
#   4. BUSCO (gene completeness)
#   5. Merqury (k-mer based QC)
#
# Usage: bash genome_qc_pipeline.sh
# Runtime: ~40-60 hours for 15 assemblies
#
# REQUIREMENTS:
#   - hifi_assembly environment (Steps 1-3, 5)
#   - busco_env environment (Step 4) - BUSCO 5.7.1 requires Python 3.11
#     Create with: conda create -n busco_env -c bioconda -c conda-forge busco=5.7.1 python=3.11
#####################################################################

set -e  # Exit on error

# Directories
WORK_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome"
RESULTS_DIR="${WORK_DIR}/results/qc"
ASSEMBLY_DIR="${WORK_DIR}/assemblies"
RAW_DATA="${WORK_DIR}/raw_data/hifi"

# Create output directories
mkdir -p ${RESULTS_DIR}/{01_stats,02_q30,03_coverage,04_busco,05_merqury,logs}

# Output summary files
STATS_SUMMARY="${RESULTS_DIR}/01_stats/assembly_statistics_summary.txt"
Q30_SUMMARY="${RESULTS_DIR}/02_q30/q30_scores_summary.txt"
COV_SUMMARY="${RESULTS_DIR}/03_coverage/coverage_summary.txt"
BUSCO_SUMMARY="${RESULTS_DIR}/04_busco/busco_summary.txt"
MERQURY_SUMMARY="${RESULTS_DIR}/05_merqury/merqury_summary.txt"

# Get all assemblies
ASSEMBLIES=()

# Phased haplotypes
for SAMPLE in Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13; do
    ASSEMBLIES+=("${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap1.fasta")
    ASSEMBLIES+=("${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap2.fasta")
done

# Scaffolded assemblies
for SAMPLE in Mo-TH-30 Mo-TH-66 Mo-US-5; do
    ASSEMBLIES+=("${ASSEMBLY_DIR}/scaffolded/${SAMPLE}/${SAMPLE}_scaffolds_scaffolds_final.fa")
done

echo "======================================================================="
echo "GENOME QC PIPELINE - Starting $(date)"
echo "======================================================================="
echo "Total assemblies to process: ${#ASSEMBLIES[@]}"
echo ""

#####################################################################
# STEP 1: ASSEMBLY STATISTICS
#####################################################################

echo "======================================================================="
echo "STEP 1/5: ASSEMBLY STATISTICS"
echo "Started: $(date)"
echo "======================================================================="

# Header
{
echo "Assembly Statistics Summary"
echo "Generated: $(date)"
echo ""
printf "%-25s | %-12s | %-10s | %-10s | %-8s | %-12s | %-8s\n" \
    "Assembly" "Total Size" "Contigs" "N50" "L50" "Longest" "GC%"
echo "------------------------------------------------------------------------------------------------------------------------"
} > ${STATS_SUMMARY}

for ASSEMBLY in "${ASSEMBLIES[@]}"; do
    NAME=$(basename ${ASSEMBLY}); NAME=${NAME%.fasta}; NAME=${NAME%.fa}; NAME=${NAME%_scaffolds_scaffolds_final}
    echo "  Processing: ${NAME}"
    
    # Run seqkit stats
    seqkit stats -a ${ASSEMBLY} &> ${RESULTS_DIR}/logs/${NAME}_stats.log
    
    # Extract key metrics
    TOTAL_SIZE=$(seqkit stats ${ASSEMBLY} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f Mb", $5/1000000}')
    NUM_SEQS=$(seqkit stats ${ASSEMBLY} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$4); print $4}')
    N50=$(seqkit stats ${ASSEMBLY} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$13); printf "%.2f Mb", $13/1000000}')
    L50=$(seqkit stats ${ASSEMBLY} 2>/dev/null | tail -1 | awk '{print $14}')
    MAX_LEN=$(seqkit stats ${ASSEMBLY} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$7); printf "%.2f Mb", $7/1000000}')
    GC=$(seqkit fx2tab -n -g ${ASSEMBLY} 2>/dev/null | awk '{sum+=$2; count++} END {printf "%.2f%%", sum/count}')
    
    printf "%-25s | %-12s | %-10s | %-10s | %-8s | %-12s | %-8s\n" \
        "${NAME}" "${TOTAL_SIZE}" "${NUM_SEQS}" "${N50}" "${L50}" "${MAX_LEN}" "${GC}" >> ${STATS_SUMMARY}
done

echo ""
echo "Step 1 complete! Summary: ${STATS_SUMMARY}"
echo ""

#####################################################################
# STEP 2: Q30 SCORES
#####################################################################

echo "======================================================================="
echo "STEP 2/5: Q30 SCORES FROM RAW READS"
echo "Started: $(date)"
echo "======================================================================="

# Header
{
echo "Q30 Scores Summary"
echo "Generated: $(date)"
echo ""
printf "%-20s | %-15s | %-15s | %-10s\n" "Sample" "Total Reads" "Total Bases" "Q30%"
echo "----------------------------------------------------------------"
} > ${Q30_SUMMARY}

# Get unique sample names
SAMPLES=(Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13 Mo-TH-30 Mo-TH-66 Mo-US-5)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Processing: ${SAMPLE}"
    
    READS="${RAW_DATA}/${SAMPLE}_hifi.fastq.gz"
    
    if [ -f "${READS}" ]; then
        seqkit stats -a ${READS} &> ${RESULTS_DIR}/logs/${SAMPLE}_q30.log
        
        TOTAL_READS=$(seqkit stats ${READS} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$4); printf "%d", $4}')
        TOTAL_BASES=$(seqkit stats ${READS} 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.2f Gb", $5/1000000000}')
        Q30=$(seqkit stats ${READS} 2>/dev/null | tail -1 | awk '{print $10}')
        
        printf "%-20s | %-15s | %-15s | %-10s\n" "${SAMPLE}" "${TOTAL_READS}" "${TOTAL_BASES}" "${Q30}" >> ${Q30_SUMMARY}
    else
        echo "    Warning: Reads not found for ${SAMPLE}"
        printf "%-20s | %-15s | %-15s | %-10s\n" "${SAMPLE}" "N/A" "N/A" "N/A" >> ${Q30_SUMMARY}
    fi
done

echo ""
echo "Step 2 complete! Summary: ${Q30_SUMMARY}"
echo ""

#####################################################################
# STEP 3: COVERAGE ANALYSIS
#####################################################################

echo "======================================================================="
echo "STEP 3/5: COVERAGE ANALYSIS WITH MOSDEPTH"
echo "Started: $(date)"
echo "======================================================================="

# Check if mosdepth is installed
if ! command -v mosdepth &> /dev/null; then
    echo "Installing mosdepth..."
    conda install -c bioconda mosdepth -y
fi

# Header
{
echo "Coverage Analysis Summary"
echo "Generated: $(date)"
echo ""
printf "%-25s | %-12s | %-12s | %-12s\n" "Assembly" "Mean Coverage" "Median Cov" "Stddev"
echo "-----------------------------------------------------------------------------------------"
} > ${COV_SUMMARY}

for ASSEMBLY in "${ASSEMBLIES[@]}"; do
    NAME=$(basename ${ASSEMBLY}); NAME=${NAME%.fasta}; NAME=${NAME%.fa}; NAME=${NAME%_scaffolds_scaffolds_final}
    SAMPLE=${NAME%.hap*}
    
    echo "  Processing: ${NAME}"
    
    # Find corresponding reads
    READS="${RAW_DATA}/${SAMPLE}_hifi.fastq.gz"
    
    if [ -f "${READS}" ]; then
        OUTDIR="${RESULTS_DIR}/03_coverage/${NAME}"
        mkdir -p ${OUTDIR}
        
        # Map reads to assembly
        echo "    Mapping reads..."
        minimap2 -ax map-hifi -t 24 ${ASSEMBLY} ${READS} 2>${RESULTS_DIR}/logs/${NAME}_mapping.log | \
            samtools sort -@ 8 -o ${OUTDIR}/${NAME}.bam 2>>${RESULTS_DIR}/logs/${NAME}_mapping.log
        
        samtools index ${OUTDIR}/${NAME}.bam
        
        # Run mosdepth
        echo "    Calculating coverage..."
        mosdepth -t 8 ${OUTDIR}/${NAME} ${OUTDIR}/${NAME}.bam 2>&1 | tee ${RESULTS_DIR}/logs/${NAME}_mosdepth.log
        
        # Extract metrics
        if [ -f "${OUTDIR}/${NAME}.mosdepth.summary.txt" ]; then
            MEAN=$(grep "total" ${OUTDIR}/${NAME}.mosdepth.summary.txt | awk '{printf "%.1fx", $4}')
            MEDIAN=$(grep "total" ${OUTDIR}/${NAME}.mosdepth.summary.txt | awk '{printf "%.1fx", $5}' || echo "N/A")
            STDDEV=$(grep "total" ${OUTDIR}/${NAME}.mosdepth.summary.txt | awk '{printf "%.1f", $6}' || echo "N/A")
            
            printf "%-25s | %-12s | %-12s | %-12s\n" "${NAME}" "${MEAN}" "${MEDIAN}" "${STDDEV}" >> ${COV_SUMMARY}
        fi
        
        # Clean up large BAM files
        rm -f ${OUTDIR}/${NAME}.bam ${OUTDIR}/${NAME}.bam.bai
        
    else
        echo "    Warning: Reads not found for ${SAMPLE}"
        printf "%-25s | %-12s | %-12s | %-12s\n" "${NAME}" "N/A" "N/A" "N/A" >> ${COV_SUMMARY}
    fi
done

echo ""
echo "Step 3 complete! Summary: ${COV_SUMMARY}"
echo ""

#####################################################################
# STEP 4: BUSCO ANALYSIS
#####################################################################

echo "======================================================================="
echo "STEP 4/5: BUSCO - GENE COMPLETENESS ASSESSMENT"
echo "Started: $(date)"
echo "======================================================================="

# Activate BUSCO environment (requires separate env due to Python 3.11 dependency)
source $(conda info --base)/etc/profile.d/conda.sh
if ! conda activate busco_env 2>/dev/null; then
    echo "ERROR: busco_env not found!"
    echo "Create it with: conda create -n busco_env -c bioconda -c conda-forge busco=5.7.1 python=3.11 -y"
    echo "Then rerun this script."
    exit 1
fi

# Download eudicots lineage if not present
LINEAGE="eudicots_odb10"
BUSCO_DOWNLOADS="${HOME}/.busco_downloads"

if [ ! -d "${BUSCO_DOWNLOADS}/lineages/${LINEAGE}" ]; then
    echo "Downloading BUSCO lineage: ${LINEAGE}..."
    busco --download ${LINEAGE}
fi

# Header
{
echo "BUSCO Analysis Summary"
echo "Lineage: ${LINEAGE}"
echo "Generated: $(date)"
echo ""
printf "%-25s | %-10s | %-10s | %-10s | %-10s | %-10s\n" "Assembly" "Complete" "Single" "Duplicated" "Fragmented" "Missing"
echo "--------------------------------------------------------------------------------------------------------------------"
} > ${BUSCO_SUMMARY}

for ASSEMBLY in "${ASSEMBLIES[@]}"; do
    NAME=$(basename ${ASSEMBLY}); NAME=${NAME%.fasta}; NAME=${NAME%.fa}; NAME=${NAME%_scaffolds_scaffolds_final}
    echo "  Processing: ${NAME} (this will take 2-4 hours)"
    
    OUTDIR="${RESULTS_DIR}/04_busco/${NAME}"
    
    # Run BUSCO
    busco -i ${ASSEMBLY} \
          -o ${NAME} \
          --out_path ${RESULTS_DIR}/04_busco \
          -l ${LINEAGE} \
          -m genome \
          -c 24 \
          --offline \
          -f 2>&1 | tee ${RESULTS_DIR}/logs/${NAME}_busco.log
    
    # Parse results
    if [ -f "${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt" ]; then
        COMPLETE=$(grep "C:" ${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt | head -1 | awk -F'C:' '{print $2}' | awk -F'%' '{print $1}')
        SINGLE=$(grep "S:" ${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt | head -1 | awk -F'S:' '{print $2}' | awk -F'%' '{print $1}')
        DUPLICATED=$(grep "D:" ${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt | head -1 | awk -F'D:' '{print $2}' | awk -F'%' '{print $1}')
        FRAGMENTED=$(grep "F:" ${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt | head -1 | awk -F'F:' '{print $2}' | awk -F'%' '{print $1}')
        MISSING=$(grep "M:" ${OUTDIR}/short_summary.specific.${LINEAGE}.${NAME}.txt | head -1 | awk -F'M:' '{print $2}' | awk -F'%' '{print $1}')
        
        printf "%-25s | %-10s | %-10s | %-10s | %-10s | %-10s\n" \
            "${NAME}" "${COMPLETE}%" "${SINGLE}%" "${DUPLICATED}%" "${FRAGMENTED}%" "${MISSING}%" >> ${BUSCO_SUMMARY}
    fi
done

echo ""
echo "Step 4 complete! Summary: ${BUSCO_SUMMARY}"
echo ""

#####################################################################
# STEP 5: MERQURY ANALYSIS
#####################################################################

echo "======================================================================="
echo "STEP 5/5: MERQURY - K-MER BASED QC"
echo "Started: $(date)"
echo "======================================================================="

# Check if merqury is available
if ! command -v merqury.sh &> /dev/null; then
    echo "WARNING: merqury.sh not found in PATH"
    echo "Skipping Merqury analysis..."
    echo "To enable Merqury, install it and ensure merqury.sh is in PATH"
else
    # Header
    {
    echo "Merqury Analysis Summary"
    echo "Generated: $(date)"
    echo ""
    printf "%-25s | %-12s | %-12s | %-12s\n" "Assembly" "QV" "Completeness" "Type"
    echo "-----------------------------------------------------------------------"
    } > ${MERQURY_SUMMARY}
    
    # Process phased samples (diploid mode)
    for SAMPLE in Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13; do
        echo "  Processing: ${SAMPLE} (diploid mode)"
        
        OUTDIR="${RESULTS_DIR}/05_merqury/${SAMPLE}"
        mkdir -p ${OUTDIR}
        cd ${OUTDIR}
        
        READS="${RAW_DATA}/${SAMPLE}_hifi.fastq.gz"
        HAP1="${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
        HAP2="${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
        
        if [ -f "${READS}" ] && [ -f "${HAP1}" ] && [ -f "${HAP2}" ]; then
            # Build k-mer database if not exists
            if [ ! -d "${SAMPLE}_reads.meryl" ]; then
                echo "    Building k-mer database..."
                meryl count k=21 memory=100 threads=24 output ${SAMPLE}_reads.meryl ${READS} 2>&1 | tee ${RESULTS_DIR}/logs/${SAMPLE}_meryl.log
            fi
            
            # Run merqury
            echo "    Running Merqury..."
            merqury.sh ${SAMPLE}_reads.meryl ${HAP1} ${HAP2} ${SAMPLE} 2>&1 | tee ${RESULTS_DIR}/logs/${SAMPLE}_merqury.log
            
            # Parse results
            if [ -f "${SAMPLE}.qv" ]; then
                QV=$(tail -1 ${SAMPLE}.qv | awk '{print $4}')
                COMPLETENESS=$(grep "both" ${SAMPLE}.completeness.stats 2>/dev/null | awk '{print $5}' || echo "N/A")
                
                printf "%-25s | %-12s | %-12s | %-12s\n" "${SAMPLE}" "${QV}" "${COMPLETENESS}%" "diploid" >> ${MERQURY_SUMMARY}
            fi
        fi
        
        cd ${WORK_DIR}
    done
    
    # Process scaffolded samples (haploid mode)
    for SAMPLE in Mo-TH-30 Mo-TH-66 Mo-US-5; do
        echo "  Processing: ${SAMPLE} (haploid mode)"
        
        OUTDIR="${RESULTS_DIR}/05_merqury/${SAMPLE}"
        mkdir -p ${OUTDIR}
        cd ${OUTDIR}
        
        READS="${RAW_DATA}/${SAMPLE}_hifi.fastq.gz"
        ASSEMBLY="${ASSEMBLY_DIR}/scaffolded/${SAMPLE}/${SAMPLE}_scaffolds_scaffolds_final.fa"
        
        if [ -f "${READS}" ] && [ -f "${ASSEMBLY}" ]; then
            # Build k-mer database if not exists
            if [ ! -d "${SAMPLE}_reads.meryl" ]; then
                echo "    Building k-mer database..."
                meryl count k=21 memory=100 threads=24 output ${SAMPLE}_reads.meryl ${READS} 2>&1 | tee ${RESULTS_DIR}/logs/${SAMPLE}_meryl.log
            fi
            
            # Run merqury
            echo "    Running Merqury..."
            merqury.sh ${SAMPLE}_reads.meryl ${ASSEMBLY} ${SAMPLE} 2>&1 | tee ${RESULTS_DIR}/logs/${SAMPLE}_merqury.log
            
            # Parse results
            if [ -f "${SAMPLE}.qv" ]; then
                QV=$(tail -1 ${SAMPLE}.qv | awk '{print $4}')
                COMPLETENESS=$(grep "all" ${SAMPLE}.completeness.stats 2>/dev/null | awk '{print $5}' || echo "N/A")
                
                printf "%-25s | %-12s | %-12s | %-12s\n" "${SAMPLE}" "${QV}" "${COMPLETENESS}%" "haploid" >> ${MERQURY_SUMMARY}
            fi
        fi
        
        cd ${WORK_DIR}
    done
    
    echo ""
    echo "Step 5 complete! Summary: ${MERQURY_SUMMARY}"
    echo ""
fi

#####################################################################
# FINAL SUMMARY
#####################################################################

echo "======================================================================="
echo "GENOME QC PIPELINE COMPLETE!"
echo "Finished: $(date)"
echo "======================================================================="
echo ""
echo "Summary files:"
echo "  1. Assembly Statistics: ${STATS_SUMMARY}"
echo "  2. Q30 Scores:          ${Q30_SUMMARY}"
echo "  3. Coverage Analysis:   ${COV_SUMMARY}"
echo "  4. BUSCO Results:       ${BUSCO_SUMMARY}"
echo "  5. Merqury Results:     ${MERQURY_SUMMARY}"
echo ""
echo "All logs saved to: ${RESULTS_DIR}/logs/"
echo "======================================================================="
