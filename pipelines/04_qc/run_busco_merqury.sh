#!/bin/bash
#####################################################################
# RUN BUSCO AND MERQURY ONLY
# Steps 1-3 already complete, this runs Steps 4-5
#####################################################################

set -e

WORK_DIR="/home/hawaii-agriculture-research-cent/moringa_pangenome"
RESULTS_DIR="${WORK_DIR}/results/qc"
ASSEMBLY_DIR="${WORK_DIR}/assemblies"
RAW_DATA="${WORK_DIR}/raw_data/hifi"

BUSCO_SUMMARY="${RESULTS_DIR}/04_busco/busco_summary.txt"
MERQURY_SUMMARY="${RESULTS_DIR}/05_merqury/merqury_summary.txt"

# Get all assemblies
ASSEMBLIES=()
for SAMPLE in Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13; do
    ASSEMBLIES+=("${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap1.fasta")
    ASSEMBLIES+=("${ASSEMBLY_DIR}/phased/${SAMPLE}/${SAMPLE}.hap2.fasta")
done
for SAMPLE in Mo-TH-30 Mo-TH-66 Mo-US-5; do
    ASSEMBLIES+=("${ASSEMBLY_DIR}/scaffolded/${SAMPLE}/${SAMPLE}_scaffolds_scaffolds_final.fa")
done

#####################################################################
# STEP 4: BUSCO ANALYSIS
#####################################################################

echo "======================================================================="
echo "STEP 4/5: BUSCO - GENE COMPLETENESS ASSESSMENT"
echo "Started: $(date)"
echo "======================================================================="

# Activate BUSCO environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate busco_env

LINEAGE="eudicots_odb10"

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
    
    # Run BUSCO
    busco -i ${ASSEMBLY} \
          -o ${NAME} \
          --out_path ${RESULTS_DIR}/04_busco \
          -l ${LINEAGE} \
          -m genome \
          -c 24 \
          -f 2>&1 | tee ${RESULTS_DIR}/logs/${NAME}_busco.log
    
    # Parse results
    SHORT_SUMMARY=$(find ${RESULTS_DIR}/04_busco/${NAME} -name "short_summary*.txt" 2>/dev/null | head -1)
    if [ -f "${SHORT_SUMMARY}" ]; then
        COMPLETE=$(grep "C:" ${SHORT_SUMMARY} | grep -oP "C:\K[0-9.]+" | head -1)
        SINGLE=$(grep "S:" ${SHORT_SUMMARY} | grep -oP "S:\K[0-9.]+" | head -1)
        DUPLICATED=$(grep "D:" ${SHORT_SUMMARY} | grep -oP "D:\K[0-9.]+" | head -1)
        FRAGMENTED=$(grep "F:" ${SHORT_SUMMARY} | grep -oP "F:\K[0-9.]+" | head -1)
        MISSING=$(grep "M:" ${SHORT_SUMMARY} | grep -oP "M:\K[0-9.]+" | head -1)
        
        printf "%-25s | %-10s | %-10s | %-10s | %-10s | %-10s\n" \
            "${NAME}" "${COMPLETE}%" "${SINGLE}%" "${DUPLICATED}%" "${FRAGMENTED}%" "${MISSING}%" >> ${BUSCO_SUMMARY}
    fi
done

echo ""
echo "Step 4 complete! Summary: ${BUSCO_SUMMARY}"

#####################################################################
# STEP 5: MERQURY ANALYSIS
#####################################################################

echo "======================================================================="
echo "STEP 5/5: MERQURY - K-MER BASED QC"
echo "Started: $(date)"
echo "======================================================================="

# Switch back to hifi_assembly for merqury
conda activate hifi_assembly

# Header
{
echo "Merqury Analysis Summary"
echo "Generated: $(date)"
echo ""
printf "%-25s | %-12s | %-15s | %-12s\n" "Assembly" "QV" "Completeness" "Type"
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
            
            printf "%-25s | %-12s | %-15s | %-12s\n" "${SAMPLE}" "${QV}" "${COMPLETENESS}%" "diploid" >> ${MERQURY_SUMMARY}
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
            
            printf "%-25s | %-12s | %-15s | %-12s\n" "${SAMPLE}" "${QV}" "${COMPLETENESS}%" "haploid" >> ${MERQURY_SUMMARY}
        fi
    fi
    
    cd ${WORK_DIR}
done

echo ""
echo "Step 5 complete! Summary: ${MERQURY_SUMMARY}"

echo "======================================================================="
echo "QC PIPELINE COMPLETE!"
echo "Finished: $(date)"
echo "======================================================================="
