#!/bin/bash

################################################################################
# Gap Estimation Analysis - ALL HETEROZYGOUS SAMPLES
# Purpose: Check for gaps in all 6 phased samples
# Runtime: <10 minutes total
################################################################################

echo "======================================================================="
echo "GAP ESTIMATION - ALL HETEROZYGOUS SAMPLES"
echo "======================================================================="
echo "Date: $(date)"
echo ""

# All heterozygous samples
SAMPLES=(Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13)

OUTDIR="../haplotype_analysis/batch_results/01_gaps"
mkdir -p "$OUTDIR"

# Summary file
SUMMARY="${OUTDIR}/gap_summary_all_samples.txt"

{
echo "Gap Analysis Summary - All Heterozygous Samples"
echo "Date: $(date)"
echo ""
echo "======================================================================="
echo ""
printf "%-10s | %-15s | %-15s | %-10s\n" "Sample" "Hap1 Gaps" "Hap2 Gaps" "Status"
echo "------------------------------------------------------------------------"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing: $SAMPLE"
    
    HAP1="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
    HAP2="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
    
    if [ -f "$HAP1" ] && [ -f "$HAP2" ]; then
        # Count Ns in Hap1
        N_H1=$(grep -v "^>" "$HAP1" | tr -d '\n' | grep -o "N" | wc -l)
        
        # Count Ns in Hap2
        N_H2=$(grep -v "^>" "$HAP2" | tr -d '\n' | grep -o "N" | wc -l)
        
        # Status
        if [ "$N_H1" -eq 0 ] && [ "$N_H2" -eq 0 ]; then
            STATUS="Gap-free ✓"
        else
            STATUS="Has gaps !"
        fi
        
        printf "%-10s | %-15s | %-15s | %-10s\n" "$SAMPLE" "$N_H1 Ns" "$N_H2 Ns" "$STATUS"
    else
        printf "%-10s | %-15s | %-15s | %-10s\n" "$SAMPLE" "NOT FOUND" "NOT FOUND" "ERROR"
    fi
done

echo ""
echo "======================================================================="
echo "CONCLUSION:"
echo "======================================================================="
echo ""
echo "All samples should be gap-free (0 Ns) for high-quality HiFi assemblies"

} | tee "$SUMMARY"

echo ""
echo "======================================================================="
echo "Summary saved to: $SUMMARY"
echo "======================================================================="
