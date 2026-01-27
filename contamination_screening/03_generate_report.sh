#!/bin/bash

# Generate Comprehensive Contamination Screening Report

BASE_DIR="$HOME/moringa_pangenome"
CONTAM_DIR="$BASE_DIR/contamination_screening"
RESULTS_DIR="$CONTAM_DIR/results"
REPORT="$CONTAM_DIR/contamination_report.txt"

echo "========================================"
echo "Generating Contamination Report"
echo "========================================"

{
echo "========================================"
echo "CONTAMINATION SCREENING REPORT"
echo "========================================"
echo "Generated: $(date)"
echo "Threshold: >80% query coverage"
echo "Based on: Dr. Jingjing Zheng's workflow section I.3.2"
echo ""
echo "Reference databases used:"
echo "  • 6 Chloroplast genomes (5 Moringa + 1 Arabidopsis)"
echo "  • 7 Mitochondrial references (2 Carica papaya genomes, 1 Brassica napus genome, 1 Raphanus sativus genome, 1 Arabidopsis genome, 2 Moringa oleifera partial genes)"
echo "  • Bacterial contaminants (E. coli, Agrobacterium)"
echo "  • Vector sequences (UniVec)"
echo ""
echo "========================================"
echo "RESULTS BY SAMPLE"
echo "========================================"
echo ""

printf "%-10s | %-6s | %-12s | %-15s | %-12s | %-10s | %-10s\n" \
    "Sample" "Contigs" "Chloroplast" "Mitochondria" "Bacteria/Vec" "Total" "% Contam"
echo "------------------------------------------------------------------------------------------------------------"

SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

TOTAL_SAMPLES=0
TOTAL_WITH_CHLORO=0
TOTAL_WITH_MITO=0
TOTAL_WITH_CONTAM=0

for SAMPLE in "${SAMPLES[@]}"; do
    SUMMARY="$RESULTS_DIR/${SAMPLE}_summary.txt"
    
    if [ -f "$SUMMARY" ]; then
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        
        # Extract values from summary
        TOTAL_CONTIGS=$(grep "Total contigs:" "$SUMMARY" | awk '{print $3}')
        CHLORO=$(grep "Chloroplast contigs:" "$SUMMARY" | awk '{print $3}')
        MITO=$(grep "Mitochondrial contigs:" "$SUMMARY" | awk '{print $3}')
        CONTAM=$(grep "Bacterial/Vector contigs:" "$SUMMARY" | awk '{print $3}')
        TOTAL=$(grep "Total contamination:" "$SUMMARY" | awk '{print $3}')
        PERCENT=$(grep "Total contamination:" "$SUMMARY" | grep -oP '\(\K[0-9.]+')
        
        # Count samples with each type
        [ "$CHLORO" -gt 0 ] && TOTAL_WITH_CHLORO=$((TOTAL_WITH_CHLORO + 1))
        [ "$MITO" -gt 0 ] && TOTAL_WITH_MITO=$((TOTAL_WITH_MITO + 1))
        [ "$CONTAM" -gt 0 ] && TOTAL_WITH_CONTAM=$((TOTAL_WITH_CONTAM + 1))
        
        printf "%-10s | %-6s | %-12s | %-15s | %-12s | %-10s | %-10s\n" \
            "$SAMPLE" "$TOTAL_CONTIGS" "$CHLORO" "$MITO" "$CONTAM" "$TOTAL" "${PERCENT}%"
    else
        printf "%-10s | %-6s | %-12s | %-15s | %-12s | %-10s | %-10s\n" \
            "$SAMPLE" "N/A" "Not screened" "Not screened" "Not screened" "-" "-"
    fi
done

echo ""
echo "========================================"
echo "SUMMARY STATISTICS"
echo "========================================"
echo ""
echo "Samples screened: $TOTAL_SAMPLES/9"
echo "Samples with chloroplast contamination: $TOTAL_WITH_CHLORO"
echo "Samples with mitochondrial contamination: $TOTAL_WITH_MITO"
echo "Samples with bacterial/vector contamination: $TOTAL_WITH_CONTAM"
echo ""
echo "========================================"
echo "INTERPRETATION & RECOMMENDATIONS"
echo "========================================"
echo ""
echo "ORGANELLAR SEQUENCES (Chloroplast/Mitochondria):"
echo "  Status: Expected and normal in plant assemblies"
echo "  Action: Save separately for organelle genome assembly"
echo "  Expected amounts:"
echo "    - Chloroplast: 1-2 contigs, ~150-165 kb total"
echo "    - Mitochondria: 1-5 contigs, ~200-500 kb total"
echo ""
echo "TRUE CONTAMINANTS (Bacteria/Vectors):"
echo "  Status: Should be minimal (<1% is acceptable)"
echo "  Action: Remove from nuclear genome assembly"
echo "  Concern levels:"
echo "    - <1%: Acceptable contamination level"
echo "    - 1-5%: Review recommended, possibly acceptable"
echo "    - >5%: Significant contamination, requires investigation"
echo ""
echo "NEXT STEPS:"
echo "  1. Review individual sample summaries in: $RESULTS_DIR"
echo "  2. Examine specific contigs if contamination >1%"
echo "  3. Run: bash 04_create_cleaned_assemblies.sh"
echo "     to create cleaned nuclear assemblies"
echo ""
echo "========================================"
echo "DETAILED SAMPLE SUMMARIES"
echo "========================================"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    SUMMARY="$RESULTS_DIR/${SAMPLE}_summary.txt"
    if [ -f "$SUMMARY" ]; then
        echo "----------------------------------------"
        cat "$SUMMARY"
        echo ""
    fi
done

echo "========================================"
echo "END OF REPORT"
echo "========================================"

} | tee "$REPORT"

echo ""
echo "========================================"
echo "Report Generation Complete"
echo "========================================"
echo ""
echo "✓ Full report saved to: $REPORT"
echo ""
echo "Quick summary:"
cat "$REPORT" | grep -A 10 "RESULTS BY SAMPLE"
echo ""
echo "To view full report:"
echo "  cat $REPORT"
echo "  less $REPORT"
echo ""
echo "========================================"
