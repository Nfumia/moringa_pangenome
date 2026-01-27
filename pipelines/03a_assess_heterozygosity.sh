#!/bin/bash

# Heterozygosity Assessment for All Samples
# Extracts peak_hom and peak_het from hifiasm logs

BASE_DIR="$HOME/moringa_pangenome"
OUTPUT_FILE="$BASE_DIR/heterozygosity_summary.txt"

echo "========================================"
echo "Heterozygosity Assessment"
echo "========================================"
echo "Analyzing hifiasm logs for all 9 samples..."
echo ""

# Create header
echo -e "Sample\tpeak_hom\tpeak_het\tStatus\tRecommendation" > "$OUTPUT_FILE"

# Sample list
SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Find the hifiasm log with specific naming pattern
    LOG="$BASE_DIR/logs/hifiasm_${SAMPLE}_*.log"
    LOG_FILE=$(ls $LOG 2>/dev/null | head -1)
    
    if [ -f "$LOG_FILE" ]; then
        # Extract heterozygosity metrics
        PEAK_HOM=$(grep "peak_hom:" "$LOG_FILE" | head -1 | sed 's/.*peak_hom: \([0-9]*\);.*/\1/')
        PEAK_HET=$(grep "peak_het:" "$LOG_FILE" | head -1 | sed 's/.*peak_het: \([-0-9]*\).*/\1/')
        
        if [ -z "$PEAK_HOM" ]; then
            PEAK_HOM="N/A"
            PEAK_HET="N/A"
            STATUS="Log parsing failed"
            RECOMMENDATION="Check log manually"
        elif [ "$PEAK_HET" = "-1" ]; then
            STATUS="Homozygous"
            RECOMMENDATION="Hi-C scaffolding only (YaHS)"
        else
            STATUS="Heterozygous"
            RECOMMENDATION="Hi-C integrated phasing (hifiasm --h1 --h2)"
        fi
        
        echo -e "$SAMPLE\t$PEAK_HOM\t$PEAK_HET\t$STATUS\t$RECOMMENDATION"
    else
        echo -e "$SAMPLE\tN/A\tN/A\tLog not found\tWait for assembly"
    fi
done | tee -a "$OUTPUT_FILE"

echo ""
echo "========================================"
echo "Summary"
echo "========================================"

# Count categories
HOMO_COUNT=$(grep "Homozygous" "$OUTPUT_FILE" | wc -l)
HETERO_COUNT=$(grep "Heterozygous" "$OUTPUT_FILE" | wc -l)
PENDING_COUNT=$(grep -E "not found|N/A|parsing failed" "$OUTPUT_FILE" | wc -l)

echo "Homozygous samples: $HOMO_COUNT"
echo "Heterozygous samples: $HETERO_COUNT"
echo "Pending/Failed: $PENDING_COUNT"
echo ""

if [ $HOMO_COUNT -gt 0 ]; then
    echo "=== HOMOZYGOUS SAMPLES (Hi-C scaffolding only) ==="
    grep "Homozygous" "$OUTPUT_FILE" | cut -f1
    echo ""
fi

if [ $HETERO_COUNT -gt 0 ]; then
    echo "=== HETEROZYGOUS SAMPLES (Hi-C phasing recommended) ==="
    grep "Heterozygous" "$OUTPUT_FILE" | cut -f1
    echo ""
fi

echo "Full report saved to: $OUTPUT_FILE"
echo "========================================"
