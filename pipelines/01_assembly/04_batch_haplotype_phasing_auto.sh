#!/bin/bash

# Batch Hifiasm Hi-C Phasing (AUTO - No Prompts)
# Processes all heterozygous samples automatically

# set -e  # Disabled

BASE_DIR="$HOME/moringa_pangenome"
HETERO_FILE="$BASE_DIR/heterozygosity_summary.txt"

if [ ! -f "$HETERO_FILE" ]; then
    echo "ERROR: Run 03a_assess_heterozygosity.sh first!"
    exit 1
fi

echo "========================================"
echo "Batch Haplotype Phasing (AUTO)"
echo "========================================"
echo "Processing heterozygous samples automatically..."
echo ""

# Extract heterozygous samples
HETERO_SAMPLES=$(grep "Heterozygous" "$HETERO_FILE" | cut -f1)

if [ -z "$HETERO_SAMPLES" ]; then
    echo "No heterozygous samples found!"
    exit 0
fi

echo "Heterozygous samples to phase:"
echo "$HETERO_SAMPLES"
echo ""

# Counters
TOTAL=$(echo "$HETERO_SAMPLES" | wc -l)
SUCCESS=0
FAILED=0
CURRENT=0

# Process each sample
for SAMPLE in $HETERO_SAMPLES; do
    ((CURRENT++))
    echo "========================================"
    echo "Processing: $SAMPLE ($CURRENT/$TOTAL)"
    echo "========================================"
    
    # Run phasing script
    if bash "$BASE_DIR/pipelines/04_hifiasm_hic_phasing.sh" "$SAMPLE"; then
        echo "✓ $SAMPLE phasing successful"
        ((SUCCESS++))
    else
        echo "✗ $SAMPLE phasing failed"
        ((FAILED++))
    fi
    
    echo ""
done

# Summary
echo "========================================"
echo "Batch Phasing Complete"
echo "========================================"
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Failed: $FAILED"
echo "========================================"
