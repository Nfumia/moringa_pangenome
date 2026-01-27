#!/bin/bash

# Batch YaHS Hi-C Scaffolding
# Processes all homozygous samples sequentially

# set -e  # Disabled

BASE_DIR="$HOME/moringa_pangenome"
HETERO_FILE="$BASE_DIR/heterozygosity_summary.txt"

if [ ! -f "$HETERO_FILE" ]; then
    echo "ERROR: Run 03a_assess_heterozygosity.sh first!"
    exit 1
fi

echo "========================================"
echo "Batch YaHS Hi-C Scaffolding"
echo "========================================"
echo "Processing homozygous samples only..."
echo ""

# Extract homozygous samples
HOMO_SAMPLES=$(grep "Homozygous" "$HETERO_FILE" | cut -f1)

if [ -z "$HOMO_SAMPLES" ]; then
    echo "No homozygous samples found!"
    exit 0
fi

echo "Homozygous samples to scaffold:"
echo "$HOMO_SAMPLES"
echo ""

# Counters
TOTAL=$(echo "$HOMO_SAMPLES" | wc -l)
SUCCESS=0
FAILED=0
CURRENT=0

# Process each sample
for SAMPLE in $HOMO_SAMPLES; do
    ((CURRENT++))
    echo "========================================"
    echo "Processing: $SAMPLE ($CURRENT/$TOTAL)"
    echo "========================================"
    
    # Run scaffolding script
    if bash "$BASE_DIR/pipelines/03_yahs_scaffolding.sh" "$SAMPLE"; then
        echo "✓ $SAMPLE scaffolding successful"
        ((SUCCESS++))
    else
        echo "✗ $SAMPLE scaffolding failed"
        ((FAILED++))
    fi
    
    echo ""
done

# Summary
echo "========================================"
echo "Batch Scaffolding Complete"
echo "========================================"
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Failed: $FAILED"
echo "========================================"
