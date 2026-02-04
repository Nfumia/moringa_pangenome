#!/bin/bash

# Batch GFA to FASTA Conversion
# Converts primary and alternate GFA files to FASTA for all samples

# set -e  # Disabled

BASE_DIR="$HOME/moringa_pangenome"
ASSEMBLY_DIR="$BASE_DIR/assemblies/production"

echo "========================================"
echo "Batch GFA to FASTA Conversion"
echo "========================================"
echo "Start time: $(date)"
echo ""

# Sample list
SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

# Counters
TOTAL=${#SAMPLES[@]}
SUCCESS=0
SKIPPED=0
FAILED=0

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "Processing: $SAMPLE ($((SUCCESS + SKIPPED + FAILED + 1))/$TOTAL)"
    echo "----------------------------------------"
    
    SAMPLE_DIR="$ASSEMBLY_DIR/$SAMPLE"
    
    # Check if assembly directory exists
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "⚠ Assembly directory not found: $SAMPLE_DIR"
        echo "Skipping $SAMPLE..."
        ((SKIPPED++))
        echo ""
        continue
    fi
    
    cd "$SAMPLE_DIR"
    
    # Define GFA and FASTA files
    PRIMARY_GFA="${SAMPLE}.p_ctg.gfa"
    ALTERNATE_GFA="${SAMPLE}.a_ctg.gfa"
    PRIMARY_FASTA="${SAMPLE}.primary.fasta"
    ALTERNATE_FASTA="${SAMPLE}.alternate.fasta"
    
    # Check if GFA files exist
    if [ ! -f "$PRIMARY_GFA" ]; then
        echo "✗ Primary GFA not found: $PRIMARY_GFA"
        ((FAILED++))
        echo ""
        continue
    fi
    
    # Check if FASTA already exists
    if [ -f "$PRIMARY_FASTA" ] && [ -f "$ALTERNATE_FASTA" ]; then
        echo "✓ FASTA files already exist, skipping..."
        ((SKIPPED++))
        echo ""
        continue
    fi
    
    # Convert primary GFA to FASTA
    echo "Converting primary GFA to FASTA..."
    if awk '/^S/{print ">"$2"\n"$3}' "$PRIMARY_GFA" > "$PRIMARY_FASTA"; then
        echo "✓ Primary FASTA created: $PRIMARY_FASTA"
    else
        echo "✗ Failed to convert primary GFA"
        ((FAILED++))
        echo ""
        continue
    fi
    
    # Convert alternate GFA to FASTA (if exists)
    if [ -f "$ALTERNATE_GFA" ]; then
        echo "Converting alternate GFA to FASTA..."
        if awk '/^S/{print ">"$2"\n"$3}' "$ALTERNATE_GFA" > "$ALTERNATE_FASTA"; then
            echo "✓ Alternate FASTA created: $ALTERNATE_FASTA"
        else
            echo "⚠ Failed to convert alternate GFA (non-critical)"
        fi
    fi
    
    # Verify FASTA files
    if [ -f "$PRIMARY_FASTA" ] && [ -s "$PRIMARY_FASTA" ]; then
        echo ""
        echo "=== Statistics ==="
        seqkit stats -a "$PRIMARY_FASTA" 2>/dev/null || echo "seqkit not available for stats"
        ((SUCCESS++))
    else
        echo "✗ Primary FASTA is empty or missing"
        ((FAILED++))
    fi
    
    echo ""
done

# Summary
echo "========================================"
echo "Batch Conversion Summary"
echo "========================================"
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Skipped (already exist or no assembly): $SKIPPED"
echo "Failed: $FAILED"
echo "End time: $(date)"
echo "========================================"

# Exit with error if any conversions failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi
