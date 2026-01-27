#!/bin/bash

# Quality Assessment of Hi-C Scaffolded Assemblies

echo "========================================"
echo "Hi-C Scaffolding Quality Assessment"
echo "========================================"
echo ""

SAMPLES=(Mo-TH-30 Mo-TH-66 Mo-US-5)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "========================================="
    echo "Sample: $SAMPLE"
    echo "========================================="
    
    PRIMARY="assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
    SCAFFOLDED="assemblies/scaffolded/$SAMPLE/${SAMPLE}_scaffolds_scaffolds_final.fa"
    
    if [ ! -f "$SCAFFOLDED" ]; then
        echo "✗ Scaffolded assembly not found"
        echo ""
        continue
    fi
    
    echo "--- BEFORE SCAFFOLDING (HiFi-only) ---"
    seqkit stats "$PRIMARY"
    
    echo ""
    echo "--- AFTER SCAFFOLDING (Hi-C) ---"
    seqkit stats "$SCAFFOLDED"
    
    echo ""
    echo "--- IMPROVEMENT SUMMARY ---"
    
    # Get key metrics
    BEFORE_N50=$(seqkit stats "$PRIMARY" | tail -1 | awk '{print $16}')
    AFTER_N50=$(seqkit stats "$SCAFFOLDED" | tail -1 | awk '{print $16}')
    BEFORE_CONTIGS=$(seqkit stats "$PRIMARY" | tail -1 | awk '{print $4}')
    AFTER_SCAFFOLDS=$(seqkit stats "$SCAFFOLDED" | tail -1 | awk '{print $4}')
    BEFORE_SIZE=$(seqkit stats "$PRIMARY" | tail -1 | awk '{print $5}')
    AFTER_SIZE=$(seqkit stats "$SCAFFOLDED" | tail -1 | awk '{print $5}')
    
    echo "  Before: $BEFORE_CONTIGS contigs, N50 = $(echo "scale=2; $BEFORE_N50/1000000" | bc) Mb"
    echo "  After:  $AFTER_SCAFFOLDS scaffolds, N50 = $(echo "scale=2; $AFTER_N50/1000000" | bc) Mb"
    
    if [ "$BEFORE_N50" -gt 0 ]; then
        N50_FOLD=$(echo "scale=2; $AFTER_N50 / $BEFORE_N50" | bc)
        echo "  N50 improvement: ${N50_FOLD}x"
    fi
    
    REDUCTION=$(echo "scale=1; (1 - $AFTER_SCAFFOLDS / $BEFORE_CONTIGS) * 100" | bc)
    echo "  Sequence count reduced: ${REDUCTION}%"
    
    echo ""
done

echo "========================================"
echo "Summary"
echo "========================================"
echo "All scaffolded assemblies should show:"
echo "  ✓ Higher N50 (ideally >20 Mb for chromosome-scale)"
echo "  ✓ Fewer sequences (contigs joined into scaffolds)"
echo "  ✓ Similar total size (no sequence loss)"
echo "========================================"
