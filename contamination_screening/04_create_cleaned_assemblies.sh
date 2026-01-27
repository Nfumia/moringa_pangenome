#!/bin/bash

# Create Cleaned Assemblies
# Remove true contaminants, save organelles separately

set -e

BASE_DIR="$HOME/moringa_pangenome"
CONTAM_DIR="$BASE_DIR/contamination_screening"
RESULTS_DIR="$CONTAM_DIR/results"

# Create output directories
mkdir -p "$BASE_DIR/assemblies/cleaned"
mkdir -p "$BASE_DIR/assemblies/organelles"

echo "========================================"
echo "Creating Cleaned Assemblies"
echo "========================================"
echo ""
echo "Strategy:"
echo "  1. Extract organelles → assemblies/organelles/"
echo "  2. Remove true contaminants → assemblies/cleaned/"
echo "  3. Generate before/after statistics"
echo ""
echo "========================================"

SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

TOTAL_PROCESSED=0
TOTAL_ORGANELLES=0
TOTAL_CONTAMINANTS=0

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "========================================="
    echo "Processing: $SAMPLE"
    echo "========================================="
    
    ORIGINAL="$BASE_DIR/assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
    CLEANED="$BASE_DIR/assemblies/cleaned/${SAMPLE}.cleaned.fasta"
    ORGANELLES="$BASE_DIR/assemblies/organelles/${SAMPLE}.organelles.fasta"
    CHLORO="$BASE_DIR/assemblies/organelles/${SAMPLE}.chloroplast.fasta"
    MITO="$BASE_DIR/assemblies/organelles/${SAMPLE}.mitochondria.fasta"
    
    if [ ! -f "$ORIGINAL" ]; then
        echo "  ✗ Original assembly not found, skipping"
        continue
    fi
    
    # Check if screening was done
    ORGANELLE_LIST="$RESULTS_DIR/${SAMPLE}_organelle_ids.txt"
    CONTAM_LIST="$RESULTS_DIR/${SAMPLE}_contaminant_ids.txt"
    
    if [ ! -f "$ORGANELLE_LIST" ]; then
        echo "  ⚠ Screening results not found"
        echo "  Please run: bash 02_blast_screen_assemblies.sh"
        continue
    fi
    
    TOTAL_PROCESSED=$((TOTAL_PROCESSED + 1))
    
    # Get counts
    BEFORE=$(grep -c "^>" "$ORIGINAL")
    ORGANELLE_COUNT=$(wc -l < "$ORGANELLE_LIST" 2>/dev/null || echo "0")
    CONTAM_COUNT=$(wc -l < "$CONTAM_LIST" 2>/dev/null || echo "0")
    
    echo "  Original assembly: $BEFORE contigs"
    echo "  Organelles found: $ORGANELLE_COUNT contigs"
    echo "  True contaminants: $CONTAM_COUNT contigs"
    
    # 1. Extract organelles (chloroplast + mitochondria combined)
    if [ -s "$ORGANELLE_LIST" ]; then
        seqkit grep -f "$ORGANELLE_LIST" "$ORIGINAL" > "$ORGANELLES"
        ORGANELLE_EXTRACTED=$(grep -c "^>" "$ORGANELLES")
        echo ""
        echo "  Step 1: Organelles extracted"
        echo "    → $ORGANELLE_EXTRACTED contigs saved to:"
        echo "      ${SAMPLE}.organelles.fasta"
        
        TOTAL_ORGANELLES=$((TOTAL_ORGANELLES + ORGANELLE_EXTRACTED))
        
        # Separate chloroplast and mitochondria
        CHLORO_LIST="$RESULTS_DIR/${SAMPLE}_chloroplast_ids.txt"
        MITO_LIST="$RESULTS_DIR/${SAMPLE}_mitochondria_ids.txt"
        
        if [ -s "$CHLORO_LIST" ]; then
            seqkit grep -f "$CHLORO_LIST" "$ORIGINAL" > "$CHLORO"
            CHLORO_COUNT=$(grep -c "^>" "$CHLORO")
            echo "    - Chloroplast: $CHLORO_COUNT contigs"
        fi
        
        if [ -s "$MITO_LIST" ]; then
            seqkit grep -f "$MITO_LIST" "$ORIGINAL" > "$MITO"
            MITO_COUNT=$(grep -c "^>" "$MITO")
            echo "    - Mitochondria: $MITO_COUNT contigs"
        fi
    else
        echo ""
        echo "  Step 1: No organelles found"
        touch "$ORGANELLES"
    fi
    
    # 2. Create cleaned assembly (remove true contaminants only)
    echo ""
    echo "  Step 2: Creating cleaned nuclear assembly"
    
    if [ -s "$CONTAM_LIST" ]; then
        seqkit grep -v -f "$CONTAM_LIST" "$ORIGINAL" > "$CLEANED"
        AFTER=$(grep -c "^>" "$CLEANED")
        REMOVED=$((BEFORE - AFTER))
        
        echo "    → Removed $REMOVED contaminant contigs"
        echo "    → Cleaned assembly: $AFTER contigs"
        
        TOTAL_CONTAMINANTS=$((TOTAL_CONTAMINANTS + REMOVED))
    else
        echo "    → No contaminants found"
        cp "$ORIGINAL" "$CLEANED"
        AFTER=$BEFORE
    fi
    
    # 3. Generate statistics
    echo ""
    echo "  Step 3: Generating statistics"
    
    # Size comparison
    BEFORE_SIZE=$(seqkit stats "$ORIGINAL" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    AFTER_SIZE=$(seqkit stats "$CLEANED" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    
    if [ -s "$ORGANELLES" ]; then
        ORG_SIZE=$(seqkit stats "$ORGANELLES" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    else
        ORG_SIZE="0.0"
    fi
    
    PERCENT_REMOVED=$(awk -v r="$REMOVED" -v b="$BEFORE" 'BEGIN {printf "%.2f", (r/b)*100}')
    
    echo "    Before: $BEFORE contigs, ${BEFORE_SIZE} Mb"
    echo "    After:  $AFTER contigs, ${AFTER_SIZE} Mb"
    echo "    Organelles: ${ORG_SIZE} Mb"
    echo "    Removed: ${PERCENT_REMOVED}% of contigs"
    
    # Save detailed summary
    {
        echo "Sample: $SAMPLE"
        echo "Date: $(date)"
        echo ""
        echo "ORIGINAL ASSEMBLY:"
        echo "  File: $ORIGINAL"
        echo "  Contigs: $BEFORE"
        echo "  Size: ${BEFORE_SIZE} Mb"
        echo ""
        echo "CONTAMINATION DETECTED:"
        echo "  Organelles: $ORGANELLE_COUNT contigs"
        echo "  True contaminants: $CONTAM_COUNT contigs"
        echo "  Total removed: $REMOVED contigs (${PERCENT_REMOVED}%)"
        echo ""
        echo "CLEANED NUCLEAR ASSEMBLY:"
        echo "  File: $CLEANED"
        echo "  Contigs: $AFTER"
        echo "  Size: ${AFTER_SIZE} Mb"
        echo ""
        echo "ORGANELLE ASSEMBLY:"
        echo "  File: $ORGANELLES"
        echo "  Size: ${ORG_SIZE} Mb"
        if [ -s "$CHLORO_LIST" ]; then
            echo "  Chloroplast contigs: $CHLORO_COUNT"
        fi
        if [ -s "$MITO_LIST" ]; then
            echo "  Mitochondrial contigs: $MITO_COUNT"
        fi
        echo ""
    } > "$BASE_DIR/assemblies/cleaned/${SAMPLE}.cleaning_report.txt"
    
    echo ""
    echo "  ✓ Complete!"
    echo "    Cleaned: assemblies/cleaned/${SAMPLE}.cleaned.fasta"
    echo "    Organelles: assemblies/organelles/${SAMPLE}.organelles.fasta"
    echo "    Report: assemblies/cleaned/${SAMPLE}.cleaning_report.txt"
    
done

# Final summary
echo ""
echo "========================================"
echo "CLEANING COMPLETE - SUMMARY"
echo "========================================"
echo ""
echo "Samples processed: $TOTAL_PROCESSED/9"
echo "Total organelle contigs extracted: $TOTAL_ORGANELLES"
echo "Total contaminant contigs removed: $TOTAL_CONTAMINANTS"
echo ""
echo "Output locations:"
echo "  Cleaned assemblies: $BASE_DIR/assemblies/cleaned/"
echo "  Organelle assemblies: $BASE_DIR/assemblies/organelles/"
echo "  Cleaning reports: $BASE_DIR/assemblies/cleaned/*_cleaning_report.txt"
echo ""
echo "========================================"
echo "RECOMMENDATIONS"
echo "========================================"
echo ""
echo "CLEANED ASSEMBLIES:"
echo "  Use these for:"
echo "    - Pangenome construction"
echo "    - Gene annotation"
echo "    - Comparative genomics"
echo "    - All nuclear genome analyses"
echo ""
echo "ORGANELLE ASSEMBLIES:"
echo "  Use these for:"
echo "    - Chloroplast genome assembly/analysis"
echo "    - Mitochondrial genome assembly/analysis"
echo "    - Phylogenetic studies"
echo "    - Organellar gene annotation"
echo ""
echo "NEXT STEPS:"
echo "  1. Review cleaning reports"
echo "  2. Run quality control on cleaned assemblies"
echo "  3. Proceed with pangenome construction using cleaned assemblies"
echo ""
echo "========================================"
