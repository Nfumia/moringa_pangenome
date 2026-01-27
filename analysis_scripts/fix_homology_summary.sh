#!/bin/bash

echo "Recalculating homology summary with correct sizes..."

OUTDIR="../haplotype_analysis/batch_results/02_homology"
mkdir -p "$OUTDIR"

# New results file
echo "Sample,Hap1_Size_Mb,Hap2_Size_Mb,Size_Diff_Mb,High_Homology,Partial_Homology,Low_Homology,Collapsed_Size_Mb" > "$OUTDIR/homology_results_fixed.csv"

SAMPLES=(Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."
    
    HAP1="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
    HAP2="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
    SAMPLE_DIR="../haplotype_analysis/${SAMPLE}/02_homology"
    
    if [ ! -f "$HAP1" ] || [ ! -f "$HAP2" ]; then
        echo "  ERROR: Files not found"
        continue
    fi
    
    # Get sizes properly
    HAP1_SIZE=$(seqkit stats "$HAP1" 2>&1 | grep -v "^file" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    HAP2_SIZE=$(seqkit stats "$HAP2" 2>&1 | grep -v "^file" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    SIZE_DIFF=$(awk -v h1="$HAP1_SIZE" -v h2="$HAP2_SIZE" 'BEGIN {printf "%.1f", h1-h2}')
    
    echo "  Hap1: ${HAP1_SIZE} Mb, Hap2: ${HAP2_SIZE} Mb, Diff: ${SIZE_DIFF} Mb"
    
    # Get homology counts from existing analysis
    if [ -f "$SAMPLE_DIR/hap1_best_matches.txt" ]; then
        HIGH=$(awk '$3 > 95 && $4 > 80' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
        PARTIAL=$(awk '$3 > 90 && $4 >= 50 && $4 <= 80' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
        LOW=$(awk '$4 < 50' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
        
        # Get collapsed size
        if [ -s "$SAMPLE_DIR/collapsed_contigs.txt" ]; then
            COLLAPSED=$(seqkit grep -f "$SAMPLE_DIR/collapsed_contigs.txt" "$HAP1" 2>&1 | \
                seqkit stats 2>&1 | grep -v "^file" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        else
            COLLAPSED=0
        fi
    else
        HIGH=0
        PARTIAL=0
        LOW=0
        COLLAPSED=0
    fi
    
    echo "$SAMPLE,$HAP1_SIZE,$HAP2_SIZE,$SIZE_DIFF,$HIGH,$PARTIAL,$LOW,$COLLAPSED" >> "$OUTDIR/homology_results_fixed.csv"
done

echo ""
echo "Fixed results saved to: $OUTDIR/homology_results_fixed.csv"
cat "$OUTDIR/homology_results_fixed.csv"

