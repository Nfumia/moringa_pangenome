#!/bin/bash

echo "=========================================="
echo "RECALCULATING SYNTENY WITH CORRECT METHOD"
echo "=========================================="
echo ""

OUTDIR="../haplotype_analysis/batch_results/03_synteny"

echo "Sample,Hap1_Size_Mb,Hap2_Size_Mb,Aligned_Mb,Alignment_Pct" > "$OUTDIR/synteny_results_CORRECTED.csv"

for SAMPLE in Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13; do
    echo "Processing: $SAMPLE"
    
    HAP1="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
    HAP2="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
    PAF="../haplotype_analysis/${SAMPLE}/03_synteny/hap1_vs_hap2.paf"
    
    if [ -f "$PAF" ] && [ -f "$HAP1" ] && [ -f "$HAP2" ]; then
        # Get sizes
        HAP1_SIZE=$(seqkit stats "$HAP1" 2>&1 | grep -v "^file" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        HAP2_SIZE=$(seqkit stats "$HAP2" 2>&1 | grep -v "^file" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        
        # CORRECT METHOD: Calculate unique covered bases in Hap1
        # For each contig, take the max end position (only count once)
        ALIGNED=$(awk '{
            contig = $1
            start = $3
            end = $4
            if (end > max_end[contig]) {
                max_end[contig] = end
            }
        }
        END {
            total = 0
            for (contig in max_end) {
                total += max_end[contig]
            }
            printf "%.1f", total/1000000
        }' "$PAF")
        
        # Calculate percentage
        PCT=$(awk -v a="$ALIGNED" -v h="$HAP1_SIZE" 'BEGIN {printf "%.1f", (a/h)*100}')
        
        echo "  Hap1: ${HAP1_SIZE} Mb, Aligned: ${ALIGNED} Mb, Coverage: ${PCT}%"
        echo "$SAMPLE,$HAP1_SIZE,$HAP2_SIZE,$ALIGNED,$PCT" >> "$OUTDIR/synteny_results_CORRECTED.csv"
    else
        echo "  ERROR: Files missing"
    fi
done

echo ""
echo "=========================================="
echo "CORRECTED RESULTS:"
echo "=========================================="
cat "$OUTDIR/synteny_results_CORRECTED.csv" | column -t -s,

echo ""
echo "Saved to: $OUTDIR/synteny_results_CORRECTED.csv"

