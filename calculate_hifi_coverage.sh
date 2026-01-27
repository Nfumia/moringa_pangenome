#!/bin/bash

# Calculate HiFi Coverage for All Samples
# Coverage = (Total bases sequenced) / (Genome size)

echo "========================================"
echo "Calculating HiFi Coverage for All Samples"
echo "========================================"
echo ""

OUTPUT="hifi_coverage_summary.txt"

{
echo "Sample Coverage Calculations"
echo "Date: $(date)"
echo ""
printf "%-10s | %-15s | %-12s | %-10s\n" "Sample" "Total Bases" "Genome Size" "Coverage"
echo "----------------------------------------------------------------"

# Sample list
SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing: $SAMPLE"
    
    # Correct path to HiFi reads
    HIFI_READS="raw_data/hifi/${SAMPLE}_hifi.fastq.gz"
    
    if [ -f "$HIFI_READS" ]; then
        # Calculate total bases from HiFi reads
        echo "  Calculating bases from: $HIFI_READS"
        
        # Using seqkit
        TOTAL_BASES=$(seqkit stats "$HIFI_READS" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); print $5}')
        
        # Get genome size from assembly
        ASSEMBLY="assemblies/production/${SAMPLE}/${SAMPLE}.primary.fasta"
        GENOME_SIZE=$(seqkit stats "$ASSEMBLY" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); print $5}')
        
        # Calculate coverage
        if [ -n "$TOTAL_BASES" ] && [ -n "$GENOME_SIZE" ] && [ "$GENOME_SIZE" -gt 0 ]; then
            COVERAGE=$(awk -v tb="$TOTAL_BASES" -v gs="$GENOME_SIZE" 'BEGIN {printf "%.1f", tb/gs}')
            
            # Convert to human readable
            TOTAL_GB=$(awk -v tb="$TOTAL_BASES" 'BEGIN {printf "%.2f", tb/1000000000}')
            GENOME_MB=$(awk -v gs="$GENOME_SIZE" 'BEGIN {printf "%.1f", gs/1000000}')
            
            printf "%-10s | %-15s | %-12s | %-10s\n" "$SAMPLE" "${TOTAL_GB} Gb" "${GENOME_MB} Mb" "${COVERAGE}x"
        else
            printf "%-10s | %-15s | %-12s | %-10s\n" "$SAMPLE" "ERROR" "ERROR" "ERROR"
            echo "  [ERROR] Could not calculate coverage"
        fi
    else
        printf "%-10s | %-15s | %-12s | %-10s\n" "$SAMPLE" "NOT FOUND" "NOT FOUND" "NOT FOUND"
        echo "  [NOT FOUND] HiFi reads not found: $HIFI_READS"
    fi
    echo ""
done

echo "================================================================"
echo ""
echo "Coverage Interpretation:"
echo "  20-30x:  Minimum for good assembly"
echo "  30-40x:  Good coverage (typical target)"
echo "  40-60x:  Excellent coverage"
echo "  >60x:    Very high coverage (diminishing returns)"
echo ""
echo "Mean Coverage Summary:"
MEAN=$(awk '/[0-9]+\.[0-9]+x/ {sum+=$NF; count++} END {if(count>0) printf "%.1fx", sum/count; else print "N/A"}' <<< "$(grep -E '[0-9]+\.[0-9]+x' hifi_coverage_summary.txt 2>/dev/null || echo '')")
echo "  Mean across all samples: $MEAN"
echo ""
echo "================================================================"
} | tee "$OUTPUT"

echo ""
echo "Results saved to: $OUTPUT"

