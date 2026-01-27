#!/bin/bash

# Comprehensive Assembly Summary - All Types

OUTPUT="assembly_summary_report.txt"

{
echo "================================================================================"
echo "                    MORINGA PANGENOME ASSEMBLY SUMMARY"
echo "================================================================================"
echo "Generated: $(date)"
echo ""
echo "========================================"
echo "1. HiFi-Only Assemblies (9 samples)"
echo "========================================"
printf "%-10s | %-10s | %-8s | %-10s | %-12s | %-8s\n" \
    "Sample" "Size(Mb)" "Contigs" "N50(Mb)" "Longest(Mb)" "Het_Level"
echo "--------------------------------------------------------------------------------"

for SAMPLE in Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5; do
    FASTA="assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
    if [ -f "$FASTA" ]; then
        STATS=$(seqkit stats -a "$FASTA" 2>/dev/null | tail -1)
        
        # Field 5=sum_len, 4=num_seqs, 8=max_len, 13=N50
        SIZE=$(echo "$STATS" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        CONTIGS=$(echo "$STATS" | awk '{gsub(/,/,"",$4); print $4}')
        LONGEST=$(echo "$STATS" | awk '{gsub(/,/,"",$8); printf "%.2f", $8/1000000}')
        N50=$(echo "$STATS" | awk '{gsub(/,/,"",$13); printf "%.2f", $13/1000000}')
        
        HET=$(grep "^$SAMPLE" heterozygosity_summary.txt 2>/dev/null | head -1 | awk '{print $3}')
        
        printf "%-10s | %-10s | %-8s | %-10s | %-12s | %-8s\n" \
            "$SAMPLE" "$SIZE" "$CONTIGS" "$N50" "$LONGEST" "$HET"
    fi
done

echo ""
echo "========================================"
echo "2. Hi-C Scaffolded Assemblies (3 homozygous)"
echo "========================================"
printf "%-10s | %-12s | %-10s | %-12s | %-14s | %-12s\n" \
    "Sample" "Size(Mb)" "Scaffolds" "N50(Mb)" "Longest(Mb)" "Improvement"
echo "--------------------------------------------------------------------------------"

for SAMPLE in Mo-TH-30 Mo-TH-66 Mo-US-5; do
    SCAFF="assemblies/scaffolded/$SAMPLE/${SAMPLE}_scaffolds_scaffolds_final.fa"
    if [ -f "$SCAFF" ]; then
        STATS=$(seqkit stats -a "$SCAFF" 2>/dev/null | tail -1)
        
        SIZE=$(echo "$STATS" | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        SCAFFOLDS=$(echo "$STATS" | awk '{gsub(/,/,"",$4); print $4}')
        LONGEST=$(echo "$STATS" | awk '{gsub(/,/,"",$8); printf "%.2f", $8/1000000}')
        N50=$(echo "$STATS" | awk '{gsub(/,/,"",$13); printf "%.2f", $13/1000000}')
        
        # Calculate improvement
        ORIG="assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
        BEFORE_N50=$(seqkit stats -a "$ORIG" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$13); print $13}')
        AFTER_N50=$(echo "$STATS" | awk '{gsub(/,/,"",$13); print $13}')
        
        if [ -n "$BEFORE_N50" ] && [ "$BEFORE_N50" != "0" ]; then
            FOLD=$(awk -v a="$AFTER_N50" -v b="$BEFORE_N50" 'BEGIN {printf "%.2f", a/b}')
        else
            FOLD="N/A"
        fi
        
        printf "%-10s | %-12s | %-10s | %-12s | %-14s | %-12s\n" \
            "$SAMPLE" "$SIZE" "$SCAFFOLDS" "$N50" "$LONGEST" "${FOLD}x N50"
    fi
done

echo ""
echo "========================================"
echo "3. Haplotype-Phased Assemblies (6 heterozygous)"
echo "========================================"
printf "%-10s | %-10s | %-10s | %-10s | %-12s | %-10s\n" \
    "Sample" "Hap1(Mb)" "Hap2(Mb)" "Total(Mb)" "Ratio" "Status"
echo "--------------------------------------------------------------------------------"

for SAMPLE in Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13; do
    HAP1="assemblies/phased/$SAMPLE/${SAMPLE}.hap1.fasta"
    HAP2="assemblies/phased/$SAMPLE/${SAMPLE}.hap2.fasta"
    
    if [ -f "$HAP1" ] && [ -f "$HAP2" ]; then
        HAP1_SIZE=$(seqkit stats "$HAP1" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        HAP2_SIZE=$(seqkit stats "$HAP2" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
        TOTAL=$(awk -v a="$HAP1_SIZE" -v b="$HAP2_SIZE" 'BEGIN {printf "%.1f", a+b}')
        
        ORIG="assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
        ORIG_SIZE=$(seqkit stats "$ORIG" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); print $5}')
        
        if [ -n "$ORIG_SIZE" ] && [ "$ORIG_SIZE" != "0" ]; then
            RATIO=$(awk -v t="$TOTAL" -v o="$ORIG_SIZE" 'BEGIN {printf "%.2f", (t*1000000)/o}')
        else
            RATIO="N/A"
        fi
        
        printf "%-10s | %-10s | %-10s | %-10s | %-12s | %-10s\n" \
            "$SAMPLE" "$HAP1_SIZE" "$HAP2_SIZE" "$TOTAL" "${RATIO}x" "✓ Complete"
    else
        printf "%-10s | %-10s | %-10s | %-10s | %-12s | %-10s\n" \
            "$SAMPLE" "-" "-" "-" "-" "⏳ Pending"
    fi
done

echo ""
echo "========================================"
echo "SUMMARY STATISTICS"
echo "========================================"
echo "Total assemblies when complete: 15"
echo "  - HiFi-only: 9 assemblies"
echo "  - Scaffolded (homozygous): 3 assemblies"
echo "  - Phased (heterozygous): 6 samples × 2 haplotypes = 12 assemblies"
echo ""
echo "Heterozygosity distribution:"
echo "  - Homozygous (peak_het = -1): 3 samples (33%)"
echo "  - Heterozygous (peak_het > 0): 6 samples (67%)"
echo ""
echo "Size range:"
echo "  - HiFi assemblies: 290-314 Mb (N50: 8.9-15.3 Mb)"
echo "  - Scaffolded: 310-312 Mb (N50: 16.5-21.6 Mb, longest: 29-42 Mb)"
echo "  - Phased total: 529-573 Mb per sample (diploid, 1.76-1.83x)"
echo "================================================================================"
} | tee "$OUTPUT"

echo ""
echo "✓ Summary saved to: $OUTPUT"
