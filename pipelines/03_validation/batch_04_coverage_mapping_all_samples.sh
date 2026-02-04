#!/bin/bash

################################################################################
# Coverage Mapping - ALL HETEROZYGOUS SAMPLES
# Purpose: Identify collapsed regions via read depth
# Runtime: 6-12 hours total (1-2 hours per sample)
################################################################################

echo "======================================================================="
echo "COVERAGE MAPPING - ALL HETEROZYGOUS SAMPLES"
echo "======================================================================="
echo "Date: $(date)"
echo ""
echo "WARNING: This will take 6-12 hours to complete!"
echo "Consider running overnight or in a screen session"
echo ""

# All heterozygous samples
SAMPLES=(Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13)

OUTDIR="../haplotype_analysis/batch_results/04_coverage"
mkdir -p "$OUTDIR"

# Results file
echo "Sample,Hap1_Cov,Hap2_Cov,High_Cov_Contigs,Collapsed_Size_Mb" > "$OUTDIR/coverage_results.csv"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "======================================================================="
    echo "Processing: $SAMPLE (this will take 1-2 hours)"
    echo "======================================================================="
    
    SAMPLE_DIR="../haplotype_analysis/${SAMPLE}/04_coverage"
    mkdir -p "$SAMPLE_DIR"
    
    HAP1="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
    HAP2="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
    READS="../raw_data/hifi/${SAMPLE}_hifi.fastq.gz"
    
    if [ ! -f "$HAP1" ] || [ ! -f "$HAP2" ] || [ ! -f "$READS" ]; then
        echo "ERROR: Files not found for $SAMPLE"
        continue
    fi
    
    echo "Step 1/5: Indexing haplotypes..."
    minimap2 -d "$SAMPLE_DIR/hap1.mmi" "$HAP1"
    minimap2 -d "$SAMPLE_DIR/hap2.mmi" "$HAP2"
    
    echo "Step 2/5: Mapping reads to Hap1..."
    minimap2 -ax map-hifi -t 24 "$SAMPLE_DIR/hap1.mmi" "$READS" | \
        samtools sort -@ 8 -o "$SAMPLE_DIR/hap1.bam"
    samtools index "$SAMPLE_DIR/hap1.bam"
    
    echo "Step 3/5: Mapping reads to Hap2..."
    minimap2 -ax map-hifi -t 24 "$SAMPLE_DIR/hap2.mmi" "$READS" | \
        samtools sort -@ 8 -o "$SAMPLE_DIR/hap2.bam"
    samtools index "$SAMPLE_DIR/hap2.bam"
    
    echo "Step 4/5: Calculating coverage..."
    samtools depth "$SAMPLE_DIR/hap1.bam" > "$SAMPLE_DIR/hap1_depth.txt"
    samtools depth "$SAMPLE_DIR/hap2.bam" > "$SAMPLE_DIR/hap2_depth.txt"
    
    echo "Step 5/5: Analyzing coverage patterns..."
    
    # Mean coverage
    H1_COV=$(awk '{sum+=$3; count++} END {printf "%.1f", sum/count}' "$SAMPLE_DIR/hap1_depth.txt")
    H2_COV=$(awk '{sum+=$3; count++} END {printf "%.1f", sum/count}' "$SAMPLE_DIR/hap2_depth.txt")
    
    # Find high coverage regions (>1.5x mean)
    awk -v mean="$H1_COV" '$3 > mean*1.5' "$SAMPLE_DIR/hap1_depth.txt" | \
        awk '{print $1}' | sort -u > "$SAMPLE_DIR/high_cov_contigs.txt"
    
    HIGH_COUNT=$(wc -l < "$SAMPLE_DIR/high_cov_contigs.txt")
    
    if [ "$HIGH_COUNT" -gt 0 ]; then
        COLLAPSED=$(seqkit grep -f "$SAMPLE_DIR/high_cov_contigs.txt" "$HAP1" 2>/dev/null | \
            seqkit stats 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    else
        COLLAPSED=0
    fi
    
    echo "$SAMPLE,$H1_COV,$H2_COV,$HIGH_COUNT,$COLLAPSED" >> "$OUTDIR/coverage_results.csv"
    
    echo "  Hap1 coverage: ${H1_COV}x"
    echo "  High-coverage contigs: $HIGH_COUNT (${COLLAPSED} Mb)"
done

echo ""
echo "======================================================================="
echo "GENERATING SUMMARY"
echo "======================================================================="

{
echo "Coverage Mapping Summary - All Heterozygous Samples"
echo "Date: $(date)"
echo ""
echo "======================================================================="
echo ""
printf "%-10s | %-10s | %-10s | %-15s | %-15s\n" "Sample" "Hap1 Cov" "Hap2 Cov" "High-Cov Ctgs" "Conserved (Mb)"
echo "--------------------------------------------------------------------------------"

while IFS=',' read -r sample h1cov h2cov highcount collapsed; do
    if [ "$sample" != "Sample" ]; then
        printf "%-10s | %-10s | %-10s | %-15s | %-15s\n" "$sample" "${h1cov}x" "${h2cov}x" "$highcount" "$collapsed"
    fi
done < "$OUTDIR/coverage_results.csv"

echo ""
echo "======================================================================="
echo "INTERPRETATION:"
echo "======================================================================="
echo ""
echo "High-coverage regions (>1.5x mean) in Hap1 indicate collapsed regions"
echo "High-coverage regions (>1.5x mean) indicate conserved sequences between haplotypes"
echo "  → These sequences are highly similar (>95% identity) between Hap1 and Hap2"
echo "  → Reads from BOTH haplotypes map here (NOT collapsed - properly phased)"
} | tee "$OUTDIR/coverage_summary_all_samples.txt"

echo ""
echo "======================================================================="
echo "All samples complete!"
echo "Summary: $OUTDIR/coverage_summary_all_samples.txt"
echo "======================================================================="
