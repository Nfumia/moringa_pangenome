#!/bin/bash

################################################################################
# Homology Analysis - ALL HETEROZYGOUS SAMPLES
# Purpose: Identify collapsed regions vs structural variants in all samples
# Runtime: 3-6 hours total (30-60 min per sample)
################################################################################

echo "======================================================================="
echo "HOMOLOGY ANALYSIS - ALL HETEROZYGOUS SAMPLES"
echo "======================================================================="
echo "Date: $(date)"
echo ""

# All heterozygous samples
SAMPLES=(Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TW-13)

OUTDIR="../haplotype_analysis/batch_results/02_homology"
mkdir -p "$OUTDIR"

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "======================================================================="
    echo "Processing: $SAMPLE"
    echo "======================================================================="
    
    SAMPLE_DIR="../haplotype_analysis/${SAMPLE}/02_homology"
    mkdir -p "$SAMPLE_DIR"
    
    HAP1="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap1.fasta"
    HAP2="../assemblies/phased/${SAMPLE}/${SAMPLE}.hap2.fasta"
    
    if [ ! -f "$HAP1" ] || [ ! -f "$HAP2" ]; then
        echo "ERROR: Files not found for $SAMPLE"
        continue
    fi
    
    echo "Step 1/5: Creating BLAST databases..."
    makeblastdb -in "$HAP1" -dbtype nucl -out "$SAMPLE_DIR/hap1_db" -logfile "$SAMPLE_DIR/makeblastdb_hap1.log"
    makeblastdb -in "$HAP2" -dbtype nucl -out "$SAMPLE_DIR/hap2_db" -logfile "$SAMPLE_DIR/makeblastdb_hap2.log"
    
    echo "Step 2/5: BLASTing Hap1 against Hap2..."
    blastn -query "$HAP1" \
           -db "$SAMPLE_DIR/hap2_db" \
           -out "$SAMPLE_DIR/hap1_vs_hap2.blast" \
           -outfmt "6 qseqid sseqid pident length qlen slen qcovs qstart qend sstart send evalue bitscore" \
           -evalue 1e-10 \
           -num_threads 24 \
           -max_target_seqs 5
    
    echo "Step 3/5: BLASTing Hap2 against Hap1..."
    blastn -query "$HAP2" \
           -db "$SAMPLE_DIR/hap1_db" \
           -out "$SAMPLE_DIR/hap2_vs_hap1.blast" \
           -outfmt "6 qseqid sseqid pident length qlen slen qcovs qstart qend sstart send evalue bitscore" \
           -evalue 1e-10 \
           -num_threads 24 \
           -max_target_seqs 5
    
    echo "Step 4/5: Analyzing results..."
    
    # Parse BLAST results
    awk '$7 > 50' "$SAMPLE_DIR/hap1_vs_hap2.blast" | \
        sort -k1,1 -k13,13nr | \
        awk '!seen[$1]++ {print $1"\t"$2"\t"$3"\t"$7"\t"$4}' \
        > "$SAMPLE_DIR/hap1_best_matches.txt"
    
    # Categorize
    HIGH=$(awk '$3 > 95 && $4 > 80' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
    PARTIAL=$(awk '$3 > 90 && $4 >= 50 && $4 <= 80' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
    LOW=$(awk '$4 < 50' "$SAMPLE_DIR/hap1_best_matches.txt" | wc -l)
    
    # Get collapsed regions
    awk '$4 < 50' "$SAMPLE_DIR/hap1_best_matches.txt" | cut -f1 > "$SAMPLE_DIR/collapsed_contigs.txt"
    
    if [ -s "$SAMPLE_DIR/collapsed_contigs.txt" ]; then
        seqkit grep -f "$SAMPLE_DIR/collapsed_contigs.txt" "$HAP1" > "$SAMPLE_DIR/collapsed_regions.fasta"
        COLLAPSED_SIZE=$(seqkit stats "$SAMPLE_DIR/collapsed_regions.fasta" 2>/dev/null | \
            tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    else
        COLLAPSED_SIZE=0
    fi
    
    # Get sizes
    HAP1_SIZE=$(seqkit stats "$HAP1" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    HAP2_SIZE=$(seqkit stats "$HAP2" 2>/dev/null | tail -1 | awk '{gsub(/,/,"",$5); printf "%.1f", $5/1000000}')
    SIZE_DIFF=$(awk -v h1="$HAP1_SIZE" -v h2="$HAP2_SIZE" 'BEGIN {printf "%.1f", h1-h2}')
    
    echo "$SAMPLE,$HAP1_SIZE,$HAP2_SIZE,$SIZE_DIFF,$HIGH,$PARTIAL,$LOW,$COLLAPSED_SIZE" >> "$OUTDIR/homology_results.csv"
    
    echo "Step 5/5: $SAMPLE complete"
    echo "  Collapsed regions: ${COLLAPSED_SIZE} Mb"
    echo "  Size difference: ${SIZE_DIFF} Mb"
done

echo ""
echo "======================================================================="
echo "GENERATING SUMMARY REPORT"
echo "======================================================================="

{
echo "Homology Analysis Summary - All Heterozygous Samples"
echo "Date: $(date)"
echo ""
echo "======================================================================="
echo ""
printf "%-10s | %-8s | %-8s | %-8s | %-12s | %-12s\n" "Sample" "Hap1(Mb)" "Hap2(Mb)" "Diff(Mb)" "Collapsed(Mb)" "% Collapsed"
echo "--------------------------------------------------------------------------------"

while IFS=',' read -r sample h1 h2 diff high partial low collapsed; do
    if [ "$diff" != "0" ] && [ "$diff" != "SIZE_DIFF" ]; then
        PCT=$(awk -v c="$collapsed" -v d="$diff" 'BEGIN {printf "%.1f", (c/d)*100}')
        printf "%-10s | %-8s | %-8s | %-8s | %-12s | %-12s\n" "$sample" "$h1" "$h2" "$diff" "$collapsed" "${PCT}%"
    fi
done < "$OUTDIR/homology_results.csv"

echo ""
echo "======================================================================="
echo "INTERPRETATION:"
echo "======================================================================="
echo ""
echo "Collapsed regions = sequences in Hap1 with no/low homology in Hap2"
echo "  → Should have been separated but collapsed due to low heterozygosity"
echo ""
echo "% Collapsed indicates what portion of size difference is technical vs biological"
echo "  >70% = Primarily phasing limitation (collapsed regions)"
echo "  <30% = Primarily biological (structural variants)"

} | tee "$OUTDIR/homology_summary_all_samples.txt"

echo ""
echo "======================================================================="
echo "All samples complete!"
echo "Summary: $OUTDIR/homology_summary_all_samples.txt"
echo "======================================================================="
