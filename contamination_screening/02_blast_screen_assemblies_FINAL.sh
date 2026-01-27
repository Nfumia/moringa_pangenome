#!/bin/bash

# BLAST-based Contamination Screening - FINAL CORRECT VERSION
# Filter: Contig is >90% chloroplast/mitochondria = remove as organellar
# Keep large contigs even if they have small organellar insertions

set -e

BASE_DIR="$HOME/moringa_pangenome"
CONTAM_DIR="$BASE_DIR/contamination_screening"
REF_DIR="$CONTAM_DIR/references"
RESULTS_DIR="$CONTAM_DIR/results"

mkdir -p "$RESULTS_DIR"

echo "========================================"
echo "CONTAMINATION SCREENING via BLASTN"
echo "========================================"
echo "STRATEGY:"
echo "  - Contig is contamination if >90% of contig length aligns"
echo "  - Prevents removing large nuclear contigs with small insertions"
echo "  - Identifies fragmented organellar sequences"
echo "========================================"
echo ""

# Check databases exist
if [ ! -f "$REF_DIR/chloroplast_db.nhr" ]; then
    echo "ERROR: BLAST databases not found!"
    exit 1
fi

# Samples
SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

THREADS=16
EVALUE="1e-50"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "========================================="
    echo "Processing: $SAMPLE"
    echo "========================================="
    
    ASSEMBLY="$BASE_DIR/assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
    
    if [ ! -f "$ASSEMBLY" ]; then
        echo "  [NOT FOUND] $ASSEMBLY"
        continue
    fi
    
    TOTAL_CONTIGS=$(grep -c "^>" "$ASSEMBLY")
    echo "  Assembly: $TOTAL_CONTIGS contigs"
    
    # 1. CHLOROPLAST screening
    echo ""
    echo "  Step 1: Chloroplast screening..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/chloroplast_db" \
        -out "$RESULTS_DIR/${SAMPLE}_chloroplast_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs 1
    
    # Filter: >90% of contig is chloroplast + >50kb alignment
    awk '$14 > 80 && $4 > 50000 && $4/$5 > 0.90' "$RESULTS_DIR/${SAMPLE}_chloroplast_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt"
    
    CHLORO_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt" | sort -u | wc -l)
    echo "     Chloroplast contigs: $CHLORO_HITS (>90% of contig is chloroplast)"
    
    # 2. MITOCHONDRIA screening
    echo ""
    echo "  Step 2: Mitochondria screening..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/mitochondria_db" \
        -out "$RESULTS_DIR/${SAMPLE}_mitochondria_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs 1
    
    # Filter: >90% of contig is mitochondria + >100kb alignment
    awk '$14 > 80 && $4 > 100000 && $4/$5 > 0.90' "$RESULTS_DIR/${SAMPLE}_mitochondria_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt"
    
    MITO_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt" | sort -u | wc -l)
    echo "     Mitochondria contigs: $MITO_HITS (>90% of contig is mitochondria)"
    
    # 3. BACTERIAL/VECTOR screening
    echo ""
    echo "  Step 3: Bacterial/Vector screening..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/contaminants_db" \
        -out "$RESULTS_DIR/${SAMPLE}_contaminants_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs 1
    
    # For bacteria: >80% of contig matches (true bacterial contamination)
    awk '$4/$5 > 0.80' "$RESULTS_DIR/${SAMPLE}_contaminants_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt"
    
    CONTAM_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt" | sort -u | wc -l)
    echo "     Bacterial/Vector contigs: $CONTAM_HITS (>80% of contig)"
    
    # 4. Create ID lists
    echo ""
    echo "  Step 4: Creating classification lists..."
    
    cut -f1 "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_chloroplast_ids.txt"
    cut -f1 "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_mitochondria_ids.txt"
    cut -f1 "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_contaminant_ids.txt"
    
    cat "$RESULTS_DIR/${SAMPLE}_chloroplast_ids.txt" \
        "$RESULTS_DIR/${SAMPLE}_mitochondria_ids.txt" | \
        sort -u > "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt"
    
    ORGANELLE_TOTAL=$(wc -l < "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt")
    
    cat "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt" \
        "$RESULTS_DIR/${SAMPLE}_contaminant_ids.txt" | \
        sort -u > "$RESULTS_DIR/${SAMPLE}_all_contam_ids.txt"
    
    TOTAL_CONTAM=$(wc -l < "$RESULTS_DIR/${SAMPLE}_all_contam_ids.txt")
    CONTAM_PERCENT=$(awk -v t="$TOTAL_CONTAM" -v c="$TOTAL_CONTIGS" 'BEGIN {printf "%.2f", (t/c)*100}')
    
    echo "     Total contamination: $TOTAL_CONTAM contigs ($CONTAM_PERCENT%)"
    echo "       - Organelles: $ORGANELLE_TOTAL"
    echo "       - Bacterial/Vector: $CONTAM_HITS"
    
    # 5. Summary
    {
        echo "Sample: $SAMPLE"
        echo "Total contigs: $TOTAL_CONTIGS"
        echo ""
        echo "CONTAMINATION:"
        echo "  Chloroplast contigs: $CHLORO_HITS"
        echo "  Mitochondrial contigs: $MITO_HITS"
        echo "  Bacterial/Vector contigs: $CONTAM_HITS"
        echo "  Total contamination: $TOTAL_CONTAM ($CONTAM_PERCENT%)"
        echo ""
        echo "CRITERIA:"
        echo "  Chloroplast: >90% of contig aligns to chloroplast"
        echo "  Mitochondria: >90% of contig aligns to mitochondria"
        echo "  Bacteria: >80% of contig aligns to bacteria"
    } > "$RESULTS_DIR/${SAMPLE}_summary.txt"
    
    echo "  [OK] Complete"
done

echo ""
echo "========================================"
echo "SCREENING COMPLETE"
echo "========================================"
