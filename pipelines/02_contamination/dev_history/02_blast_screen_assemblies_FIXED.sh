#!/bin/bash

# BLAST-based Contamination Screening - FIXED VERSION
# Strategy: Filter based on SUBJECT coverage (contig covers >80% of reference)
# NOT query coverage (which would flag nuclear contigs with small chloroplast insertions)

set -e

BASE_DIR="$HOME/moringa_pangenome"
CONTAM_DIR="$BASE_DIR/contamination_screening"
REF_DIR="$CONTAM_DIR/references"
RESULTS_DIR="$CONTAM_DIR/results"

mkdir -p "$RESULTS_DIR"

echo "========================================"
echo "CONTAMINATION SCREENING via BLASTN"
echo "========================================"
echo "FIXED STRATEGY:"
echo "  - Filter by SUBJECT coverage (scovs), not query coverage"
echo "  - Contig is contamination if it covers >80% of reference genome"
echo "  - Example: Contig covering >80% of chloroplast = chloroplast contig"
echo ""
echo "This prevents false positives from:"
echo "  - Small chloroplast insertions in nuclear genome"
echo "  - NUMTs (nuclear mitochondrial DNA)"
echo "========================================"
echo ""

# Check databases exist
if [ ! -f "$REF_DIR/chloroplast_db.nhr" ]; then
    echo "ERROR: BLAST databases not found!"
    echo "Please run: bash 01_download_contam_references.sh"
    exit 1
fi

# Samples to screen
SAMPLES=(Mo-TH-30 Mo-TH-55 Mo-TH-16 Mo-TH-43 Mo-TH-6 Mo-TH-63 Mo-TH-66 Mo-TW-13 Mo-US-5)

# BLAST parameters
THREADS=16
EVALUE="1e-50"  # More stringent
MAX_TARGET_SEQS=5

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "========================================="
    echo "Processing: $SAMPLE"
    echo "========================================="
    
    ASSEMBLY="$BASE_DIR/assemblies/production/$SAMPLE/${SAMPLE}.primary.fasta"
    
    if [ ! -f "$ASSEMBLY" ]; then
        echo "  [NOT FOUND] Assembly: $ASSEMBLY"
        echo "  Skipping..."
        continue
    fi
    
    # Get assembly stats
    TOTAL_CONTIGS=$(grep -c "^>" "$ASSEMBLY")
    echo "  Assembly: $TOTAL_CONTIGS contigs"
    
    # 1. Screen against CHLOROPLAST database
    echo ""
    echo "  Step 1: Screening against chloroplast database..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/chloroplast_db" \
        -out "$RESULTS_DIR/${SAMPLE}_chloroplast_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs "$MAX_TARGET_SEQS"
    
    # Filter: Subject coverage >80% (contig covers >80% of chloroplast reference)
    # Column 14 = scovhsp (subject coverage per HSP)
    # Also require minimum alignment length of 50kb
    awk '$14 > 80 && $4 > 50000' "$RESULTS_DIR/${SAMPLE}_chloroplast_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt"
    
    CHLORO_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt" | sort -u | wc -l)
    echo "     Chloroplast contigs: $CHLORO_HITS (>80% subject coverage, >50kb alignment)"
    
    # 2. Screen against MITOCHONDRIA database
    echo ""
    echo "  Step 2: Screening against mitochondria database..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/mitochondria_db" \
        -out "$RESULTS_DIR/${SAMPLE}_mitochondria_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs "$MAX_TARGET_SEQS"
    
    # Filter: Subject coverage >80%, minimum 100kb alignment (mitochondria are larger)
    awk '$14 > 80 && $4 > 100000' "$RESULTS_DIR/${SAMPLE}_mitochondria_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt"
    
    MITO_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt" | sort -u | wc -l)
    echo "     Mitochondria contigs: $MITO_HITS (>80% subject coverage, >100kb alignment)"
    
    # 3. Screen against CONTAMINANTS database (bacteria + vectors)
    echo ""
    echo "  Step 3: Screening against contaminants database..."
    blastn \
        -query "$ASSEMBLY" \
        -db "$REF_DIR/contaminants_db" \
        -out "$RESULTS_DIR/${SAMPLE}_contaminants_blast.txt" \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -max_target_seqs "$MAX_TARGET_SEQS"
    
    # For bacteria/vectors: use query coverage (contig is mostly bacterial)
    # OR high subject coverage with large alignment
    awk '($13 > 80) || ($14 > 80 && $4 > 10000)' "$RESULTS_DIR/${SAMPLE}_contaminants_blast.txt" | \
        sort -k1,1 -k12,12nr | \
        awk '!seen[$1]++' > "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt"
    
    CONTAM_HITS=$(cut -f1 "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt" | sort -u | wc -l)
    echo "     Bacterial/Vector contigs: $CONTAM_HITS"
    
    # 4. Create master lists of contaminated contig IDs
    echo ""
    echo "  Step 4: Creating classification lists..."
    
    # Chloroplast IDs
    cut -f1 "$RESULTS_DIR/${SAMPLE}_chloroplast_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_chloroplast_ids.txt"
    
    # Mitochondria IDs
    cut -f1 "$RESULTS_DIR/${SAMPLE}_mitochondria_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_mitochondria_ids.txt"
    
    # True contaminants IDs
    cut -f1 "$RESULTS_DIR/${SAMPLE}_contaminants_hits.txt" | sort -u > "$RESULTS_DIR/${SAMPLE}_contaminant_ids.txt"
    
    # Combined organelles
    cat "$RESULTS_DIR/${SAMPLE}_chloroplast_ids.txt" \
        "$RESULTS_DIR/${SAMPLE}_mitochondria_ids.txt" | \
        sort -u > "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt"
    
    ORGANELLE_TOTAL=$(wc -l < "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt")
    
    # All contamination
    cat "$RESULTS_DIR/${SAMPLE}_organelle_ids.txt" \
        "$RESULTS_DIR/${SAMPLE}_contaminant_ids.txt" | \
        sort -u > "$RESULTS_DIR/${SAMPLE}_all_contam_ids.txt"
    
    TOTAL_CONTAM=$(wc -l < "$RESULTS_DIR/${SAMPLE}_all_contam_ids.txt")
    
    # Calculate percentages
    CONTAM_PERCENT=$(awk -v t="$TOTAL_CONTAM" -v c="$TOTAL_CONTIGS" 'BEGIN {printf "%.2f", (t/c)*100}')
    
    echo "     Total contamination: $TOTAL_CONTAM contigs ($CONTAM_PERCENT%)"
    echo "       - Organelles: $ORGANELLE_TOTAL contigs"
    echo "       - True contaminants: $CONTAM_HITS contigs"
    
    # 5. Create detailed summary
    {
        echo "Sample: $SAMPLE"
        echo "Assembly: $ASSEMBLY"
        echo "Total contigs: $TOTAL_CONTIGS"
        echo ""
        echo "CONTAMINATION SUMMARY:"
        echo "  Chloroplast contigs: $CHLORO_HITS"
        echo "  Mitochondrial contigs: $MITO_HITS"
        echo "  Bacterial/Vector contigs: $CONTAM_HITS"
        echo "  Total organelles: $ORGANELLE_TOTAL"
        echo "  Total contamination: $TOTAL_CONTAM ($CONTAM_PERCENT%)"
        echo ""
        echo "FILTERING CRITERIA:"
        echo "  Chloroplast: >80% subject coverage + >50kb alignment"
        echo "  Mitochondria: >80% subject coverage + >100kb alignment"
        echo "  Bacteria/Vectors: >80% query coverage OR (>80% subject + >10kb)"
        echo ""
    } > "$RESULTS_DIR/${SAMPLE}_summary.txt"
    
    echo ""
    echo "  [OK] Summary saved"
    
done

echo ""
echo "========================================"
echo "BLAST SCREENING COMPLETE!"
echo "========================================"
echo ""
echo "Results location: $RESULTS_DIR"
echo ""
echo "NEXT STEPS:"
echo "  1. Review results: bash 03_generate_report.sh"
echo "  2. Create cleaned assemblies: bash 04_create_cleaned_assemblies.sh"
echo ""
echo "========================================"
