#!/bin/bash

################################################################################
# Master Script - Run All Analyses on All Samples
# WARNING: Will take 24-48 hours total
################################################################################

echo "======================================================================="
echo "COMPLETE HAPLOTYPE ANALYSIS SUITE - ALL 6 SAMPLES"
echo "======================================================================="
echo ""
echo "This will run all 5 analyses on all 6 heterozygous samples:"
echo "  1. Gap estimation       (~10 min)"
echo "  2. Homology analysis    (~6 hours)"
echo "  3. Synteny analysis     (~3 hours)"
echo "  4. Coverage mapping     (~12 hours)"
echo "  5. K-mer analysis       (~24 hours)"
echo ""
echo "Total estimated time: 24-48 hours"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
fi

START_TIME=$(date +%s)

echo ""
echo "Starting analysis suite at $(date)"
echo ""

# Run each analysis
bash batch_01_gap_estimation_all_samples.sh
bash batch_02_homology_analysis_all_samples.sh
bash batch_03_synteny_analysis_all_samples.sh
bash batch_04_coverage_mapping_all_samples.sh
bash batch_05_kmer_analysis_all_samples.sh

END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 3600 ))

echo ""
echo "======================================================================="
echo "ALL ANALYSES COMPLETE!"
echo "======================================================================="
echo "Total time: ${ELAPSED} hours"
echo "Results saved in: ../haplotype_analysis/batch_results/"
echo ""
echo "Summary files:"
echo "  01_gaps/gap_summary_all_samples.txt"
echo "  02_homology/homology_summary_all_samples.txt"
echo "  03_synteny/synteny_summary_all_samples.txt"
echo "  04_coverage/coverage_summary_all_samples.txt"
echo "  05_kmers/kmer_summary_all_samples.txt"
echo "======================================================================="
