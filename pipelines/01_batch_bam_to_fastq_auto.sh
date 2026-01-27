#!/bin/bash

# Batch BAM to FASTQ conversion script (AUTO VERSION - no prompts)
# Automatically skips existing files without asking

# set -e  # Disabled to prevent exit on counter increment

# Directories
BAM_DIR="/home/hawaii-agriculture-research-cent/10193_78f4e76953b74a288de0f6e0bbd93029/central_memory2/pm_storage/caolinzhi/shujuxiazai/zhikong/163_20251121/SG251015HK01S67N1_福建农林大学Maringa_HiFi、HiC及转录组建库测序/BC2025100554-PB-hifi-9samples"
OUTPUT_DIR="$HOME/moringa_pangenome/raw_data/hifi"
LOG_DIR="$HOME/moringa_pangenome/logs"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Sample list (all 9 samples)
SAMPLES=(
    "Mo-TH-30"
    "Mo-TH-55"
    "Mo-TH-16"
    "Mo-TH-43"
    "Mo-TH-6"
    "Mo-TH-63"
    "Mo-TH-66"
    "Mo-TW-13"
    "Mo-US-5"
)

echo "========================================"
echo "Batch BAM to FASTQ Conversion (AUTO)"
echo "========================================"
echo "Start time: $(date)"
echo "BAM directory: $BAM_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOG_DIR"
echo ""

# Counters
TOTAL=${#SAMPLES[@]}
SUCCESS=0
SKIPPED=0
FAILED=0

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "Processing: $SAMPLE ($((SUCCESS + SKIPPED + FAILED + 1))/$TOTAL)"
    echo "----------------------------------------"
    
    # Define file paths
    BAM_FILE="$BAM_DIR/$SAMPLE/${SAMPLE}.bam"
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}_hifi.fastq.gz"
    LOG_FILE="$LOG_DIR/bam_to_fastq_${SAMPLE}_$(date +%Y%m%d_%H%M%S).log"
    
    # Check if output file already exists
    if [ -f "$OUTPUT_FILE" ]; then
        echo "Output file already exists: $OUTPUT_FILE"
        echo "Skipping $SAMPLE..."
        ((SKIPPED++))
        echo ""
        continue
    fi
    
    # Check if BAM file exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "ERROR: BAM file not found: $BAM_FILE"
        ((FAILED++))
        echo ""
        continue
    fi
    
    echo "Input BAM: $BAM_FILE"
    echo "Output FASTQ: $OUTPUT_FILE"
    echo "Log file: $LOG_FILE"
    echo ""
    
    # Start conversion
    echo "Starting conversion at $(date)..."
    {
        echo "=== BAM to FASTQ Conversion Log ==="
        echo "Sample: $SAMPLE"
        echo "Start time: $(date)"
        echo "Command: samtools view $BAM_FILE | awk '{print \"@\"\$1\"\\n\"\$10\"\\n+\\n\"\$11}' | gzip > $OUTPUT_FILE"
        echo ""
        
        # Run conversion
        samtools view "$BAM_FILE" | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > "$OUTPUT_FILE"
        
        echo ""
        echo "End time: $(date)"
        echo ""
        
        # Validate output
        echo "=== Validation ==="
        if [ -f "$OUTPUT_FILE" ]; then
            # Count reads in BAM
            BAM_READS=$(samtools view -c "$BAM_FILE")
            echo "Reads in BAM: $BAM_READS"
            
            # Count reads in FASTQ
            FASTQ_READS=$(zcat "$OUTPUT_FILE" | grep -c "^@" || true)
            echo "Reads in FASTQ: $FASTQ_READS"
            
            # Check if counts match
            if [ "$BAM_READS" -eq "$FASTQ_READS" ]; then
                echo "✓ Read counts match!"
                echo "✓ Conversion successful!"
            else
                echo "✗ Read count mismatch!"
                echo "✗ Conversion may have failed!"
            fi
        else
            echo "✗ Output file not created!"
        fi
        
    } > "$LOG_FILE" 2>&1
    
    # Check if conversion was successful
    if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
        echo "✓ Conversion completed successfully!"
        ((SUCCESS++))
    else
        echo "✗ Conversion failed!"
        ((FAILED++))
    fi
    
    echo ""
done

# Summary
echo "========================================"
echo "Batch Conversion Summary"
echo "========================================"
echo "Total samples: $TOTAL"
echo "Successful: $SUCCESS"
echo "Skipped (already exist): $SKIPPED"
echo "Failed: $FAILED"
echo "End time: $(date)"
echo "========================================"

# Exit with error if any conversions failed
if [ $FAILED -gt 0 ]; then
    exit 1
fi
