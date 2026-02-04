#!/bin/bash

# Download Contamination Reference Databases
# Based on colleague's workflow section I.3.2

set -e

BASE_DIR="$HOME/moringa_pangenome/contamination_screening"
REF_DIR="$BASE_DIR/references"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "========================================"
echo "Downloading Contamination References"
echo "========================================"
echo "Based on: RMcontamination workflow"
echo "Priority: Chloroplast > Mitochondria > Contaminants"
echo ""

# Check if makeblastdb is available
if ! command -v makeblastdb &> /dev/null; then
    echo "ERROR: makeblastdb not found!"
    echo "Please install BLAST+ or activate the correct conda environment"
    echo "  conda activate hifi_assembly"
    exit 1
fi

echo "BLAST+ version:"
makeblastdb -version | head -1
echo ""

# 1. PRIORITY: All 5 Moringa oleifera chloroplast genomes
echo "========================================="
echo "1. CHLOROPLASTS (Highest Priority)"
echo "========================================="
echo "Downloading 5 Moringa oleifera chloroplast genomes..."
echo ""

CHLORO_ACC=(NC_041432.1 MH939149.1 MK713965.1 ON855354.1 MK726020.1)
for ACC in "${CHLORO_ACC[@]}"; do
    echo "   Downloading $ACC..."
    wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ACC&rettype=fasta&retmode=text" -O "moringa_chloroplast_${ACC}.fasta"
    if [ -s "moringa_chloroplast_${ACC}.fasta" ]; then
        SIZE=$(grep -v "^>" "moringa_chloroplast_${ACC}.fasta" | tr -d '\n' | wc -c)
        SIZE_KB=$((SIZE / 1000))
        echo "   [OK] $ACC (~${SIZE_KB} kb)"
    else
        echo "   [FAILED] $ACC"
        exit 1
    fi
done

# Combine all Moringa chloroplasts
cat moringa_chloroplast_*.fasta > moringa_chloroplast_all.fasta
echo ""
echo "   [OK] 5 Moringa chloroplasts combined"

# Arabidopsis chloroplast for validation
echo "   Downloading Arabidopsis chloroplast..."
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000932.1&rettype=fasta&retmode=text" -O arabidopsis_chloroplast.fasta
echo "   [OK] Arabidopsis chloroplast (NC_000932.1, ~155 kb)"

# 2. PRIORITY: Mitochondrial genomes from related species
echo ""
echo "========================================="
echo "2. MITOCHONDRIA (High Priority)"
echo "========================================="
echo "Downloading mitochondrial genomes from Brassicales..."
echo ""

# Papaya (closest relative to Moringa in Brassicales)
echo "   Carica papaya (closest relative):"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012116.1&rettype=fasta&retmode=text" -O papaya_mito_1.fasta
echo "   [OK] NC_012116.1"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=EU431224.1&rettype=fasta&retmode=text" -O papaya_mito_2.fasta
echo "   [OK] EU431224.1"

# Brassica napus (Brassicales)
echo "   Brassica napus:"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_008285.1&rettype=fasta&retmode=text" -O brassica_mito.fasta
echo "   [OK] NC_008285.1"

# Raphanus sativus (Brassicales)
echo "   Raphanus sativus:"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_018551.1&rettype=fasta&retmode=text" -O raphanus_mito.fasta
echo "   [OK] NC_018551.1"

# Arabidopsis (model plant)
echo "   Arabidopsis thaliana:"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_037304.1&rettype=fasta&retmode=text" -O arabidopsis_mito.fasta
echo "   [OK] NC_037304.1"

# Moringa partial mitochondrial genes (supplementary)
echo "   Moringa oleifera (partial genes):"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KU739782.1&rettype=fasta&retmode=text" -O moringa_mito_rps3.fasta
echo "   [OK] KU739782.1 (rps3 partial)"
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AF451593.1&rettype=fasta&retmode=text" -O moringa_mito_nd4.fasta
echo "   [OK] AF451593.1 (ND4 partial)"

# Combine all mitochondrial references
cat papaya_mito_*.fasta brassica_mito.fasta raphanus_mito.fasta arabidopsis_mito.fasta moringa_mito_*.fasta > all_mitochondria.fasta
echo ""
echo "   [OK] 7 mitochondrial references combined"

# 3. Bacterial contaminants
echo ""
echo "========================================="
echo "3. BACTERIAL CONTAMINANTS"
echo "========================================="

echo "   E. coli K-12..."
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" -O ecoli_k12.fasta
echo "   [OK] NC_000913.3 (~4.6 Mb)"

echo "   Agrobacterium tumefaciens..."
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_003063.2&rettype=fasta&retmode=text" -O agrobacterium.fasta
echo "   [OK] NC_003063.2 (~2.8 Mb)"

# 4. Vector sequences
echo ""
echo "========================================="
echo "4. VECTOR SEQUENCES"
echo "========================================="

echo "   UniVec Core database..."
wget -q "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" -O univec.fasta
if [ -s univec.fasta ]; then
    UNIVEC_COUNT=$(grep -c "^>" univec.fasta 2>/dev/null || echo "0")
    echo "   [OK] UniVec downloaded ($UNIVEC_COUNT sequences)"
else
    echo "   [WARNING] UniVec download failed (non-critical)"
    touch univec.fasta
fi

# 5. Create categorized BLAST databases
echo ""
echo "========================================="
echo "5. CREATING BLAST DATABASES"
echo "========================================="
echo ""

# Function to create database and verify
create_db() {
    local INPUT=$1
    local OUTPUT=$2
    local TITLE=$3
    
    echo "Creating $OUTPUT..."
    if makeblastdb -in "$INPUT" -dbtype nucl -out "$OUTPUT" -title "$TITLE" -parse_seqids > /dev/null 2>&1; then
        # Verify database files were created
        if [ -f "${OUTPUT}.nhr" ] && [ -f "${OUTPUT}.nin" ] && [ -f "${OUTPUT}.nsq" ]; then
            SEQ_COUNT=$(grep -c "^>" "$INPUT" 2>/dev/null || echo "?")
            echo "  [OK] $OUTPUT created ($SEQ_COUNT sequences)"
            return 0
        else
            echo "  [ERROR] $OUTPUT database files not created!"
            return 1
        fi
    else
        echo "  [ERROR] makeblastdb failed for $OUTPUT"
        return 1
    fi
}

# Chloroplast database (5 Moringa + 1 Arabidopsis)
cat moringa_chloroplast_all.fasta arabidopsis_chloroplast.fasta > all_chloroplasts.fasta
create_db "all_chloroplasts.fasta" "chloroplast_db" "Chloroplasts" || exit 1

# Mitochondria database (5 complete + 2 partial)
create_db "all_mitochondria.fasta" "mitochondria_db" "Mitochondria" || exit 1

# Organelles combined
cat all_chloroplasts.fasta all_mitochondria.fasta > organelles.fasta
create_db "organelles.fasta" "organelles_db" "Organelles" || exit 1

# Contaminants (bacteria + vectors)
cat ecoli_k12.fasta agrobacterium.fasta univec.fasta > contaminants.fasta
create_db "contaminants.fasta" "contaminants_db" "Contaminants" || exit 1

# All combined
cat organelles.fasta contaminants.fasta > all_contam.fasta
create_db "all_contam.fasta" "all_contam_db" "All_Contamination" || exit 1

# Final summary
echo ""
echo "========================================"
echo "DOWNLOAD COMPLETE - SUMMARY"
echo "========================================"
echo ""
echo "CHLOROPLAST REFERENCES (6 genomes):"
echo "  - 5 Moringa oleifera complete genomes"
echo "    NC_041432.1, MH939149.1, MK713965.1"
echo "    ON855354.1, MK726020.1"
echo "  - 1 Arabidopsis thaliana (NC_000932.1)"
echo "  Expected size: ~150-165 kb"
echo ""
echo "MITOCHONDRIAL REFERENCES (7 sequences):"
echo "  - 2 Carica papaya complete (NC_012116.1, EU431224.1)"
echo "  - 1 Brassica napus complete (NC_008285.1)"
echo "  - 1 Raphanus sativus complete (NC_018551.1)"
echo "  - 1 Arabidopsis thaliana complete (NC_037304.1)"
echo "  - 2 Moringa oleifera partial genes (KU739782.1, AF451593.1)"
echo "  Expected size: ~200-700 kb (varies by species)"
echo ""
echo "CONTAMINANT REFERENCES:"
echo "  - E. coli K-12 (NC_000913.3)"
echo "  - Agrobacterium tumefaciens (NC_003063.2)"
echo "  - UniVec database"
echo ""
echo "BLAST DATABASES CREATED AND VERIFIED:"
echo "  1. chloroplast_db - 6 genomes"
echo "  2. mitochondria_db - 7 references"
echo "  3. organelles_db - Combined (13 total)"
echo "  4. contaminants_db - Bacteria + vectors"
echo "  5. all_contam_db - Everything"
echo ""
echo "Database location: $REF_DIR"
echo ""
echo "All database files (.nhr, .nin, .nsq) verified present"
echo ""
echo "NEXT STEP:"
echo "  bash 02_blast_screen_assemblies.sh"
echo "========================================"
