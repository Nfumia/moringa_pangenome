# Moringa Pangenome Assembly Pipeline

Complete workflow for generating chromosome-scale genome assemblies from PacBio HiFi and Hi-C sequencing data.

---

## Pipeline Overview
```
Raw Reads (HiFi + Hi-C)
    ↓
[1] HiFi Assembly (hifiasm)
    ↓
[2] Heterozygosity Assessment (GenomeScope2)
    ↓
    ├─→ Homozygous samples → [3] Hi-C Scaffolding (YaHS)
    │
    └─→ Heterozygous samples → [4] Haplotype Phasing (hifiasm + Hi-C)
    ↓
[5] Contamination Screening (BLAST)
    ↓
[6] Quality Control (BUSCO, Quast, Merqury)
    ↓
[7] Pangenome Construction (Minigraph-Cactus)
```

---

## Step 1: HiFi Assembly

**Input:** PacBio HiFi reads (~30-60x coverage)  
**Tool:** Hifiasm v0.19.8  
**Runtime:** ~30 minutes per sample
```bash
# For initial assembly (all samples)
hifiasm -o sample.asm \
    -t 48 \
    hifi_reads.fastq.gz

# Convert GFA to FASTA
awk '/^S/{print ">"$2; print $3}' sample.asm.bp.p_ctg.gfa > sample.primary.fasta
```

**Output:**
- Primary assembly (gap-free, high contiguity)
- N50: 8-15 Mb
- Size: 290-314 Mb

---

## Step 2: Heterozygosity Assessment

**Input:** HiFi reads  
**Tool:** GenomeScope2 v2.0  
**Purpose:** Classify samples as homozygous vs heterozygous
```bash
# Count k-mers
jellyfish count -C -m 21 -s 3G -t 24 hifi_reads.fastq.gz -o reads.jf

# Export histogram
jellyfish histo -t 24 reads.jf > reads.histo

# Run GenomeScope2
genomescope2 -i reads.histo -o genomescope_output -k 21
```

**Classification:**
- `peak_het = -1` → Homozygous (use Hi-C scaffolding)
- `peak_het > 0` → Heterozygous (use haplotype phasing)

---

## Step 3: Hi-C Scaffolding (Homozygous Samples)

**Input:** Primary assembly + Hi-C reads  
**Tool:** YaHS v1.2a.2  
**Runtime:** ~2-4 hours per sample
```bash
# Index assembly
bwa-mem2 index assembly.fasta
samtools faidx assembly.fasta

# Align Hi-C reads
bwa-mem2 mem -5SP -t 48 assembly.fasta hic_R1.fastq.gz hic_R2.fastq.gz | \
    samtools view -@ 8 -Sb - | \
    samtools sort -@ 8 -o aligned.bam -
samtools index aligned.bam

# Scaffold with YaHS
yahs assembly.fasta aligned.bam -o scaffolds

# Convert scaffold output
awk '/^S/{print ">"$2; print $3}' scaffolds_scaffolds_final.agp > scaffolds_final.fasta
```

**Results:**
- Chromosome-scale scaffolds (29-42 Mb longest)
- N50 improvement: 1.2-2.0x
- Minimal sequence loss (<0.1%)

---

## Step 4: Haplotype Phasing (Heterozygous Samples)

**Input:** HiFi reads + Hi-C reads  
**Tool:** Hifiasm v0.19.8 with Hi-C integration  
**Runtime:** ~20 minutes per sample
```bash
# Run hifiasm with Hi-C phasing
hifiasm -o sample.asm \
    --h1 hic_R1.fastq.gz \
    --h2 hic_R2.fastq.gz \
    --n-hap 2 \
    -t 48 \
    hifi_reads.fastq.gz

# Extract haplotypes
awk '/^S/{print ">"$2; print $3}' sample.asm.hic.hap1.p_ctg.gfa > sample.hap1.fasta
awk '/^S/{print ">"$2; print $3}' sample.asm.hic.hap2.p_ctg.gfa > sample.hap2.fasta
```

**Output:**
- Two phased haplotypes per sample
- Phasing ratio: 1.76-1.83x (expected ~2.0x)
- Size differences represent biological variation

---

## Step 5: Contamination Screening

**Input:** Assemblies  
**Tool:** BLAST+ v2.14.0  
**Runtime:** ~4-6 hours per sample
```bash
# Download and index references
# (chloroplast, mitochondria, bacterial/vector)

# BLAST against chloroplast
blastn -query assembly.fasta \
       -db chloroplast_db \
       -outfmt "6 qseqid sseqid pident length qlen slen qcovs" \
       -evalue 1e-10 \
       -num_threads 24 \
       -out chloroplast_hits.txt

# Filter: Remove contigs where >90% aligns to organelles
awk '$7 > 90 && $4 > 50000' chloroplast_hits.txt > remove_list.txt

# Create cleaned assembly
seqkit grep -v -f remove_list.txt assembly.fasta > assembly.cleaned.fasta
```

**Results:**
- Remove ~6-8 Mb chloroplast
- Remove ~0.2-0.5 Mb mitochondria
- Final size: ~302-304 Mb per sample

---

## Step 6: Quality Control

### 6a. BUSCO (Gene Completeness)
```bash
busco -i assembly.fasta \
      -o sample_busco \
      -m genome \
      -l embryophyta_odb10 \
      -c 24
```

**Expected:** >95% complete BUSCOs

### 6b. QUAST (Assembly Statistics)
```bash
quast assembly.fasta \
      -r reference.fasta \
      -o quast_output \
      --threads 24 \
      --large
```

### 6c. Merqury (K-mer based QC)
```bash
meryl k=21 count output reads.meryl reads.fastq.gz
merqury reads.meryl assembly.fasta output_prefix
```

**Expected QV:** >40 (>99.99% accuracy)

---

## Step 7: Pangenome Construction

**Input:** All cleaned assemblies (24 total: 9 primary + 3 scaffolded + 12 haplotypes)  
**Tool:** Minigraph-Cactus  
**Runtime:** Several days
```bash
# Create sample list
cat > samples.txt << EOF
sample1 path/to/assembly1.fasta
sample2 path/to/assembly2.fasta
...
