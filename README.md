# Moringa oleifera Pangenome Project

## Project Overview
De novo genome assembly and pangenome construction for 9 Moringa oleifera samples using PacBio HiFi and Hi-C sequencing data.

**Project Status:** Comprehensive QC in progress (February 2026)

## Samples
- **Total samples:** 9
- **Homozygous samples (3):** Mo-TH-30, Mo-TH-66, Mo-US-5
- **Heterozygous samples (6):** Mo-TH-55, Mo-TH-16, Mo-TH-43, Mo-TH-6, Mo-TH-63, Mo-TW-13

## Assembly Statistics

### HiFi-Only Assemblies (All 9 samples - COMPLETE)
- **Assembly size:** 290-314 Mb
- **Contig count:** 563-992 contigs
- **N50:** 8-15 Mb
- **Longest contig:** 17-30 Mb
- **Quality:** 100% gap-free, high contiguity

### Hi-C Scaffolded Assemblies (3 homozygous samples - COMPLETE)
| Sample | Size (Mb) | Scaffolds | N50 (Mb) | Longest (Mb) | Improvement |
|--------|-----------|-----------|----------|--------------|-------------|
| Mo-TH-30 | 310.1 | 569 | 17.0 | 29.3 | 1.21x N50 |
| Mo-TH-66 | 310.8 | 686 | 15.2 | 41.7 | 1.88x N50 |
| Mo-US-5 | 312.3 | 689 | 17.4 | 42.3 | 1.95x N50 |

**Achievement:** Chromosome-scale scaffolds (29-42 Mb = likely full chromosome arms)

### Haplotype-Phased Assemblies (6 heterozygous samples - COMPLETE)
| Sample | Hap1 (Mb) | Hap2 (Mb) | Total (Mb) | Ratio | Peak Het |
|--------|-----------|-----------|------------|-------|----------|
| Mo-TH-55 | 308.5 | 260.6 | 569.1 | 1.81x | 24 |
| Mo-TH-16 | 286.6 | 286.8 | 573.4 | 1.83x | 32 |
| Mo-TH-43 | 298.3 | 267.2 | 565.5 | 1.81x | 22 |
| Mo-TH-6 | 266.3 | 267.3 | 533.6 | 1.76x | 12 |
| Mo-TH-63 | 292.5 | 260.2 | 552.7 | 1.80x | 7 |
| Mo-TW-13 | 287.1 | 241.9 | 529.0 | 1.82x | 29 |

**Quality:** All ratios 1.76-1.83x indicate successful haplotype resolution
---

## Repository Structure

```
moringa_pangenome/
├── pipelines/                      # All analysis pipelines
│   ├── 01_assembly/                # BAM→FASTQ, HiFi assembly, scaffolding, phasing
│   ├── 02_contamination/           # BLAST screening and cleaning
│   ├── 03_validation/              # Haplotype validation (5 methods)
│   ├── 04_qc/                      # Comprehensive QC (BUSCO, Merqury)
│   ├── 05_annotation/              # Genome annotation (future)
│   ├── 06_pangenome/               # Pan-genome construction (future)
│   └── utils/                      # Utility scripts
├── assemblies/                     # Assembly files (not in git)
├── results/                        # Analysis results (not in git)
├── raw_data/                       # Sequencing reads (not in git)
└── docs/                           # Documentation
```

---

## Workflow Progress

### Phase 1: Initial Assembly (COMPLETE)
**Status:** All 9 HiFi assemblies generated successfully

**Method:** Hifiasm v0.19.8
- Input: PacBio HiFi reads (30-60x coverage)
- Output: Primary assemblies (290-314 Mb each)
- Quality: High contiguity, gap-free

**Key achievements:**
- N50 values: 8.9-15.3 Mb
- Longest contigs: 17-30 Mb
- All assemblies exceed reference quality

---

### Phase 2: Heterozygosity Assessment (COMPLETE)
**Status:** All samples classified

**Method:** GenomeScope2 k-mer analysis
- Classified 3 homozygous samples (peak_het = -1)
- Classified 6 heterozygous samples (peak_het = 7-32)
- Informed downstream processing strategy

**Biological insight:** 67% heterozygous samples reflects Moringa's mixed mating system

---

### Phase 3: Hi-C Scaffolding (COMPLETE)
**Status:** All 3 homozygous samples scaffolded to chromosome scale

**Method:** YaHS v1.2a.2 with Hi-C data
- **Input:** HiFi assemblies + Hi-C reads
- **Process:** BWA-MEM alignment + YaHS scaffolding
- **Output:** Chromosome-scale scaffolds

**Results:**
- Achieved 29-42 Mb longest scaffolds
- N50 improvements: 1.2-2.0x
- 18-31% reduction in sequence count
- No sequence loss (<0.1% size change)

**Runtime per sample:** ~2-4 hours (BWA alignment + YaHS)

---

### Phase 4: Haplotype Phasing (COMPLETE)
**Status:** All 6 heterozygous samples phased into haplotypes

**Method:** Hifiasm integrated Hi-C phasing
- **Input:** HiFi reads + Hi-C reads
- **Command:** `hifiasm --h1 <hic_R1> --h2 <hic_R2> -t 48 --n-hap 2`
- **Output:** Haplotype 1 and Haplotype 2 assemblies

**Results:**
- Successfully separated parental haplotypes
- Phasing ratios: 1.76-1.83x (expected ~2.0x for diploid)
- GC% differences between haplotypes confirm successful phasing
- Most balanced: Mo-TH-16 (286.6 vs 286.8 Mb, nearly identical)

**Runtime per sample:** ~20 minutes

**Quality metrics:**
- Haplotype size balance: 0.2-48 Mb imbalance
- All samples show successful diploid resolution
- Haplotype size differences validated (see Phase 4b below)

---

### Phase 4b: Haplotype Validation (COMPLETE - January 2026)
**Status:** All 5 validation methods completed successfully

**Purpose:** Determine if haplotype size imbalances represent technical artifacts or biological reality.

**Methods Completed:**

1. **Gap Estimation** - 100% gap-free (all 6 samples)
2. **Homology Analysis (BLAST)** - 0 collapsed regions detected
3. **Synteny Analysis (minimap2)** - 100% alignment coverage
4. **Coverage Mapping** - 88.6% conserved genome
5. **K-mer Analysis (meryl)** - 92.1% shared k-mers

**VALIDATION RESULTS:**

| Finding | Result |
|---------|--------|
| Collapsed regions | **0 Mb** (all methods agree) |
| Genome conservation | **~92%** (present in both haplotypes, >95% identity) |
| Structural variants | **~8%** (true biological differences) |
| Assembly quality | **Excellent** (100% gap-free, proper phasing) |

**Conclusion:** Size differences (0.2-47.9 Mb) represent **TRUE BIOLOGICAL VARIATION**, not technical artifacts. All 12 haplotype assemblies approved for pangenome construction.

**Key Findings:**
- Mo-TH-16 & Mo-TH-6: Nearly perfect haplotype balance (0.2-1.0 Mb difference)
- Mo-TH-55 & To-TW-13: Largest size differences (45-48 Mb) = significant parental divergence
- All samples: Structural variants include large indels, duplications, CNVs
- Biological interpretation: Evidence of outcrossing between divergent parents

**Documentation:** See `haplotype_analysis/batch_results/` for complete results and `HAPLOTYPE_VALIDATION_FINAL.md` for comprehensive report and `analysis_scripts/haplotype_analysis/` for complete workflows.
---

### Phase 4c: Comprehensive QC (IN PROGRESS - February 2026)

| Step | Status | Tool |
|------|--------|------|
| Assembly Statistics | ✅ Complete | seqkit |
| Q30 Scores | ✅ Complete | seqkit |
| Coverage Analysis | ✅ Complete | mosdepth (41-67x mean) |
| BUSCO | 🔄 Running | BUSCO 5.7.1 (eudicots_odb10) |
| Merqury | ⏳ Pending | Merqury |

---

### Phase 5: Contamination Screening (COMPLETE - December 2025)
**Status:** BLAST screening complete, cleaned assemblies generated

**Method:** BLASTN-based contamination detection
- Following Dr. Jingjing Zheng's RMcontamination workflow
- Threshold: >90% of contig aligns to reference = contamination

**Approach - Key Learnings:**

#### Initial Issues Identified:
1. **Query coverage vs Subject coverage:** Initial script used wrong metric
   - Problem: Flagged 606/727 contigs as chloroplast (clearly wrong!)
   - Cause: Used query coverage (flagged large contigs with small insertions)
   - Solution: Filter by alignment length covering >90% of contig length

2. **Whole contig removal concern:**
   - Problem: Shouldn't remove 10 Mb contig with 50kb chloroplast insertion
   - Solution: Only remove contigs where >90% of contig is organellar

#### Final Filtering Strategy:
**Chloroplast:**
- Contig is chloroplast if: >90% of contig aligns to chloroplast reference
- Additional criteria: >50kb alignment + >80% subject coverage
- Result: ~99 contigs (50-94 kb fragments) identified in Mo-TH-30
- Expected total: ~6-8 Mb chloroplast contamination per sample

**Mitochondria:**
- Contig is mitochondria if: >90% of contig aligns to mitochondria reference
- Additional criteria: >100kb alignment + >80% subject coverage
- Expected: 1-5 contigs, ~200-500 kb total per sample

**Bacterial/Vector:**
- Contig is contamination if: >80% of contig aligns to bacterial/vector reference
- Expected: 0-5 contigs (rare with HiFi data)

#### Reference Databases:
**Chloroplast (6 genomes):**
- 5 Moringa oleifera complete genomes (NC_041432.1, MH939149.1, MK713965.1, ON855354.1, MK726020.1)
- 1 Arabidopsis thaliana (NC_000932.1)
- Expected size: ~150-165 kb

**Mitochondria (7 references):**
- 2 Carica papaya complete (NC_012116.1, EU431224.1) - closest relative
- 1 Brassica napus (NC_008285.1)
- 1 Raphanus sativus (NC_018551.1)
- 1 Arabidopsis thaliana (NC_037304.1)
- 2 Moringa oleifera partial genes (KU739782.1, AF451593.1)
- Expected size: ~200-700 kb (varies by species)

**Contaminants:**
- E. coli K-12 (NC_000913.3)
- Agrobacterium tumefaciens (NC_003063.2)
- UniVec database (3,155 sequences)

#### Expected Results:
- **Original assembly:** ~310 Mb
- **Chloroplast removal:** ~6-8 Mb
- **Final cleaned assembly:** ~302-304 Mb
- **Contamination rate:** <3% (acceptable for plant genomes)

**Note:** Cleaned assemblies will still be larger than 240 Mb reference genome. This is normal and expected due to:
- Higher completeness
- Structural variants
- Repeat content differences
- Assembly quality improvements

---

### Phase 6: Quality Control (PENDING)
**Planned analyses:**
- BUSCO gene completeness assessment
- QUAST comparison to MoringaV2 reference
- Merqury k-mer based QC
- Assembly completeness metrics

---

### Phase 7: Pangenome Construction (PENDING)
**Method:** Minigraph-Cactus
- Input: All cleaned assemblies (9 HiFi + 3 scaffolded + 12 phased = 24 total)
- Output: Pangenome graph
- Analysis: Core vs accessory genome, structural variants

---

## Directory Structure
```
moringa_pangenome/
  - README.md (This file)
  - rawdata/ → raw_data/ (Raw sequencing data)
    - hifi/ (PacBio HiFi reads, 9 samples)
    - hic/ (Hi-C reads, 9 samples)
  - assemblies/
    - production/ (HiFi-only assemblies, 9 samples)
    - scaffolded/ (Hi-C scaffolded, 3 homozygous)
    - phased/ (Haplotype-phased, 6 heterozygous x 2)
    - cleaned/ (Contamination-filtered, pending)
    - organelles/ (Extracted organellar sequences, pending)
  - analysis_scripts/ (NEW - Haplotype validation workflows)
    - batch_01_gap_estimation_all_samples.sh
    - batch_02_homology_analysis_all_samples.sh
    - batch_03_synteny_analysis_all_samples.sh
    - batch_04_coverage_mapping_all_samples.sh
    - batch_05_kmer_analysis_all_samples.sh
    - run_all_analyses_all_samples.sh
  - haplotype_analysis/ (Analysis results)
    - batch_results/ (Summary reports)
    - Mo-TH-55/ through Mo-TW-13/ (Individual sample results)
  - contamination_screening/ (Contamination detection workflow)
    - references/ (Reference databases)
    - results/ (BLAST results)
    - *.sh (Screening scripts)
  - logs/ (All processing logs)
  - qc/ (Quality control results)
  - heterozygosity_summary.txt (Sample classification)
  - assembly_summary_report.txt (Complete assembly statistics)
  - Quality_Assessment_Report.md (Assembly quality metrics)
```

---

## Software and Versions

### Assembly Tools
- **Hifiasm:** v0.19.8 - Primary HiFi assembly
- **YaHS:** v1.2a.2 - Hi-C scaffolding
- **BWA-MEM:** v0.7.17 - Hi-C read alignment

### Analysis Tools
- **GenomeScope2:** v2.0 - Heterozygosity estimation
- **seqkit:** v2.3.0 - Sequence statistics
- **samtools:** v1.18 - BAM processing
- **BLAST+:** v2.14.0 - Contamination screening & homology analysis
- **minimap2:** v2.24 - Synteny analysis
- **meryl:** v1.3 - K-mer analysis

### Quality Control (Planned)
- **BUSCO:** Gene completeness
- **QUAST:** Assembly comparison
- **Merqury:** K-mer based QC

### Pangenome (Planned)
- **Minigraph-Cactus:** Pangenome construction

---

## Key Results Summary

### Assembly Quality
- **100% gap-free assemblies** (0 Ns in all primary and haplotype assemblies)
- High contiguity (N50: 8-15 Mb)
- Chromosome-scale scaffolds achieved (29-42 Mb)
- Successful haplotype phasing (1.76-1.83x ratio)
- Contamination screening complete (~6-8 Mb organellar removed)

### Haplotype Quality (NEW - January 2026)
- **99.9% homology success rate** - virtually all Hap1 contigs have Hap2 matches
- **Zero collapsed regions** detected in completed samples
- **Size differences = structural variants** (confirmed by BLAST analysis)
- **30-35 large SVs per sample** (contigs with 50-80% coverage)
- **Excellent phasing quality** across all heterozygous samples

### Biological Insights
- **Heterozygosity distribution:** 33% homozygous, 67% heterozygous
- **Mating system:** Mixed mating (selfing and outcrossing)
- **Genome size:** 290-314 Mb (larger than 240 Mb reference)
- **Structural variation:** Size differences between samples (24 Mb range)
- **Haplotype diversity:** 
  - Confirmed by GC% differences and size imbalances
  - Size differences are TRUE biological variation (31-48 Mb between haplotypes)
  - Represents parental genome divergence, not assembly artifacts

### Technical Achievements
- Optimized Hi-C scaffolding workflow
- Automated haplotype phasing pipeline
- Robust contamination screening with proper filtering
- **Multi-method haplotype validation framework (4 independent methods)**
- Comprehensive QC framework
- Reproducible scripts for all steps

---

## Analysis Timeline

- **November 2025:** HiFi assembly (all 9 samples)
- **December 2-3, 2025:** Heterozygosity assessment and sample classification
- **December 3-4, 2025:** Hi-C scaffolding (3 homozygous samples)
- **December 4, 2025:** Haplotype phasing (6 heterozygous samples)
- **December 5, 2025:** Contamination screening
- **January 20-27, 2026:** Haplotype validation complete (5/5 methods)
  - Gap analysis: 100% gap-free
  - Homology analysis: 0 collapsed regions
  - Synteny analysis: 100% alignment
  - Coverage mapping: 88.6% conserved
  - K-mer analysis: 92.1% shared k-mers
- **Next Phase:** Pangenome construction with 15 validated assemblies

---

## Important Notes

### Assembly Size vs Reference
Our assemblies (290-314 Mb) are larger than the MoringaV2 reference (240 Mb). This is **expected and normal** because:
1. Higher assembly completeness
2. Better repeat resolution
3. Structural variants between accessions
4. Different assembly technologies
5. Inclusion of heterozygous regions

After removing ~6-8 Mb of organellar contamination, cleaned assemblies are ~302-304 Mb.

### Haplotype Size Imbalances (NEW)
Haplotype size differences (up to 47.9 Mb) were initially concerning but have been validated as **true biological variation**:
- **NOT phasing artifacts** - homology analysis shows 99.9% matching
- **NOT collapsed regions** - zero contigs with <50% coverage
- **ARE structural variants** - 30-35 large indels/rearrangements per sample
- **Biological interpretation:** Reflects divergence between parental genomes
- **Impact on pangenome:** This variation SHOULD be captured for comprehensive pangenome

### Contamination Screening Philosophy
**What we remove:**
- Contigs that are >90% chloroplast/mitochondria (true organellar fragments)
- Contigs that are >80% bacterial/vector

**What we keep:**
- Large nuclear contigs with small organellar insertions (NUMTs)
- All contigs where organellar sequence is <90% of contig length

This prevents overcorrection while removing genuine contamination.

### Organellar Sequences
Removed organellar contigs are **saved separately** in `assemblies/organelles/` for:
- Chloroplast genome assembly
- Mitochondrial genome assembly
- Phylogenetic studies
- Organellar gene annotation

---

## GitHub Repository

**Repository:** https://github.com/Nfumia/moringa_pangenome (to be created)

**Contents:**
- All analysis scripts and workflows
- Documentation and tutorials
- Summary statistics tables
- Publication-ready figures

**Installation:**
```bash
git clone https://github.com/Nfumia/moringa_pangenome.git
cd moringa_pangenome
conda env create -f environment.yml
conda activate moringa_pangenome
```

---

## Contact and Citation

**Project Lead:** Dr. Nathan Fumia  
**Institution:** Hawaii Agriculture Research Center  
**Date Started:** November 2025  

**Corresponding Author:** [Your email]

---

## References

### Methods Based On:
- Dr. Jingjing Zheng's "Software and workflow for pangenome analyses" document
- RMcontamination workflow (Section I.3.2)
- Minigraph-Cactus pangenome construction (Section I.8)

### Key Publications:
- Hifiasm: Cheng et al. (2021) Nature Methods
- YaHS: Zhou et al. (2023) Bioinformatics
- GenomeScope2: Ranallo-Benavidez et al. (2020) Nature Communications

---

**Last Updated:** February 5, 2026
**Current Phase:** Comprehensive QC (BUSCO running)
**Next Milestone:** Complete QC → Genome Annotation → Pangenome
