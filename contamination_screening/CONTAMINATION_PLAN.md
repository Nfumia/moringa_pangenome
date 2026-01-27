# Moringa Pangenome Contamination Screening Plan

## Priority Focus (Dr. Jingjing Zheng's Recommendations)
1. **Chloroplast sequences** (highest priority)
2. **Mitochondrial sequences** (highest priority)
3. Bacterial contamination
4. Vector sequences
5. Fungal sequences (lower priority)
6. 45S rRNA sequences

## Databases to Download

### 1. Moringa-Specific Organelles
- **Chloroplast:** Moringa oleifera chloroplast complete genome
  - NCBI: NC_056961.1 (164,805 bp)
- **Mitochondria:** Search for Moringa or use closest relative
  - Fallback: Brassicales mitochondrial genomes

### 2. Related Species Organelles (Brassicales Order)
- Arabidopsis thaliana chloroplast (NC_000932.1)
- Arabidopsis thaliana mitochondria (NC_037304.1)
- Carica papaya organelles (if available)

### 3. Bacterial Contaminants
- E. coli K-12 (NC_000913.3)
- Common lab bacteria
- Agrobacterium (plant transformation vector)

### 4. Vector Sequences
- UniVec Core database (NCBI)
- Common cloning vectors

### 5. Fungal Sequences (Optional)
- Common plant fungal pathogens
- Saccharomyces cerevisiae (common lab contaminant)

### 6. 45S rRNA Sequences
- Plant 45S rDNA sequences from closely related species

## Workflow Steps

### Step 1: Download Reference Databases
- Automated download script
- Create combined BLAST database
- Index with makeblastdb

### Step 2: BLAST Screening (Per Colleague: >80% coverage threshold)
- Screen all 9 HiFi primary assemblies
- Output format: tabular with qcovs (query coverage)
- Filter: qcovs > 80 = contamination

### Step 3: Classification and Reporting
- Separate contigs by contamination type:
  - Chloroplast contigs → Save separately (useful!)
  - Mitochondrial contigs → Save separately (useful!)
  - Bacterial/Vector/Fungal → Flag for removal
  
### Step 4: Decision Making
- **Organelles:** Keep separate but don't include in nuclear assembly
- **True contaminants:** Remove entirely
- Generate before/after statistics

### Step 5: Create Cleaned Assemblies
- Filter contaminated contigs
- Save to assemblies/cleaned/
- Update assembly statistics

## Expected Findings

### Likely Contaminants
1. **Chloroplast:** 1-2 contigs, ~150-165 kb total
2. **Mitochondria:** 1-3 contigs, ~300-500 kb total
3. **Bacterial:** Possibly 0-5 contigs (if any)
4. **Vectors:** Unlikely with HiFi data (already high quality)

### Action Plan
- Chloroplast/Mitochondria: **Keep as separate organelle assemblies**
- True contaminants: **Remove from nuclear assembly**
- Report contamination rate per sample

## Success Criteria
- Clear identification of organellar sequences
- <1% contamination rate (excluding organelles)
- Cleaned assemblies ready for pangenome construction
