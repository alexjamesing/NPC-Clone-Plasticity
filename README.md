# NPC-Clone-Plasticity
Adaptive Clonal Evolution in the Developing Mouse Brain 




## Related data
 * The raw data reported in this work can be found in [GEO]() GSE281055
 * The related data is available in **~/data**:


## Script Overview

### 1. Data Preprocessing and Clustering
- **Description**:  
  The raw sequencing data were processed using the CellRanger pipeline (v7.0.0, 10x Genomics), aligning reads to a custom genome reference. Afterward, raw count matrices were processed to remove doublet cells using Scrublet. Seurat (v4.3.0) was used for downstream analysis, including normalization, identification of highly variable genes, PCA, and batch integration using Harmony. The integrated data were clustered using the Leiden algorithm.

- **Scripts**:  
  - `0_Doublets.R`: Identifies and filters out doublet cells.  
  - `0_doublet_prediction.py`: Python script for doublet predictions.  
  - `1_1_Preprocess.R`: Performs preprocessing steps such as normalization and gene filtering.  
  - `1_2_Integration.R`: Integrates datasets using Harmony.

---

### 2. Meta-Cells and Cell Type Annotation
- **Description**:  
  High-resolution clustering was performed to generate meta-cells, improving the accuracy of cell type annotations. Label transfer using Canonical Correlation Analysis (CCA) was employed to annotate cells, with the final annotations validated by canonical marker genes.

- **Scripts**:  
  - `1_3_Metacells_Hclu.R`: Generates meta-cells for cell type annotation.  
  - `1_4_Annotation.R`: Annotates cells using CCA and marker genes.

---

### 3. Differential Expression Analysis
- **Description**:  
  Differentially expressed genes (DEGs) were identified between conditions. This analysis was used to pinpoint genes and pathways driving cellular differences.

- **Script**:  
  - `2_1_DE.R`: Identifies DEGs between conditions.

---

### 4. Cell Type Proportion Statistical Analysis
- **Description**:  
  Cell type compositions across biological replicates were compared using a Generalized Linear Model (GLM) with a binomial distribution, followed by ANOVA for statistical significance.

- **Script**:  
  - `2_2_cell_prop_Fig5A_S4.R`: Analyzes and plots cell type proportions.

---

### 5. S Phase Score Analysis
- **Description**:  
  Cell cycle scores were assigned using Seurat’s `CellCycleScoring` function, and the differences between conditions were tested using the Kolmogorov-Smirnov test.

- **Script**:  
  - `2_3_CellCycle_Fig5B.R`: Performs cell cycle scoring and plots results.

---

### 6. Pathway Enrichment Analysis (GSEA)
- **Description**:  
  Gene Set Enrichment Analysis (GSEA) was performed using `fgsea` (v1.30.0) with the Hallmark gene sets for *Mus musculus*. The analysis compared conditions to identify enriched pathways, with input derived from a Wilcoxon rank-sum test.

- **Script**:  
  - `2_4_GSEA_Fig5C.R`: Performs GSEA and visualizes pathway enrichment results.

---

### 7. Regulon Analysis
- **Description**:  
  SCENIC was used to infer gene regulatory networks (GRNs) and identify conserved regulons using the mm10 motif database and transcription factors. Regulon activity differences between conditions were analyzed. Regulons with significant differences (p ≤ 0.01) between NesCre::R26-DTA and R26-DTA conditions were identified via permutation tests. Associated pathways were analyzed using Hallmark gene sets and Fisher’s exact test (Odds Ratio ≥1, p ≤ 0.01). Regulons with low gene set overlap (<10th percentile) were excluded. Kolmogorov-Smirnov tests on module scores validated pathway contributions, retaining only significant results (p ≤ 0.01). **Use pySCENIC output**

- **Script**:  
  - `2_5_RegulAssociPathways_FigS4.R`: Analyzes regulon activity and explores associated pathways.



    
