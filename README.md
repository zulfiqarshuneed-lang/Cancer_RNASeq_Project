# Identification of Candidate Biomarkers Through Differential Gene Expression Analysis of Breast Cancer Using R and Bioconductor

## Dataset

This project used the publicly available breast cancer transcriptomic dataset **GSE42568** obtained from the NCBI Gene Expression Omnibus (GEO).

- **Repository:** NCBI Gene Expression Omnibus (GEO)
- **Accession:** GSE42568
- **Species:** *Homo sapiens*
- **Disease:** Breast Cancer
- **Data Type:** Gene expression microarray
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array

The dataset contains transcriptomic profiles from normal breast tissue and breast cancer samples, enabling comparative differential gene expression analysis to identify genes associated with breast cancer progression.

## Biological Question

Which genes are significantly differentially expressed between normal breast tissue and breast cancer samples, and what biological insights can be obtained from these transcriptomic alterations?

This project aims to identify dysregulated transcripts associated with breast cancer and provide computational evidence for potential biomarkers and therapeutic targets.

## Methodology

The analysis was performed using R and Bioconductor packages following a standard transcriptomic differential expression workflow.

### Workflow

1. Downloaded the GSE42568 dataset using GEOquery.
2. Extracted gene expression and phenotype information.
3. Constructed the experimental design matrix for normal and breast cancer samples.
4. Performed differential gene expression analysis using the limma package.
5. Applied empirical Bayes moderation (eBayes) to improve statistical estimation.
6. Corrected for multiple hypothesis testing using the Benjamini–Hochberg False Discovery Rate (FDR).
7. Identified significantly differentially expressed transcripts.
8. Visualized results using:
   - Volcano Plot
   - Principal Component Analysis (PCA)
   - Hierarchical Clustering Heatmap
9. Investigated highly dysregulated genes associated with breast cancer biology.

## Software and Packages

The analysis was performed using R (Bioconductor) with the following packages:

- GEOquery
- limma
- edgeR
- EnhancedVolcano
- pheatmap
- clusterProfiler
- org.Hs.eg.db

Additional packages were used for visualization and statistical analysis.

## Dataset Summary

| Metric | Value |
|--------|------:|
| Total transcripts analysed | **54,675** |
| Significant transcripts (FDR < 0.05) | **16,803** |
| Upregulated transcripts | **8,240** |
| Downregulated transcripts | **8,563** |

The results indicate widespread transcriptomic alterations between normal breast tissue and breast cancer samples, highlighting extensive molecular dysregulation associated with tumor development.

## Top Differentially Expressed Transcripts

| Rank | Transcript (Probe ID) | logFC | Adjusted p-value (FDR) |
| ---: | --------------------- | ----: | ---------------------: |
| 1 | 237351_at | -4.16 | 7.46 × 10⁻³⁸ |
| 2 | 214505_s_at | -5.37 | 6.08 × 10⁻³⁷ |
| 3 | 89977_at | -4.37 | 1.33 × 10⁻³⁶ |
| 4 | 230595_at | -3.58 | 1.43 × 10⁻³⁶ |
| 5 | 214582_at | -4.41 | 1.43 × 10⁻³⁶ |
| 6 | 223764_x_at | -3.51 | 4.86 × 10⁻³⁶ |
| 7 | 219059_s_at | -4.44 | 1.06 × 10⁻³⁵ |
| 8 | 235670_at | -3.59 | 1.93 × 10⁻³⁵ |
| 9 | 206930_at | -4.49 | 2.00 × 10⁻³⁵ |
| 10 | 221295_at | -3.93 | 2.98 × 10⁻³⁵ |

## Key Findings

- Analysed 54,675 transcript probes from public breast cancer transcriptomic data.
- Identified 16,803 significantly differentially expressed transcripts after false discovery rate correction.
- Detected 8,240 upregulated and 8,563 downregulated transcripts in breast cancer samples.
- Generated volcano plots, PCA visualizations, and clustered heatmaps to characterize transcriptomic differences.
- Observed substantial molecular differences between normal breast tissue and breast cancer, supporting the role of transcriptomic profiling in identifying disease-associated genes and candidate biomarkers.

## Figures

### Volcano Plot

![Enh Volcano Plot] [<img width="800" height="600" alt="Enhanced_Volcano_Plot" src="https://github.com/user-attachments/assets/9adf5604-08c4-4315-9b47-7e66712c8bf4" />]
![Volcano Plot] [<img width="306" height="523" alt="Volcano_Plot" src="https://github.com/user-attachments/assets/6848fe35-9d8f-46e9-83fd-2344243c085c" />]

---

### Principal Component Analysis (PCA)

![PCA] [<img width="534" height="523" alt="PCA_Plot" src="https://github.com/user-attachments/assets/8b862434-ab1b-4e3d-8057-8fad7fa57818" />]

---

### Heatmap of Differentially Expressed Transcripts

![Heatmap] [<img width="306" height="523" alt="Heat_Map" src="https://github.com/user-attachments/assets/822d1aa2-6c9d-47db-aa99-cf8405edfa2b" />]
