# Differential Gene Expression between HER2 Sub-type Breast Cancer

This project involves differential gene expression analysis of HER2 sub-type breast cancer using RNA sequencing and copy number alteration (CNA) data from the TCGA dataset. The workflow includes data pre-processing, statistical analysis, pathway enrichment analysis, and survival analysis.

## Steps:

### 1. **Data Loading**
   - **Clinical Data**: Patient data is loaded from `data_clinical_patient.txt`.
   - **RNA-Seq Data**: mRNA expression data is loaded from `data_mrna_seq_v2_rsem.txt`.
   - **CNA Data**: Copy number alteration data is loaded from `data_cna.txt`.

### 2. **Data Pre-processing**
   - **ERBB2 Identification**: The `ERBB2` gene is identified from the CNA data and its counts across patients are extracted.
   - **Patient-ID Matching**: RNA-seq data and CNA data are matched with patient IDs from the clinical data to ensure consistency across datasets.

### 3. **HER2 Status Determination**
   - HER2 status is assigned based on the ERBB2 gene copy number: Patients with a positive ERBB2 count are considered HER2 amplified (assigned value 1), while others are non-amplified (assigned value 0).

### 4. **Differential Gene Expression Analysis**
   - **DESeq2**: The `DESeq2` package is used to analyze differential gene expression between HER2 amplified and non-amplified groups.
   - **PCA Plot**: Principal Component Analysis (PCA) is performed to visualize the distribution of gene expression data between the two HER2 status groups.
   - **Top Differential Genes**: The top 10 most differentially expressed genes based on adjusted p-value are displayed.

### 5. **Pathway Enrichment Analysis**
   - **Gene Ontology (GO)**: Over-represented biological processes are identified for both upregulated and downregulated genes using the `clusterProfiler` package.
   - **KEGG Pathways**: Pathway enrichment for both over-expressed and under-expressed genes is performed using KEGG pathways.
   - **Reactome Pathways**: Enrichment analysis is also performed for Reactome pathways.

### 6. **Heatmap Visualization**
   - A heatmap of the top 20 most differentially expressed genes is created to visualize gene expression patterns across samples, grouped by HER2 status.

### 7. **Survival Analysis**
   - **Cox Proportional Hazard Model**: The survival data (Overall Survival Status and Months) are used to perform a Cox proportional hazard model to investigate the correlation between HER2 status and survival outcomes.
   - **Survival Curves**: Survival curves are plotted for HER2 amplified and non-amplified groups.

### 8. **Output**
   - The differentially expressed genes (DEGs) with adjusted p-values < 0.05 are saved in an Excel file (`DEG.xlsx`).

## R Packages Used:
- **DESeq2**: Differential expression analysis.
- **clusterProfiler**: Gene Ontology and pathway enrichment analysis.
- **ReactomePA**: Reactome pathway enrichment.
- **survival** and **glmnet**: Survival analysis using Cox proportional hazard models.
- **pheatmap**: Visualizing gene expression patterns through heatmaps.

## Data and File Requirements:
1. **Clinical Data**: `data_clinical_patient.txt`
2. **RNA-Seq Data**: `data_mrna_seq_v2_rsem.txt`
3. **CNA Data**: `data_cna.txt`
4. **R Packages**: `DESeq2`, `clusterProfiler`, `ReactomePA`, `survival`, `glmnet`, `pheatmap`

## Summary:
This analysis provides a comprehensive view of the gene expression landscape in HER2 positive vs. negative breast cancer subtypes, identifies potential biomarkers, and evaluates the survival impact of HER2 status. The project combines differential expression analysis, pathway enrichment, and survival analysis to gain insights into the biological mechanisms underlying HER2-driven breast cancer.
