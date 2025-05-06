# Differential Gene Expression and Functional Enrichment Analysis Basic Workflow

This repository provides a workflow for performing differential gene expression analysis using RNA-seq data, followed by functional enrichment analysis, including Gene Ontology (GO) and KEGG pathway analysis. The workflow is implemented in R using several bioinformatics packages.

## Overview

The workflow involves the following steps:
1. **Data Import and Preprocessing**: Importing RNA-seq count data and sample metadata, subsetting required columns, and preparing the data for analysis.
2. **Differential Gene Expression Analysis**: Using the DESeq2 package to identify differentially expressed genes.
3. **Functional Enrichment Analysis**: Conducting GO and KEGG pathway enrichment analysis on the significant genes.
4. **Visualization**: Visualizing the results using various plots, including MA plots, barplots, dotplots, and cnetplots.

## Workflow Details

### 1. Data obtain and Preprocessing

- Raw RNA-Seq counts data was downloaded from NCBI [https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE52778]
- Import RNA-seq count data from local disk using the `rio` package.
- Subset the data to include only relevant columns.
- Replace row names with gene IDs and remove the gene ID column.
- Import sample metadata and ensure the order of columns matches the count data.

### 2. Differential Gene Expression Analysis

- Create a DESeq2 dataset from the count data and metadata.
- Pre-filter low-count genes[Not compulsory].
- Run the DESeq2 pipeline to identify differentially expressed genes[You can run DESeq2 on raw RNA-Seq count data].
- Extract significant results based on an adjusted p-value threshold.

### 3. Functional Enrichment Analysis

- **Gene Ontology (GO) Analysis**:
  - Perform GO enrichment analysis for Biological Processes (BP) and Cellular Components (CC) using `clusterProfiler`.
  - Visualize the results with barplots.

- **KEGG Pathway Analysis**:
  - Conduct KEGG pathway enrichment analysis using `clusterProfiler`.
  - Visualize enriched pathways with dotplots and cnetplots.

### 4. Visualization

- MA plots to visualize differentially expressed genes.
- ![MA plot](https://github.com/user-attachments/assets/50f5ccb4-4534-49e4-a17b-79f20553e7c9)

  **Barplots and dotplots for enriched GO terms and KEGG pathways**.
  
- Bar plot of enriched terms [BP-Biological Process enrichment]:
- ![BP enrichment](https://github.com/user-attachments/assets/3a59946b-4860-46e5-a25f-1fc62e627ab1)
  
  
- Bar plot of enriched terms [CC- Cellular component enrichment]:
- ![CC GO enrichment](https://github.com/user-attachments/assets/8ced73ed-c55f-4a4c-ad3f-80f49aabfc0c)
  
  
- Dot Plot of pathway enriched terms:
  
  ![KEGG](https://github.com/user-attachments/assets/367a674e-cf73-4f34-9dd1-4518f340bcf9)
  


- Cnetplots to explore the relationship between genes and enriched pathways.
**NOTE:** Both `barplot()` and `dotplot()` display only the most significant or selected enriched terms, but one may want to see which genes are involved in these terms. To address this complexity where a gene may belong to multiple annotation categories, the `cnetplot()` function was developed to capture and display these complex associations. 

- Narrowing down to 5 highly significantly enriched pathway terms
  
  ![network of terms_1](https://github.com/user-attachments/assets/75a62258-1322-43ef-baeb-be47cf73ba50)

  Differentially expressed genes associated with the 5 enriched pathway terms
  ![genes involved in enriched network](https://github.com/user-attachments/assets/24f1b8ea-ad5d-4856-aaad-c77f636c4b2a)

  

  ![circular plot of enriched pathway genes](https://github.com/user-attachments/assets/74fa1ef0-4920-4bf8-b3b1-63615465066b)


## Installation

To run this workflow, you need to install the following R packages:

```r
install.packages(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "DOSE", "rio"))
