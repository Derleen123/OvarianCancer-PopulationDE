# Early-Stage Ovarian Cancer RNA-Seq Analysis in Diverse Populations

## Overview

This repository contains R scripts designed for the analysis of RNA-Seq data from early-stage Ovarian Serous Cystadenocarcinoma (OV) samples, sourced from The Cancer Genome Atlas (TCGA). The primary objective is to uncover distinctive gene expression patterns by performing differential expression analysis on two specific populations: Black/African-American and Asian. This analysis contrasts these populations with normal ovarian tissue samples, offering insights into population-specific molecular characteristics.

## Key Features

* **Data Acquisition**: Retrieval of RNA-Seq data from TCGA and normal tissue data from GTEx.
* **Data Preprocessing**: Normalization and filtration of low-expression genes.
* **Population-Specific Differential Expression Analysis**: Identification of differentially expressed genes (DEGs) in Black/African-American and Asian populations.
* **Pathway Enrichment Analysis**: GO and KEGG enrichment analyses to understand the biological pathways associated with DEGs.
* **Comparative Analysis**: Cross-population comparison of enriched pathways and GO terms.
* **Visualization**: Creation of heatmaps, enrichment plots, and Venn diagrams to illustrate findings.
* **Export of Results**: Generation of CSV files and visual outputs summarizing the analysis.

## Prerequisites

Before running the scripts, ensure you have R installed along with the necessary packages:

* **Core Packages**: 
  * `TCGAbiolinks`
  * `SummarizedExperiment`
  * `DESeq2`
  * `clusterProfiler`
  * `org.Hs.eg.db`
  * `AnnotationDbi`
  * `pheatmap`
  * `recount`
  * `ggplot2`
  * `dplyr`
  * `GOSemSim`
  * `VennDiagram`

Install these packages using the following R code:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
                       "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", 
                       "pheatmap", "recount", "ggplot2", "dplyr", 
                       "GOSemSim", "VennDiagram"))
```

## Usage

To perform the analysis, run the main script available in this repository:

* Run the main script: https://github.com/Derleen123/OvarianCancer-PopulationDE/blob/main/ovarian_cancer_analysis.R

## Workflow

1. **Data Acquisition**:
   * Retrieve RNA-Seq data for early-stage OV samples from Black/African-American and Asian populations using TCGA.
   * Download normal ovarian tissue data from the GTEx database.

2. **Data Preprocessing**:
   * Normalize the RNA-Seq data.
   * Filter out genes with low expression levels.

3. **Differential Expression Analysis**:
   * Identify DEGs by comparing each population's cancer samples to normal ovarian tissue samples.

4. **Pathway Enrichment Analysis**:
   * Conduct Gene Ontology (GO) and KEGG pathway enrichment analyses on the identified DEGs for each population.

5. **Comparative Analysis**:
   * Compare the enriched pathways and GO terms between Black/African-American and Asian populations.

6. **Visualization**:
   * Generate visual representations, including heatmaps, enrichment plots, and Venn diagrams, to illustrate the findings.

7. **Results Export**:
   * Save the analysis results and visualizations for further interpretation.

## Output

The script generates the following output files:

* `black_african_american_de_results.csv`: DEGs for the Black/African-American population.
* `asian_de_results.csv`: DEGs for the Asian population.
* `black_african_american_go_enrichment.csv`: GO enrichment results for the Black/African-American population.
* `black_african_american_kegg_enrichment.csv`: KEGG enrichment results for the Black/African-American population.
* `asian_go_enrichment.csv`: GO enrichment results for the Asian population.
* `asian_kegg_enrichment.csv`: KEGG enrichment results for the Asian population.
* Various PNG files for visualizations such as heatmaps, enrichment plots, and Venn diagrams.

## Acknowledgments

This project seeks to advance our understanding of early-stage ovarian cancer across different populations. We acknowledge the following resources:

* **Data Providers**:
   * The Cancer Genome Atlas (TCGA) https://www.cancer.gov/ccg/research/genome-sequencing/tcga
   * Genotype-Tissue Expression (GTEx) Project https://gtexportal.org/home/
* **Course**: This project was developed as part of the *Advanced Genomics Course: Bioinformatics for Cancer Biology* from HackBio. For more details, visit HackBio's course page: https://course.thehackbio.com
