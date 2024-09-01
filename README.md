# Early-Stage Ovarian Cancer RNA-Seq Analysis in Diverse Populations

## Overview

This repository contains R scripts for analyzing early-stage Ovarian Serous Cystadenocarcinoma (OV) RNA-Seq data from The Cancer Genome Atlas (TCGA). The project aims to identify unique gene expression patterns in early-stage ovarian cancers by conducting differential expression analysis for two human populations: Black/African-American and Asian. The analysis compares these populations to normal ovarian tissue samples.

## Key Features

- Data acquisition from TCGA and GTEx
- Data preprocessing and normalization
- Population-specific differential expression analysis
- Pathway enrichment analysis (GO and KEGG)
- Comparative analysis between populations
- Visualization of results
- Export of analysis results

## Prerequisites

Ensure you have R (version 4.0.0 or higher) installed along with the following packages:

- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- AnnotationDbi
- pheatmap
- recount
- ggplot2
- dplyr
- GOSemSim
- VennDiagram

Install the required packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

BiocManager::install(c(
    "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "clusterProfiler",
    "org.Hs.eg.db", "AnnotationDbi", "pheatmap", "recount", "ggplot2",
    "dplyr", "GOSemSim", "VennDiagram"
))

## Usage

To perform the analysis, run the main script available in this repository:

* Run the main script: `ovarian_cancer_analysis.R`

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
   *The Cancer Genome Atlas (TCGA))
   * Genotype-Tissue Expression (GTEx) Project
* **Course**: This project was developed as part of the *Advanced Genomics Course: Bioinformatics for Cancer Biology* from HackBio. For more details, visit HackBio's course page.


