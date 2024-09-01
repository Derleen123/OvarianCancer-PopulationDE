# Load necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(recount)
library(ggplot2)
library(dplyr)
library(GOSemSim)
library(VennDiagram)

# Function to query and download TCGA data
query_tcga_data <- function() {
  query <- GDCquery(
    project = "TCGA-OV",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor",
    access = "open"
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  return(data)
}

# Download the TCGA data
tcga_data <- query_tcga_data()

# Define the relevant early stages
early_stages <- c("stage ic", "stage iia", "stage iib", "stage iic")

# Filter data for Black/African-American and Asian patients
data_black <- tcga_data[, 
                        !is.na(colData(tcga_data)$race) & 
                          !is.na(colData(tcga_data)$figo_stage) & 
                          colData(tcga_data)$race == "black or african american" & 
                          tolower(colData(tcga_data)$figo_stage) %in% early_stages]

data_asian <- tcga_data[, 
                        !is.na(colData(tcga_data)$race) & 
                          !is.na(colData(tcga_data)$figo_stage) & 
                          colData(tcga_data)$race == "asian" & 
                          tolower(colData(tcga_data)$figo_stage) %in% early_stages]

# Download and prepare GTEx data for normal ovarian tissue
gtex_study <- "SRP012682"
ovary_gtex <- recount::download_study(gtex_study, type = "rse-gene")
load(file.path(gtex_study, "rse_gene.Rdata"))

ovary_samples <- rse_gene[, colData(rse_gene)$smts == "Ovary"]
ovary_counts <- assay(ovary_samples)

# Randomly select 5 normal samples
set.seed(42)
normal_data <- SummarizedExperiment(assays = list(counts = ovary_counts[, sample(ncol(ovary_counts), 5)]))

# Get common genes between tumor and normal data
common_genes <- intersect(rownames(assay(data_black)), rownames(assay(normal_data)))

# Subset both datasets to include only common genes
counts_black <- assay(data_black)[common_genes, ]
counts_asian <- assay(data_asian)[common_genes, ]
counts_normal <- assay(normal_data)[common_genes, ]

# Combine the datasets
combined_counts_black <- cbind(counts_black, counts_normal)
combined_counts_asian <- cbind(counts_asian, counts_normal)

# Function for data preprocessing using DESeq2
preprocess_data <- function(counts, sample_labels) {
  sample_type <- factor(sample_labels)
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = data.frame(sample_type = sample_type),
                                design = ~ sample_type)
  
  dds <- dds[rowSums(counts(dds)) >= 10,]  # Filter low-count genes
  dds <- estimateSizeFactors(dds)  # Normalize data
  
  return(dds)
}

# Preprocess data for both populations
dds_black <- preprocess_data(combined_counts_black, c(rep("Tumor", ncol(counts_black)), rep("Normal", ncol(counts_normal))))
dds_asian <- preprocess_data(combined_counts_asian, c(rep("Tumor", ncol(counts_asian)), rep("Normal", ncol(counts_normal))))

# Function to perform differential expression analysis
perform_de_analysis <- function(dds) {
  dds <- DESeq(dds)
  results <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))
  return(results)
}

# Perform the analysis
de_results_black <- perform_de_analysis(dds_black)
de_results_asian <- perform_de_analysis(dds_asian)

# Process Black/African-American data
clean_ensembl_ids_black <- sub("\\..*", "", rownames(de_results_black))
entrez_ids_black <- mapIds(org.Hs.eg.db, keys = clean_ensembl_ids_black, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
de_results_black$EntrezID <- entrez_ids_black
de_results_black_df <- as.data.frame(de_results_black)

# Filter significant genes for Black/African-American
sig_genes_black <- de_results_black_df %>%
  filter(!is.na(EntrezID)) %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) >= 1) %>%
  pull(EntrezID) %>%
  unique()

# Perform enrichment analysis for Black/African-American
if (length(sig_genes_black) > 0) {
  go_enrich_black <- enrichGO(gene = sig_genes_black, 
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", 
                              keyType = "ENTREZID", 
                              pAdjustMethod = "BH", 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.2, 
                              readable = TRUE)
  
  kegg_enrich_black <- enrichKEGG(gene = sig_genes_black, 
                                  organism = 'hsa', 
                                  keyType = "ncbi-geneid", 
                                  pAdjustMethod = "BH", 
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.2)
} else {
  warning("No significant genes found for enrichment analysis with the specified thresholds.")
  go_enrich_black <- NULL
  kegg_enrich_black <- NULL
}

# Process Asian data
clean_ensembl_ids_asian <- sub("\\..*", "", rownames(de_results_asian))
entrez_ids_asian <- mapIds(org.Hs.eg.db, keys = clean_ensembl_ids_asian, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
de_results_asian$EntrezID <- entrez_ids_asian
de_results_asian_df <- as.data.frame(de_results_asian)

# Filter significant genes for Asian
sig_genes_asian <- de_results_asian_df %>%
  filter(!is.na(EntrezID)) %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) >= 1) %>%
  pull(EntrezID) %>%
  unique()

# Perform enrichment analysis for Asian
if (length(sig_genes_asian) > 0) {
  go_enrich_asian <- enrichGO(gene = sig_genes_asian, 
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", 
                              keyType = "ENTREZID", 
                              pAdjustMethod = "BH", 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.2, 
                              readable = TRUE)
  
  kegg_enrich_asian <- enrichKEGG(gene = sig_genes_asian, 
                                  organism = 'hsa', 
                                  keyType = "ncbi-geneid", 
                                  pAdjustMethod = "BH", 
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.2)
} else {
  warning("No significant genes found for enrichment analysis with the specified thresholds.")
  go_enrich_asian <- NULL
  kegg_enrich_asian <- NULL
}

# Visualization and analysis
# GO and KEGG enrichment plots for Black/African-American
if (!is.null(go_enrich_black)) {
  barplot(go_enrich_black, showCategory = 10, title = "GO Enrichment - Black/African-American")
  dotplot(go_enrich_black, showCategory = 10, title = "GO Enrichment - Black/African-American")
  emapplot(go_enrich_black)
  cnetplot(go_enrich_black)
}

if (!is.null(kegg_enrich_black)) {
  barplot(kegg_enrich_black, showCategory = 10, title = "KEGG Enrichment - Black/African-American")
  dotplot(kegg_enrich_black, showCategory = 10, title = "KEGG Enrichment - Black/African-American")
}

# GO and KEGG enrichment plots for Asian
if (!is.null(go_enrich_asian)) {
  barplot(go_enrich_asian, showCategory = 10, title = "GO Enrichment - Asian")
  dotplot(go_enrich_asian, showCategory = 10, title = "GO Enrichment - Asian")
  emapplot(go_enrich_asian)
  cnetplot(go_enrich_asian)
}

if (!is.null(kegg_enrich_asian)) {
  barplot(kegg_enrich_asian, showCategory = 10, title = "KEGG Enrichment - Asian")
  dotplot(kegg_enrich_asian, showCategory = 10, title = "KEGG Enrichment - Asian")
}

# Save results
write.csv(de_results_black_df, file = "black_african_american_de_results.csv")
write.csv(de_results_asian_df, file = "asian_de_results.csv")

if (!is.null(go_enrich_black)) write.csv(as.data.frame(go_enrich_black), file = "black_african_american_go_enrichment.csv")
if (!is.null(kegg_enrich_black)) write.csv(as.data.frame(kegg_enrich_black), file = "black_african_american_kegg_enrichment.csv")
if (!is.null(go_enrich_asian)) write.csv(as.data.frame(go_enrich_asian), file = "asian_go_enrichment.csv")
if (!is.null(kegg_enrich_asian)) write.csv(as.data.frame(kegg_enrich_asian), file = "asian_kegg_enrichment.csv")

# Heatmaps
if (length(sig_genes_black) > 0) {
  sig_gene_expression_black <- assay(dds_black)[rownames(dds_black) %in% sig_genes_black, ]
  pheatmap(sig_gene_expression_black, 
           scale = "row", 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           show_rownames = FALSE, 
           show_colnames = TRUE, 
           main = "Heatmap - Significant Genes (Black/African-American)")
}

if (length(sig_genes_asian) > 0) {
  sig_gene_expression_asian <- assay(dds_asian)[rownames(dds_asian) %in% sig_genes_asian, ]
  pheatmap(sig_gene_expression_asian, 
           scale = "row", 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           show_rownames = FALSE, 
           show_colnames = TRUE, 
           main = "Heatmap - Significant Genes (Asian)")
}

# Comparative analysis
if (!is.null(go_enrich_black) && !is.null(go_enrich_asian)) {
  sig_go_terms_black <- go_enrich_black@result$Description
  sig_go_terms_asian <- go_enrich_asian@result$Description
  
  unique_go_black <- setdiff(sig_go_terms_black, sig_go_terms_asian)
  unique_go_asian <- setdiff(sig_go_terms_asian, sig_go_terms_black)
  shared_go <- intersect(sig_go_terms_black, sig_go_terms_asian)
  
  venn.diagram(
    x = list(
      Black_African_American = sig_go_terms_black,
      Asian = sig_go_terms_asian
    ),
    category.names = c("Black/African-American", "Asian"),
    filename = "venn_go_terms.png",
    output = TRUE
  )
  
  go_sim <- mgoSim(sig_go_terms_black, sig_go_terms_asian, semData = godata("org.Hs.eg.db", ont = "BP"), measure = "Wang")
}

if (!is.null(kegg_enrich_black) && !is.null(kegg_enrich_asian)) {
  sig_kegg_paths_black <- kegg_enrich_black@result$Description
  sig_kegg_paths_asian <- kegg_enrich_asian@result$Description
  
  unique_kegg_black <- setdiff(sig_kegg_paths_black, sig_kegg_paths_asian)
  unique_kegg_asian <- setdiff(sig_kegg_paths_asian, sig_kegg_paths_black)
  shared_kegg <- intersect(sig_kegg_paths_black, sig_kegg_paths_asian)
  
  venn.diagram(
    x = list(
      Black_African_American = sig_kegg_paths_black,
      Asian = sig_kegg_paths_asian
    ),
    category.names = c("Black/African-American", "Asian"),
    filename = "venn_kegg_paths.png",
    output = TRUE
  )
}

# Combined visualizations
if (!is.null(go_enrich_black) && !is.null(go_enrich_asian)) {
  go_results_combined <- rbind(
    data.frame(Term = go_enrich_black@result$Description,
               pvalue = go_enrich_black@result$pvalue,
               Population = "Black/African-American"),
    data.frame(Term = go_enrich_asian@result$Description,
               pvalue = go_enrich_asian@result$pvalue,
               Population = "Asian")
  )
  
  ggplot(go_results_combined, aes(x = reorder(Term, -pvalue), y = -log10(pvalue), fill = Population)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "GO Enrichment - Top Terms by Population",
         x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
}

if (!is.null(kegg_enrich_black) && !is.null(kegg_enrich_asian)) {
  kegg_results_combined <- rbind(
    data.frame(Pathway = kegg_enrich_black@result$Description,
               pvalue = kegg_enrich_black@result$pvalue,
               Population = "Black/African-American"),
    data.frame(Pathway = kegg_enrich_asian@result$Description,
               pvalue = kegg_enrich_asian@result$pvalue,
               Population = "Asian")
  )
  
  ggplot(kegg_results_combined, aes(x = reorder