suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(limma)
})
#1.Configuration & Cohort Dedinition
cancer_project    <- "TCGA-BRCA"
cancer_label      <- "Breast Cancer (BRCA)"
n_samples_per_grp <- 20 
message(sprintf("Initiating data query for project: %s", cancer_project))
# 2. Query GDC Portal & Download Data
query_tcga <- GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
results_df <- getResults(query_tcga)
tumor_barcodes <- results_df$cases[results_df$sample_type == "Primary Tumor"][1:n_samples_per_grp]
normal_barcodes <- results_df$cases[results_df$sample_type == "Solid Tissue Normal"][1:n_samples_per_grp]
query_sub <- GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(tumor_barcodes, normal_barcodes)
)
GDCdownload(query_sub, method = "api", files.per.chunk = 10)
tcga_se <- GDCprepare(query_sub)
# 3. Extract Assay Data & Format Metadata
counts_mat <- assay(tcga_se, "unstranded")
gene_info <- as.data.frame(rowData(tcga_se))
sample_metadata <- as.data.frame(colData(tcga_se)) %>%
  transmute(
    barcode = barcode,
    sample_type = sample_type,
    Condition = factor(
      ifelse(definition == "Primary solid Tumor", "Tumour", "Normal"),
      levels = c("Normal", "Tumour")
    )
  )
#4. Low-Count Filtering & Voom Normalization
dge <- DGEList(counts = counts_mat, group = sample_metadata$Condition)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, plot = FALSE)
log_expr <- v$E
#5 save Processed artifacts
write.csv(sample_metadata, "data/sample_metadata.csv", row.names = FALSE)
write.csv(gene_info, "data/gene_info.csv", row.names = TRUE)
saveRDS(log_expr, "data/normalised_log_expression.rds")
saveRDS(v, "data/voom_object.rds")
message("Preprocessing complete! Cleaned data saved to dara/")
