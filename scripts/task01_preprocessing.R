suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(limma)
})
project_brca <- "TCGA-BRCA"
project_luad <- "TCGA-LUAD"
label_brca <- "Breast Cancer(BRCA)"
label_luad <- "Lung Adenocarcinoma(LUAD)"
n_samples_per_grp <- 20
query_brca <- GDCquery(
  project= project_brca,
  data.category = "Transcriptome Profiling",
  data.type ="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)
query_luad <- GDCquery(
  project= project_luad,
  data.category = "Transcriptome Profiling",
  data.type ="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)
results_brca <- getResults(query_brca)
results_luad <- getResults(query_luad)
barcodes_brca <- results_brca$cases[results_brca$sample_type == "Primary Tumor"][1:n_samples_per_grp]
barcodes_luad <- results_luad$cases[results_luad$sample_type == "Primary Tumor"][1:n_samples_per_grp]
query_sub_brca <- GDCquery(
  project= project_brca,
  data.category = "Transcriptome Profiling",
  data.type ="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = barcodes_brca
)
query_sub_luad <- GDCquery(
  project= project_luad,
  data.category = "Transcriptome Profiling",
  data.type ="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = barcodes_luad
)
GDCdownload(query_sub_brca, method="api", files.per.chunk= 10)
GDCdownload(query_sub_luad, method="api", files.per.chunk= 10)
se_brca <- GDCprepare(query_sub_brca)
se_luad <-GDCprepare(query_sub_luad)
dim(se_brca)
dim(se_luad)
brca_counts <- assay(se_brca,"unstranded")
luad_counts <- assay(se_luad,"unstranded")
common_genes <- intersect(
  rownames(brca_counts),
  rownames(luad_counts)
)
brca_counts_sub <- brca_counts[common_genes, ]
luad_counts_sub <- luad_counts[common_genes, ]
combined_counts <- cbind(brca_counts_sub, luad_counts_sub)
gene_info <- as.data.frame(rowData(se_brca))
gene_info <- gene_info[common_genes, , drop = FALSE]
metadata_brca <- as.data.frame(colData(se_brca)) %>%
  transmute(
    barcode= barcode,
    sample_type= sample_type,
    Cancer_Type ="BRCA"
  )
metadata_luad <- as.data.frame(colData(se_luad)) %>%
  transmute(
    barcode=barcode,
    sample_type=sample_type,
    Cancer_Type="LUAD"
  )
merged_metadata <- bind_rows(metadata_brca,metadata_luad)
merged_metadata$Cancer_Type <- factor(
  merged_metadata$Cancer_Type,
  levels = c("BRCA", "LUAD")
)
identical(
  colnames(combined_counts),
  merged_metadata$barcode
)
dge <- DGEList(
  counts = combined_counts,
  group = merged_metadata$Cancer_Type
)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ Cancer_Type, data = merged_metadata)
v <- voom(dge, design = design, plot = FALSE)
log_expr <- v$E
dir.create("data/task01",recursive=TRUE, showWarnings=FALSE)
write.csv(
  merged_metadata,
  "data/task01/smaple_metadata.csv",
  row.names= FALSE
)
write.csv(
  gene_info,
  "data/task01/gene_info.csv",
  row.names=FALSE
)
saveRDS(log_expr,"data/task01/normalised_log_expression.rds")
saveRDS(v,"data/task01/voom_object.rds")
message("Preprocessing complete!")
