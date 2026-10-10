suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(matrixStats)
})
sample_metadata <- read.csv("data/task01/sample_metadata.csv")
gene_info <- read.csv("data/task01/gene_info.csv")
colnames(gene_info)
head(gene_info[, 1:3])
log_expr <- readRDS("data/task01/normalised_log_expression.rds")
v <- readRDS("data/task01/voom_object.rds")
sample_metadata$Cancer_Type <-factor(
  sample_metadata$Cancer_Type,
  levels=c("BRCA","LUAD")
)
table(sample_metadata$Cancer_Type, useNA = "ifany")
design <- model.matrix(~ Cancer_Type, data = sample_metadata)
fit <- lmFit(v,design)
fit <- eBayes(fit)
de_results <- topTable(fit, coef = "Cancer_TypeLUAD", number = Inf) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by = "gene_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 1.5  ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < -1.5 ~ "Down-regulated",
      TRUE                           ~ "Not Significant"
    )
  )
var_genes <- head(
  order(rowVars(log_expr),decreasing = TRUE),
  500
)
pca_res <- prcomp(
  t(log_expr[var_genes, ]),
  scale. = TRUE
)
var_explained <- summary(pca_res)$importance[2, ] * 100
var_explained[1:2]
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("barcode") %>%
  inner_join(sample_metadata, by= "barcode")
head(pca_df[, c("barcode", "PC1", "PC2", "Cancer_Type")])
attr(pca_df, "var_explained") <- var_explained
dir.create("output/task01",recursive = TRUE, showWarnings = FALSE)
write.csv(de_results, "output/task01/de_results.csv", row.names= FALSE)
write.csv(pca_df, "output/task01/pca_results.csv",row.names = FALSE)
saveRDS(var_explained, "output/task01/variance_explained.rds")
message("Statistical modeling and PCA complete!")
