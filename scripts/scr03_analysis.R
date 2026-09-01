suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(matrixStats)
})
sample_metadata<- read.csv("data/sample_metadata.csv")
gene_info <- read.csv("data/gene_info.csv", row.names = 1)
log_expr <- readRDS("data/normalised_log_expression.rds")
v <-readRDS("data/voom_object.rds")
sample_metadata$Condition <-factor(sample_metadata$Condition, levels=c("Normal","Tumour"))
#1.Differential Expression Analysis (limma)
design <- model.matrix(~Condition, data = sample_metadata)
fit <- lmFit(v,design)
fit <- eBayes(fit)
de_results <- topTable(fit, coef = "ConditionTumour", number = Inf) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by ="gene_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
    Significance = case_when(
      adj.P.Val <0.05 & logFC >1.5 ~ "Up-regulated",
      adj.P.Val <0.05 & logFC < -1.5 ~ "Down-regulated",
      TRUE                          ~ "Not Significant"
    )
  )
#2.Dimensionality Reduction (PCA)
var_genes <- head(order(rowVars(log_expr),decreasing = TRUE),500)
pca_res <- prcomp (t(log_expr[var_genes, ]),scale= TRUE)
var_explained <- summary(pca_res)$importance[2, ]* 100
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("barcode")%>%
  inner_join(sample_metadata, by = "barcode")
attr(pca_df,"var_explained") <-var_explained
#3.Save Analytical Outputs
write.csv(de_results,"output/de_results.csv", row.names = FALSE)
write.csv(pca_df, "output/pca_results.csv", row.names = FALSE)
saveRDS(var_explained, "output/variance_explained.rds")
message("Statistical modeling and PCA complete! Results exported to output/")
