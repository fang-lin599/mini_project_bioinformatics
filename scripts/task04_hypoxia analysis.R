if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})
sample_metadata <- read.csv("data/task01/sample_metadata.csv")
gene_info <- read.csv("data/task01/gene_info.csv")
log_expr <- readRDS("data/task01/normalised_log_expression.rds")
sample_metadata$Cancer_Type<-factor(
  sample_metadata$Cancer_Type,
  levels=c("BRCA","LUAD")
)
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}
hallmark_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
target_genes <- hallmark_sets %>%
  filter(gs_name == "HALLMARK_HYPOXIA") %>%
  pull(gene_symbol) %>%
  unique()
target_info <- gene_info %>%
  filter(gene_name %in% target_genes)
mat_subset<- log_expr[
  rownames(log_expr) %in% target_info$gene_id,
  ,
  drop =FALSE
]
rownames(mat_subset) <- target_info$gene_name[
  match(rownames(mat_subset), target_info$gene_id)
]
long_df <- as.data.frame(mat_subset) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "barcode",
    values_to = "Expression"
  ) %>%
  left_join(sample_metadata, by = "barcode")
plot_genes <- head(sort(unique(long_df$Gene)), 12)
plot_df <- long_df %>%
  filter(Gene %in% plot_genes)
p_box <- ggplot(
  plot_df,
  aes(x = Cancer_Type, y = Expression, fill = Cancer_Type)
) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.2) +
  facet_wrap(~ Gene, ncol = 4) +
  scale_fill_manual(values = c(
    "BRCA" = "#2b5c8f",
    "LUAD" = "#d95f02"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title = "Hallmark Hypoxia: BRCA vs LUAD",
    subtitle = "First 12 available genes in alphabetical order",
    x = "Cancer Type",
    y = "Expression (log2 CPM)"
  ) +
  theme(legend.position = "none")
print(p_box)
ggsave(
  "output/task01/hypoxia_boxplot_first12.png",
  plot = p_box,
  width = 10,
  height = 8,
  dpi = 300
)
anno_df <- data.frame(
  Cancer_Type=sample_metadata$Cancer_Type
)
rownames(anno_df) <-sample_metadata$barcode
anno_colors <- list(
  Cancer_Type = c(
    "BRCA" = "#2b5c8f",
    "LUAD" = "#d95f02"
  )
)
pheatmap(
  mat_subset,
  annotation_col = anno_df,
  annotation_colors = anno_colors,
  scale = "row",
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Hallmark Hypoxia: BRCA vs LUAD",
  filename = "output/task01/hypoxia_heatmap.png",
  width = 8,
  height = 10
)
message("Hypoxia gene analysis complete! Results saved to output/task01/")
