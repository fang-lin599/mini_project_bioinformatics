suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
cancer_label <- "BRCA vs LUAD Tumours"
de_results <- read.csv("output/task01/de_results.csv")
pca_df <- read.csv("output/task01/pca_results.csv")
var_explained <- readRDS("output/task01/variance_explained.rds")
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cancer_Type)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c(
    "BRCA" = "#2b5c8f",
    "LUAD" = "#d95f02"
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("PCA: %s", cancer_label),
    x = sprintf("PC1 (%.1f%% Variance)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% Variance)", var_explained[2]),
    color = "Cancer Type"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )
print(p_pca)
ggsave(
  "output/task01/pca_plot.png",
  plot = p_pca,
  width = 6,
  height = 5,
  dpi = 300
)
top_genes <- de_results %>%
  filter(Significance != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(10)
p_volcano <- ggplot(
  de_results,
  aes(x = logFC, y = -log10(adj.P.Val), color = Significance)
) +
  geom_point(alpha = 0.5, size = 1.8) +
  scale_color_manual(values = c(
    "Down-regulated" = "#3182bd",
    "Not Significant" = "#969696",
    "Up-regulated" = "#de2d26"
  )) +
  geom_vline(
    xintercept = c(-1.5, 1.5),linetype = "dashed", alpha = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),linetype = "dashed", alpha = 0.5
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3.5,
    max.overlaps = 15,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15) +
  labs(
    title = sprintf("Differential Expression: %s", cancer_label),
    x = "Log2 Fold Change (LUAD relative to BRCA)",
    y = "-Log10 Adjusted P-value"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )
print(p_volcano)
ggsave(
  "output/task01/volcano_plot.png",
  plot = p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)
ggsave( "output/task01/volcano_plot.png", plot = p_volcano, width = 7,height = 6, dpi = 300)
message("Visualisations complete! PNG figures saved to output/task01/")