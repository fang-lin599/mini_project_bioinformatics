if(!requireNamespace("pheatmap",quietly=TRUE)){
  install.packages("pheatmap")
}
supressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})
sample_metadata <- read.csv("data/sample_metadata.csv")
gene_info <-read.csv("data/gene_info.csv",row.name=1)
log_expr <-readRDS("data/normalised_log_expression.rds")
sample_metadata$Condition <- factor(sample_metadata$Condition,levels=c("Normal","Tumour"))
#2 Gene set extraction
target_genes <- c("CD274", "CTLA4", "PDCD1", "LAG3", "CDKN1A", "CDKN2A", "IL6", "CXCL8")
target_info<- gene_info %>%
  filter(gene_name %in% target_genes)
mat_subset <-log_expr[rownames(log_expr)%in% target_info$gene_id, ]
rownames(mat_subset) <- target_info$gene_name[match(rownames(mat_subset),target_info$gene_id)]
#3.Grouped Bosplots(ggplot2)
long_df <-as.data.frame(mat_subset)%>%
  tibble::rownames_to_column("Gene")%>%
  pivot_longer(cols = -Gene,names_to = "barcode",values_to = "Expression")%>%
  left_join(sample_metadata,by="barcode")
library(tidyr)

p_box <-ggplot(long_df,aes(x=Gene,y=Expression,fill=Condition))+
  geom_boxplot(alpha=0.8,outlier.size=1.2)+
  scale_fill_manual(values=c("Normal"="#2b5c8f","Tumour"="#d95f02"))+
  theme_classic(base_size=14)+
  labs(
    title = "Expression of Immune & Senescence Markers",
    x = "Gene Symbool",
    y= expression("Log"[2] ~ "Expression (CPM)")
  )+
  theme(plot.title=element_text(face="bold",hjust = 0.5))
print(p_box)
ggsave("output/geneset_boxplot.png", p_box, width = 8, height = 5, dpi = 300)

#4.Clustered Heatmap (pheatmap)
anno_df <- data.frame(Condition = sample_metadata$Condition)
rownames(anno_df) <- sample_metadata$barcode
anno_colors <- list(
  Condition = c("Normal" = "#2b5c8f", "Tumour" = "#d95f02")
)
pheatmap(
  mat_subset,
  annotation_col = anno_df,
  annotation_colors = anno_colors,
  scale = "row",
  show_colnames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Targeted Gene Expression Profile",
  filename = "output/geneset_heatmap.png",
  width = 7, 
  height = 5
)
library(pheatmap)
message("Gene set visualisations complete! PNG figures saved to output/")
