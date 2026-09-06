suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
cancer_label <-"Breast Cancer (BRCA)"
de_results <-read.csv("output/de_results.csv")
pca_df <-read.csv("output/pca_results.csv")
var_explained <- readRDS("output/variance_explained.rds")
pca_df$Condition <- factor(pca_df$Condition,levels=c("Normal","Tumour"))
#1.Principal Componet Analysis (PCA) Plot
p_pca <-ggplot(pca_df,aes(x= PC1,y= PC2, color = Condition)) +
  geom_point(size = 3.5,alpha =0.85)+
  scale_color_manual(values = c("Normal"="#2b5c8f","Tumour"="#d95f02"))+
  theme_minimal(base_size=12)+
  labs(
    title=sprintf("PCA:%s (Tumour vs Nomral)",cancer_label),
    x= sprintf("PC1(%.1f%% Variance",var_explained[1]),
    y= sprintf("PC2(%.1f%% Variance",var_explained[2]),
    color="Group"
  )+
  theme(
    plot.title = element_text(face="bold",hjust = 0.5),
    legend.position = "top"
  )
