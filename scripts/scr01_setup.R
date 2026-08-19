if (!requireNamespace("BiocManager",quietly= TRUE)) {
  install.packages("BiocManager")
}
cran_packages <- c("dplyr", "tidyr", "ggplot2", "ggrepel")
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment", "limma", "edgeR")
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(new_cran) > 0) {
  install.packages(new_cran)
}
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if(length(new_bioc)>0) {
  BiocManager::install(new_bioc)
}
message("Environment successfully setup! All dependencies are installed.")
