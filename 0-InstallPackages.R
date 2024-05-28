# Packages from CRAN

install.packages("dplyr") 
install.packages("Seurat") 
install.packages("patchwork") 
install.packages("metap") 
install.packages("ggplot2") 
install.packages("Matrix") 
install.packages("tidyr") 
install.packages("tidyverse") 
install.packages("cowplot") 
install.packages("RColorBrewer") 
install.packages("readxl") 
install.packages("reshape2") 
install.packages("SeuratObject") 
install.packages("circlize") 
install.packages("reactable") 
install.packages("sctransform") 
install.packages("shiny") 
install.packages("shinyWidgets") 
install.packages("shinyFeedback") 
install.packages("rclipboard") 
install.packages("future") 
install.packages("ggthemes") 
install.packages("DT") 
install.packages("hdf5r") 
install.packages("scales") 
install.packages("utils") 
install.packages("vroom") 
install.packages("svglite") 
install.packages("WGCNA") 
install.packages("hdWGCNA") 


# Packages from Bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder") #Bioconductor
BiocManager::install("scater") #Bioconductor
BiocManager::install("DropletUtils") #Bioconductor
BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))

# Packages from GitHub

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github("NightingaleHealth/ggforestplot")
devtools::install_github('smorabit/hdWGCNA', ref='dev')
