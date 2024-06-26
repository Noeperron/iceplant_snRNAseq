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
install.packages("doRNG")

# Packages from Bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder") #Bioconductor
BiocManager::install("scater") #Bioconductor
BiocManager::install("DropletUtils") #Bioconductor
BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))
BiocManager::install("GENIE3")
BiocManager::install("edgeR")

# Packages from GitHub

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github("NightingaleHealth/ggforestplot")
devtools::install_github('smorabit/hdWGCNA', ref='dev')

# Load Packages

library(dplyr) 
library(Seurat) 
library(patchwork) 
library(metap) 
library(ggplot2) 
library(Matrix) 
library(tidyr) 
library(tidyverse) 
library(cowplot) 
library(RColorBrewer) 
library(readxl) 
library(reshape2) 
library(SeuratObject) 
library(circlize) 
library(reactable) 
library(sctransform) 
library(shiny) 
library(shinyWidgets) 
library(shinyFeedback) 
library(rclipboard) 
library(future) 
library(ggthemes) 
library(DT) 
library(hdf5r) 
library(scales) 
library(utils) 
library(vroom) 
library(svglite) 
library(doRNG)
library(scDblFinder)
library(scater)
library(DropletUtils)
library(WGCNA)
library(igraph)
library(devtools)
library(GeneOverlap)
library(ggrepel)
library(UCell)
library(GENIE3)
library(edgeR)
library(DoubletFinder)
library(hdWGCNA)
library(ggforestplot)

