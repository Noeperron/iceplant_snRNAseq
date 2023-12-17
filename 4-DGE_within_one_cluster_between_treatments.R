set.seed(3112)
library(dplyr)
library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)
library(Matrix)
library(DoubletFinder)
library(tidyr)
library(tidyverse)

seurat_object <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

# Split integrated object into salt-treated and control subsets
# Create logical vectors for treatments within the metadata
seurat_object@meta.data$salt_treatment <- seurat_object@meta.data$orig.ident %in% c("D8 - Dark Salt", "D8 - Light Salt")
seurat_object@meta.data$control_treatment <- seurat_object@meta.data$orig.ident %in% c("D8 - Dark Control", "D8 - Light Control")


# Now loop through the clusters and perform differential expression
for (cluster_id in 0:16) {
  # Find cells in the cluster and that are salt-treated
  salt_cells <- WhichCells(seurat_object, 
                           expression = salt_treatment & seurat_object@meta.data$seurat_clusters == cluster_id)
  
  # Find cells in the cluster and that are control
  control_cells <- WhichCells(seurat_object, 
                              expression = control_treatment & seurat_object@meta.data$seurat_clusters == cluster_id)
  
  # Perform differential expression analysis
  DGE <- FindMarkers(seurat_object, ident.1 = salt_cells, ident.2 = control_cells)
  
  # Check if DGE is not empty then write to file
  if (!is.null(DGE) && nrow(DGE) > 0) {
    write.csv(DGE, 
              file = paste0("DGE/Salt_vs_control/DEG_saltvscontrol_cluster", cluster_id, ".csv"),
              row.names = TRUE)
  }
}
