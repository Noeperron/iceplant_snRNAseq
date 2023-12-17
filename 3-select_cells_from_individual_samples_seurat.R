set.seed(3112)
library(dplyr)
library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)

datasets.combined <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

# Define the sample name you want to keep
samples_to_keep <- c("D8 - Dark Control", "D8 - Dark Salt")

# Filter the cells from the specified sample
datasets.combined_filtered <- subset(datasets.combined, orig.ident == samples_to_keep)

# Run clustering on the filtered object (optional)
#datasets.combined_filtered <- FindNeighbors(datasets.combined_filtered)
#datasets.combined_filtered <- FindClusters(datasets.combined_filtered)

# Access the clustering results
clusters <- Idents(datasets.combined_filtered)

DimPlot(datasets.combined_filtered, reduction = "umap", label = TRUE, repel = TRUE)

saveRDS(datasets.combined_filtered, file="rds_files/OnlyControls_Split.rds")

