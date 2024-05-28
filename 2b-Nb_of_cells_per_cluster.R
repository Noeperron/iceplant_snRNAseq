set.seed(3112)

datasets.combined <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

#Get number of cells per cluster and per sample of origin
Nb_cells_per_cluster <- table(datasets.combined$orig.ident, datasets.combined@meta.data$seurat_clusters)

write.csv(Nb_cells_per_cluster, file = "Nb_cells_per_cluster.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

