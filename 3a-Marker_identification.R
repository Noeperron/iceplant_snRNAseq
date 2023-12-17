set.seed(3112)

#find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(datasets.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(all.markers, file = "All_conserved_markers.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

#Top 50 most specific markers in each cluster
top50 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

write.csv(top50, file = "top50_gene_per_cluster_specificity.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

