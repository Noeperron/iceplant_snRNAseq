set.seed(3112)

#DEG between two specific clusters

DEG <- FindMarkers(datasets.combined, ident.1 = "5", ident.2 = "6") #DEG analysis between clusters 5 and 6

write.csv(DEG, file = "DEG_5vs6.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)
