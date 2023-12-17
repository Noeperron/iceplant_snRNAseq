set.seed(3112)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

###################### Plot gene expression along lineages ######################

######## For a single gene ######## 

sigGeneStart <- "Mcr-001290"
smoother <- plotSmoothers(sce_fitted, counts, gene = sigGeneStart)
feature <- FeaturePlot(cds, features = sigGeneStart, pt.size=0.7) 
UMAP <- plotGeneCount(sds2, counts, gene = sigGeneStart)


ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/sigGeneStart/smoother.png"),
                smoother,
                height=15,
                width=40,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/sigGeneStart/feature.png"),
                feature,
                height=15,
                width=20,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/sigGeneStart/UMAP.png"),
                UMAP,
                height=15,
                width=20,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

######## For a list of genes ######## 

# Read the gene list and corresponding names from the CSV file
gene_data <- read.csv("ListsKnownMarkers.csv")
gene_list <- gene_data[,1]
gene_names <- gene_data[,2]

# Loop over each gene and its corresponding name
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  gene_name <- gene_names[i]
  
  # Perform your analyses
  smoother <- plotSmoothers(sce_fitted, counts, gene = gene)
  feature <- FeaturePlot(cds, features = gene, pt.size=1, order=TRUE) 
  UMAP <- plotGeneCount(sds2, counts, gene = gene)
  
  # Save the plots
  ggsave(
    filename = paste0("Trajectory_C6/Expression_along_trajectory/", gene_name, "/smoother.png"),
    plot = smoother,
    height = 15,
    width = 15,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  
  ggsave(
    filename = paste0("Trajectory_C6/Expression_along_trajectory/", gene_name, "/feature.png"),
    plot = feature,
    height = 15,
    width = 20,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  
  ggsave(
    filename = paste0("Trajectory_C6/Expression_along_trajectory/", gene_name, "/UMAP.png"),
    plot = UMAP,
    height = 15,
    width = 20,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  
}
