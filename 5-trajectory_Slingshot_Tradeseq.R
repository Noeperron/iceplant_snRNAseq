set.seed(3112)
library(Seurat)
library(cowplot)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggplot2)

# CRAN
suppressMessages( require(tidyverse) )
suppressMessages( require(SeuratObject) )
suppressMessages( require(patchwork) )
suppressMessages( require(circlize) )
suppressMessages( require(reactable) )
suppressMessages( require(sctransform) )
suppressMessages( require(shiny) )
suppressMessages( require(shinyWidgets) )
suppressMessages( require(shinyFeedback) )
suppressMessages( require(rclipboard) )
suppressMessages( require(future) )
suppressMessages( require(ggthemes) )
suppressMessages( require(shinycssloaders) )
suppressMessages( require(DT) )
suppressMessages( require(dplyr) )
suppressMessages( require(hdf5r) )
suppressMessages( require(scales) )
suppressMessages( require(utils) )
suppressMessages( require(vroom) )
suppressMessages( require(svglite) )

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

cds <- readRDS("rds_files/OnlyMesophyllClusters.rds")

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 12

# Save the objects as separate matrices for input in slingshot
counts <- as.matrix(cds@assays$RNA@counts)
counts <- counts[rowSums(counts) > 0, ]
sce <- SingleCellExperiment(assays = List(counts = counts ) )

dimred <- cds@reductions$umap@cell.embeddings #Matrix dimensionality reduction
clustering <- cds@meta.data[["seurat_clusters"]] #Matrix with cluster assignment
STARTING_CLUSTER = 6

sds_trade <- slingshot::slingshot(
  sce,
  clusterLabels = clustering,
  start.clus = STARTING_CLUSTER,
  reducedDim = dimred
)

slingshot::slingLineages(sds_trade) #Get lineages from slingshot

sds2 <- SlingshotDataSet(sds_trade)
counts <- as.matrix(sds_trade@assays@data@listData[["counts"]])

# To evaluate the number of knots to be used on the modeling phase
icMat <- evaluateK(counts = counts,
                   sds = sds2,
                   k = 3:20,
                   BPPARAM = BPPARAM,
                   nGenes = 200,
                   verbose = T)

ggplot2::ggsave("Reclustering_Mesophyll/Trajectory/images/Tradeseq_knots.png",
                icMat,
                height=30,
                width=20,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

sce_fitted <- fitGAM(counts = counts,
                     sds = sds2,
                     nknots = 15,
                     verbose = T,
                     parallel = TRUE,
                     BPPARAM = BPPARAM)


saveRDS(sce_fitted,
        file = paste0("TRADEseq.rds")
)

sce_fitted <- readRDS("Trajectory_C6/rds_files/TRADEseq_C6.rds")



feat_importances <- tradeSeq::associationTest(sce_fitted, lineages=TRUE)

feat_importances$fdr <- stats::p.adjust(feat_importances$pvalue,
                                        method = "fdr",
                                        n = length(feat_importances$pvalue))

feat_importances$genes <- rownames(feat_importances)
feat_importances2 <- feat_importances[ ,c(ncol(feat_importances), 1:( ncol(feat_importances) - 1)) ]

foldername <- "Reclustering_Mesophyll/Trajectory/Markers/"
system( paste0("mkdir -p ", foldername) )
write.table( feat_importances2,
             paste0(foldername, "Association_Test_Markers.tsv"),
             col.names = T,
             row.names = F,
             sep = "\t",
             quote = F
)

#Discovering progenitor marker genes
startRes <- startVsEndTest(sce_fitted, lineages=TRUE)
write.table( startRes,
             paste0(foldername, "Progenitor_marker_genes.tsv"),
             col.names = T,
             row.names = T,
             sep = "\t",
             quote = F
)


#Discovering differentiated cell type markers
endRes <- diffEndTest(sce_fitted, pairwise=TRUE)
write.table( endRes,
             paste0(foldername, "Differentiated_celltypes_markers.tsv"),
             col.names = T,
             row.names = T,
             sep = "\t",
             quote = F
)

#Discovering genes with different expression patterns across lineages
patternRes <- patternTest(sce_fitted)
write.table( patternRes,
             paste0(foldername, "Different_exp_patterns_across_lineages.tsv"),
             col.names = T,
             row.names = T,
             sep = "\t",
             quote = F
)


oPat <- order(patternRes$waldStat, decreasing = TRUE)
plotSmoothers(sce_fitted, counts, gene = rownames(patternRes)[oPat][1]) #Plot expression of gene with most variation across lineages

#Early drivers of differentiation
plotGeneCount(curve = sds2, counts = counts,
              clusters = apply(slingClusterLabels(sds2), 1, which.max),
              models = sce_fitted)

earlyDERes <- earlyDETest(sce_fitted, knots = c(1, 5)) #DE analysis between two knots defined here. 
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)

write.table( earlyDERes,
             paste0(foldername, "Early_drivers_of_differentation.tsv"),
             col.names = T,
             row.names = T,
             sep = "\t",
             quote = F
)

plotSmoothers(sce_fitted, counts, gene = rownames(earlyDERes)[oEarly][1])
plotGeneCount(sds2, counts, gene = rownames(earlyDERes)[oEarly][1])

###################### Plot gene expression along lineages ######################

#oStart <- order(startRes$waldStat, decreasing = TRUE) #These two lines are for if i want to plot the expression of the most variable gene according to given test. Change startRes if want ot plot for a different test.
#sigGeneStart <- names(sce_fitted)[oStart[1]]
sigGeneStart <- "Mcr-001290"
smoother <- plotSmoothers(sce_fitted, counts, gene = sigGeneStart)
feature <- FeaturePlot(cds, features = sigGeneStart, pt.size=0.7) 
UMAP <- plotGeneCount(sds2, counts, gene = sigGeneStart)


ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/LHY_wide/smoother.png"),
                smoother,
                height=15,
                width=40,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/PPCK1b/feature.png"),
                feature,
                height=15,
                width=20,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


ggplot2::ggsave(filename = paste0("Trajectory_C6/Expression_along_trajectory/PPCK1b/UMAP.png"),
                UMAP,
                height=15,
                width=20,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

######## For a list of genes ######## 

# Read the gene list and corresponding names from the CSV file
gene_data <- read.csv("ListsKnownMarkers/DEcircardian_C3vsCAM_trajectories.csv")
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
    filename = paste0("Trajectory_C6/Expression_along_trajectory/DECircadian/", gene_name, "/smoother.png"),
    plot = smoother,
    height = 15,
    width = 15,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  
  ggsave(
    filename = paste0("Trajectory_C6/Expression_along_trajectory/DECircadian/", gene_name, "/feature.png"),
    plot = feature,
    height = 15,
    width = 20,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  
  ggsave(
    filename = paste0("Trajectory_C6/Expression_along_trajectory/DECircadian/", gene_name, "/UMAP.png"),
    plot = UMAP,
    height = 15,
    width = 20,
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )

}
  
  ### Save all smoother plots in one big image ###
  
  # Load necessary library
  library(gridExtra)
  library(ggplot2)
  
  # Read the gene list and corresponding names from the CSV file
  gene_data <- read.csv("ListsKnownMarkers/DEcircardian_C3vsCAM_trajectories.csv")
  gene_list <- gene_data[,1]
  gene_names <- gene_data[,2]
  
  # Initialize a list to store plots
  smoother_plots <- list()
  
  # Loop over each gene and its corresponding name
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    gene_name <- gene_names[i]
    
    # Perform your analyses
    smoother <- plotSmoothers(sce_fitted, counts, gene = gene)
    
    # Add a title to the plot with the gene ID
    smoother_with_title <- smoother + ggtitle(gene)
    
    # Add the smoother plot with title to the list
    smoother_plots[[i]] <- smoother_with_title
  }
  
  # Define the layout of the grid
  ncol <- 3 # Number of columns (adjust as needed)
  nrow <- 2 # Number of rows (adjust as needed)
  
  # Combine all smoother plots into one big plot
  combined_plot <- do.call(gridExtra::grid.arrange, c(smoother_plots, ncol = ncol, nrow = nrow))
  
  # Save the combined plot
  ggsave(
    filename = "Trajectory_C6/Expression_along_trajectory/0-DECircadian/0-combined_smoother_plots.png",
    plot = combined_plot,
    height = 15, # specify height
    width = 35,  # specify width
    dpi = 300,
    units = "cm",
    bg = "#FFFFFF"
  )
  