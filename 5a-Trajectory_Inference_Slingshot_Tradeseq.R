set.seed(3112)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

# Isolate mesophyll clusters to perform cell trajectory analysis
to_select <- c(0,1,3,5,6,15)
cds <- subset(x = datasets.combined, idents = to_select)

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

ggplot2::ggsave("/Trajectory/images/Tradeseq_knots.png",
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