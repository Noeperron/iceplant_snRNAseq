set.seed(3112)

theme_set(theme_cowplot())
allowWGCNAThreads(nThreads = 12)

# load the integrated snRNA-seq dataset

seurat_obj <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds") 

p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)
p

#Set up Seurat object for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "IcePlant" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
  reduction = "pca", # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = "seurat_clusters" # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix
seurat_obj <- NormalizeMetacells(seurat_obj)

#Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name=c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,14,15,16),
  group.by=c("seurat_clusters"), # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = "RNA", # using RNA assay
  slot = "data" # using normalized data
)

##Select soft-power threshold##
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed" # you can also use "unsigned" or "signed hybrid"
)

# plot the results
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

#Construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=7, 
  setDatExpr=FALSE,
  networkType = "signed",
  TOMType = "signed",
  tom_name = "AllClusters", # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE 
)

PlotDendrogram(seurat_obj, main="hdWGCNA Dendrogram")


#Save the dendrogram
png(filename = "dendrogram.png", width=15, height=10, units="cm", res=600, bg="#FFFFFF")
p1 <- PlotDendrogram(seurat_obj)
dev.off()


TOM <- GetTOM(seurat_obj)

##############Compute harmonized module eigengenes##############

# need to run ScaleData first or else harmony throws an error
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "seurat_clusters",
  group_name=c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,14,15,16)
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Module "
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=2)
p

ggplot2::ggsave(filename = paste0("top_gene_per_module.png"),
                p,
                height=23,
                width=15,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

##Getting the module assignment table

modules <- GetModules(seurat_obj)
head(modules[,1:6]) #show the first 6 columns

write.csv(modules, file = "WGCNA_modules_ALL.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

write.csv(hub_df, file = "Hub_genes_ALL.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

saveRDS(seurat_obj, file="hdWGCNA_object.rds")
seurat_obj <- readRDS("hdWGCNA_object.rds")

# compute gene scoring for the top 25 hub genes by kME for each module

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method="UCell"
)

##############Visualization##############

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=5)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features="scores", # plot the hub gene scores
  order="shuffle", # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
feature <- wrap_plots(plot_list, ncol=5)

ggplot2::ggsave(filename = paste0("FeaturePlot_modules.png"),
                feature,
                height=20,
                width=60,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


# plot module correlogram
ModuleCorrelogram(seurat_obj)

#Save the correlogram
png(filename = "Correlogram.png", width=15, height=10, units="cm", res=600, bg="#FFFFFF")
p1 <- ModuleCorrelogram(seurat_obj)
dev.off()


# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
dotplot <- DotPlot(seurat_obj, features=mods, group.by = NULL, dot.scale = 9, scale = T, scale.by = "radius")

# flip the x/y axes, rotate the axis labels, and change color scheme:
dotplot <- dotplot +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

ggplot2::ggsave(filename = paste0("DotPlot_modules.png"),
                dotplot,
                height=15,
                width=17,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")
