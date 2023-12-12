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

#Data without doublets
D8DS <- Read10X(data.dir = "/data_cleaned/D8DS")
D8DC <- Read10X(data.dir = "/data_cleaned/D8DC")
D8LS <- Read10X(data.dir = "/data_cleaned/D8LS")
D8LC <- Read10X(data.dir = "/data_cleaned/D8LC")

#########Create seurat objects##########

D8DS <- CreateSeuratObject(counts = D8DS, project = "Dusk Salt", min.cells = 3)
D8DC <- CreateSeuratObject(counts = D8DC, project = "Dusk Control", min.cells = 3)
D8LS <- CreateSeuratObject(counts = D8LS, project = "Dawn Salt", min.cells = 3)
D8LC <- CreateSeuratObject(counts = D8LC, project = "Dawn Control", min.cells = 3)

VlnPlot(D8DS, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(D8DC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(D8LS, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(D8LC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

####################################

D8DS <- subset(D8DS, subset = nFeature_RNA > 230)
D8DC <- subset(D8DC, subset = nFeature_RNA > 460)
D8LS <- subset(D8LS, subset = nFeature_RNA > 330)
D8LC <- subset(D8LC, subset = nFeature_RNA > 500)

D8DS <- NormalizeData(D8DS)
D8DC <- NormalizeData(D8DC)
D8LS <- NormalizeData(D8LS)
D8LC <- NormalizeData(D8LC)

D8DS <- FindVariableFeatures(D8DS, selection.method = "vst", nfeatures = 2000)
D8DC <- FindVariableFeatures(D8DC, selection.method = "vst", nfeatures = 2000)
D8LS <- FindVariableFeatures(D8LS, selection.method = "vst", nfeatures = 2000)
D8LC <- FindVariableFeatures(D8LC, selection.method = "vst", nfeatures = 2000)

ifnb.list <- c(D8DS, D8DC, D8LC, D8LS)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

#Perform integration
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
datasets.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(datasets.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
datasets.combined <- ScaleData(datasets.combined, verbose = FALSE)
datasets.combined <- RunPCA(datasets.combined, npcs = 30, verbose = FALSE)
datasets.combined <- RunUMAP(datasets.combined, reduction = "pca", dims = 1:30)
datasets.combined <- FindNeighbors(datasets.combined, reduction = "pca", dims = 1:30)
datasets.combined <- FindClusters(datasets.combined, resolution = 0.55) #0.55 gives 16 clusters
DefaultAssay(datasets.combined) <- "RNA"

#saveRDS(datasets.combined, file="rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

datasets.combined <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")
datasets.combined <- readRDS("rds_files/OnlyControls_Split.rds")
datasets.combined <- readRDS("rds_files/OnlySalts_Split.rds")
datasets.combined <- readRDS("rds_files/OnlyMesophyllClusters.rds")


#(Optional) Select only some clusters for graphs 
to_select <- c(0,1,3,5,6,15)
datasets.combined <- subset(x = datasets.combined, idents = to_select)

# Visualization

p1 <- DimPlot(datasets.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(datasets.combined, reduction = "umap", label = T, repel = TRUE,
              pt.size = 0.8,
              label.size=6)
p1 + p2


p3 <- DimPlot(datasets.combined, 
        reduction = "umap", 
        split.by = "orig.ident", 
        label=T,
        label.size = 4, 
        label.color = "black", 
        pt.size=0.5)

system("mkdir -p images/")
ggplot2::ggsave(filename = paste0("images/Sample_proportion.png"),
                p1,
                height=10,
                width=15,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/UMAP_OnlyMesophyll.png"),
                p2,
                height=20,
                width=20,
                dpi=500,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/UMAP_OnlyMesophyll_per_sample.png"),
                p3,
                height=10,
                width=35,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

#Get number of cells per cluster and per sample of origin
Nb_cells_per_cluster <- table(datasets.combined$orig.ident, datasets.combined@meta.data$seurat_clusters)

write.csv(Nb_cells_per_cluster, file = "Nb_cells_per_cluster_NoDoublets.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

#Identify conserved cell type markers for one cluster

C6_markers <- FindConservedMarkers(datasets.combined, 
                                   ident.1 = 4, #Change to cluster number here
                                   grouping.var = "orig.ident", 
                                   verbose = FALSE)


write.csv(C6_markers, file = "C4_conserved_markers.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
            col.names = TRUE)

#find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(datasets.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(all.markers, file = "All_conserved_markers_res_0.5.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

#Top 50 markers of each cluster
top50 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

write.csv(top50, file = "SplitSamples/D8DC/top50_per_cluster.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

#Visualize expression markers
####Heatmaps####

Markers <- read.csv("ListsKnownMarkers/Annotation.csv")
Markers <- Markers[,1]

datasets.combined_scaled <- ScaleData(datasets.combined, features = rownames(datasets.combined), assay="RNA", verbose = FALSE)
heatmap <- DoHeatmap(datasets.combined_scaled, slot="scale.data", features=Markers, assay = "RNA")

ggplot2::ggsave(filename = paste0("images/Heatmaps/C3+CAM.png"),
                heatmap,
                height=20,
                width=50,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

####Dot plots######
# Define an order of cluster identities
my_levels <- c(0,1,6,5,3,15,7,8,9,14,2,11,16,10,12,13,4)
my_levels <- c(0,1,6,5,3,15)


datasets.combined$seurat_clusters <- factor(
  x = datasets.combined$seurat_clusters, 
  levels = as.character(my_levels)
)

d1 <- DotPlot(object = datasets.combined, features = Markers, 
        scale.by = "size", #size or radius
        dot.scale = 21, 
        scale = TRUE, 
        group.by = "seurat_clusters") + 
  RotatedAxis() + 
  coord_flip() + 
  theme(legend.position = "top") # Moves the legend to the top

ggplot2::ggsave(filename = paste0("images/DotPlots/Dotplot_Annotation.png"),
                d1,
                height=45, #For CAM and circadian: 25 #For annotation:45
                width=28,  #For CAM and circadian: 20 #For annotation:34
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


####Violin plots####

# Define an order of cluster identities
my_levels <- c(3,0,1,5,6,15)

datasets.combined$seurat_clusters <- factor(
  x = datasets.combined$seurat_clusters, 
  levels = as.character(my_levels)
)

v1 <- VlnPlot(
  datasets.combined, 
  features = Markers, 
  group.by = "seurat_clusters", 
  ncol = 5,
  adjust = 1, 
  same.y.lims = TRUE
)

ggplot2::ggsave(filename = paste0("OnlySalts/images/ViolinPlots/Violin_CAM_salts.png"),
                v1,
                height=8,
                width=50,
                dpi=500,
                units="cm",
                bg = "#FFFFFF")


####Feature plots####

FeaturePlot(datasets.combined, features = Markers) #Feature plots

f0 <- FeaturePlot(datasets.combined, features = c("Mcr-004318")) # Feature plots combined

ggplot2::ggsave(filename = paste0("images/FeaturePlots/CAM.png"),
                f0,
                height=20,
                width=40,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

f1 <- FeaturePlot(datasets.combined, coord.fixed=T, features = c("Mcr-009574"), split.by = "orig.ident", pt.size=1, order=T) # Feature plots by sample

ggplot2::ggsave(filename = paste0("images/FeaturePlots/PPCK1_by_sample.png"),
                f1,
                height=10,
                width=40,
                dpi=400,
                units="cm",
                bg = "#FFFFFF")

# Visualize co-expression of two features simultaneously
f2 <- FeaturePlot(datasets.combined,
                  features = c("Mcr-009574", "Mcr-002592"), 
                  pt.size = 2, 
                  blend=TRUE, 
                  order=TRUE, 
                  cols = c("lightgrey", "red", "blue"),
                  blend.threshold=0.5,
                  interactive=FALSE)

ggplot2::ggsave(filename = paste0("images/FeaturePlotsCoExpression/PPCK1_and_FULL.png"),
                f2,
                height=15,
                width=40,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


#DEG between two specific clusters

DEG <- FindMarkers(datasets.combined, ident.1 = "0", ident.2 = "1")

write.csv(DEG, file = "DEG_0vs1.csv", append = FALSE, quote = TRUE,
          row.names = TRUE,
          col.names = TRUE)

