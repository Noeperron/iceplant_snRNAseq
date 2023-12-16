set.seed(3112)
library(dplyr)
library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)
library(Matrix)
library(scDblFinder)
library(tidyr)
library(tidyverse)
library(scater)
library(DropletUtils)

#Import raw data
D8DS <- Read10X(data.dir = "/Raw_sc_data/D8DS/")
D8DC <- Read10X(data.dir = "/Raw_sc_data/D8DC/")
D8LS <- Read10X(data.dir = "/Raw_sc_data/D8LS/")
D8LC <- Read10X(data.dir = "/Raw_sc_data/D8LC/")

#Create a new dgCMatrix without the genes to remove (chloroplastic genes)
genes_to_remove <- read.csv('List_of_genes_to_remove.csv', header = FALSE)
genes_to_remove <- genes_to_remove$V1

filtered_genes <- rownames(D8DS)[!rownames(D8DS) %in% genes_to_remove]
D8DS <- D8DS[filtered_genes, ]

filtered_genes <- rownames(D8DC)[!rownames(D8DC) %in% genes_to_remove]
D8DC<- D8DC[filtered_genes, ]

filtered_genes <- rownames(D8LS)[!rownames(D8LS) %in% genes_to_remove]
D8LS <- D8LS[filtered_genes, ]

filtered_genes <- rownames(D8LC)[!rownames(D8LC) %in% genes_to_remove]
D8LC <- D8LC[filtered_genes, ]

#########Create seurat objects##########

D8DS <- CreateSeuratObject(counts = D8DS, project = "Dark - Salt", min.cells = 3)
D8DC <- CreateSeuratObject(counts = D8DC, project = "Dark - Control", min.cells = 3)
D8LS <- CreateSeuratObject(counts = D8LS, project = "Light - Salt", min.cells = 3)
D8LC <- CreateSeuratObject(counts = D8LC, project = "Light - Control", min.cells = 3)

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

all.genes_D8DS <- rownames(D8DS)
all.genes_D8DC <- rownames(D8DC)
all.genes_D8LS <- rownames(D8LS)
all.genes_D8LC <- rownames(D8LC)

D8DS <- ScaleData(D8DS, features = all.genes_D8DS)
D8DC <- ScaleData(D8DC, features = all.genes_D8DC)
D8LS <- ScaleData(D8LS, features = all.genes_D8LS)
D8LC <- ScaleData(D8LC, features = all.genes_D8LC)

D8DS <- RunPCA(D8DS, features = VariableFeatures(object = D8DS), npcs = 30)
D8DC <- RunPCA(D8DC, features = VariableFeatures(object = D8DC), npcs = 30)
D8LS <- RunPCA(D8LS, features = VariableFeatures(object = D8LS), npcs = 30)
D8LC <- RunPCA(D8LC, features = VariableFeatures(object = D8LC), npcs = 30)

D8DS <- RunUMAP(D8DS, reduction = "pca", dims = 1:30)
D8DC <- RunUMAP(D8DC, reduction = "pca", dims = 1:30)
D8LS <- RunUMAP(D8LS, reduction = "pca", dims = 1:30)
D8LC <- RunUMAP(D8LC, reduction = "pca", dims = 1:30)

D8DS <- FindNeighbors(D8DS, reduction = "pca", dims = 1:30)
D8DC <- FindNeighbors(D8DC, reduction = "pca", dims = 1:30)
D8LS <- FindNeighbors(D8LS, reduction = "pca", dims = 1:30)
D8LC <- FindNeighbors(D8LC, reduction = "pca", dims = 1:30)

D8DS <- FindClusters(D8DS, resolution = 0.3)
D8DC <- FindClusters(D8DC, resolution = 0.3)
D8LS <- FindClusters(D8LS, resolution = 0.3)
D8LC <- FindClusters(D8LC, resolution = 0.3)

counts_D8DS <- as.matrix(D8DS@assays[["RNA"]]@counts)
counts_D8DC <- as.matrix(D8DC@assays[["RNA"]]@counts)
counts_D8LS <- as.matrix(D8LS@assays[["RNA"]]@counts)
counts_D8LC <- as.matrix(D8LC@assays[["RNA"]]@counts)

####Identifying empty droplets using emptyDrops####
#D8DS
limit <- 230 #All cells with less than this number of UMIs are considered as "empty"
D8DS.out <- emptyDrops(counts_D8DS, lower=limit, test.ambient=TRUE)
summary(D8DS.out$FDR <= 0.001)
sum(D8DS.out$Total <= limit & D8DS.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected
 
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8DS.out$PValue[D8DS.out$Total <= limit & D8DS.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.
 
#D8DC
limit <- 460 #All cells with less than this number of UMIs are considered as "empty"
D8DC.out <- emptyDrops(counts_D8DC, lower=limit, test.ambient=TRUE)
summary(D8DC.out$FDR <= 0.001)
sum(D8DC.out$Total <= limit & D8DC.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8DC.out$PValue[D8DC.out$Total <= limit & D8DC.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.
 
#D8LS
limit <- 330 #All cells with less than this number of UMIs are considered as "empty"
D8LS.out <- emptyDrops(counts_D8LS, lower=limit, test.ambient=TRUE)
summary(D8LS.out$FDR <= 0.001)
sum(D8LS.out$Total <= limit & D8LS.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8LS.out$PValue[D8LS.out$Total <= limit & D8LS.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.

#D8LC
limit <- 500 #All cells with less than this number of UMIs are considered as "empty"
D8LC.out <- emptyDrops(counts_D8LC, lower=limit, test.ambient=TRUE)
summary(D8LC.out$FDR <= 0.001)
sum(D8LC.out$Total <= limit & D8LC.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8LC.out$PValue[D8LC.out$Total <= limit & D8LC.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.
############################################

#Convert to single-cell experiment object
D8DS_sce <- SingleCellExperiment(assays = List(counts = counts_D8DS ) )
D8DC_sce <- SingleCellExperiment(assays = List(counts = counts_D8DC ) )
D8LS_sce <- SingleCellExperiment(assays = List(counts = counts_D8LS ) )
D8LC_sce <- SingleCellExperiment(assays = List(counts = counts_D8LC ) )

#Doublet Identification
D8DS_doublets <- scDblFinder(D8DS_sce, clusters=D8DS@meta.data[["seurat_clusters"]])
D8DC_doublets <- scDblFinder(D8DC_sce, clusters=D8DC@meta.data[["seurat_clusters"]])
D8LS_doublets <- scDblFinder(D8LS_sce, clusters=D8LS@meta.data[["seurat_clusters"]])
D8LC_doublets <- scDblFinder(D8LC_sce, clusters=D8LC@meta.data[["seurat_clusters"]])

#Create a logical vector where TRUE indicates singlets
is_singlet_D8DS <- D8DS_doublets$scDblFinder.class == "singlet"
is_singlet_D8DC <- D8DC_doublets$scDblFinder.class == "singlet"
is_singlet_D8LS <- D8LS_doublets$scDblFinder.class == "singlet"
is_singlet_D8LC <- D8LC_doublets$scDblFinder.class == "singlet"

#Subset the original SCE object to keep only singlets (keep only columns with singlet label)
D8DS_no_doublets <- D8DS_sce[, is_singlet_D8DS]
D8DC_no_doublets <- D8DC_sce[, is_singlet_D8DC]
D8LS_no_doublets <- D8LS_sce[, is_singlet_D8LS]
D8LC_no_doublets <- D8LC_sce[, is_singlet_D8LC]

#Converting to a sparse matrix
sparse_counts_D8DS <- as(D8DS_no_doublets@assays@data@listData[["counts"]], "CsparseMatrix")
sparse_counts_D8DC <- as(D8DC_no_doublets@assays@data@listData[["counts"]], "CsparseMatrix")
sparse_counts_D8LS <- as(D8LS_no_doublets@assays@data@listData[["counts"]], "CsparseMatrix")
sparse_counts_D8LC <- as(D8LC_no_doublets@assays@data@listData[["counts"]], "CsparseMatrix")

#Write new 10X files from the cleaned objects
write10xCounts("/Users/herve/Desktop/Seurat/data_cleaned/D8DS/", sparse_counts_D8DS, type="sparse")
write10xCounts("/Users/herve/Desktop/Seurat/data_cleaned/D8DC/", sparse_counts_D8DC, type="sparse")
write10xCounts("/Users/herve/Desktop/Seurat/data_cleaned/D8LS/", sparse_counts_D8LS, type="sparse")
write10xCounts("/Users/herve/Desktop/Seurat/data_cleaned/D8LC/", sparse_counts_D8LC, type="sparse")
