set.seed(3112)

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

#######################################

D8DS <- subset(D8DS, subset = nFeature_RNA > 230) # Based on the results of the emtpyDrops() function specific to this sample
D8DC <- subset(D8DC, subset = nFeature_RNA > 460) # Based on the results of the emtpyDrops() function specific to this sample
D8LS <- subset(D8LS, subset = nFeature_RNA > 330) # Based on the results of the emtpyDrops() function specific to this sample
D8LC <- subset(D8LC, subset = nFeature_RNA > 500) # Based on the results of the emtpyDrops() function specific to this sample

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

saveRDS(datasets.combined, file="rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

#(Optional) Select a subset of clusters for graphs 
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

ggplot2::ggsave(filename = paste0("images/UMAP.png"),
                p2,
                height=20,
                width=20,
                dpi=500,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/UMAP_per_sample.png"),
                p3,
                height=10,
                width=35,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")
