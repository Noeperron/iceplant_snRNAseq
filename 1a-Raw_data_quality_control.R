set.seed(3112)

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

#Convert to single-cell experiment object
D8DS_sce <- SingleCellExperiment(assays = List(counts = counts_D8DS ) )
D8DC_sce <- SingleCellExperiment(assays = List(counts = counts_D8DC ) )
D8LS_sce <- SingleCellExperiment(assays = List(counts = counts_D8LS ) )
D8LC_sce <- SingleCellExperiment(assays = List(counts = counts_D8LC ) )