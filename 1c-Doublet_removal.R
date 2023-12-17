set.seed(3112)

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
write10xCounts("/data_cleaned/D8DS/", sparse_counts_D8DS, type="sparse")
write10xCounts("/data_cleaned/D8DC/", sparse_counts_D8DC, type="sparse")
write10xCounts("/data_cleaned/D8LS/", sparse_counts_D8LS, type="sparse")
write10xCounts("/data_cleaned/D8LC/", sparse_counts_D8LC, type="sparse")
