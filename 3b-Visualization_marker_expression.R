set.seed(3112)

datasets.combined <- readRDS("rds_files/D8_Nodoublets_WithoutChloroplasts_0.55.rds")

#Visualize expression markers
####Heatmaps####

Markers <- read.csv("ListsKnownMarkers.csv") # Input a list of known markers with the gene ID in the first column
Markers <- Markers[,1]

datasets.combined_scaled <- ScaleData(datasets.combined, features = rownames(datasets.combined), assay="RNA", verbose = FALSE)
heatmap <- DoHeatmap(datasets.combined_scaled, slot="scale.data", features=Markers, assay = "RNA")

ggplot2::ggsave(filename = paste0("images/Heatmaps.png"),
                heatmap,
                height=20,
                width=50,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

####Dot plots######
# Define an order of cluster identities
my_levels <- c(0,1,6,5,3,15,7,8,9,14,2,11,16,10,12,13,4) #Organized by cell-type

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

ggplot2::ggsave(filename = paste0("images/DotPlots.png"),
                d1,
                height=45, 
                width=28,  
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

ggplot2::ggsave(filename = paste0("/images/Violin.png"),
                v1,
                height=8,
                width=50,
                dpi=500,
                units="cm",
                bg = "#FFFFFF")


####Feature plots####

FeaturePlot(datasets.combined, features = Markers) #Feature plots

f0 <- FeaturePlot(datasets.combined, features = c("Mcr-004318")) # Feature plots combined

ggplot2::ggsave(filename = paste0("images/FeaturePlots/PPC1.png"),
                f0,
                height=20,
                width=40,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

f1 <- FeaturePlot(datasets.combined, coord.fixed=T, features = c("Mcr-004318"), split.by = "orig.ident", pt.size=1, order=T) # Feature plots by sample

ggplot2::ggsave(filename = paste0("images/FeaturePlots/PPC1_by_sample.png"),
                f1,
                height=10,
                width=40,
                dpi=400,
                units="cm",
                bg = "#FFFFFF")

