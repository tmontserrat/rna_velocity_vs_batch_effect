# Check integration of a simulated batch effect on an spliced matrix

### LIBRARIES, DATA AND SEURAT OBJECTS ###

# Libraries and functions
library(Seurat)
library(tidyverse)

# Load the data
spliced_matrix <- read.csv("../../initial_data/spliced.csv", row.names = 1)
cell_types <- read.csv("../../initial_data/cell_types.csv", row.names = 1)

# Create Seurat objects
spliced <- CreateSeuratObject(counts = t(spliced_matrix))

# Add annotations
metadata <- merge(spliced@meta.data, cell_types, 
                  by.x = "row.names",
                  by.y = "row.names", 
                  all = FALSE,
                  sort = FALSE)

spliced@meta.data$cell_type <- metadata$clusters


### STANDARD WORKFLOW ###

# Normalize the data
spliced <- SCTransform(spliced, verbose = TRUE,
                       return.only.var.genes = TRUE,
                       variable.features.n = 4000,
                       min_cells = 3)

# Dimensionality reduction
spliced <- RunPCA(spliced)
ElbowPlot(spliced, ndims = 30)

# Non-linear dimensionality reduction
spliced <- RunUMAP(spliced, dims = 1:30,
                   n.neighbors = 45L)

# Visualize the dataset
plot_cell_type <- DimPlot(spliced, reduction = 'umap', group.by = c('cell_type'), label = TRUE)
plot_cell_type

ggsave("../results/plots/seurat_umap.pdf",
       device = "pdf",
       width = 2100,
       height = 1200,
       units = "px")

# Extract embeddings
seurat_umap <- as.data.frame(spliced@reductions$umap@cell.embeddings)
seurat_pca <- as.data.frame(spliced@reductions$pca@cell.embeddings)
sct_x <- t(as.data.frame(spliced@assays$SCT@scale.data))

write.csv(seurat_umap,
          file = "../data/seurat_umap.csv",
          quote = FALSE,
          row.names = TRUE)
write.csv(seurat_pca,
          file = "../data/seurat_pca.csv",
          quote = FALSE,
          row.names = TRUE)
write.csv(sct_x,
          file = "../data/sct_x.csv",
          quote = FALSE,
          row.names = TRUE)
