# Check integration of a simulated batch effect on an spliced matrix

### LIBRARIES, DATA AND SEURAT OBJECTS ###

# Libraries and functions
library(Seurat)
library(tidyverse)

# Load the data
load(file = "../data/spliced_b1.RData")
load(file = "../data/spliced_b2.RData")
cell_types <- read.csv("../../initial_data/cell_types.csv", row.names = 1)

# Create Seurat objects
spliced_b1 <- CreateSeuratObject(counts = t(spliced_b1))
spliced_b2 <- CreateSeuratObject(counts = t(spliced_b2))

# Merge datasets
merged_seurat <- merge(x = spliced_b1, y = spliced_b2,
                       add.cell.ids = c("b1", "b2"),
                       project = "batch_velocity")

# Create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split the sample column to create barcode and batch columns
merged_seurat@meta.data <- separate(merged_seurat@meta.data,
                                    col = "sample",
                                    into = c("batch", "barcode"),
                                    sep = "_")

# Add annotations
metadata <- merge(merged_seurat@meta.data, cell_types, 
                                by.x = "barcode",
                                by.y = "row.names", 
                  all = TRUE,
                  sort = FALSE)

merged_seurat@meta.data$cell_type <- metadata$clusters


### STANDARD WORKFLOW ###

# Normalize the data
merged_seurat <- SCTransform(merged_seurat, verbose = TRUE,
                             return.only.var.genes = TRUE,
                             min_cells = 3)

# Dimensionality reduction
merged_seurat <- RunPCA(merged_seurat)

# Non-linear dimensionality reduction
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)
ElbowPlot(merged_seurat)

# Visualize the dataset
plot_cell_type <- DimPlot(merged_seurat, reduction = 'umap', 
                          group.by = c('cell_type'), label = TRUE) + 
  NoLegend()

plot_batch <- DimPlot(merged_seurat, reduction = 'umap', group.by = "batch")

plots <- plot_cell_type | plot_batch

plots

ggsave("../results/combined_no_integrated/plots/seurat_umap_no_integrated.pdf",
       device = "pdf",
       width = 2100,
       height = 1200,
       units = "px")


# Extract embeddings
seurat_umap <- as.data.frame(merged_seurat@reductions$umap@cell.embeddings)
seurat_pca <- as.data.frame(merged_seurat@reductions$pca@cell.embeddings)
pearson_residuals_x <- t(as.data.frame(merged_seurat@assays$SCT@scale.data))

write.csv(pearson_residuals_x,
          file = "../data/sct_x_no_anchors.csv",
          quote = FALSE,
          row.names = TRUE)
write.csv(seurat_pca,
          file = "../data/seurat_pca_no_anchors.csv",
          quote = FALSE,
          row.names = TRUE)
write.csv(seurat_umap,
          file = "../data/seurat_umap_no_anchors.csv",
          quote = FALSE,
          row.names = TRUE)

### INTEGRATION WITH ANCHORS ###

# Load the data for a fresh start
load(file = "../data/spliced_b1.RData")
load(file = "../data/spliced_b2.RData")
cell_types <- read.csv("../data/cell_types.csv", row.names = 1)

# Create Seurat objects
spliced_b1 <- CreateSeuratObject(counts = t(spliced_b1))
spliced_b2 <- CreateSeuratObject(counts = t(spliced_b2))

# Merge datasets
merged_seurat <- merge(x = spliced_b1, y = spliced_b2,
                       add.cell.ids = c("b1", "b2"),
                       project = "batch_velocity")

# Create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split the sample column to create barcode and batch columns
merged_seurat@meta.data <- separate(merged_seurat@meta.data,
                                    col = "sample",
                                    into = c("batch", "barcode"),
                                    sep = "_")

# Add annotations
metadata <- merge(merged_seurat@meta.data, cell_types, 
                  by.x = "barcode",
                  by.y = "row.names", 
                  all = TRUE,
                  sort = FALSE)

merged_seurat@meta.data$cell_type <- metadata$clusters


# Split object by batch
obj.list <- SplitObject(merged_seurat, split.by = "batch")

obj.list <- lapply(X = obj.list,  FUN = SCTransform)


# Seelct integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 4500)

# Preparation step
obj.list <- PrepSCTIntegration(object.list = obj.list,
                               anchor.features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

# Perform the integration
seurat_integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT")


### DIMENSIONALITY REDUCTION ###

# PCA
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
ElbowPlot(seurat_integrated, ndims = 30)

# UMAP
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30,
                             n.neighbors = 45L)

plot_cell_type_int <- DimPlot(seurat_integrated, reduction = 'umap', 
                              group.by = c('cell_type'), label = TRUE) +
  # cols = c("purple", "yellow", "orange", "blue", "black", "pink", "green", "cyan")
  NoLegend()
plot_batch_int <- DimPlot(seurat_integrated, reduction = 'umap', group.by = "batch")

plots <- plot_cell_type_int | plot_batch_int

plots
ggsave("../results/combined_integrated/plots/seurat_umap_integrated.pdf",
       device = "pdf",
       width = 2100,
       height = 1200,
       units = "px")

# Extract embeddings
seurat_umap_anchors <- as.data.frame(seurat_integrated@reductions$umap@cell.embeddings)
seurat_pca_anchors <- as.data.frame(seurat_integrated@reductions$pca@cell.embeddings)
corrected_x <- t(as.data.frame(seurat_integrated@assays$integrated@scale.data))

write.csv(corrected_x,
          file = "../data/corrected_x.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(seurat_pca_anchors,
          file = "../data/seurat_pca_anchors.csv",
          quote = FALSE,
          row.names = TRUE)
write.csv(seurat_umap_anchors,
          file = "../data/seurat_umap_anchors.csv",
          quote = FALSE,
          row.names = TRUE)
