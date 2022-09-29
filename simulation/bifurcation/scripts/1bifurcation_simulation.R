### Simulations 5000 cells per batch ###

### Simulating batch effect with dyngen ###

# Libraries
library(tidyverse)
library(dyngen)
library(dynwrap)
library(dynplot2)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(anndata)

# Pseudorandom seed
set.seed(2022)

# Type of backbone
backbone <- backbone_bifurcating()

# Number of cells and genes
num_cells <- 1500
# num_features <- 100

# Genes are divided up such that there are nrow(module_info) TFs
# and the rest is split between targets and HKs evenly
num_tfs <- nrow(backbone$module_info)
# num_targets <- round((num_features - num_tfs) / 2)
# num_hks <- num_features - num_targets - num_tfs

# Main simulations parameters. Same parameters as "Inferring single-cell 
# dynamics with structured dynamical representations of RNA velocity"
config <- 
  initialise_model(
    num_tfs = num_tfs,
    num_targets = 70,
    num_hks = 0,
    backbone = backbone,
    num_cells = num_cells,
    verbose = interactive(),
    # download_cache_dir = tools::R_user_dir("dyngen", "data"),
    # gold_standard_params = gold_standard_default(tau = 0.01),
    simulation_params = simulation_default(
      # total_time = 600,
      simtime_from_backbone(backbone, burn = TRUE) * 1.5,
      compute_rna_velocity = TRUE,
      store_reaction_propensities = TRUE,
      census_interval = 1,
      ssa_algorithm = ssa_etl(tau = 0.01),
      experiment_params = simulation_type_wild_type(num_simulations = 100)
    )
  )


# 'Common' part of the dyngen simulation
model_common <- 
  config %>% 
  generate_tf_network() %>% 
  generate_feature_network()

# Simulate two samples (a batch effect will be present)
model_a <- model_common %>% 
  generate_kinetics() %>% 
  generate_gold_standard() %>% 
  generate_cells() %>% 
  generate_experiment()

model_b <- model_common %>% 
  generate_kinetics() %>% 
  generate_gold_standard() %>% 
  generate_cells() %>% 
  generate_experiment()

save(model_a, file = "../models/bifurcated_model_a.RData")
save(model_b, file = "../models/bifurcated_model_b.RData")

# Convert to Seurat objects
seurat_a <- as_seurat(model_a)
seurat_b <- as_seurat(model_b)


# Merge both objects into one bigger Seurat object
merged_seurat <- merge(x = seurat_a, y = c(seurat_b),
                       add.cell.ids = c("b1", "b2"),
                       project = "batch_velocity")

DefaultAssay(merged_seurat) <- "spliced"

# Create a column in metadata called sample
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split the sample column to create cells and batch columns
merged_seurat@meta.data <- separate(merged_seurat@meta.data,
                                    col = "sample",
                                    into = c("batch", "barcode"),
                                    sep = "_")

# All genes names
features <- rownames(merged_seurat@assays$RNA@counts)

# Normalize the data
merged_seurat <- SCTransform(merged_seurat, verbose = TRUE,
                             return.only.var.genes = FALSE,
                             min_cells = 3)

# Dimensionality reduction
merged_seurat <- RunPCA(merged_seurat, features = features)
ElbowPlot(merged_seurat)

# Non-linear dimensionality reduction
merged_seurat <- RunUMAP(merged_seurat, dims = 1:6,
                         n.neighbors = 30)

# Visualize the dataset
merged_seurat@meta.data$sim_time <- merged_seurat@meta.data$sim_time
plot_batch <- DimPlot(merged_seurat, 
                      reduction = "umap", 
                      group.by = "batch")
plot_batch_time <- FeaturePlot(object = merged_seurat, 
                               reduction = "umap",
                               features = "sim_time")

plot_batch | plot_batch_time

seurat_umap <- as.data.frame(merged_seurat@reductions$umap@cell.embeddings)
seurat_pca <- as.data.frame(merged_seurat@reductions$pca@cell.embeddings)

spliced_matrix <- merged_seurat@assays$spliced@counts
unspliced_matrix <- merged_seurat@assays$unspliced@counts
velocity_matrix <- merged_seurat@assays$rnavelocity@counts

seurat_integrated_metadata <- merged_seurat@meta.data[, c("batch", "sim_time")]

write.csv(t(as.matrix(spliced_matrix)),
          file = "../data/spliced_matrix.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(t(as.matrix(unspliced_matrix)),
          file = "../data/unspliced_matrix.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(t(as.matrix(velocity_matrix)),
          file = "../data/velocity_matrix.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(seurat_umap,
          file = "../data/seurat_umap_no_anchors.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(seurat_pca,
          file = "../data/seurat_pca_batch_no_anchors.csv",
          quote = FALSE,
          row.names = TRUE)


# SaveH5Seurat(merged_seurat, filename = "../data/merged_seurat.h5Seurat")
# Convert("../data/merged_seurat.h5Seurat", dest = "h5ad")



################## ANCHORS ####################

# Split object by batch
obj.list <- SplitObject(merged_seurat, split.by = "batch")

obj.list <- lapply(X = obj.list,  FUN = SCTransform)


# Seelct integration features
# features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)

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


# PCA
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE, aprox = FALSE)
ElbowPlot(seurat_integrated)



# UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:6, n.neighbors = 20L)

plot_batch_int <- DimPlot(seurat_integrated, reduction = 'umap', group.by = "batch")

plot_batch_time <- FeaturePlot(object = seurat_integrated, features = "sim_time")

plot_batch_int | plot_batch_time


# Extract embeddings
seurat_umap_anchors <- as.data.frame(seurat_integrated@reductions$umap@cell.embeddings)
seurat_pca_anchors <- as.data.frame(seurat_integrated@reductions$pca@cell.embeddings)
corrected_x <- t(as.data.frame(seurat_integrated@assays$integrated@scale.data))

# Extract metadata
seurat_integrated_metadata <- seurat_integrated@meta.data[, c("batch", "sim_time")]


write.csv(seurat_umap_anchors,
          file = "../data/seurat_umap_anchors.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(seurat_pca_anchors,
          file = "../data/seurat_pca_anchors.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(seurat_integrated_metadata,
          file = "../data/seurat_integrated_metadata.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(corrected_x,
          file = "../data/corrected_x.csv",
          quote = FALSE,
          row.names = TRUE)
