### Create two batches with batch effect for the first experiment

source("../../functions/select_populations.R")
source("../../functions/is_nan_data_frame.R")

set.seed(111)

# Paths
SPLICED <- "../../initial_data/spliced.csv"
UNSPLICED <- "../../initial_data/unspliced.csv"
CELL_TYPES <- "../../initial_data/cell_types.csv"
DATA <- "../data/"

# Load the data
spliced <- read.csv(SPLICED, row.names = 1)
unspliced <- read.csv(UNSPLICED, row.names = 1)
cell_types <- read.csv(CELL_TYPES, row.names = 1)

# Add the barcodes as a new column in the 'cell_types' object
cell_types$barcode <- rownames(cell_types)

# Proportion of each cell type: Pre-endocrine, Ductal, Alpha, Ngn3 high EP, Delta, Beta, Ngn3 low EP, Epsilon 
# proportions <- runif(8, min = 0.45, max = 0.55) # Add some random noise

### Dividing the dataset in two halves ###

# Select the indices for one half
indices_one_half <- sample(rownames(spliced), 
                           round(nrow(spliced)/2),
                           replace = FALSE)

# Compute the proportion of each cell in the first batch
# cells_b1 <- select_populations(cell_types, proportions)
# indices_cells_b1 <- which(rownames(spliced) %in% cells_b1)

# Split the datasets in two batches
spliced_b1 <- spliced[indices_one_half, ]
unspliced_b1 <- unspliced[indices_one_half, ]
spliced_b2 <- spliced[-which(rownames(spliced) %in% indices_one_half), ]
unspliced_b2 <- unspliced[-which(rownames(spliced) %in% indices_one_half), ]

# Compute the ratio matrix for each batch
# ratio_b1 <- (unspliced_b1) / (spliced_b1)
ratio_b2 <- (unspliced_b2) / (spliced_b2)

# Simulate a simple batch effect
batch_effect <- round(runif(nrow(spliced_b2),
                            min = 2,
                            max = 5))

# Save the vector as a matrix
batch_effect_matrix <- matrix(batch_effect, 
                              length(batch_effect), 
                              ncol(spliced_b2))

# Matrices of 0 and 1
zeros_spliced_b2 <- ifelse(spliced_b2 == 0, 0, 1)
zeros_unspliced_b2 <- ifelse(unspliced_b2 == 0, 0, 1)

# BATCH SPLICED MATRIX
# Add the vector to the gene counts of batch 2
spliced_b2 <- as.data.frame(apply(spliced_b2, 
                                  2,
                                  function(x) x + batch_effect))

# Keep zeros
spliced_b2 <- spliced_b2 * zeros_spliced_b2

# BATCH UNSPLICED MATRIX
# Apply is.nan to a dataframe NaN to 0
ratio_b2[is.nan(ratio_b2)] <- 0

# Regenerate the unspliced matrix keeping the proportions
unspliced_b2 <- as.data.frame(ratio_b2 * spliced_b2)

# Unspliced matrix without batch effect
unspliced_b2_no_batch <- unspliced[-which(rownames(spliced) %in% indices_one_half), ]

# Change the infinite values (now codified as NaN) by the original number plus the batch effect
unspliced_b2[is.nan(unspliced_b2)] <- unspliced_b2_no_batch[is.nan(unspliced_b2)] + batch_effect_matrix[is.nan(unspliced_b2)]

# Round the counts from 'unspliced_b2'
unspliced_b2 <- round(unspliced_b2)

# Convert to dataframe and give the appropriate structure to export the matrices
spliced_b1_df <- as.data.frame(spliced_b1)
unspliced_b1_df <- as.data.frame(unspliced_b1)

spliced_b2_df <- as.data.frame(spliced_b2)
unspliced_b2_df <- as.data.frame(unspliced_b2)

spliced_b1_df <- cbind(barcodes = rownames(spliced_b1_df), spliced_b1_df)
rownames(spliced_b1_df) <- NULL

unspliced_b1_df <- cbind(barcodes = rownames(unspliced_b1_df), unspliced_b1_df)
rownames(unspliced_b1_df) <- NULL

spliced_b2_df <- cbind(barcodes = rownames(spliced_b2_df), spliced_b2_df)
rownames(spliced_b2_df) <- NULL

unspliced_b2_df <- cbind(barcodes = rownames(unspliced_b2_df), unspliced_b2_df)
rownames(unspliced_b2_df) <- NULL

# Export the matrices to downstream analysis
write.csv(spliced_b1_df,
          file = paste0(DATA, "spliced_b1.csv"),
          quote = FALSE,
          row.names = FALSE)

write.csv(unspliced_b1_df,
          file = paste0(DATA, "unspliced_b1.csv"),
          quote = FALSE,
          row.names = FALSE)

write.csv(spliced_b2_df,
          file = paste0(DATA, "spliced_b2.csv"),
          quote = FALSE,
          row.names = FALSE)

write.csv(unspliced_b2_df,
          file = paste0(DATA,"unspliced_b2.csv"),
          quote = FALSE,
          row.names = FALSE)

# Save RData for use in downstream analysis
save(spliced_b1, file = paste0(DATA, "spliced_b1.RData"))
save(spliced_b2, file = paste0(DATA, "spliced_b2.RData"))