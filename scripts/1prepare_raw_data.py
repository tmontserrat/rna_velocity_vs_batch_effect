import scanpy as sc

# Load the dataset
adata = sc.read_h5ad('../data/endocrinogenesis_day15.h5ad')

# Get the spliced and unspliced count matrices
spliced = adata.to_df(layer='spliced')
unspliced = adata.to_df(layer='unspliced')

# Cell types information
cell_types = adata.obs[['clusters']]

# Save the data to use in Seurat
spliced.to_csv('../data/spliced.csv', header=True)
unspliced.to_csv('../data/unspliced.csv', header=True)
cell_types.to_csv('../data/cell_types.csv', header=True)