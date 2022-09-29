import simulation_velocity_analysis as sva

# Integrated with principal components from Pearson residuals corrected
sva.perform_velocity_analysis(path_spliced="../data/spliced_matrix.csv",
                            path_unspliced="../data/unspliced_matrix.csv",
                            path_umap_embeddings="../data/seurat_umap_anchors.csv",
                            path_metadata="../data/seurat_integrated_metadata.csv",
                            velocity_plot_file_name="integrated/seurat_pca/velocity_anchors.png",
                            latent_time_plot_file_name="integrated/seurat_pca/latent_time_anchors.pdf",
                            pseudotime_plot_file_name="integrated/seurat_pca/pseudotime_anchors.pdf",
                            adata_processed_object_file_name="gt_anchors.h5ad",
                            use_velocity_ground_truth=True,
                            path_velocity_ground_truth="../data/velocity_matrix.csv",
                            corrected_pca=True,
                            path_pca_embeddings="../data/seurat_pca_anchors.csv",
                            n_jobs=6,
                            vel_mode="dynamical",
                            n_recurse_neighbors=2)

# Integrated with principal components from raw counts
sva.perform_velocity_analysis(path_spliced="../data/spliced_matrix.csv",
                            path_unspliced="../data/unspliced_matrix.csv",
                            path_umap_embeddings="../data/seurat_umap_anchors.csv",
                            path_metadata="../data/seurat_integrated_metadata.csv",
                            velocity_plot_file_name="integrated/raw_pca/velocity_anchors.png",
                            latent_time_plot_file_name="integrated/raw_pca/latent_time_anchors.pdf",
                            pseudotime_plot_file_name="integrated/raw_pca/pseudotime_anchors.pdf",
                            adata_processed_object_file_name="gt_anchors_raw_pca.h5ad",
                            use_velocity_ground_truth=True,
                            path_velocity_ground_truth="../data/velocity_matrix.csv",
                            corrected_pca=False,
                            path_pca_embeddings=None,
                            n_jobs=6,
                            vel_mode="dynamical",
                            n_recurse_neighbors=2)


# Unintegrated
sva.perform_velocity_analysis(path_spliced="../data/spliced_matrix.csv",
                            path_unspliced="../data/unspliced_matrix.csv",
                            path_umap_embeddings="../data/seurat_umap_no_anchors.csv",
                            path_metadata="../data/seurat_integrated_metadata.csv",
                            velocity_plot_file_name="unintegrated/velocity_no_anchors.png",
                            latent_time_plot_file_name="unintegrated/latent_time_no_anchors.pdf",
                            pseudotime_plot_file_name="unintegrated/pseudotime_no_anchors.pdf",
                            adata_processed_object_file_name="gt_unintegrated.h5ad",
                            use_velocity_ground_truth=True,
                            path_velocity_ground_truth="../data/velocity_matrix.csv",
                            corrected_pca=False,
                            path_pca_embeddings=None,
                            n_jobs=6,
                            vel_mode="dynamical",
                            n_recurse_neighbors=2)