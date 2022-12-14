import velocity_analysis as va

# Integrated with principal components from Pearson residuals corrected
va.perform_velocity_analysis_2_samples(path_spliced1="../data/spliced_b1.csv",
                                    path_spliced2="../data/spliced_b2.csv",
                                    path_unspliced1="../data/unspliced_b1.csv",
                                    path_unspliced2="../data/unspliced_b2.csv",
                                    integrated=True,
                                    corrected_pca = True,
                                    path_cell_types="../../initial_data/cell_types.csv",
                                    path_umap_embeddings="../data/seurat_umap_anchors.csv",
                                    velocity_plot_file_name="combined_integrated/plots/velocity_anchors.png",
                                    latent_time_plot_file_name="combined_integrated/plots/latent_time_anchors.pdf",
                                    pseudotime_plot_file_name="combined_integrated/plots/pseudotime_anchors.pdf",
                                    velocity_confidence_plot_file_name="combined_integrated/tables/confidence_anchors.pdf",
                                    velocity_confidence_data_file_name="combined_integrated/tables/confidence_table_anchors.csv",
                                    velocity_graph_file_name="combined_integrated/plots/graph_anchors.pdf",
                                    paga_plot_file_name="combined_integrated/plots/paga_anchors.pdf",
                                    adata_processed_object_file_name="../processed_objects/anchors_processed.h5ad",
                                    path_pca_embeddings="../data/seurat_pca_anchors.csv",
                                    path_corrected_x="../data/corrected_x.csv",
                                    n_recurse_neighbors=4,
                                    n_jobs=6)

# Integrated with principals components from normalized raw counts
va.perform_velocity_analysis_2_samples(path_spliced1="../data/spliced_b1.csv",
                                    path_spliced2="../data/spliced_b2.csv",
                                    path_unspliced1="../data/unspliced_b1.csv",
                                    path_unspliced2="../data/unspliced_b2.csv",
                                    integrated=True,
                                    corrected_pca = False,
                                    path_cell_types="../../initial_data/cell_types.csv",
                                    path_umap_embeddings="../data/seurat_umap_anchors.csv",
                                    velocity_plot_file_name="combined_integrated_raw_pca/plots/velocity_anchors.png",
                                    latent_time_plot_file_name="combined_integrated_raw_pca/plots/latent_time_anchors.pdf",
                                    pseudotime_plot_file_name="combined_integrated_raw_pca/plots/pseudotime_anchors.pdf",
                                    velocity_confidence_plot_file_name="combined_integrated_raw_pca/tables/confidence_anchors.pdf",
                                    velocity_confidence_data_file_name="combined_integrated_raw_pca/tables/confidence_table_anchors.csv",
                                    velocity_graph_file_name="combined_integrated_raw_pca/plots/graph_anchors.pdf",
                                    paga_plot_file_name="combined_integrated_raw_pca/plots/paga_anchors.pdf",
                                    adata_processed_object_file_name="../processed_objects/anchors_processed_raw_pca.h5ad",
                                    path_pca_embeddings=None,
                                    path_corrected_x=None,
                                    n_recurse_neighbors=4,
                                    n_jobs=6)

# Unintegrated
va.perform_velocity_analysis_2_samples(path_spliced1="../data/spliced_b1.csv",
                                    path_spliced2="../data/spliced_b2.csv",
                                    path_unspliced1="../data/unspliced_b1.csv",
                                    path_unspliced2="../data/unspliced_b2.csv",
                                    integrated=False,
                                    corrected_pca=False,
                                    path_cell_types="../../initial_data/cell_types.csv",
                                    path_umap_embeddings="../data/seurat_umap_no_anchors.csv",
                                    velocity_plot_file_name="combined_no_integrated/plots/velocity_no_anchors.png",
                                    latent_time_plot_file_name="combined_no_integrated/plots/latent_time_no_anchors.pdf",
                                    pseudotime_plot_file_name="combined_no_integrated/plots/pseudotime_no_anchors.pdf",
                                    velocity_graph_file_name="combined_no_integrated/plots/graph_anchors.pdf",
                                    paga_plot_file_name="combined_no_integrated/plots/paga_anchors.pdf",
                                    velocity_confidence_plot_file_name="combined_no_integrated/tables/confidence_no_anchors.pdf",
                                    velocity_confidence_data_file_name="combined_no_integrated/tables/confidence_table_no_anchors.csv",
                                    adata_processed_object_file_name="../processed_objects/no_anchors_processed.h5ad",
                                    n_recurse_neighbors=4,
                                    n_jobs=6)