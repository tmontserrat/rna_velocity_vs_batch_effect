import velocity_analysis as va

va.perform_velocity_analysis(path_spliced="../../initial_data/spliced.csv", 
                            path_unspliced="../../initial_data/unspliced.csv",
                            path_cell_types="../../initial_data/cell_types.csv",
                            path_umap_embeddings="../data/seurat_umap.csv",
                            path_pca_embeddings="../data/seurat_pca.csv",
                            velocity_plot_file_name="plots/velocity.png",
                            latent_time_plot_file_name="plots/latent_time.pdf",
                            pseudotime_plot_file_name="plots/pseudotime.pdf",
                            velocity_confidence_plot_file_name="plots/confidence.pdf",
                            velocity_confidence_data_file_name="tables/confidence_table.csv",
                            adata_processed_object_file_name="../processed_objects/control_processed.h5ad",
                            velocity_graph_file_name="plots/velocity_graph.pdf",
                            paga_plot_file_name="plots/paga_plot.pdf",
                            path_corrected_x="../data/sct_x.csv",
                            n_jobs = 6,
                            vel_mode='dynamical',
                            n_recurse_neighbors=4
                            )