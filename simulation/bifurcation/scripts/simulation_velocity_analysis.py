import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd

def perform_velocity_analysis(path_spliced, 
                            path_unspliced,
                            path_umap_embeddings,
                            path_metadata,
                            velocity_plot_file_name,
                            latent_time_plot_file_name,
                            pseudotime_plot_file_name,
                            adata_processed_object_file_name,
                            use_velocity_ground_truth=True,
                            path_velocity_ground_truth=None,
                            corrected_pca = True,
                            path_pca_embeddings=None,
                            n_jobs = 1,
                            vel_mode='dynamical',
                            n_recurse_neighbors=2
                            ):
    
    """Perform the velocity analysis for the spliced and unspliced RNA matrices supplied."""


    ### PREPARING THE AnnData OBJECT ###

    print("Preparing the AnnData object.")

    # Load the spliced and unspliced matrices as AnnData objects
    with open(path_spliced) as spliced:
        adata = anndata.read_csv(spliced)
    with open(path_unspliced) as unspliced:
        adata_unspliced = anndata.read_csv(unspliced)
    
    # Build the AnnData object with the two matrices
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = adata_unspliced.X

    if use_velocity_ground_truth:
        with open(path_velocity_ground_truth) as velocity:
            adata_velocity = anndata.read_csv(velocity)
        
        adata.layers['velocity'] = adata_velocity.X

    # Add the batch and simulation metadata
    obs_metadata = pd.read_csv(path_metadata,
                            index_col=0)

    adata.obs['batch'] = obs_metadata['batch'].astype('category')
    adata.obs['sim_time'] = obs_metadata['sim_time']

    # Load the cell UMAP embeddings
    umap_embeddings = pd.read_csv(path_umap_embeddings,
                                index_col=0)

    # Add the embeddings to the AnnData object
    adata.obsm['X_umap'] = umap_embeddings.loc[adata.obs_names].values


    if corrected_pca:
        # Load the cell PCA embeddings
        pca_embeddings = pd.read_csv(path_pca_embeddings,
                                    index_col=0)
        adata.obsm['X_pca'] = pca_embeddings.loc[adata.obs_names].values


    ### RNA VELOCITY ANALYSIS ###

    print("Starting RNA velocity analysis.")

    # Standard workflow
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    
    if corrected_pca:
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    else:
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print("Recovering dynamics. This can take a while.")

    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)

    if use_velocity_ground_truth:
        scv.tl.velocity_graph(adata, n_jobs=n_jobs, n_recurse_neighbors=n_recurse_neighbors)
    else:
        scv.tl.velocity(adata, mode=vel_mode)
    

    # Velocity embedding stream plot
    scv.pl.velocity_embedding_stream(adata, 
                                    color='batch',
                                    basis='umap', 
                                    dpi=200, 
                                    size=15,
                                    alpha=0.2,
                                    smooth  = 1,
                                    palette = ['blue', 'red'],
                                    save=f'../results/{velocity_plot_file_name}',
                                    legend_loc="right margin")
    
    # Compute latent time
    scv.tl.latent_time(adata)

    # Latent time and pseudotime plots
    scv.pl.scatter(adata, 
                color='latent_time', 
                color_map='gnuplot',
                dpi=200, 
                size=80,
                alpha=0.8,
                save=f'../results/{latent_time_plot_file_name}')

    scv.pl.scatter(adata, 
            color='sim_time', 
            color_map='gnuplot', 
            dpi=200,
            size=80,
            alpha=0.8,
            save=f'../results/{pseudotime_plot_file_name}')

    scv.pl.scatter(adata, 
                color='velocity_pseudotime', 
                color_map='gnuplot', 
                dpi=200,
                size=80,
                alpha=0.8,
                save=f'../results/{pseudotime_plot_file_name}')

    adata.write(f"../processed_objects/{adata_processed_object_file_name}")

    print("Your analysis has finished. Your results are in the 'results' directory.")