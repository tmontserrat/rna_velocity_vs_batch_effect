import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd

def perform_velocity_analysis(path_spliced, 
                            path_unspliced,
                            path_cell_types,
                            path_umap_embeddings,
                            path_pca_embeddings,
                            velocity_plot_file_name,
                            latent_time_plot_file_name,
                            pseudotime_plot_file_name,
                            velocity_confidence_plot_file_name,
                            velocity_confidence_data_file_name,
                            adata_processed_object_file_name,
                            velocity_graph_file_name,
                            paga_plot_file_name,
                            path_corrected_x=None,
                            n_jobs = 1,
                            vel_mode='dynamical',
                            n_recurse_neighbors=2
                            ):
    
    """Perform the velocity analysis for the spliced and unspliced RNA matrices supplied."""


    ### PREPARING THE AnnData OBJECT ###

    print("Preparing the AnnData object.")

    # Load the spliced, unspliced and x matrices as AnnData objects
    with open(path_spliced) as spliced:
        adata = anndata.read_csv(spliced)
    with open(path_unspliced) as unspliced:
        adata_unspliced = anndata.read_csv(unspliced)
    with open(path_corrected_x) as corrected_x:
        adata_x = anndata.read_csv(corrected_x)
    
    # Build the AnnData object with the two matrices
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = adata_unspliced.X

    # Load the cell type information
    cell_types = pd.read_csv(path_cell_types,
                            index_col='index')

    # Filter the barcodes to keep the cells in the current batch
    cell_types = cell_types.loc[adata.obs_names, :]

    # Add the information to the metadata
    adata.obs['cell_types'] = cell_types['clusters']

    # Load the cell PCA embeddings
    pca_embeddings = pd.read_csv(path_pca_embeddings,
                                index_col=0)

    # Add the UMAP embeddings
    adata.obsm['X_pca'] = pca_embeddings.loc[adata.obs_names].values

    # Load the cell UMAP embeddings
    umap_embeddings = pd.read_csv(path_umap_embeddings,
                            index_col=0)

    # Add the embeddings to the AnnData object
    adata.obsm['X_umap'] = umap_embeddings.loc[adata.obs_names].values

    ### RNA VELOCITY ANALYSIS ###

    print("Starting RNA velocity analysis.")

    # Standard workflow
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

    # Take the variables genes from scvelo
    gene_names = list(adata.var_names)

    # Keep the intersect from variables genes from
    # scvelo and those from anchors
    genes_to_keep = []
    for gene in adata_x.var_names:
        if gene in gene_names:
            genes_to_keep.append(gene)
    # Make the two AnnData objects equals in size
    adata_x = adata_x[:, genes_to_keep]
    adata = adata[:, genes_to_keep]
    adata.X = adata_x.X

    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print("Recovering dynamics. This can take a while.")

    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode=vel_mode)
    scv.tl.velocity_graph(adata, n_jobs=n_jobs, n_recurse_neighbors=n_recurse_neighbors)

    # Velocity embedding stream plot
    scv.pl.velocity_embedding_stream(adata, 
                                    color='cell_types',
                                    basis='umap', 
                                    dpi=200, 
                                    size=80,
                                    alpha=0.8,
                                    save=f'../results/{velocity_plot_file_name}')
    
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
                color='velocity_pseudotime', 
                color_map='gnuplot', 
                dpi=200,
                size=80,
                alpha=0.8,
                save=f'../results/{pseudotime_plot_file_name}')

    # Velocity length and velocity confidence plots and metrics
    scv.tl.velocity_confidence(adata)
    keys = ('velocity_length', 'velocity_confidence')

    scv.pl.scatter(adata, 
                c=keys, 
                cmap='coolwarm', 
                perc=[5, 95],
                save=f'../results/{velocity_confidence_plot_file_name}')
    
    df = adata.obs.groupby('cell_types')[keys].mean().T
    df.style.background_gradient(cmap='coolwarm', axis=1)
    df.to_csv(f"../results/{velocity_confidence_data_file_name}")

    

    # Graph plot
    scv.pl.velocity_graph(adata, 
                        threshold=.1, 
                        color="cell_types",
                        save=f"../results/{velocity_graph_file_name}")

    # PAGA plot
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    scv.tl.paga(adata, groups='cell_types')

    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,
            save=f"../results/{paga_plot_file_name}")

    adata.write(f"../data/{adata_processed_object_file_name}")


def perform_velocity_analysis_2_samples(path_spliced1,
                                    path_spliced2,
                                    path_unspliced1,
                                    path_unspliced2,
                                    integrated,
                                    path_cell_types,
                                    path_umap_embeddings,
                                    velocity_plot_file_name,
                                    latent_time_plot_file_name,
                                    pseudotime_plot_file_name,
                                    velocity_confidence_plot_file_name,
                                    velocity_confidence_data_file_name,
                                    adata_processed_object_file_name,
                                    velocity_graph_file_name,
                                    paga_plot_file_name,
                                    n_recurse_neighbors=2,
                                    path_pca_embeddings=None,
                                    path_corrected_x=None,
                                    n_jobs = 1,
                                    vel_mode='dynamical',
                                    ):
    """
    Perform the velocity analysis for the spliced 
    and unspliced RNA matrices supplied comming from two
    different samples.
    """

    ### PREPARING THE AnnData OBJECT ###

    print("Preparing the AnnData object.")

    # Load the spliced and unspliced matrices as AnnData objects
    with open(path_spliced1) as spliced1:
        adata = anndata.read_csv(spliced1)

    with open(path_unspliced1) as unspliced1:
        adata_u = anndata.read_csv(unspliced1)
    
    with open(path_spliced2) as spliced2:
        adata_spliced2 = anndata.read_csv(spliced2)

    with open(path_unspliced2) as unspliced2:
        adata_unspliced2 = anndata.read_csv(unspliced2)

    if integrated:
        print("integrated")
        # Load the corrected X matrix if the dataset has been integrated
        with open(path_corrected_x) as corrected_x:
            adata_x = anndata.read_csv(corrected_x)
        
        # Remove the batch label fromo the cell's names
        cell_names_x = []
        for cell in adata_x.obs_names:
            cell_names_x.append(cell.split("_")[1])

        # Update the cell's names of the x matrix
        adata_x.obs_names = cell_names_x

    # Keep track the cells from each batch
    cells_b1 = adata.obs_names
    cells_b2 = adata_spliced2.obs_names

    # Concatenate both batches
    adata = anndata.concat([adata, adata_spliced2])
    adata_u = anndata.concat([adata_u, adata_unspliced2])

    # Build the spliced and unspliced layers
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = adata_u.X

    # Add batch information in the cell's metadata
    adata.obs['batch'] = "NaN"
    adata.obs['batch'].loc[cells_b1] = "b1"
    adata.obs['batch'].loc[cells_b2] = "b2"

    # Load the cell types
    cell_types = pd.read_csv(path_cell_types,
                            index_col="index")

    # Cell taypes annotation
    adata.obs['cell_types'] = cell_types['clusters']

    if integrated:
        print("integrated")
        # Load the PCA embeddings from Seurat
        pca_embeddings = pd.read_csv(path_pca_embeddings,
                                    index_col=0)
        pca_embeddings_cell_names = list(pca_embeddings.index)
        
        # Remove the batch label from the cell names
        cell_names = []
        for cell in pca_embeddings_cell_names:
            cell_names.append(cell.split("_")[1])

        # Correct the names
        pca_embeddings.index = cell_names

        # Add the UMAP embeddings
        adata.obsm['X_pca'] = pca_embeddings.loc[adata.obs_names].values

    
    if integrated:
        print("integrated")
        # Load UMAP embeddings from Seurat
        umap_embeddings = pd.read_csv(path_umap_embeddings,
                                    index_col=0)
        umap_embeddings_cell_names = list(umap_embeddings.index)
        
        # Remove the batch label from the cell names
        cell_names = []
        for cell in umap_embeddings_cell_names:
            cell_names.append(cell.split("_")[1])

        # Correct the names
        umap_embeddings.index = cell_names

        # Add the UMAP embeddings
        adata.obsm['X_umap'] = umap_embeddings.loc[adata.obs_names].values
    else:
        # Load UMAP embeddings from Seurat
        umap_embeddings = pd.read_csv(path_umap_embeddings,
                                    index_col=0)
        umap_embeddings_cell_names = list(umap_embeddings.index)
        
        # Remove the batch label from the cell names
        cell_names = []
        for cell in umap_embeddings_cell_names:
            cell_names.append(cell.split("_")[1])

        # Correct the names
        umap_embeddings.index = cell_names

        # Add the UMAP embeddings
        adata.obsm['X_umap'] = umap_embeddings.loc[adata.obs_names].values

    ### RNA VELOCITY ANALYSIS ###

    print("Starting RNA velocity analysis.")

    # Standard workflow
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

    if integrated:
        print("integrated")
        # Take the variables genes from scvelo
        gene_names = list(adata.var_names)

        # Keep the intersect from variables genes from
        # scvelo and those from anchors
        keep_genes = []
        for gene in adata_x.var_names:
            if gene in gene_names:
                keep_genes.append(gene)
        # Make the two AnnData objects equals in size
        adata_x = adata_x[:, keep_genes]
        adata = adata[:, keep_genes]
        adata.X = adata_x.X

        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    else:
        print("no integrated")
        # Perform standard workflow
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    
    print("Recovering dynamics. This can take a while.")

    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode=vel_mode)
    scv.tl.velocity_graph(adata, n_jobs=n_jobs, n_recurse_neighbors=n_recurse_neighbors)

    # Velocity embedding stream plot
    scv.pl.velocity_embedding_stream(adata, 
                                    color='cell_types',
                                    basis='umap', 
                                    dpi=200, 
                                    size=80,
                                    alpha=0.8,
                                    save=f'../results/{velocity_plot_file_name}')

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
                color='velocity_pseudotime', 
                color_map='gnuplot', 
                dpi=200,
                size=80,
                alpha=0.8,
                save=f'../results/{pseudotime_plot_file_name}')

    scv.tl.velocity_confidence(adata)
    keys = ('velocity_length', 'velocity_confidence')

    scv.pl.scatter(adata, 
                c=keys, 
                cmap='coolwarm', 
                perc=[5, 95],
                save=f'../results/{velocity_confidence_plot_file_name}')
    
    df = adata.obs.groupby('cell_types')[keys].mean().T
    df.style.background_gradient(cmap='coolwarm', axis=1)
    df.to_csv(f"../results/{velocity_confidence_data_file_name}")

    adata.write(f"../data/{adata_processed_object_file_name}")

    # Graph plot
    scv.pl.velocity_graph(adata, 
                        threshold=.1, 
                        color="cell_types",
                        save=f"../results/{velocity_graph_file_name}")

    # PAGA plot
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

    scv.tl.paga(adata, groups='cell_types')

    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,
            save=f"../results/{paga_plot_file_name}")

    print("Your analysis has finished. Your results are in the 'results' directory.")