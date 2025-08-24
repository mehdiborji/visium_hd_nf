/*
 * compute KNN using pynndescent followed by UMAP and Leiden clustering
 */
process doKNN_UMAP_Leiden {

    memory '32 GB'

    container 'visium_hd_env'

    publishDir params.results_folder, mode: 'symlink'

    input:
        path adata_SVD
        val bin_size
        val n_components
        val metric
        val n_neighbors
        val min_dist
        val spread
        val n_epochs
        val resolution
        

    output:
        path "adata_${bin_size}_clustered.h5ad", emit: adata_clustered
        path "${bin_size}_leiden_umap.pdf"
        path "${bin_size}_leiden_spatial.pdf"

    script:
    """
    #!/usr/bin/env python

    import umap
    import scanpy as sc
    import matplotlib.pyplot as plt

    adata = sc.read("${adata_SVD}")
    select_pcs = adata.obsm["X_pca"][:, :$n_components]

    reducer = umap.UMAP(
        metric="${metric}",
        n_neighbors=$n_neighbors,
        min_dist=$min_dist,
        spread=$spread,
        low_memory=False,
        n_components=2,
        verbose=True,
        n_epochs=$n_epochs,
    )

    embedding = reducer.fit_transform(select_pcs)

    adata.obsm["X_umap"] = embedding
    adata.obsp["connectivities"] = reducer.graph_

    sc.tl.leiden(
        adata,
        resolution=$resolution,
        flavor="igraph",
        n_iterations=1,
        obsp="connectivities"
    )

    plt.rcParams["figure.figsize"] = (15, 15)
    sc.pl.spatial(
        adata,
        color="leiden",
        spot_size=20,
        show=False,
        frameon=False
    )
    plt.savefig("${bin_size}_leiden_spatial.pdf", bbox_inches="tight")
    plt.show()
    plt.close()

    plt.rcParams["figure.figsize"] = (5, 5)
    sc.pl.umap(
        adata,
        color=['leiden'],
        legend_loc='on data',
        s=1,
        show=False,
        frameon=False
    )
    plt.savefig("${bin_size}_leiden_umap.pdf", bbox_inches="tight")
    plt.show()
    plt.close()

    adata.write_h5ad(f"adata_${bin_size}_clustered.h5ad", compression="gzip")
    """
}