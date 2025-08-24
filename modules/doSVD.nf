/*
 * do SVD on variable genes after scaling
 */
process doSVD {

    memory '48 GB'

    container 'visium_hd_env'

    publishDir params.results_folder, mode: 'symlink'

    input:
        path adata_filtered
        val bin_size
        val n_variable_genes
        val n_components
        val method

    output:
        path "adata_${bin_size}_SVD.h5ad", emit: adata_SVD
        path "${bin_size}_pca_variance.png"

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    import matplotlib.pyplot as plt
    from sklearn.decomposition import TruncatedSVD

    adata = sc.read("${adata_filtered}")
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable_rank < $n_variable_genes].copy()
    sc.pp.scale(adata, max_value=10)

    svd = TruncatedSVD(n_components=$n_components, algorithm="${method}")
    X_svd = svd.fit_transform(adata.X)
    adata.obsm["X_pca"] = X_svd
    adata.varm["PCs"] = svd.components_.T
    adata.uns["pca"] = dict(
        variance=svd.explained_variance_,
        variance_ratio=svd.explained_variance_ratio_,
    )

    adata = adata.raw.to_adata()

    plt.rcParams["figure.figsize"] = (5, 3)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=$n_components)
    plt.savefig("${bin_size}_pca_variance.png", bbox_inches="tight")
    plt.close()

    adata.write_h5ad(f"adata_${bin_size}_SVD.h5ad", compression="gzip")
    """
}