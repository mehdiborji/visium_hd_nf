/*
 * filter AnnData object based on quality metrics
 */

process filterAdata {

    container 'visium_hd_env'

    publishDir params.results_folder, mode: 'symlink'

    input:
        path adata_raw
        val bin_size
        val min_counts_per_gene
        val min_cells_per_gene
        val min_counts_per_cell
        val min_genes_per_cell

    output:
        path "adata_${bin_size}_filtered.h5ad", emit: adata_filtered
        path "${bin_size}_log10_n_counts_per_gene.png"
        path "${bin_size}_log10_n_cells_per_gene.png"
        path "${bin_size}_log10_n_counts_per_cell.png"
        path "${bin_size}_log10_n_genes_per_cell.png"
        path "${bin_size}_log10_n_counts_spatial.pdf"

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    adata = sc.read("${adata_raw}")

    ### filter genes ###
    sc.pp.filter_genes(adata, min_counts=$min_counts_per_gene)
    sc.pp.filter_genes(adata, min_cells=$min_cells_per_gene)

    adata.var["log10_n_counts"] = np.log10(adata.var["n_counts"])
    adata.var["log10_n_cells"] = np.log10(adata.var["n_cells"])

    adata.var = adata.var[["log10_n_counts", "log10_n_cells"]].copy()

    plt.rcParams["figure.figsize"] = (5, 3)
    sns.histplot(adata.var["log10_n_counts"], bins=50)
    plt.savefig("${bin_size}_log10_n_counts_per_gene.png", bbox_inches="tight")
    plt.close()

    plt.rcParams["figure.figsize"] = (5, 3)
    sns.histplot(adata.var["log10_n_cells"], bins=50)
    plt.savefig("${bin_size}_log10_n_cells_per_gene.png", bbox_inches="tight")
    plt.close()

    ### filter cells ###
    sc.pp.filter_cells(adata, min_genes=$min_genes_per_cell)
    sc.pp.filter_cells(adata, min_counts=$min_counts_per_cell)

    adata.obs["log10_n_counts"] = np.log10(adata.obs["n_counts"])
    adata.obs["log10_n_genes"] = np.log10(adata.obs["n_genes"])

    adata.obs = adata.obs[["log10_n_counts", "log10_n_genes"]].copy()

    plt.rcParams["figure.figsize"] = (5, 3)
    sns.histplot(adata.obs["log10_n_counts"], bins=50)
    plt.savefig("${bin_size}_log10_n_counts_per_cell.png", bbox_inches="tight")
    plt.close()

    plt.rcParams["figure.figsize"] = (5, 3)
    sns.histplot(adata.obs["log10_n_genes"], bins=50)
    plt.savefig("${bin_size}_log10_n_genes_per_cell.png", bbox_inches="tight")
    plt.close()

    plt.rcParams["figure.figsize"] = (15, 15)
    sc.pl.spatial(
        adata, 
        color="log10_n_counts", 
        spot_size=20, 
        vmax="p99.5", 
        show=False, 
        frameon=False
    )
    plt.savefig("${bin_size}_log10_n_counts_spatial.pdf", bbox_inches="tight")
    plt.show()

    sc.pp.highly_variable_genes(
        adata, 
        flavor="seurat_v3", 
        n_top_genes=10000, 
        inplace=True
    )

    mean_counts = (10**adata.var.log10_n_counts).mean()
    sc.pp.normalize_total(adata, target_sum=mean_counts)
    sc.pp.log1p(adata)

    adata.write_h5ad("adata_${bin_size}_filtered.h5ad", compression="gzip")
    """
}
