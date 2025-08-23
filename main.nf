#!/usr/bin/env nextflow

params.dataset_folder = "/Users/mb/visium/Visium_HD_Human_Lymph_Node_FFPE"
input_folder = Channel.of(params.dataset_folder)

params.results_folder = "results"

process MAKE_ADATA_BINNED {

    container 'visium_hd_env'

    publishDir params.results_folder, mode: 'symlink'

    input:
        path input_folder
        val bin_size

    output:
        path "adata_${bin_size}_raw.h5ad", emit: adata_raw

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    import pandas as pd

    binned_folder = "${input_folder}/binned_outputs/square_${bin_size}"

    raw_h5 = f"{binned_folder}/raw_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(raw_h5)
    adata.var_names_make_unique()

    spatial_pos = f"{binned_folder}/spatial/tissue_positions.parquet"
    positions = pd.read_parquet(spatial_pos)
    positions.set_index("barcode", inplace=True)

    y = "pxl_row_in_fullres"
    x = "pxl_col_in_fullres"

    coords = adata.obs.join(positions[[x, y]])[[x, y]].values

    adata.obsm["spatial"] = coords

    adata = adata[:, ~adata.var.gene_ids.str.contains("DEPRECATED_")].copy()

    adata.write_h5ad("adata_${bin_size}_raw.h5ad", compression="gzip")
    """
}


process FILTER_ADATA {

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

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=10000, inplace=True)

    mean_counts = (10**adata.var.log10_n_counts).mean()
    sc.pp.normalize_total(adata, target_sum=mean_counts)
    sc.pp.log1p(adata)

    adata.write_h5ad("adata_${bin_size}_filtered.h5ad", compression="gzip")
    """
}


workflow {
    
    MAKE_ADATA_BINNED(input_folder, "008um")

    FILTER_ADATA(MAKE_ADATA_BINNED.out.adata_raw, "008um", 10, 10, 15, 20)

    
}

