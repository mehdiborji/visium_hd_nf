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
        path "adata_${bin_size}_raw.h5ad" emit: adata_raw

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
    positions.set_index('barcode',inplace=True)

    y = 'pxl_row_in_fullres'
    x = 'pxl_col_in_fullres'

    coords = adata.obs.join(positions[[x,y]])[[x,y]].values

    adata.obsm['spatial'] = coords

    adata = adata[:,~adata.var.gene_ids.str.contains('DEPRECATED_')].copy()

    adata.write_h5ad("adata_${bin_size}_raw.h5ad", compression="gzip")

    """
}

workflow {
    
    MAKE_ADATA_BINNED(input_folder, "008um")
}

