#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { makeAdata } from './modules/makeAdata.nf'
include { filterAdata } from './modules/filterAdata.nf'
include { doSVD } from './modules/doSVD.nf'
include { doKNN_UMAP_Leiden } from './modules/doKNN_UMAP_Leiden.nf'

workflow {
    
    input_folder = Channel.of(params.dataset_folder)

    makeAdata(input_folder, params.bin_size)

    filterAdata(
        makeAdata.out.adata_raw,
        params.bin_size,
        params.min_counts_per_gene,
        params.min_cells_per_gene,
        params.min_counts_per_cell,
        params.min_genes_per_cell
    )

    doSVD(
        filterAdata.out.adata_filtered,
        params.bin_size,
        params.n_variable_genes,
        params.n_components,
        params.svd_method
    )

    doKNN_UMAP_Leiden(
        doSVD.out.adata_SVD,
        params.bin_size,
        params.umap_n_components,
        params.umap_metric,
        params.n_neighbors,
        params.min_dist,
        params.spread,
        params.n_epochs,
        params.leiden_resolution
    )
}

