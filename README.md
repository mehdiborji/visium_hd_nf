# Nextflow Pipeline For Processing Visium HD from 10x Genomics


## Install Prerequisites

Docker and Nextflow are required prior to running the pipeline
```
https://docs.docker.com/get-started/get-docker/
```
```
https://nextflow.io/docs/latest/install.html
```
## Set Up the Pipeline

Download Nextflow Pipeline
```
cd $HOME
mkdir nextflow
https://github.com/mehdiborji/visium_hd_nf
cd visium_hd_nf
```
Using `tree` you should be able to see the following:
```
.
├── conda.yml
├── Dockerfile
├── main.nf
├── modules
│   ├── doKNN_UMAP_Leiden.nf
│   ├── doSVD.nf
│   ├── filterAdata.nf
│   └── makeAdata.nf
└── nextflow.config
```

Within `visium_hd_nf` build Docker image:

```
docker build -t visium_hd_env .
```


## Download Input Files

The dataset can be obtained from a variety of publicly available 10x Genomics datasets at [10x Genomics](https://www.10xgenomics.com/datasets).

Here, we use a recent Visium HD run on an FFPE Human Lymph Node. This dataset is particularly challenging for cell calling due to the small size and dense cellularity of lymphoid tissue. Moreover, in terms of total transcript counts, this dataset is not among the best and serves as a representative example of a suboptimal run.


```
cd $HOME
mkdir visium_hd_lymph_node
cd visium_hd_lymph_node
curl -O https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Lymph_Node_FFPE/Visium_HD_Human_Lymph_Node_FFPE_binned_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Lymph_Node_FFPE/Visium_HD_Human_Lymph_Node_FFPE_segmented_outputs.tar.gz
tar -xzf Visium_HD_Human_Lymph_Node_FFPE_binned_outputs.tar.gz
tar -xzf Visium_HD_Human_Lymph_Node_FFPE_segmented_outputs.tar.gz
```
Using tree -L 2, you should see the following inside the visium_hd_lymph_node folder:
```
.
├── binned_outputs
│   ├── square_002um
│   ├── square_008um
│   └── square_016um
├── segmented_outputs
│   ├── analysis
│   ├── cell_segmentations.geojson
│   ├── cloupe.cloupe
│   ├── filtered_feature_cell_matrix
│   ├── filtered_feature_cell_matrix.h5
│   ├── graphclust_annotated_cell_segmentations.geojson
│   ├── graphclust_annotated_nucleus_segmentations.geojson
│   ├── nucleus_segmentations.geojson
│   ├── raw_feature_cell_matrix
│   ├── raw_feature_cell_matrix.h5
│   └── spatial
├── Visium_HD_Human_Lymph_Node_FFPE_binned_outputs.tar.gz
└── Visium_HD_Human_Lymph_Node_FFPE_segmented_outputs.tar.gz
```
The `segmented_outputs` correspond to results from cell segmentation of H&E image data using StarDist and Python scripts from:
`https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd`


## Run the Pipeline

Execute the pipeline using the following sample commands for different binning strategies:

For binned outputs: 
```
nextflow run $HOME/nextflow/visium_hd_nf \
    --dataset_folder $HOME/visium_hd_lymph_node \
    --bin_size "016um" \
    --min_genes_per_cell 250
```

```
nextflow run $HOME/nextflow/visium_hd_nf \
    --dataset_folder $HOME/visium_hd_lymph_node \
    --bin_size "008um" \
    --min_genes_per_cell 60
```

While the pipeline can run with a `002um` bin size, it is generally not advisable due to the sparse nature of the data and the very large number of spots.

For segmented outputs: 
```
nextflow run $HOME/nextflow/visium_hd_nf \
    --dataset_folder $HOME/visium_hd_lymph_node \
    --bin_size "cell" \
    --min_genes_per_cell 20
```


## System Requirements and Notes

With default pipeline parameters, the analysis at `016um` runs in about two minutes on an Intel Core i9 2019 MacBook Pro and uses up to ~8 GB of RAM.  
At `008um`, the analysis takes around six minutes and requires up to ~25 GB of RAM.  
For segmented output (`cell` mode), the analysis takes about ten minutes and uses up to ~40 GB of RAM.