# Nextflow Pipeline For Processing Visium HD from 10x Genomics


## Setting Up the Pipeline

Install Docker
```
https://docs.docker.com/get-started/get-docker/
```

Install Nextflow
```
https://nextflow.io/docs/latest/install.html
```

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

The dataset can be from a variety of publicly available 10x genomics data from `https://www.10xgenomics.com/datasets`

```
cd $HOME
mkdir visium_hd_lymph_node
cd visium_hd_lymph_node
curl -O https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Lymph_Node_FFPE/Visium_HD_Human_Lymph_Node_FFPE_binned_outputs.tar.gz
curl -O https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Lymph_Node_FFPE/Visium_HD_Human_Lymph_Node_FFPE_segmented_outputs.tar.gz
tar -xzf Visium_HD_Human_Lymph_Node_FFPE_binned_outputs.tar.gz
tar -xzf Visium_HD_Human_Lymph_Node_FFPE_segmented_outputs.tar.gz
```

using `tree -L 2` the following should be shown within `visium_hd_lymph_node` folder:
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
Segmented outpus correspond to out from cell segmentation of H&E image data using StarDist and the corresponding python scipts from 
`https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd`



## Run the Pipeline

Execute the pipeline with:

```
nextflow run $HOME/nextflow/visium_hd_nf  --dataset_folder $HOME/visium_hd_lymph_node --bin_size "016um" --min_genes_per_cell 250
```