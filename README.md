# GAGE-seq

__Concurrent profiling of multiscale 3D genome organization and gene expression in single mammalian cells__

Our preprint is available at [bioRxiv](https://doi.org/10.1038/s41588-022-01256-z).

The organization of mammalian genomes within the nucleus features a complex, multiscale three-dimensional (3D) architecture. The functional significance of these 3D genome features, however, remains largely elusive due to limited single-cell technologies that can concurrently profile genome organization and transcriptional activities. Here, we report GAGE-seq, a highly scalable, robust single-cell co-assay that simultaneously measures 3D genome structure and transcriptome within the same cell. Employing GAGE-seq on mouse brain cortex and human bone marrow CD34+ cells, we comprehensively characterized the intricate relationships between 3D genome and gene expression. We found that these multiscale 3D genome features collectively inform cell type-specific gene expressions, hence contributing to defining cell identity at the single-cell level. Integration of GAGE-seq data with spatial transcriptomic data revealed in situ variations of the 3D genome in mouse cortex. Moreover, our observations of lineage commitment in normal human hematopoiesis unveiled notable discordant changes between 3D genome organization and gene expression, underscoring a complex, temporal interplay at the single-cell level that is more nuanced than previously appreciated. Together, GAGE-seq provides a powerful, cost-effective approach for interrogating genome structure and gene expression relationships at the single-cell level across diverse biological contexts.

This repository contains the code for processing GAGE-seq raw data and integrating GAGE-seq data with other data.

# Running the integration notebooks

## Dependencies

We strongly recommend installing everything via Conda. Below is the list of required packages.

- r-seurat=4.1.1
- r-ggplot2=3.3.6
- jupyter>=1.0.0 or jupyterlab>=3.0.0
- tqdm=4.62.0
- numpy=1.20.3
- scikit-learn=0.24.1
- scipy=1.9.0
- statsmodels=0.12.2
- matplotlib_venn
- anndata=0.9.1
- scanpy=1.9.3
- pandas=2.0.0
- umap-learn=0.5.1
- matplotlib=3.6.1
- seaborn=0.11.1

## Input data

Our integration framework takes two or more co-assayed single-cell datasets as input. The co-assayed single-cell datasets are required to share one modality. For each dataset, 1) a processed file for each modality and 2) metadata of cells are required. The processed file will be internally generalized to a cell-by-feature matrix. For example, the matrix is the expression count matrix for scRNA-seq, where each gene is a feature; the scATAC profiles of multiple chromosomes will be concatenated along chromosomes and each feature is one genomic locus at a particular chromosome; the scHi-C contact map for each cell will be flattened into a vector and concatenated just as scATAC profiles, where each feature is a locus pair.

## Overview of the framework

The framework has two major steps
1. Integrate the datasets based on the shared modality using Seurat. Seurat can be replaced with any other integration algorithm.
1. Train a regression model that predicts a non-shared modality from the embedding space.

## Examples

As an example, we integrate two datasets: the GAGE-seq data and the Paired-seq data that share the scRNA modality. We use Seurat and the K-NN regression model in this example. The two major steps in the framework are implemented in the notebooks [integrate-PairedTag-Seurat-mBC.ipynb](./scripts_analysis/integrate-PairedTag-Seurat-mBC.ipynb) and [integrate-PairedTag-Seurat-mBC-post.ipynb](./scripts_analysis/integrate-PairedTag-Seurat-mBC-post.ipynb), respectively.

In addition to the framework, we also included example analyses enabled by the integration in [integrate-PairedTag-Seurat-mBC-post.ipynb](./scripts_analysis/integrate-PairedTag-Seurat-mBC-post.ipynb). Specifically, in section 5, we investigate the joint influence of 3D genome structure and accessibility on expression, and divide genes into 4 groups based on the influence patterns. In section 6, we visualize the differential accessibility at a user-selected region between inhibitory subtypes.

# Contact
Please email [tianming@andrew.cmu.edu](tianming@andrew.cmu.edu) or raise an issue in the github repository with any questions about installation or usage.
