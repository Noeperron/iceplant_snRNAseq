<p align="center">
  <!-- <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <!-- <h2 align="center">Asc-Seurat</h2> -->

  <p align="center">
    <h3 align="center"> Scripts to reproduce the snRNA-seq analysis of <i> Mesembryanthemum crystallinum </i> leaves under various conditions</h3>
  </p>
</p>

<!-- ABOUT THE PROJECT -->
## About the _M. crystallinum_ snRNA-seq analysis

This repository contains a series of R scripts used for comprehensive analysis of the single-nuclei RNA sequencing data and reproduction of the results published by Perron et al., 2024 (add doi once accepted). The analyses progress through various stages of data preprocessing, exploration, and advanced statistical testing. The scripts are organized in a specific order, starting from initial data quality control to trajectory inference.

## Downloading the data and installing R packages.

### Installing R packages used in the analyses.

All analyses were conducted using R v.4.2.2 in macOS 13.2.1 A list of necessary R packages, and the commands to install them, can be found in the script *_0-InstallPackages.R_*.

Execute the command below to run the script:

```sh
Rscript 0-InstallPackages.R
```
### Downloading the data

The four single-nuclei RNA seq datasets (Dawn Salt, Dawn Control, Dusk Salt, Dusk Control) can be downloaded using the following link:

## Running the analysis

To reproduce the results from the article, it is recommended to run the scripts in the order listed below. Each script builds upon the results of the previous one, leading to a comprehensive analysis of the single-cell data. This can be achieved by copying and pasting the commands below in your computer terminal. Description of the tasks performed using each script is provided in the following sections.
The R files contain detailed comments on the actions performed by the various functions used.

```sh
Rscript 1a-Raw_data_quality_control.R
Rscript 1b-Identify_empty_droplets.R
Rscript 1c-Doublet_removal.R
Rscript 2a-Seurat_data_processing_integration.R
Rscript 2b-Nb_of_cells_per_cluster.R
Rscript 3a-Marker_identification.R
Rscript 3b-Visualization_marker_expression.R
Rscript 4a-DGE_between_clusters.R
Rscript 4b-DGE_within_one_cluster_between_treatments.R
Rscript 5a-Trajectory_Inference_Slingshot_Tradeseq.R
Rscript 5b-Different_expression_patterns_across_lineages.R
Rscript 5c-Plot_expression_along_trajectories.R
```

### 1-Data Pre-processing

1a. Raw Data Quality Control (1a-Raw_data_quality_control.R)

Purpose: Import and perform initial quality control on raw single-cell data. This includes the removal of chloroplast genes from the single-cell expression matrix to prevent contamination (need to download supplementary table 8 to use as input .csv file).

1b. Identify Empty Droplets (1b-Identify_empty_droplets.R)

Purpose: Detect and exclude empty droplets from the datasets, an essential step in single-cell data preprocessing.

1c. Doublet Removal (1c-Doublet_removal.R)

Purpose: Identify and remove doublets (multiple cells captured in a single droplet) from each sample, ensuring data quality and accuracy. Writes new 10X files from the cleaned objects to input in the downstream analyses.

### 2-Data Processing and Integration

2a. Seurat Data Processing Integration (2a-Seurat_data_processing_integration.R)

Purpose: Further process the cleaned data, including integration and normalization, using Seurat objects.

2b. Calculate Number of Cells Per Cluster and Sample (2b-Nb_of_cells_per_cluster.R)

Purpose: Calculate and analyze the number of cells present in each identified cluster and respective sample. Output can be used to make the barplot presented in Figure 2C.


### 3-Marker Identification and Visualization

3a. Marker Identification (3a-Marker_identification.R)

Purpose: Identify marker genes for various cell clusters to characterize different cell types within the data.

3b. Visualization of Marker Expression (3b-Visualization_marker_expression.R)

Purpose: Visualize the expression of marker genes across different clusters or cell types, often using heatmaps.

### 4-Differential Gene Expression Analysis

4a. DGE Between Clusters (4a-DGE_between_clusters.R)

Purpose: Perform differential gene expression analysis between different cell clusters. By default, reproduces the differential gene expression analysis between cluster 5 and 6.

4b. DGE Within One Cluster Between Treatments (4b-DGE_within_one_cluster_between_treatments.R)

Purpose: Conduct differential expression analysis between the cells within a single cluster from the different experimental treatments (i.e. gene expression in cells originating from the control samples vs. salt-treated samples)

### 5-Cell Trajectory Inference

5a. Trajectory Inference Using Slingshot and tradeSeq (5a-Trajectory_Inference_Slingshot_Tradeseq.R)

Purpose: Infer cellular trajectories, tracing transitions between different cell states. By default, reproduces the trajectory analysis with cluster 6 set as the starting point of the trajectory.

5b. Different expression patterns across trajectories (5b-Different_expression_patterns_across_lineages.R)

Purpose: Identification of genes presenting drastically different expression patterns between the two trajectories.

5c. Plot Expression Along Trajectories (5c-Plot_expression_along_trajectories.R)

Purpose: Visualize gene expression along identified cell trajectories. Inputs can be a single gene or a list of genes provided in the first column of a .csv file.
