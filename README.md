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

Run the scripts in the order presented below to reproduce the analyses from the article. This can be achieved by copying and pasting the commands below in your computer terminal. Description of the tasks performed using each script is provided in the following sections.

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

