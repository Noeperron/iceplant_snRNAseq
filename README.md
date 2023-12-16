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
In this repository, you will find the source code and instructions to reproduce the results published by Perron et al., 2024 (add doi once accepted).

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

### Nuclei quality control and filtering in each sample separately

The command below:

1. Imports each sample's raw data into R.
2. Removes chloroplast genes from the single-cell expression matrix to prevent contamination (need to download supplementary table 8)
3. Creates a Seurat object for each sample
4. Processes the data (normalization, dimensionality reduction, clustering).
5. Uses emptyDrops() from the DropletUtils package to determine a low-quality nuclei threshold for each sample
6. Removes doublets using ScDblFinder
7. Writes new 10X files from the cleaned objects

```sh
Rscript 1-raw_data_quality_control.R
```

### Integration of the datasets and clustering

