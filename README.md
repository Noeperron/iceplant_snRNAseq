  <p align="center">
    <h3 align="center"> Part I: <i> Mesembryanthemum crystallinum </i> (Ice Plant) Genome Analysis Pipeline
</h3>
  </p>

All single-nuclei RNA sequencing data analyses presented in the "Mesophyll-Specific Circadian Dynamics of CAM Induction in the Ice Plant Unveiled by Single-Cell Transcriptomics" manuscript were performed using the ice plant genome sequence and annotation of the Kirst lab at the University of Florida.
The IcePlantGenome.sh script included in this repository allows for reproduction of the workflow used for assembling, polishing and annotating the Ice Plant genome.

Before running the script, ensure you have downloaded the raw genome files from the NCBI SRA database using the provided accession numbers (to be added once the paper is published). Also, install the following software packages as per their documentation:

### Software pacckages and versions used for assembly, polishing, and quality control

nanoFilt/2.7.1

flye/2.9.3

pilon/1.24 

bowtie2/2.4.5 

samtools/1.18 

busco/5.3.0 

assembly-stats/1.0.1

### Software pacckages used for annotation

maker/2.31.6

repeatmodeler/2.0

ltr_finder/1.07

ltrretriever/2.5

pacbio/8.0.0

isoseq3/3.2.2

lima/2.7.1

augustus/3.4.0

gff3toolkit/2.0.3

iprscan/5.60

ncbi_blast/2.14.1

snap/2.0.3


## Workflow Overview

1. Read Filtering: The script begins by filtering raw nanopore reads using NanoFilt. This step includes removing reads shorter than 1000 bp, with an average quality score < Q10, and trimming the first 500 bp of each read.

2. Genome Assembly: Flye is used to assemble the filtered nanopore reads, followed by creation of a draft assembly file (Assembly.fasta).

3. Assembly Polishing: The draft assembly is polished using Pilon, which involves mapping Illumina short reads to the draft assembly using Bowtie2, and subsequent processing with Samtools.

4. Assembly Metrics and Completeness: The polished assembly (IcePlantAssembly_Polished.fasta) is analyzed for assembly metrics using assembly-stats and genome completeness is assessed using BUSCO.

5. Genome Annotation: The final assembly file is ready for annotation, performed using the MAKER2 pipeline.

## Running the script

```sh
bash IcePlantGenome.sh
```



  <p align="center">
    <h3 align="center"> Part II: Scripts to reproduce the snRNA-seq analysis of <i> Mesembryanthemum crystallinum </i> leaves under various conditions</h3>
  </p>

<!-- ABOUT THE PROJECT -->
## About the _M. crystallinum_ snRNA-seq analysis

This repository contains a series of R scripts used for comprehensive analysis of the single-nuclei RNA sequencing data and reproduction of the results published by Perron et al., 2024 (add doi once accepted). The analyses progress through various stages of data preprocessing, exploration, and advanced statistical testing. The scripts are organized in a specific order, starting from initial data quality control to trajectory inference.

## Downloading the data and installing R packages.

### Installing R packages used in the analyses.

All analyses were conducted using R v.4.2.2 and Python3 in macOS 13.2.1 A list of necessary R packages, and the commands to install them, can be found in the script *_0-InstallPackages.R_*.

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

### 1 - Data Pre-processing

1a. Raw Data Quality Control (1a-Raw_data_quality_control.R)

Purpose: Import and perform initial quality control on raw single-cell data. This includes the removal of chloroplast genes from the single-cell expression matrix to prevent contamination (need to download supplementary table 8 to use as input .csv file).

1b. Identify Empty Droplets (1b-Identify_empty_droplets.R)

Purpose: Detect and exclude empty droplets from the datasets, an essential step in single-cell data preprocessing.

1c. Doublet Removal (1c-Doublet_removal.R)

Purpose: Identify and remove doublets (multiple cells captured in a single droplet) from each sample, ensuring data quality and accuracy. Writes new 10X files from the cleaned objects to input in the downstream analyses.

### 2 - Data Processing and Integration

2a. Seurat Data Processing Integration (2a-Seurat_data_processing_integration.R)

Purpose: Further process the cleaned data, including integration and normalization, using Seurat objects.

2b. Calculate Number of Cells Per Cluster and Sample (2b-Nb_of_cells_per_cluster.R)

Purpose: Calculate and analyze the number of cells present in each identified cluster and respective sample. Output can be used to make the barplot presented in Figure 2C.


### 3 - Marker Identification and Visualization

3a. Marker Identification (3a-Marker_identification.R)

Purpose: Identify marker genes for various cell clusters to characterize different cell types within the data.

3b. Visualization of Marker Expression (3b-Visualization_marker_expression.R)

Purpose: Visualize the expression of marker genes across different clusters or cell types, often using heatmaps.

### 4 - Differential Gene Expression Analysis

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

## Add gene names to the output files

Some of the scripts mentioned above generate lists of genes with the gene IDs in the first column as output files. To add a column containing the gene names to these files, the "AddGeneName.py" script can be used. 

### Usage

Replace "file.csv" by the name of the input .csv file you would like to append gene names to. 
The script must be ran in a directory containing the protein sequence file of the ice plant genome presented in the Perron et al. article.
To run the script, execute the command below:

```sh
Python3 AddGeneName.py
```

