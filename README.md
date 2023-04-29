## Understand the cell-cell interaction in human tissue computationally

Understanding how cells communicate and interact with each other is critical for detangling the complexities of biological systems. Single-cell RNA sequencing (scRNA-seq) has become a powerful tool for deciphering high-resolution cell-cell communication networks.

ScRNA-seq is a technique that allows researchers to measure gene expression in individual cells. By analyzing the expression patterns of thousands of genes across many individual cells, scRNA-seq can reveal which cells communicate with each other, what types of signals are being exchanged, and how these signals are being transmitted.


<img width="571" alt="Screenshot 2023-04-29 at 1 34 32 AM" src="https://user-images.githubusercontent.com/117299113/235289729-32f80a37-5673-43ba-95bc-489cb08982f1.png">

## Overview of the Project

<img width="1007" alt="Screenshot 2023-04-29 at 1 54 43 AM" src="https://user-images.githubusercontent.com/117299113/235289265-29dc277d-6278-407f-9a2d-2fd5abb5c180.png">

The code R is created to perform data processing of scRNA seq data and use curated data for CellChat package[1]

## Installation

The script requires R 4.2.3 version comtatible packages. I have used below packages:
  1. anndata
  2. Seurat
  3. SingleCellSignalR
  4. SingleCellExperiment)
  5. Matrix
  6. methods
  7. CellChat
  8. patchwork
  9. dplyr
  10. topGO
  11. tibble
  12. biomaRt
  13. tidyr

The first part of the script is to create a folder structure.
1. Signaling
2. Siganling/Data
3. Signaling/Images
4. Signaling/Ouput

Data contains the input file i.e .rds and .h5ad file. Images folder contains all the images plotted in the process and Output folder gives the list of ligands and Receptor generated from the CellChat R package.

Note: While cloning the repo and running the code please make sure to have the directory is changed to the directory in which the repository is cloned.










