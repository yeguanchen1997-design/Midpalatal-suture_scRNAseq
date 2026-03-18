# Midpalatal Suture Single-Cell RNA Sequencing Analysis Pipeline

================================================================================
OVERVIEW
================================================================================
This repository contains the R code used for single-cell RNA sequencing (scRNA-seq) analysis of
mouse midpalatal tissue, as described in the associated manuscript.

NOTE ON DATA PREPROCESSING:
Raw sequencing data were processed using Cell Ranger 5.0.0 (alignment to GRCm38)
and quality-controlled using Seurat v5.3.0 (gene filtering: >3 cells; cell
filtering: 200 < nFeature_RNA < 8,000, nCount_RNA < 60,000, percent.mt < 12%).
Cell cycle effects were regressed using CellCycleScoring. Clustering was
performed with FindNeighbors (dims=1:20) and FindClusters (resolution=0.8),
followed by UMAP visualization. These steps are described in full in the Methods
section of the manuscript. The resulting annotated Seurat objects are provided
directly as .rds files for reproducibility.

The analysis script (scRNAseq_analysis.R) covers:
  - Cell type marker visualization (palate + mesenchymal subclusters)
  - GO Biological Process enrichment analysis (clusterProfiler)
  - Pseudotime trajectory inference (Monocle 3)
  - Ligand-receptor interaction analysis (CellChat)

================================================================================
SYSTEM REQUIREMENTS
================================================================================
Operating System:
  - Tested on: Windows 10 x64 (build 19045)
  - Compatible with: macOS 13.x, Ubuntu 20.04 LTS

R version:
  - R 4.4.2 (2024-10-31 ucrt)

Required R packages and versions:
  - Seurat           v5.3.0
  - SeuratWrappers   v0.4.0
  - CellChat         v1.6.1
  - monocle3         v1.4.26
  - clusterProfiler  v4.14.6
  - org.Mm.eg.db     v3.20.0 (Bioconductor)
  - ggplot2          v4.0.0
  - dplyr            v1.1.4
  - patchwork        v1.3.2
  - ggalluvial       v0.12.5
  - igraph           v2.2.0
  - gridExtra        v2.3
  - magrittr         v2.0.3

Non-standard hardware:
  - None required
  - Minimum 16 GB RAM recommended (32 GB preferred for CellChat analysis)

================================================================================
INSTALLATION
================================================================================
Step 1 - Install R from https://cran.r-project.org/

Step 2 - Install CRAN packages:
  install.packages(c("Seurat","ggplot2","dplyr","patchwork",
                     "ggalluvial","igraph","gridExtra","magrittr"))

Step 3 - Install SeuratWrappers:
  devtools::install_github("satijalab/seurat-wrappers")

Step 4 - Install Bioconductor packages:
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))

Step 5 - Install monocle3:
  BiocManager::install(c("BiocGenerics","DelayedArray","DelayedMatrixStats",
                         "limma","lme4","S4Vectors","SingleCellExperiment",
                         "SummarizedExperiment","batchelor","HDF5Array",
                         "terra","ggrastr"))
  devtools::install_github("cole-trapnell-lab/monocle3")

Step 6 - Install CellChat:
  devtools::install_github("sqjin/CellChat")


================================================================================
INPUT DATA
================================================================================
All input data files are available from NCBI GEO (accession: GSE325000).

Download the following files from GEO and place them in the indicated
directories before running the script:

  ./data/raw/
    barcodes.tsv    Cell barcodes produced by Cell Ranger 5.0.0
    features.tsv    Gene feature list (aligned to GRCm38)
    matrix.mtx      Sparse UMI count matrix

  ./data/rds/
    palate.rds
      Full annotated Seurat object (13 cell types: Mesenchymal cells,
      Epithelial cells, Pericytes, Endothelial cells, Chondrocytes,
      Neuronal cells, Glial cells, B cells, Neutrophils, Macrophages,
      T cells, Erythrocytes, Cycling cells).

    mesenchyme.rds
      Mesenchymal subcluster Seurat object (7 subclusters, dims=1:20,
      resolution=0.5).

    cellchat.rds
      Fully computed CellChat object.

To load the raw count matrix:
  Seurat::Read10X(data.dir = "./data/raw/")

================================================================================
DEMO
================================================================================
The demo runs the two most representative steps of the pipeline — a UMAP
visualization and a CellChat interaction plot — to verify correct installation
and data loading.

  1. Download palate.rds and cellchat.rds from GEO (accession: GSE325000)
     and place them in ./data/rds/
  2. Open R and run the following:

     palate   <- readRDS("./data/rds/palate.rds")
     cellchat <- readRDS("./data/rds/cellchat.rds")

     # Demo 1: UMAP of all annotated cell types
     library(Seurat)
     DimPlot(palate, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

     # Demo 2: FGF signaling network (CellChat)
     library(CellChat)
     netVisual_heatmap(cellchat, signaling = "FGF",
                       color.heatmap = "Reds", measure = "weight")

Expected output:
  - A UMAP plot showing 13 annotated cell populations
  - A heatmap showing FGF signaling communication weights between cell types

================================================================================
INSTRUCTIONS FOR USE
================================================================================
Load the three provided .rds files and run scRNAseq_analysis.R
in order (Sections 0-7).

All parameters in the script match exactly those reported in the manuscript
Methods section. No modifications are required for reproduction.


================================================================================
FILE CONTENTS 
================================================================================
  scRNAseq_analysis.R              -- Main annotated analysis script
  README(Github).txt               -- This file
  LICENSE                          -- MIT License

  Note: Input data files (barcodes.tsv, features.tsv, matrix.mtx,
  palate.rds, mesenchyme.rds, cellchat.rds) are not included in this
  repository. Download them from GEO (accession: GSE325000) and place
  them in ./data/raw/ and ./data/rds/ respectively.

================================================================================
LICENSE
================================================================================
This code is released under the MIT License (OSI approved).
See LICENSE file for full terms.

================================================================================
CITATION
================================================================================
If you use this code or data, please cite the associated manuscript:
Restoring Midpalatal Suture Homeostasis via Local Targeting of Gli1+ Cells
Rescues Maxillary Hypoplasia.

================================================================================
CONTACT
================================================================================
For questions, please contact: yeguanchen@zju.edu.cn

