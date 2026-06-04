# Single-Cell RNA Sequencing Analysis Pipeline

================================================================================
OVERVIEW
================================================================================
This repository contains the R code, raw count matrix data, and pre-processed
Seurat objects used for single-cell RNA sequencing (scRNA-seq) analysis of
mouse midpalatal suture tissue, as described in the associated manuscript.

The analysis script (scRNAseq_analysis.R) covers:
  - Data loading, quality control, and preprocessing (Sections 1-2)
  - Cell type annotation based on canonical marker gene expression (Section 2)
  - Cell type marker visualization — full palate and mesenchymal subclusters
    (Sections 4-5)
  - Ligand-receptor interaction analysis (CellChat) (Section 6)

NOTE: Sections 1-2 reproduce the full QC and annotation pipeline from raw
count matrices. Users who wish to reproduce downstream analysis only may
start from Section 3 using the pre-processed RDS objects provided on GEO.

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
  - CellChat         v1.6.1
  - clusterProfiler  v4.14.6
  - org.Mm.eg.db     v3.20.0 (Bioconductor)
  - ggplot2          v4.0.0
  - dplyr            v1.1.4
  - stringr          v1.5.1
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
  install.packages(c("Seurat","ggplot2","dplyr","stringr","patchwork",
                     "ggalluvial","igraph","gridExtra","magrittr"))

Step 3 - Install Bioconductor packages:
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))

Step 4 - Install CellChat:
  devtools::install_github("sqjin/CellChat")

Typical installation time: 20-40 minutes

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
To reproduce the full analysis from raw data:
  Run Sections 1-2 for QC, preprocessing, and cell type annotation,
  then continue with Sections 3-7 for downstream analysis.

To reproduce downstream analysis only (from pre-processed objects):
  Download RDS files from GEO (accession: GSE325000), place them in
  ./data/rds/, and run from Section 3 onward.

All parameters in the script match exactly those reported in the manuscript
Methods section.

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
