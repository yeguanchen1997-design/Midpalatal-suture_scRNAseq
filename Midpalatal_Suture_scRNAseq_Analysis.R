################################################################################
# Single-Cell RNA Sequencing Analysis Pipeline
# Manuscript: Restoring Midpalatal Suture Homeostasis via Local Targeting of
#             Gli1+ Cells Rescues Maxillary Hypoplasia
# Description: Analysis of mouse midpalatal suture scRNA-seq data
#
# QC parameters applied: 200 < nFeature_RNA < 8,000; nCount_RNA < 60,000;
#                        percent.mt < 12%; FindNeighbors(dims=1:20);
#                        FindClusters(resolution=0.8)
#
# Required input files:
#   Raw data (place in ./data/raw/):
#     - barcodes.tsv, features.tsv, matrix.mtx  (Cell Ranger output)
#   Pre-processed RDS objects (place in ./data/rds/):
#     - palate.rds     : Full annotated Seurat object (all cell types)
#     - mesenchyme.rds : Mesenchymal subcluster Seurat object
#     - cellchat.rds   : CellChat object (ligand-receptor analysis)
#
# R version: 4.4.2 (2024-10-31 ucrt)
# Contact:   yeguanchen@zju.edu.cn
################################################################################


# ==============================================================================
# 0. INSTALL AND LOAD PACKAGES
# ==============================================================================

suppressMessages({
  if (!require(ggplot2))    install.packages("ggplot2")
  if (!require(patchwork))  install.packages("patchwork")
  if (!require(ggalluvial)) install.packages("ggalluvial")
  if (!require(igraph))     install.packages("igraph")
  if (!require(dplyr))      install.packages("dplyr")
  if (!require(gridExtra))  install.packages("gridExtra")
  if (!require(stringr))    install.packages("stringr")
  options(stringsAsFactors = FALSE)
  options(future.globals.maxSize = 2 * 1024^3)
})

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(gridExtra)
library(SingleCellExperiment)
library(magrittr)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require(org.Mm.eg.db))
  BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(clusterProfiler)


# ==============================================================================
# 1. DATA LOADING, QC, AND PREPROCESSING
# ==============================================================================

# Load raw count matrix (Cell Ranger output)
palate.data <- Read10X(data.dir = "./data/raw/")

# Create Seurat object; retain genes expressed in >= 3 cells and
# cells with >= 200 detected genes
palate <- CreateSeuratObject(
  counts       = palate.data,
  project      = "palate",
  min.cells    = 3,
  min.features = 200
)

# Calculate mitochondrial gene percentage (mouse mt genes: lowercase "mt-")
palate[["percent.mt"]] <- PercentageFeatureSet(palate, pattern = "^mt-")

# Show QC metrics for the first 10 cells
head(palate@meta.data, 10)

# Visualize QC metrics before filtering
VlnPlot(palate,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol     = 3,
        pt.size  = 0.5)

# Scatter plots to inspect relationships between QC metrics
plot1 <- FeatureScatter(palate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(palate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells: 200 < nFeature_RNA < 8,000; nCount_RNA < 60,000; percent.mt < 12%
palate <- subset(palate,
                 subset = nFeature_RNA > 200  &
                          nFeature_RNA < 8000 &
                          nCount_RNA   < 60000 &
                          percent.mt   < 12)

# Normalize data using log normalization
palate <- NormalizeData(palate,
                        normalization.method = "LogNormalize",
                        scale.factor         = 10000)

# Identify highly variable features
palate <- FindVariableFeatures(palate,
                               selection.method = "vst",
                               nfeatures        = 2000)

# Visualize top 10 most variable genes
top10 <- head(VariableFeatures(palate), 10)
plot1 <- VariableFeaturePlot(palate)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Cell cycle scoring using built-in Seurat marker lists (Tirosh et al., 2015)
# Convert human gene symbols to mouse format (title case)
s.genes   <- cc.genes$s.genes   %>% str_to_lower() %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_lower() %>% str_to_title()
palate <- CellCycleScoring(palate,
                           s.features   = s.genes,
                           g2m.features = g2m.genes,
                           set.ident    = TRUE)

# Scale data and regress out cell cycle scores
# features = rownames(palate) scales all genes (required for heatmap)
palate <- ScaleData(palate,
                    vars.to.regress = c("S.Score", "G2M.Score"),
                    features        = rownames(palate))

# Run PCA
palate <- RunPCA(palate,
                 features = VariableFeatures(palate),
                 npcs     = 50)

# Inspect PCA results and determine dimensionality
print(palate[["pca"]], dims = 1:5, nfeatures = 5)
ElbowPlot(palate, ndims = 40)

# Optional: PCA heatmap to assess signal in higher PCs
# DimHeatmap(palate, dims = 20:26, cells = 500, balanced = TRUE)

# Optional: JackStraw test for statistical determination of significant PCs
# NOTE: This process can take a long time for large datasets
# palate <- JackStraw(palate, num.replicate = 100)
# palate <- ScoreJackStraw(palate, dims = 1:20)
# JackStrawPlot(palate, dims = 1:20)

# Cluster cells (dims = 1:20; resolution = 0.8)
palate <- FindNeighbors(palate, dims = 1:20)
palate <- FindClusters(palate,  resolution = 0.8)
head(Idents(palate), 5)

# Run UMAP for visualization
palate <- RunUMAP(palate, dims = 1:20)
DimPlot(palate, label = TRUE, label.size = 4, pt.size = 0.5)

# Find marker genes for all clusters
palate.markers <- FindAllMarkers(palate,
                                 only.pos        = TRUE,
                                 min.pct         = 0.25,
                                 logfc.threshold = 0.25)

# Top 10 markers per cluster for heatmap
top10 <- palate.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(palate, features = top10$gene, label = TRUE, size = 3) + NoLegend()


# ==============================================================================
# 2. CELL TYPE ANNOTATION
# ==============================================================================

# Cell types were annotated based on top marker genes identified by
# FindAllMarkers() combined with canonical markers from published literature.
# The following marker genes were used to define each cell type:
#
#   Mesenchymal cells  : Col1a1, Mfap5, Pdgfra
#   Epithelial cells   : Krt8, Krt18, Krt14
#   Pericytes          : Acta2, Myh11, Tagln
#   Endothelial cells  : Cdh5, Ptprb, Cyyr1
#   Chondrocytes       : Col2a1, Col9a1, Acan
#   Neuronal cells     : Fstl5, Tubb3, Scn9a
#   Glial cells        : Plp1, Gjc3, Sox10
#   B cells            : Rag1, Gm8369, Ms4a1
#   Neutrophils        : Stfa2, Stfa2l1, Retnlg
#   Macrophages        : Ms4a6c, Ccl9, Ccr2
#   T cells            : Il2rb, Trbc1, Ctsw
#   Erythrocytes       : Gypa, Slc4a1, Rhag
#   Cycling cells      : Top2a, Hist1h1a, Hist1h2ab
#
# NOTE: Cluster-to-cell-type assignments depend on clustering results which
# may vary across computing environments. Users should visualize the marker
# genes below and apply RenameIdents() according to their own clustering output.

# Visualize marker genes to confirm cluster identity
FeaturePlot(palate,
            features = c("Col1a1", "Mfap5",    "Pdgfra",
                         "Krt8",   "Krt18",    "Krt14",
                         "Acta2",  "Myh11",    "Tagln",
                         "Cdh5",   "Ptprb",    "Cyyr1",
                         "Col2a1", "Col9a1",   "Acan",
                         "Fstl5",  "Tubb3",    "Scn9a",
                         "Plp1",   "Gjc3",     "Sox10",
                         "Rag1",   "Gm8369",   "Ms4a1",
                         "Stfa2",  "Stfa2l1",  "Retnlg",
                         "Ms4a6c", "Ccl9",     "Ccr2",
                         "Il2rb",  "Trbc1",    "Ctsw",
                         "Gypa",   "Slc4a1",   "Rhag",
                         "Top2a",  "Hist1h1a", "Hist1h2ab"),
            ncol = 4)

# Apply cell type labels using RenameIdents() based on marker gene visualization
# Example (update cluster numbers to match your clustering results):
# palate <- RenameIdents(palate,
#   "X"  = "Mesenchymal cells",
#   "X"  = "Epithelial cells",
#   "X"  = "Pericytes",
#   "X"  = "Endothelial cells",
#   "X"  = "Chondrocytes",
#   "X"  = "Neuronal cells",
#   "X"  = "Glial cells",
#   "X"  = "B cells",
#   "X"  = "Neutrophils",
#   "X"  = "Macrophages",
#   "X"  = "T cells",
#   "X"  = "Erythrocytes",
#   "X"  = "Cycling cells"
# )
palate$celltype <- Idents(palate)

# UMAP with annotated cell type labels
DimPlot(palate, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Subset mesenchymal cells for subclustering
mesenchyme <- subset(palate, idents = "Mesenchymal cells")

# Re-cluster mesenchymal subpopulations
mesenchyme <- FindNeighbors(mesenchyme, dims = 1:20)
mesenchyme <- FindClusters(mesenchyme,  resolution = 0.5)
mesenchyme <- RunUMAP(mesenchyme, dims = 1:20)
DimPlot(mesenchyme, label = TRUE, label.size = 4, pt.size = 0.5)

# Save processed objects
saveRDS(palate,     "./data/rds/palate.rds")
saveRDS(mesenchyme, "./data/rds/mesenchyme.rds")


# ==============================================================================
# 3. LOAD PRE-PROCESSED RDS OBJECTS
# ==============================================================================
# NOTE: If running from pre-processed objects, start here and skip Sections 1-2.

palate     <- readRDS("./data/rds/palate.rds")       # Full annotated Seurat object
mesenchyme <- readRDS("./data/rds/mesenchyme.rds")   # Mesenchymal subcluster object
cellchat   <- readRDS("./data/rds/cellchat.rds")     # CellChat object


# ==============================================================================
# 4. FULL PALATE — CELL TYPE OVERVIEW
# ==============================================================================

# Check cell type levels and counts
levels(palate)
table(Idents(palate))

# UMAP of all annotated cell types
DimPlot(palate, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Find marker genes for all cell types
palate.markers <- FindAllMarkers(
  palate,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Top 15 markers per cluster by avg_log2FC
top15palate <- palate.markers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

# Feature plots for key cell-type marker genes
marker_genes_palate <- c("Col1a1",  # Mesenchymal marker
                         "Mfap5",   # Progenitor cell marker
                         "Runx2",   # Osteogenic cell marker
                         "Krt14",   # Epithelial marker
                         "Acta2",   # Pericyte/smooth muscle marker
                         "Tubb3")   # Neuronal marker

for (gene in marker_genes_palate) {
  print(FeaturePlot(palate, features = gene, reduction = "umap"))
}


# ==============================================================================
# 5. MESENCHYMAL SUBCLUSTERS — OVERVIEW
# ==============================================================================

levels(mesenchyme)
DimPlot(mesenchyme, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(mesenchyme,
            features = c("Mfap5", "Runx2"),
            pt.size  = 0,
            ncol     = 2)

mesenchyme <- RenameIdents(
  mesenchyme,
  "0" = "Progenitor cells",
  "1" = "Progenitor cells",
  "3" = "Progenitor cells",
  "2" = "Osteogenic cells",
  "4" = "Osteogenic cells",
  "5" = "Osteogenic cells",
  "6" = "Osteogenic cells"
)

# Save new cell type labels to metadata
mesenchyme$celltype <- Idents(mesenchyme)

# Verify
levels(mesenchyme)

DimPlot(mesenchyme, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5) + NoLegend()

# Violin plots for key mesenchymal marker genes
VlnPlot(mesenchyme,
        features = c("Mfap5", "Runx2"),
        pt.size  = 0,
        ncol     = 2)

# ==============================================================================
# 6. CELLCHAT LIGAND–RECEPTOR INTERACTION ANALYSIS
# ==============================================================================

# --- 6a. Create CellChat Object ---
selected_clusters <- c("Mesenchymal cells", "Epithelial cells", "Pericytes",
                       "Endothelial cells", "Chondrocytes",
                       "Neuronal cells",    "Glial cells")

palate_subset <- subset(palate, idents = selected_clusters)
cellchat <- createCellChat(palate_subset, group.by = "ident")

# --- 6b. Load and set mouse ligand–receptor database ---
CellChatDB     <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling",
                           key    = "annotation")
cellchat@DB    <- CellChatDB.use

# --- 6c. Preprocessing ---
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# --- 6d. Infer communication probabilities ---
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract communication table
df.net <- subsetCommunication(cellchat)
head(df.net)

# Summarize at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

# --- 6e. Global interaction overview ---
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Number of interactions"
)

# --- 6f. FGF signaling pathway visualization ---
netVisual_aggregate(cellchat, signaling = "FGF", layout = "circle")
netVisual_aggregate(cellchat, signaling = "FGF", layout = "chord",
                    width = 8, height = 6)
netVisual_heatmap(cellchat, signaling = "FGF",
                  color.heatmap = "Reds", measure = "weight")

# --- 6g. Network centrality analysis ---
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

netAnalysis_signalingRole_network(
  cellchat, signaling = "FGF",
  width = 8, height = 4,
  font.size.title = 20, font.size = 14
)

# Contribution of individual LR pairs to FGF signaling
netAnalysis_contribution(cellchat, signaling = "FGF")

# --- 6h. Ligand–receptor bubble plot ---
pairLR.use <- extractEnrichedLR(cellchat, signaling = "FGF")

netVisual_bubble(
  cellchat,
  pairLR.use      = pairLR.use,
  font.size       = 20,
  font.size.title = 20,
  show.legend     = TRUE,
  grid.on         = TRUE,
  angle.x         = 45
) +
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 80))

# --- 6i. FGF ligand/receptor expression violin plots ---
VlnPlot(palate,
        features = c("Fgf2", "Fgf1", "Fgf7", "Fgf18"),
        pt.size  = 0,
        ncol     = 4) &
  theme(axis.text.x = element_text(size = 9))

VlnPlot(mesenchyme,
        features = c("Fgfr1", "Fgfr2", "Fgfr3"),
        pt.size  = 0,
        ncol     = 3)

# Col2a1 expression across all cell types
VlnPlot(palate,
        features = "Col2a1",
        pt.size  = 0) +
  theme(
    axis.text.x     = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x    = element_text(size = 14),
    legend.position = "none",
    plot.margin     = margin(20, 10, 10, 10)
  )

# Save CellChat object
saveRDS(cellchat, "./data/rds/cellchat.rds")


# ==============================================================================
# 7. SAVE SESSION
# ==============================================================================

save.image(file = "./output/session_downstream_analysis.RData")
sessionInfo()

################################################################################
# END OF SCRIPT
################################################################################
