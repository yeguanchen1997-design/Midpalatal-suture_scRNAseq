################################################################################
# Single-Cell RNA Sequencing Analysis Pipeline
# Manuscript: Restoring Midpalatal Suture Homeostasis via Local Targeting of
#             Gli1+ Cells Rescues Maxillary Hypoplasia
# Description: Analysis of mouse midpalatal suture scRNA-seq data
#
# NOTE: Raw data preprocessing and QC were performed using Seurat v5.3.0.
# QC parameters applied: 200 < nFeature_RNA < 8,000; nCount_RNA < 60,000;
#                        percent.mt < 12%; FindNeighbors(dims=1:20);
#                        FindClusters(resolution=0.8)
#
# Required input RDS files (place in ./data/rds/):
#   - palate.rds        : Full annotated Seurat object (all cell types)
#   - mesenchyme.rds    : Mesenchymal subcluster Seurat object
#   - cellchat.rds      : CellChat object (ligand-receptor analysis)
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
  options(stringsAsFactors = FALSE)
  options(future.globals.maxSize = 2 * 1024^3)
})

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(monocle3)
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
# 1. LOAD PRE-PROCESSED RDS OBJECTS
# ==============================================================================

palate     <- readRDS("./data/rds/palate.rds")       # Full annotated Seurat object
mesenchyme <- readRDS("./data/rds/mesenchyme.rds")   # Mesenchymal subcluster object
cellchat   <- readRDS("./data/rds/cellchat.rds")     # CellChat object


# ==============================================================================
# 2. FULL PALATE — CELL TYPE OVERVIEW
# ==============================================================================

# Check cell type levels and counts
levels(palate)
table(Idents(palate))

# UMAP of all annotated cell types
DimPlot(palate, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Find marker genes for all cell types
palate.markers <- FindAllMarkers(
  palate,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25
)

# Top 15 markers per cluster by avg_log2FC
top15palate <- palate.markers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

# Feature plots for key cell-type marker genes
marker_genes_palate <- c("Col1a1",  # Mesenchymal marker
                         "Mfap5",   # Mesenchymal marker
                         "Runx2",   # Osteogenic marker
                         "Krt14",   # Epithelial marker
                         "Acta2",   # Pericyte/smooth muscle marker
                         "Tubb3")   # Neuronal marker

for (gene in marker_genes_palate) {
  print(FeaturePlot(palate, features = gene, reduction = "umap"))
}


# ==============================================================================
# 3. MESENCHYMAL SUBCLUSTERS — OVERVIEW
# ==============================================================================

levels(mesenchyme)
DimPlot(mesenchyme, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Find marker genes for mesenchymal subclusters
mesenchyme.markers <- FindAllMarkers(
  mesenchyme,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Top 15 markers per subcluster
top15MSC <- mesenchyme.markers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

# Violin plots for key mesenchymal marker genes across subclusters
VlnPlot(
  mesenchyme,
  features = c("Igfbp2", "Clec3b", "Ccnd1", "Grem1",
               "Sox6",   "Sp7",    "Ibsp",  "Bglap3"),
  pt.size = 0,
  ncol    = 4
)

# Feature plots for mesenchymal subcluster markers
marker_genes_mesen <- c("Igfbp2", "Clec3b", "Ccnd1", "Sox6",
                        "Grem1",  "Sp7",    "Ibsp",
                        "Prrx1",  "Ctsk",   "Axin2", "Bmpr1a", "Gli1")

for (gene in marker_genes_mesen) {
  print(FeaturePlot(mesenchyme, features = gene, reduction = "umap"))
}

# Twist1 in full palate
FeaturePlot(palate, features = "Twist1", reduction = "umap")

# Twist1 violin in mesenchyme
VlnPlot(mesenchyme, features = "Twist1")


# ==============================================================================
# 4. GO ENRICHMENT ANALYSIS — MESENCHYMAL SUBCLUSTERS
# ==============================================================================

# Helper function: find top 100 DEGs and convert to Entrez/Ensembl IDs
get_cluster_genes <- function(obj, cluster_id) {
  markers <- FindMarkers(obj, ident.1 = cluster_id, min.pct = 0.25)
  top100  <- markers %>% top_n(n = -100, wt = p_val_adj) %>% rownames()
  bitr(top100,
       fromType = "SYMBOL",
       toType   = c("ENTREZID", "ENSEMBL"),
       OrgDb    = org.Mm.eg.db)
}

# Run for clusters 0–6
cluster_ids  <- 0:6
gene_dfs     <- lapply(cluster_ids, function(i) get_cluster_genes(mesenchyme, i))
names(gene_dfs) <- paste0("cluster", cluster_ids)

# GO Biological Process enrichment for each cluster
run_enrichGO <- function(gene_df) {
  enrichGO(
    gene          = gene_df$SYMBOL,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    minGSSize     = 10,
    maxGSSize     = 500,
    readable      = FALSE,
    pool          = FALSE
  )
}

go_results  <- lapply(gene_dfs, run_enrichGO)
go_tables   <- lapply(go_results, function(x) x@result)

# Bubble plot function: top 5 GO terms by q-value, x-axis = Rich Factor
erichbgratioplot <- function(data4plot) {
  # Select top 5 terms by q-value
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = FALSE)[1:5], ]

  # Calculate Rich Factor = (DEGs in term) / (total genes in term)
  data4plot$BgRatio <- apply(data4plot, 1, function(x) {
    as.numeric(strsplit(x["GeneRatio"], "/")[[1]][1]) /
    as.numeric(strsplit(x["BgRatio"],   "/")[[1]][1])
  })

  ggplot(data4plot, aes(BgRatio, reorder(Description, BgRatio))) +
    geom_point(aes(size = Count, color = qvalue)) +
    scale_colour_gradient(low = "red", high = "blue") +
    scale_size_continuous(range = c(4, 8)) +
    labs(color = "q-value", size = "Count",
         x = "Rich Factor", y = "") +
    theme_bw() +
    theme(
      axis.text.y  = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_blank()
    ) +
    guides(size  = guide_legend(order = 1),
           color = guide_colorbar(order = 2))
}

# Generate bubble plots for all clusters
for (i in seq_along(go_tables)) {
  cat("Cluster", cluster_ids[i], "\n")
  print(erichbgratioplot(go_tables[[i]]))
}


# ==============================================================================
# 5. PSEUDOTIME TRAJECTORY ANALYSIS (Monocle 3)
# ==============================================================================

library(SeuratWrappers)
library(SingleCellExperiment)

# Build CellDataSet from raw counts to avoid log-normalized data being used
counts_matrix <- GetAssayData(mesenchyme, assay = "RNA", slot = "counts")
cell_metadata <- mesenchyme@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix),
                             row.names = rownames(counts_matrix))

cds <- new_cell_data_set(
  counts_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Transfer Seurat UMAP coordinates to preserve original cell layout
reducedDims(cds)[["UMAP"]] <- Embeddings(mesenchyme, reduction = "umap")

# Preprocess and dimensionality reduction
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP",
                        umap.min_dist = 0.1,
                        preprocess_method = "PCA")

# Re-assign Seurat UMAP coordinates after reduce_dimension
reducedDims(cds)[["UMAP"]] <- Embeddings(mesenchyme, reduction = "umap")

# Cluster cells
cds <- cluster_cells(cds)

# Inspect partitions
plot_cells(cds, color_cells_by = "partition")

# Learn trajectory graph
cds <- learn_graph(cds)

# Visualize trajectory on UMAP
plot_cells(
  cds,
  color_cells_by          = "partition",
  label_groups_by_cluster = FALSE,
  label_leaves            = FALSE,
  label_branch_points     = FALSE
)

# Order cells by pseudotime (interactive: select root node in plot)
cds <- order_cells(cds)

# Pseudotime UMAP
plot_cells(
  cds,
  color_cells_by      = "pseudotime",
  label_cell_groups   = FALSE,
  label_leaves        = FALSE,
  label_branch_points = FALSE,
  graph_label_size    = 1.5
)

# Gene expression dynamics along pseudotime for selected genes
trajectory_genes <- c("Mfap5", "Osr2", "Tnmd", "Dkk2", "Tnc", "Alpl")

# Add gene_short_name to rowData (required for plot_genes_in_pseudotime)
rowData(cds)$gene_short_name <- rownames(rowData(cds))

# Plot in 2 rows x 3 columns without legend
plots <- lapply(trajectory_genes, function(gene) {
  plot_genes_in_pseudotime(
    cds[gene, ],
    color_cells_by = "pseudotime",
    min_expr = 0.5
  ) + theme(legend.position = "none")
})

wrap_plots(plots, nrow = 2, ncol = 3)


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
VlnPlot(
  palate,
  features = c("Fgf2", "Fgf1", "Fgf7", "Fgf18"),
  pt.size  = 0,
  ncol     = 4
) & theme(axis.text.x = element_text(size = 9))

VlnPlot(
  mesenchyme,
  features = c("Fgfr1", "Fgfr2", "Fgfr3"),
  pt.size  = 0,
  ncol     = 3
)

# Col2a1 expression across all cell types
VlnPlot(
  palate,
  features = "Col2a1",
  pt.size  = 0
) +
  theme(
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.x     = element_text(size = 14),
    legend.position  = "none",
    plot.margin      = margin(20, 10, 10, 10)
  )


# ==============================================================================
# 7. SAVE SESSION
# ==============================================================================

save.image(file = "./output/session_downstream_analysis.RData")
sessionInfo()

################################################################################
# END OF SCRIPT
################################################################################
