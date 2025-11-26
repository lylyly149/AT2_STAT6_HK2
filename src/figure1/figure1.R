#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1: scRNA-seq analysis of STAT6 pathway in IPF and healthy lungs
# ==============================================================================
#
# This script generates all panels for Figure 1:
#   - Figure 1A: UMAP showing cell types in IPF vs healthy lung
#   - Figure 1B: tSNE of epithelial cells colored by IPF/NOM and STAT6 GSVA score
#   - Figure 1C: Violin plot of STAT6 GSVA in epithelial cells (IPF vs NOM)
#   - Figure 1D: Bar plot showing STAT6 high/low ratio in epithelial cells
#   - Figure 1E: tSNE of IPF epithelial subclusters
#   - Figure 1F: Dot plot of marker genes in IPF epithelial clusters
#   - Figure 1G: Violin plot of STAT6 GSVA across IPF epithelial clusters
#   - Figure 1H: tSNE colored by AUC score of cluster 0 signature
#   - Figure 1I: Monocle pseudotime trajectory analysis of AT2 cells
#
# Runtime environment: conda activate r_4_3 (R 4.3+)
# ==============================================================================

set.seed(123)

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(GSVA)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(ggpubr)
  library(pals)
  library(AUCell)
  library(viridis)
  library(ggrepel)
  library(monocle)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Get script directory and project root
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

# Load configuration
config <- yaml::read_yaml(file.path(project_root, "config/paths.yaml"))

# Set working directory to project root
setwd(project_root)

# Helper functions for data summary
source(config$raw_data$violin_plot_functions)

# Create output directory
fig_dir <- "figures/figure1"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 1: Load processed Seurat object (required)
# ==============================================================================
sce_all_path <- config$raw_data$seurat_all %||% stop("config$raw_data$seurat_all is not set")

if (!file.exists(sce_all_path)) {
  stop(sprintf("Seurat object not found at %s. Please provide a processed sce_all.rds via config$raw_data$seurat_all.", sce_all_path))
}

cat(sprintf("=== Loading processed Seurat object from %s ===\n", sce_all_path))
sce <- readRDS(sce_all_path)
cat(sprintf("Loaded Seurat object with %d cells\n", ncol(sce)))

required_reductions <- c("umap")
missing_reductions <- setdiff(required_reductions, names(sce@reductions))
if (length(missing_reductions) > 0) {
  stop("Seurat object is missing required reductions: ", paste(missing_reductions, collapse = ", "),
       ". Please ensure the processed object includes these embeddings.")
}

# ==============================================================================
# Step 4: Cell type annotation
# ==============================================================================
cat("=== Cell type annotation ===\n")

sce@meta.data[["nomipf"]] <- ifelse(substr(sce@meta.data[["orig.ident"]], 1, 3) == "IPF", "IPF", "NOM")

# Manual cell type annotation based on cluster markers
celltype <- data.frame(ClusterID = 0:11, celltype = 0:11)
celltype[celltype$ClusterID %in% c("0","3","7"), 2] <- "Epithelial"
celltype[celltype$ClusterID %in% c("1","2","9"), 2] <- "Macro/mono"
celltype[celltype$ClusterID %in% c("4","5"), 2] <- "B cells"
celltype[celltype$ClusterID %in% "6", 2] <- "T cells"
celltype[celltype$ClusterID %in% "8", 2] <- "Endothelial"
celltype[celltype$ClusterID %in% "10", 2] <- "Mast cells"
celltype[celltype$ClusterID %in% "11", 2] <- "Fibroblasts"

sce@meta.data$celltype <- "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]), "celltype"] <- celltype$celltype[i]
}

cat("Cell type distribution:\n")
print(table(sce@meta.data$celltype))

# ==============================================================================
# Figure 1A: UMAP showing cell types (IPF vs NOM)
# ==============================================================================
cat("\n=== Generating Figure 1A: UMAP cell types ===\n")

umap_theme <- theme(panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(size = 1, colour = "black"))

p_nomipf <- DimPlot(sce, reduction = "umap", group.by = "nomipf") +
  scale_color_manual(values = c(IPF = "#ea5514", NOM = "#5aa4ae")) +
  ggtitle("UMAP: IPF vs NOM") +
  umap_theme

p_celltype <- DimPlot(sce, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Cell types") +
  umap_theme

combined_umap <- p_nomipf | p_celltype
ggsave(file.path(fig_dir, "figure1a_umap_nomipf_celltype.png"), combined_umap, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "figure1a_umap_nomipf_celltype.pdf"), combined_umap, width = 12, height = 6)

# ==============================================================================
# Step 5: Calculate STAT6 GSVA scores
# ==============================================================================
cat("\n=== Calculating STAT6 GSVA scores ===\n")

# Handle Assay5 if needed
if (inherits(sce[["RNA"]], "Assay5")) {
  cat("Joining multi-layer assay into single counts layer for GSVA\n")
  sce <- JoinLayers(sce)
}

# Load STAT6 gene set
geneset_file <- config$raw_data$stat6_geneset
sel_gmt <- read.csv(geneset_file, header = T, na.strings = c(" "), check.names = F)
sets <- as.list(sel_gmt)
sets <- lapply(sets, function(x) x[!is.na(x)])

# Calculate GSVA scores
expr <- GetAssayData(sce, layer = "counts")
gsva_matrix <- gsva(as.matrix(expr), sets, kcdf='Poisson', method="ssgsea", min.sz = 10)
gsva_matrix <- t(scale(t(gsva_matrix)))
STAT6GSVA <- gsva_matrix %>% t() %>% as.data.frame()

sce@meta.data$STAT6NGSVA <- STAT6GSVA[,1]
sce@meta.data$STAT6PGSVA <- STAT6GSVA[,2]
sce@meta.data$STAT6GSVA <- STAT6GSVA[,2] - STAT6GSVA[,1]

# Save integrated Seurat object
processed_dir <- "data/processed"
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, file.path(processed_dir, "sce_all.rds"))

# ==============================================================================
# Step 6: Epithelial subset analysis - use existing or create new
# ==============================================================================
cat("\n=== Epithelial subset analysis ===\n")

epi_sce_config <- config$raw_data$seurat_epi
epi_sce_processed <- file.path(processed_dir, "sce_epi.rds")
epi_sce_path <- if (!is.null(epi_sce_config) && file.exists(epi_sce_config)) {
  epi_sce_config
} else if (file.exists(epi_sce_processed)) {
  epi_sce_processed
} else {
  NULL
}

if (!is.null(epi_sce_path)) {
  cat(sprintf("Loading existing epithelial subset from %s\n", epi_sce_path))
  epi_sce <- readRDS(epi_sce_path)
} else {
  cat("Creating epithelial subset from scratch\n")
  epi_sce <- subset(sce, seurat_clusters %in% c("0","3"))
  epi_sce <- NormalizeData(epi_sce, normalization.method = "LogNormalize", scale.factor = 1e4)
  epi_sce <- FindVariableFeatures(epi_sce, selection.method = "vst", nfeatures = 2000)
  epi_sce <- ScaleData(epi_sce, vars.to.regress = "percent.mt")
  epi_sce <- RunPCA(epi_sce, features = VariableFeatures(object = epi_sce))
  epi_sce <- FindNeighbors(epi_sce, dims = 1:11)
  epi_sce <- FindClusters(epi_sce, resolution = 0.9)
  epi_sce <- RunTSNE(epi_sce, dims = 1:11, perplexity = 60)
  saveRDS(epi_sce, epi_sce_processed)
}

if (is.null(epi_sce@meta.data$nomipf)) {
  epi_sce@meta.data[["nomipf"]] <- ifelse(substr(epi_sce@meta.data[["orig.ident"]], 1, 3) == "IPF", "IPF", "NOM")
}

# Calculate GSVA for epithelial subset
if (inherits(epi_sce[["RNA"]], "Assay5")) {
  epi_sce <- JoinLayers(epi_sce)
}
expr_epi <- GetAssayData(epi_sce, layer = "counts")
gsva_epi <- gsva(as.matrix(expr_epi), sets, kcdf='Poisson', method="ssgsea", min.sz = 10)
gsva_epi <- t(scale(t(gsva_epi)))
STAT6GSVA_epi <- gsva_epi %>% t() %>% as.data.frame()

epi_sce@meta.data$STAT6NGSVA <- STAT6GSVA_epi[,1]
epi_sce@meta.data$STAT6PGSVA <- STAT6GSVA_epi[,2]
epi_sce@meta.data$STAT6group <- ifelse(STAT6GSVA_epi[,1] < median(STAT6GSVA_epi[,1]) &
                                         STAT6GSVA_epi[,2] > median(STAT6GSVA_epi[,2]),
                                       "STAT6 high", "STAT6 low")
epi_sce@meta.data$STAT6GSVA <- STAT6GSVA_epi[,2] - STAT6GSVA_epi[,1]

# Save epithelial subset
saveRDS(epi_sce, file.path(processed_dir, "sce_epi.rds"))

# ==============================================================================
# Figure 1B: tSNE of epithelial cells (IPF/NOM and STAT6 GSVA)
# ==============================================================================
cat("\n=== Generating Figure 1B: Epithelial tSNE ===\n")

tsne_pos <- Embeddings(epi_sce, 'tsne')
tsne_theme <- theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"))

# tSNE colored by IPF/NOM
nomdat <- cbind(tsne_pos, data.frame(cell = rownames(epi_sce@meta.data),
                                     phe = epi_sce@meta.data$nomipf))
p_tsne_nomipf <- ggplot(nomdat, aes(x=tSNE_1, y=tSNE_2, color=phe)) +
  geom_point(size=1, alpha=1) +
  scale_color_manual(values = c(IPF = "#ea5514", NOM = "#5aa4ae")) +
  tsne_theme + guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle("tSNE: IPF vs NOM (epithelial)")
ggsave(file.path(fig_dir, "figure1b_tsne_epi_nomipf.png"), p_tsne_nomipf, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1b_tsne_epi_nomipf.pdf"), p_tsne_nomipf, width = 7, height = 5)

# tSNE colored by STAT6 GSVA score
STAT6_2 <- data.frame(cell = rownames(epi_sce@meta.data),
                      score = epi_sce@meta.data$STAT6GSVA)
TSNSTAT6GSVAdata <- cbind(tsne_pos, STAT6_2)
TSNSTAT6GSVAdata$score <- scale(TSNSTAT6GSVAdata$score)
p_tsne_stat6 <- ggplot(TSNSTAT6GSVAdata, aes(x=tSNE_1, y=tSNE_2, color=score)) +
  geom_point(size=0.5, alpha=1) +
  scale_color_gradient2(low = "blue", mid="white", high = "red") +
  tsne_theme +
  ggtitle("tSNE: STAT6 GSVA (epithelial)")
ggsave(file.path(fig_dir, "figure1b_tsne_epi_stat6gsva.png"), p_tsne_stat6, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1b_tsne_epi_stat6gsva.pdf"), p_tsne_stat6, width = 7, height = 5)

# ==============================================================================
# Figure 1C: Violin plot of STAT6 GSVA (epithelial IPF vs NOM)
# ==============================================================================
cat("\n=== Generating Figure 1C: STAT6 GSVA violin plot ===\n")

vio0 <- data.frame(cell = rownames(epi_sce@meta.data),
                   group = epi_sce@meta.data$nomipf)
design <- vio0
colnames(design) <- c('Sample','Group')
rownames(design) <- design$Sample
vio1 <- STAT6_2
vio1$sample <- vio1$cell
vio2 <- merge(vio1, design, by.x = "sample", by.y="Sample")
vio2 <- vio2[,-2]
vio2$Group <- as.character(vio2$Group)
vio3 <- melt(vio2)
colnames(vio3) <- c("sample","group","gene","expression")
vio3$group <- factor(vio3$group, levels = c("IPF","NOM"))
Data_summary <- summarySE(vio3, measurevar="expression", groupvars=c("group"))

mytheme <- theme(plot.title = element_text(size = 12, color="black", hjust = 0.5),
                 axis.title = element_text(size = 12, color ="black"),
                 axis.text = element_text(size= 12, color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

p_violin <- ggplot(vio3, aes(x=group, y=expression, fill=group)) +
  geom_violin(trim= FALSE, color=NA, scale = "width", width=0.5, adjust=1) +
  geom_errorbar(data = Data_summary, aes(ymin = expression-ci, ymax= expression+ci),
                width= 0.15,
                position= position_dodge(0.9),
                color="black",
                alpha = 0.8,
                size= 0.8) +
  geom_point(data = Data_summary, aes(x= group, y= expression), pch=19, color="white",
             position=position_dodge(0.2), size= 1.3) +
  geom_boxplot(width=0.15, outlier.colour = NA, fill="white") +
  labs(y=("Log2 expression"), x=NULL, title = "Violin graph STAT6") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     plot.title = element_text(size = 12, color="black", hjust = 0.5)) +
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "kruskal.test",
                     label.y = max(vio3$expression),
                     hide.ns = TRUE) +
  mytheme
ggsave(file.path(fig_dir, "figure1c_violin_epi_stat6_nomipf.png"), p_violin, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1c_violin_epi_stat6_nomipf.pdf"), p_violin, width = 6, height = 5)

# ==============================================================================
# Figure 1D: Bar plot showing STAT6 high/low ratio
# ==============================================================================
cat("\n=== Generating Figure 1D: STAT6 group ratio bar plot ===\n")

Cellratio <- prop.table(table(epi_sce$STAT6group, epi_sce$nomipf), margin = 2)
Cellratio <- as.data.frame(Cellratio)
p_bar <- ggplot(Cellratio) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = NA) +
  labs(x='Sample', y = 'Ratio') +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"))
ggsave(file.path(fig_dir, "figure1d_bar_epi_stat6group_nomipf.png"), p_bar, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1d_bar_epi_stat6group_nomipf.pdf"), p_bar, width = 6, height = 5)

# ==============================================================================
# Step 7: IPF epithelial re-clustering - use existing or create new
# ==============================================================================
cat("\n=== IPF epithelial re-clustering ===\n")

epi_ipf_config <- config$raw_data$seurat_ipf
epi_ipf_processed <- file.path(processed_dir, "sce_ipf.rds")
epi_ipf_path <- if (!is.null(epi_ipf_config) && file.exists(epi_ipf_config)) {
  epi_ipf_config
} else if (file.exists(epi_ipf_processed)) {
  epi_ipf_processed
} else {
  NULL
}

if (!is.null(epi_ipf_path)) {
  cat(sprintf("Loading existing IPF epithelial subset from %s\n", epi_ipf_path))
  epi_ipf <- readRDS(epi_ipf_path)
  # Ensure tSNE reduction exists
  if (is.null(epi_ipf@reductions$tsne)) {
    cat("Running tSNE for IPF epithelial subset\n")
    epi_ipf <- RunTSNE(epi_ipf, dims=1:10, perplexity=80)
  }
} else {
  cat("Creating IPF epithelial subset from scratch\n")
  epi_ipf <- subset(epi_sce, nomipf == "IPF")
  epi_ipf <- NormalizeData(epi_ipf, normalization.method = "LogNormalize", scale.factor = 1e4)
  epi_ipf <- FindVariableFeatures(epi_ipf, selection.method = 'vst', nfeatures = 2000)
  epi_ipf <- ScaleData(epi_ipf, vars.to.regress = "percent.mt")
  epi_ipf <- RunPCA(epi_ipf, features = VariableFeatures(object = epi_ipf))
  epi_ipf <- FindNeighbors(epi_ipf, dims = 1:10)
  epi_ipf <- FindClusters(epi_ipf, resolution = 0.9)
  epi_ipf <- RunTSNE(epi_ipf, dims=1:10, perplexity=80)
  saveRDS(epi_ipf, epi_ipf_processed)
}

# ==============================================================================
# Figure 1E: tSNE of IPF epithelial clusters
# ==============================================================================
cat("\n=== Generating Figure 1E: IPF epithelial clusters ===\n")

tsne_pos_ipf <- Embeddings(epi_ipf, 'tsne')
phe <- data.frame(cell = rownames(epi_ipf@meta.data),
                  cluster = epi_ipf@meta.data$seurat_clusters)
dat_ipf <- cbind(tsne_pos_ipf, phe)
cluster_levels <- levels(epi_ipf$seurat_clusters)
base_cols <- pals::glasbey()[c(1:3,5,7,11,13,16,14)]
if (length(cluster_levels) > length(base_cols)) {
  base_cols <- pals::glasbey()[seq_len(length(cluster_levels))]
}
define_color1 <- stats::setNames(base_cols[seq_len(length(cluster_levels))], cluster_levels)
centers_ipf <- aggregate(tsne_pos_ipf, by=list(cluster=epi_ipf$seurat_clusters), FUN=median)
p_tsne_ipf <- ggplot(dat_ipf, aes(x=tSNE_1, y=tSNE_2, color=cluster)) +
  geom_point(size=1.5, alpha=1) +
  scale_color_manual(values = define_color1, drop = FALSE) +
  geom_text_repel(data=centers_ipf, aes(x=tSNE_1, y=tSNE_2, label=cluster),
                  color="black", size=3, inherit.aes = FALSE, max.overlaps = Inf) +
  tsne_theme + guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle("tSNE: IPF epithelial clusters")
ggsave(file.path(fig_dir, "figure1e_tsne_ipf_epi_clusters.png"), p_tsne_ipf, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1e_tsne_ipf_epi_clusters.pdf"), p_tsne_ipf, width = 7, height = 5)

# ==============================================================================
# Figure 1F: Dot plot of marker genes
# ==============================================================================
cat("\n=== Generating Figure 1F: Marker genes dot plot ===\n")

genes_to_check_ipf <- c("SFTPC","SFTPA1","SFTPD","AGER","CLDN18","CLIC5",
                        "SCGB1A1","MUC5B","PIGR","KRT17","KRT5","TP63",
                        "VIM","FN1","ACTG1","CDH1","TM4SF1")
p_dot_ipf <- DotPlot(epi_ipf, features=unique(unlist(genes_to_check_ipf)), assay="RNA") + RotatedAxis()
ggsave(file.path(fig_dir, "figure1f_dot_ipf_epi_markers.png"), p_dot_ipf, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "figure1f_dot_ipf_epi_markers.pdf"), p_dot_ipf, width = 8, height = 6)

# ==============================================================================
# Figure 1G: Violin plot of STAT6 GSVA across IPF clusters
# ==============================================================================
cat("\n=== Generating Figure 1G: STAT6 GSVA by IPF cluster ===\n")

vio0_ipf <- data.frame(cell = rownames(epi_ipf@meta.data),
                       group = epi_ipf@meta.data$seurat_clusters)
design_ipf <- vio0_ipf
colnames(design_ipf) <- c('Sample','Group')
rownames(design_ipf) <- design_ipf$Sample
STAT6_2_ipf <- data.frame(cell = rownames(epi_ipf@meta.data),
                          score = epi_ipf@meta.data$STAT6GSVA)
vio1_ipf <- STAT6_2_ipf
vio1_ipf$sample <- vio1_ipf$cell
vio2_ipf <- merge(vio1_ipf, design_ipf, by.x = "sample", by.y="Sample")
vio2_ipf <- vio2_ipf[,-2]
vio2_ipf$Group <- as.character(vio2_ipf$Group)
vio3_ipf <- melt(vio2_ipf)
colnames(vio3_ipf) <- c("sample","group","gene","expression")
vio3_ipf$group <- factor(vio3_ipf$group)
Data_summary_ipf <- summarySE(vio3_ipf, measurevar="expression", groupvars=c("group"))
p_violin_ipf <- ggplot(vio3_ipf, aes(x=group, y=expression, fill=group)) +
  geom_violin(trim= FALSE, color=NA, scale = "width", width=0.5, adjust=1) +
  geom_errorbar(data = Data_summary_ipf, aes(ymin = expression-ci, ymax= expression+ci),
                width= 0.15,
                position= position_dodge(0.9),
                color="black",
                alpha = 0.8,
                size= 0.8) +
  geom_point(data = Data_summary_ipf, aes(x= group, y= expression), pch=19, color="white",
             position=position_dodge(0.2), size= 1.3) +
  geom_boxplot(width=0.15, outlier.colour = NA, fill="white") +
  labs(y=("Log2 expression"), x=NULL, title = "Violin graph STAT6 by cluster") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     plot.title = element_text(size = 12, color="black", hjust = 0.5)) +
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "kruskal.test",
                     label.y = max(vio3_ipf$expression),
                     hide.ns = TRUE) +
  mytheme
ggsave(file.path(fig_dir, "figure1g_violin_ipf_epi_stat6_cluster.png"), p_violin_ipf, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1g_violin_ipf_epi_stat6_cluster.pdf"), p_violin_ipf, width = 7, height = 5)

# ==============================================================================
# Figure 1H: tSNE colored by AUC score of cluster 0 signature
# ==============================================================================
cat("\n=== Generating Figure 1H: AUC score visualization ===\n")

if (inherits(epi_ipf[["RNA"]], "Assay5")) {
  epi_ipf <- JoinLayers(epi_ipf)
}
expr_auc <- GetAssayData(epi_ipf, layer = "counts")
# Disable the default QC histogram to avoid generating Rplots.pdf when running non-interactively
cells_rankings <- AUCell_buildRankings(expr_auc, plotStats = FALSE)
markers_ipf <- FindAllMarkers(epi_ipf)
marker0 <- markers_ipf[markers_ipf$cluster == 0 & markers_ipf$p_val_adj < 0.05, ]
gene_set <- unique(marker0$gene)
cells_AUC <- AUCell_calcAUC(list(geneSet = gene_set), cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
aucs <- as.numeric(getAUC(cells_AUC)["geneSet", ])
epi_ipf$AUC <- aucs
tsne_auc_df <- data.frame(epi_ipf@meta.data, Embeddings(epi_ipf, 'tsne'))
p_auc <- ggplot(tsne_auc_df, aes(tSNE_1, tSNE_2, color=AUC)) +
  geom_point(size=1.5) +
  scale_color_viridis(option="A") +
  tsne_theme +
  ggtitle("tSNE: Cluster0 signature AUC")
ggsave(file.path(fig_dir, "figure1h_tsne_ipf_epi_auc.png"), p_auc, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "figure1h_tsne_ipf_epi_auc.pdf"), p_auc, width = 7, height = 5)

# ==============================================================================
# Figure 1I: Monocle pseudotime trajectory analysis
# ==============================================================================
cat("\n=== Generating Figure 1I: Monocle pseudotime ===\n")

# Load existing AT2 subset or create new
at2_config <- config$raw_data$seurat_at2
at2_processed <- file.path(processed_dir, "sce_at2.rds")
at2_path <- if (!is.null(at2_config) && file.exists(at2_config)) {
  at2_config
} else if (file.exists(at2_processed)) {
  at2_processed
} else {
  NULL
}

if (!is.null(at2_path)) {
  cat(sprintf("Loading existing AT2 subset from %s\n", at2_path))
  sce_at2 <- readRDS(at2_path)
} else {
  cat("Creating AT2 subset from epithelial cells\n")
  sce_at2 <- subset(epi_sce, seurat_clusters == "0")
  saveRDS(sce_at2, at2_processed)
}

tryCatch({
  # Handle Assay5 format for monocle
  if (inherits(sce_at2[["RNA"]], "Assay5")) {
    cat("Converting Assay5 to Assay for monocle compatibility\n")
    sce_at2[["RNA"]] <- as(object = sce_at2[["RNA"]], Class = "Assay")
  }

  # Extract counts matrix
  if (!is.null(sce_at2@assays$RNA@counts)) {
    data_mat <- as(as.matrix(sce_at2@assays$RNA@counts), "sparseMatrix")
  } else if (!is.null(sce_at2@assays$RNA@data)) {
    data_mat <- as(as.matrix(sce_at2@assays$RNA@data), "sparseMatrix")
  } else {
    stop("Cannot find counts or data slot in RNA assay")
  }

  pd <- new("AnnotatedDataFrame", data = sce_at2@meta.data)
  fData <- data.frame(gene_short_name = row.names(data_mat), row.names = row.names(data_mat))
  fd <- new("AnnotatedDataFrame", data = fData)
  monocle_cds <- newCellDataSet(data_mat,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.01,
                                expressionFamily = negbinomial.size())
  Mono.cds <- estimateSizeFactors(monocle_cds)
  Mono.cds <- estimateDispersions(Mono.cds)
  disp_table <- dispersionTable(Mono.cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)
  Mono.cds <- reduceDimension(Mono.cds, max_components = 2, method = 'DDRTree')
  Mono.cds <- orderCells(Mono.cds)
  Mono.cds <- orderCells(Mono.cds, root_state = 6)

  p_pseudo <- plot_cell_trajectory(Mono.cds, color_by = "Pseudotime")
  ggsave(file.path(fig_dir, "figure1i_monocle_pseudotime.png"), p_pseudo, width = 6, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "figure1i_monocle_pseudotime.pdf"), p_pseudo, width = 6, height = 5)

  traj_plot <- plot_cell_trajectory(
    Mono.cds,
    color_by = "STAT6GSVA",
    cell_size = 1.5,
    alpha = 0.8
  ) +
    facet_wrap(~nomipf, nrow = 1) +
    scale_color_gradientn(
      colours = c("blue","white","red"),
      breaks = c(-4,2,8),
      name = "STAT6 Score"
    )
  ggsave(file.path(fig_dir, "figure1i_monocle_stat6gsva.png"), traj_plot, width = 8, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "figure1i_monocle_stat6gsva.pdf"), traj_plot, width = 8, height = 4)
}, error = function(e) {
  cat("Monocle pseudotime failed:", conditionMessage(e), "\n")
})

cat("\n=== Figure 1 analysis complete ===\n")
cat("All figures saved to:", fig_dir, "\n")
