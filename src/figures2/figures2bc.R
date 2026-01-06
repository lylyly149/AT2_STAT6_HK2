#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 2B/2C: Epithelial subset tSNE and marker DotPlot
# - Load full Seurat object (sce_all.rds from config)
# - Subset epithelial cells (same clusters as reference pipeline)
# - Re-run normalization, HVG selection, scaling, PCA, neighbors, clusters, tSNE
# - Plot tSNE colored by clusters (prefix figures2b_)
# - Plot DotPlot of selected markers (prefix figures2c_)
# ------------------------------------------------------------------------------

set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Project root and config
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

# Paths
sce_all_path <- cfg$raw_data$seurat_all %||% stop("config raw_data$seurat_all not set")
fig_dir <- file.path(project_root, "figures", "figures2")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object from:", sce_all_path, "\n")
sce <- readRDS(sce_all_path)

if (inherits(sce[["RNA"]], "Assay5")) {
  cat("Joining multi-layer assay into single counts layer for Seurat operations\n")
  sce <- JoinLayers(sce)
}

if (!"seurat_clusters" %in% colnames(sce@meta.data)) {
  stop("seurat_clusters column is missing in sce_all.rds meta.data")
}

# Subset epithelial clusters (as in reference code)
epi_clusters <- c("0", "3")
epi <- subset(sce, seurat_clusters %in% epi_clusters)
cat(sprintf("Epithelial subset contains %d cells across %d genes\n", ncol(epi), nrow(epi)))

# Re-run preprocessing and clustering on epithelial subset
epi <- NormalizeData(epi, normalization.method = "LogNormalize", scale.factor = 1e4)
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)

vars_to_regress <- intersect("percent.mt", colnames(epi@meta.data))
if (length(vars_to_regress) == 0) {
  vars_to_regress <- NULL
}
epi <- ScaleData(epi, vars.to.regress = vars_to_regress)
epi <- RunPCA(epi, features = VariableFeatures(object = epi))
epi <- FindNeighbors(epi, dims = 1:11)
epi <- FindClusters(epi, resolution = 1)
epi <- RunTSNE(epi, dims = 1:11, perplexity = 60)

epi$epi_clusters <- Idents(epi)

# Figure 2B: tSNE colored by epithelial clusters
p_tsne <- DimPlot(epi, reduction = "tsne", group.by = "epi_clusters",
                  label = TRUE, repel = TRUE) +
  ggtitle("tSNE of epithelial subset")

png_tsne <- file.path(fig_dir, "figures2b_tsne_epi_clusters.png")
pdf_tsne <- file.path(fig_dir, "figures2b_tsne_epi_clusters.pdf")
ggsave(png_tsne, p_tsne, width = 8, height = 6, dpi = 300)
ggsave(pdf_tsne, p_tsne, width = 8, height = 6)

# Figure 2C: DotPlot of selected markers
genes_to_check <- c(
  "SFTPC","SFTPA1","SFTPD","AGER","CLDN18","CLIC5","SCGB1A1","MUC5B","PIGR",
  "KRT17","KRT5","TP63","VIM","FN1","ACTG1","CDH1","TM4SF1","FOXJ1"
)

cat("Markers used for DotPlot:\n", paste(genes_to_check, collapse = ", "), "\n")

p_dot <- DotPlot(epi, features = genes_to_check, assay = "RNA", group.by = "epi_clusters") +
  RotatedAxis() +
  scale_x_discrete(limits = genes_to_check) +
  theme_bw() +
  ggtitle("Epithelial cluster marker expression")

png_dot <- file.path(fig_dir, "figures2c_dotplot_epi_markers.png")
pdf_dot <- file.path(fig_dir, "figures2c_dotplot_epi_markers.pdf")
ggsave(png_dot, p_dot, width = 10, height = 6, dpi = 300)
ggsave(pdf_dot, p_dot, width = 10, height = 6)

cat("Saved outputs:\n  ", png_tsne, "\n  ", pdf_tsne, "\n  ", png_dot, "\n  ", pdf_dot, "\n", sep = "")
