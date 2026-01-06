#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 2A: DotPlot of major cell types and marker genes
# Reads processed Seurat object (sce_all.rds from config), plots DotPlot with
# celltype on y-axis and marker genes on x-axis.
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

if (is.null(sce@meta.data$celltype)) {
  stop("celltype column is missing in sce_all.rds meta.data")
}

# Marker genes (from scRNA seq code.R)
markers <- c(
  "EPCAM", "KRT19",
  "C1QA", "C1QB", "CD68", "CD163",   # macrophage/monocyte
  "CD79A", "CD79B",                  # B cells
  "CD3D", "CD2",                     # T cells
  "VWF", "PECAM1",                   # endothelial
  "TPSAB1", "TPSB2",                 # mast cells
  "COL1A1", "PDGFRA"                 # fibroblasts
)

cat("Markers used:\n", paste(markers, collapse = ", "), "\n")

p_dot <- DotPlot(sce, features = unique(markers), assay = "RNA", group.by = "celltype") +
  RotatedAxis() +
  scale_x_discrete(limits = markers) +
  theme_bw() +
  ggtitle("Cell type marker expression")

png_out <- file.path(fig_dir, "figures2a_dotplot_celltype_markers.png")
pdf_out <- file.path(fig_dir, "figures2a_dotplot_celltype_markers.pdf")

ggsave(png_out, p_dot, width = 9, height = 6, dpi = 300)
ggsave(pdf_out, p_dot, width = 9, height = 6)

cat("Saved:\n  ", png_out, "\n  ", pdf_out, "\n", sep = "")
