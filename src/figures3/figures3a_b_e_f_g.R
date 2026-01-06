#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figures 3A/B/E/F/G: Macrophage annotations and STAT6 activity
# Replicates wd/macrophage_annotated_analysis.R but starts from upstream
# sce_all.rds (config raw_data$seurat_all) and only outputs:
#   - Macrophage_UMAP_Annotated.pdf
#   - Macrophage_DotPlot_Annotated.pdf
#   - Macrophage_STAT6_<CellType>_violin_plot.pdf (one per annotated subtype)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggsci)
  library(reshape2)
  library(yaml)
  library(limma)
  library(edgeR)
  library(ggpubr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ------------------------------------------------------------------------------
# Paths and configuration
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$parameters$seed %||% 123)

sce_all_path <- cfg$raw_data$seurat_all %||% stop("config raw_data$seurat_all not set")
stat6_geneset_path <- cfg$raw_data$stat6_geneset %||% stop("config raw_data$stat6_geneset not set")
violin_fun_path <- cfg$raw_data$violin_plot_functions %||% NA

fig_dir <- file.path(project_root, "figures", "figures3")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object from:", sce_all_path, "\n")
sce_all <- readRDS(sce_all_path)

if (is.null(sce_all@meta.data$celltype)) {
  stop("celltype column missing in sce_all; run figure1/processing to annotate cell types.")
}
if (is.null(sce_all@meta.data$nomipf)) {
  if ("orig.ident" %in% colnames(sce_all@meta.data)) {
    sce_all@meta.data$nomipf <- ifelse(substr(sce_all@meta.data$orig.ident, 1, 3) == "IPF", "IPF", "NOM")
  } else {
    stop("nomipf not found and orig.ident unavailable to derive IPF/NOM status.")
  }
}

# ------------------------------------------------------------------------------
# Subset macrophages from upstream sce_all
# ------------------------------------------------------------------------------
macro_cells <- subset(sce_all, subset = celltype == "Macro/mono")
cat("Total macrophage/monocyte cells:", ncol(macro_cells), "\n")
cat("By disease:\n")
print(table(macro_cells$nomipf))

# ------------------------------------------------------------------------------
# Reprocess macrophage subset (same settings as wd/macrophage_STAT6_analysis_v2.R)
# ------------------------------------------------------------------------------
macro_cells <- NormalizeData(macro_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
macro_cells <- FindVariableFeatures(macro_cells, selection.method = "vst", nfeatures = 2000)
if ("percent.mt" %in% colnames(macro_cells@meta.data)) {
  macro_cells <- ScaleData(macro_cells, vars.to.regress = "percent.mt")
} else {
  macro_cells <- ScaleData(macro_cells)
}
macro_cells <- RunPCA(macro_cells, features = VariableFeatures(object = macro_cells))
macro_cells <- FindNeighbors(macro_cells, dims = 1:10)
macro_cells <- FindClusters(macro_cells, resolution = 0.3)
macro_cells <- RunUMAP(macro_cells, dims = 1:10)

cat("Cluster counts:\n")
print(table(macro_cells$seurat_clusters))

# ------------------------------------------------------------------------------
# Manual annotation (copied from wd/macrophage_annotated_analysis.R)
# ------------------------------------------------------------------------------
cluster_annotation <- data.frame(
  cluster = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
  cell_type = c("Alveolar Macrophages",         # 0
                "Alveolar Macrophages",         # 1
                "Interstitial Macrophages",     # 2
                "Alveolar Macrophages",         # 3
                "Alveolar Macrophages",         # 4
                "Monocyte-derived Macrophages", # 5
                "Monocytes",                    # 6
                "Alveolar Macrophages",         # 7
                "Alveolar Macrophages",         # 8
                "Monocyte-derived Macrophages"), # 9
  stringsAsFactors = FALSE
)

macro_cells$annotated_celltype <- cluster_annotation$cell_type[
  match(as.character(macro_cells$seurat_clusters), cluster_annotation$cluster)
]

cat("\nAnnotated cell type distribution:\n")
print(table(macro_cells$annotated_celltype))
cat("\nBy disease status:\n")
print(table(macro_cells$annotated_celltype, macro_cells$nomipf))

# ------------------------------------------------------------------------------
# Marker definitions (copied from wd/macrophage_annotated_analysis.R)
# ------------------------------------------------------------------------------
monocyte_markers <- c("CD14", "FCGR3A", "CCR2", "LYZ", "CSF1R",
                      "S100A8", "S100A9", "CD68")
alveolar_markers <- c("ITGAX", "MARCO", "PPARG", "FABP4", "APOE",
                      "APOC1", "C1QA", "GPNMB", "SIGLEC1")
interstitial_markers <- c("LYVE1", "CD36", "HLA-DRA", "FOLR2", "CCL4",
                          "CXCL10", "CXCL11", "ISG15", "ISG20")
monocyte_derived_markers <- c("CD14", "CSF1R", "LYZ", "CD69", "LGMN")

macro_cells$annotated_celltype <- factor(
  macro_cells$annotated_celltype,
  levels = c("Monocytes",
             "Monocyte-derived Macrophages",
             "Alveolar Macrophages",
             "Interstitial Macrophages")
)

monocyte_unique <- setdiff(monocyte_markers, monocyte_derived_markers)
monocyte_shared <- intersect(monocyte_markers, monocyte_derived_markers)
monocyte_derived_unique <- setdiff(monocyte_derived_markers, monocyte_markers)
ordered_markers <- c(monocyte_unique, monocyte_shared,
                     monocyte_derived_unique, alveolar_markers, interstitial_markers)

celltype_colors <- c(
  "Alveolar Macrophages" = "#FF7F00",
  "Monocytes" = "#E31A1C",
  "Interstitial Macrophages" = "#6A3D9A",
  "Monocyte-derived Macrophages" = "#33A02C"
)

# ------------------------------------------------------------------------------
# Output 1: UMAP annotated
# ------------------------------------------------------------------------------
p_umap_annotated <- DimPlot(macro_cells,
                            reduction = "umap",
                            group.by = "annotated_celltype",
                            cols = celltype_colors,
                            pt.size = 1) +
  ggtitle("Macrophage Subtypes (Annotated)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

umap_out <- file.path(fig_dir, "figures3a_Macrophage_UMAP_Annotated.pdf")
ggsave(umap_out, p_umap_annotated, width = 10, height = 7, dpi = 300)
cat("Saved:", umap_out, "\n")

# ------------------------------------------------------------------------------
# Output 2: DotPlot annotated
# ------------------------------------------------------------------------------
p_dotplot_annotated <- DotPlot(macro_cells,
                               features = ordered_markers,
                               group.by = "annotated_celltype",
                               cols = c("lightgrey", "red"),
                               dot.scale = 10) +
  RotatedAxis() +
  labs(title = "Marker Expression by Annotated Cell Type",
       x = "Marker Genes",
       y = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11),
        legend.position = "right")

dotplot_out <- file.path(fig_dir, "figures3b_Macrophage_DotPlot_Annotated.pdf")
ggsave(dotplot_out, p_dotplot_annotated, width = 14, height = 6, dpi = 300)
cat("Saved:", dotplot_out, "\n")

# ------------------------------------------------------------------------------
# STAT6 violin plots per annotated subtype
# ------------------------------------------------------------------------------
if (file.exists(violin_fun_path)) {
  source(violin_fun_path)
}

metadata <- macro_cells@meta.data
if (!("STAT6GSVA" %in% colnames(metadata))) {
  if (all(c("STAT6PGSVA", "STAT6NGSVA") %in% colnames(metadata))) {
    metadata$STAT6GSVA <- metadata$STAT6PGSVA - metadata$STAT6NGSVA
    macro_cells@meta.data$STAT6GSVA <- metadata$STAT6GSVA
    cat("Computed STAT6GSVA = STAT6PGSVA - STAT6NGSVA\n")
  } else {
    stop("STAT6GSVA not found in metadata and STAT6PGSVA/STAT6NGSVA missing; please run upstream STAT6 scoring (see src/figure1/figure1.R).")
  }
}

prefix_map <- c(
  "Monocytes" = "figures3d",
  "Monocyte-derived Macrophages" = "figures3e",
  "Interstitial Macrophages" = "figures3f",
  "Alveolar Macrophages" = "figures3g"
)

for (ct in levels(macro_cells$annotated_celltype)) {
  cat("\nCreating STAT6 violin for:", ct, "\n")

  ct_cells <- which(macro_cells$annotated_celltype == ct)
  ct_metadata <- metadata[ct_cells, ]

  n_IPF <- sum(ct_metadata$nomipf == "IPF")
  n_NOM <- sum(ct_metadata$nomipf == "NOM")

  if (n_IPF < 10 | n_NOM < 10) {
    cat("  - Skipped (insufficient cells: IPF =", n_IPF, ", NOM =", n_NOM, ")\n")
    next
  }

  STAT6_2 <- data.frame(
    cell = rownames(ct_metadata),
    score = ct_metadata$STAT6GSVA
  )
  STAT6_2$score <- scale(STAT6_2$score)[, 1]

  vio0 <- data.frame(
    cell = rownames(ct_metadata),
    group = ct_metadata$nomipf
  )

  design <- vio0
  colnames(design) <- c("Sample", "Group")
  rownames(design) <- design$Sample
  vio1 <- STAT6_2
  vio1$sample <- vio1$cell
  vio2 <- merge(vio1, design, by.x = "sample", by.y = "Sample")
  vio2 <- vio2[, -2]
  vio2$Group <- as.character(vio2$Group)

  vio3 <- melt(vio2)
  colnames(vio3) <- c("sample", "group", "gene", "expression")

  if (exists("summarySE")) {
    Data_summary <- summarySE(vio3, measurevar = "expression", groupvars = c("group"))
  } else {
    Data_summary <- vio3 %>%
      group_by(group) %>%
      summarise(
        expression = mean(expression),
        ci = sd(expression) / sqrt(n()) * 1.96,
        .groups = "drop"
      )
  }

  p <- ggplot(vio3, aes(x = group, y = expression, fill = group)) +
    geom_violin(trim = FALSE, color = NA, scale = "width", width = 0.5, adjust = 1) +
    geom_errorbar(data = Data_summary, aes(ymin = expression - ci, ymax = expression + ci),
                  width = 0.15,
                  position = position_dodge(0.9),
                  color = "black",
                  alpha = 0.8,
                  size = 0.8) +
    geom_point(data = Data_summary, aes(x = group, y = expression), pch = 19, color = "white",
               position = position_dodge(0.2), size = 1.3) +
    geom_boxplot(width = 0.15, outlier.colour = NA, fill = "white") +
    labs(y = "Standardized STAT6 GSVA Score", x = NULL,
         title = paste(ct, "STAT6 Activity")) +
    scale_fill_manual(values = c("IPF" = "#E31A1C", "NOM" = "#1F78B4")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 12, color = "black", hjust = 0.5)) +
    stat_compare_means(aes(group = group),
                       label = "p.format",
                       method = "wilcox.test",
                       label.y = max(vio3$expression),
                       hide.ns = FALSE)

  ct_safe <- gsub(" ", "_", ct)
  prefix <- prefix_map[[ct]]
  prefix <- if (is.null(prefix)) "figures3x" else prefix
  filename <- file.path(fig_dir, paste0(prefix, "_Macrophage_STAT6_", ct_safe, "_violin_plot.pdf"))
  ggsave(filename, p, width = 6, height = 6, dpi = 300)
  cat("  Saved:", filename, "\n")
}

cat("\nDone. Outputs written to:", fig_dir, "\n")
