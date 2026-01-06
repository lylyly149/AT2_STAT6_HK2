#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 3C: STAT6 violin plot for fibroblasts
# Uses the upstream processed Seurat object (config raw_data$seurat_all),
# subsets fibroblasts, and outputs only:
#   - STAT6_Fibroblasts_violin_plot.pdf
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(yaml)
  library(limma)
  library(edgeR)
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
violin_fun_path <- cfg$raw_data$violin_plot_functions %||% NA

fig_dir <- file.path(project_root, "figures", "figures3")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object from:", sce_all_path, "\n")
sce_all <- readRDS(sce_all_path)

if (is.null(sce_all@meta.data$celltype)) {
  stop("celltype column missing in sce_all; run upstream annotation first.")
}
if (is.null(sce_all@meta.data$nomipf)) {
  if ("orig.ident" %in% colnames(sce_all@meta.data)) {
    sce_all@meta.data$nomipf <- ifelse(substr(sce_all@meta.data$orig.ident, 1, 3) == "IPF", "IPF", "NOM")
  } else {
    stop("nomipf not found and orig.ident unavailable to derive IPF/NOM status.")
  }
}

# ------------------------------------------------------------------------------
# Subset fibroblasts
# ------------------------------------------------------------------------------
fibro_sce <- subset(sce_all, subset = celltype == "Fibroblasts")
cat("Fibroblast cells:", ncol(fibro_sce), "\n")
cat("By disease:\n")
print(table(fibro_sce$nomipf))

# ------------------------------------------------------------------------------
# STAT6 score preparation
# ------------------------------------------------------------------------------
metadata <- fibro_sce@meta.data
if (!("STAT6GSVA" %in% colnames(metadata))) {
  if (all(c("STAT6PGSVA", "STAT6NGSVA") %in% colnames(metadata))) {
    metadata$STAT6GSVA <- metadata$STAT6PGSVA - metadata$STAT6NGSVA
    fibro_sce@meta.data$STAT6GSVA <- metadata$STAT6GSVA
    cat("Computed STAT6GSVA = STAT6PGSVA - STAT6NGSVA\n")
  } else {
    stop("STAT6GSVA not found and STAT6PGSVA/STAT6NGSVA missing; run upstream STAT6 scoring.")
  }
}

if (file.exists(violin_fun_path)) {
  source(violin_fun_path)
}

# ------------------------------------------------------------------------------
# Violin plot
# ------------------------------------------------------------------------------
n_IPF <- sum(metadata$nomipf == "IPF")
n_NOM <- sum(metadata$nomipf == "NOM")
if (n_IPF < 10 | n_NOM < 10) {
  stop(paste("Insufficient cells for violin plot (IPF =", n_IPF, ", NOM =", n_NOM, ")."))
}

STAT6_2 <- data.frame(
  cell = rownames(metadata),
  score = metadata$STAT6GSVA
)
STAT6_2$score <- scale(STAT6_2$score)[, 1]

vio0 <- data.frame(
  cell = rownames(metadata),
  group = metadata$nomipf
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
       title = "Fibroblasts STAT6 Activity") +
  scale_fill_manual(values = c("IPF" = "#E31A1C", "NOM" = "#1F78B4")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, color = "black", hjust = 0.5)) +
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "wilcox.test",
                     label.y = max(vio3$expression),
                     hide.ns = FALSE)

outfile <- file.path(fig_dir, "figure3c_STAT6_Fibroblasts_violin_plot.pdf")
ggsave(outfile, p, width = 8, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")

cat("\nDone. Output written to", outfile, "\n")
